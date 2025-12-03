! This file is part of xtb.
!
! Lightweight runtime control for threading across SCC and Hessian
! calculations. It keeps a process-wide preferred thread count that is
! derived from simple throughput measurements and reapplied to OpenMP
! and threaded BLAS backends.
module xtb_parallel_threads
   use xtb_mctc_accuracy, only : wp
   use, intrinsic :: iso_c_binding, only : c_char, c_int, c_null_char
   use omp_lib
   implicit none
   private

   integer, save :: preferred_scc_threads = 0
   integer, save :: baseline_threads = 0

   public :: xtb_get_preferred_threads
    ! track total thread budget even when per-SCC preference lowers omp defaults
   public :: xtb_get_max_available_threads
   public :: xtb_apply_thread_preference
   public :: xtb_record_thread_preference

   interface
      function c_setenv(name, value, overwrite) bind(C, name="setenv")
         import :: c_char, c_int
         character(c_char), dimension(*), intent(in) :: name
         character(c_char), dimension(*), intent(in) :: value
         integer(c_int), value :: overwrite
         integer(c_int) :: c_setenv
      end function c_setenv
   end interface

contains

   !> Return the preferred SCC thread count, clamped to what is available
   integer function xtb_get_preferred_threads(max_threads) result(n_threads)
      integer, intent(in) :: max_threads

      if (baseline_threads == 0) baseline_threads = omp_get_max_threads()
      n_threads = max_threads
      if (preferred_scc_threads > 0) &
         & n_threads = min(max_threads, preferred_scc_threads)
      n_threads = max(1, n_threads)
   end function xtb_get_preferred_threads


   !> Return best guess of the total threads available (before any throttling)
   integer function xtb_get_max_available_threads() result(n_threads)
      if (baseline_threads == 0) baseline_threads = omp_get_max_threads()
      n_threads = max(1, max(baseline_threads, omp_get_max_threads()))
   end function xtb_get_max_available_threads


   !> Apply the persisted preference to OpenMP and BLAS env vars
   subroutine xtb_apply_thread_preference(max_threads)
      integer, intent(in) :: max_threads
      call xtb_set_threads(xtb_get_preferred_threads(max_threads))
   end subroutine xtb_apply_thread_preference


   !> Persist a newly measured preference when it improved throughput
   subroutine xtb_record_thread_preference(max_threads, threads_per_task, &
         & baseline_tp, best_tp, verbose, unit)
      integer, intent(in) :: max_threads
      integer, intent(in) :: threads_per_task
      real(wp), intent(in) :: baseline_tp, best_tp
      logical, intent(in) :: verbose
      integer, intent(in) :: unit

      real(wp) :: gain
      integer  :: chosen

      if (threads_per_task <= 0) return
      if (best_tp <= 0.0_wp .or. baseline_tp <= 0.0_wp) return

      gain = best_tp / baseline_tp
      chosen = min(max_threads, threads_per_task)

      if (gain > 1.05_wp .or. preferred_scc_threads == 0) then
         preferred_scc_threads = max(1, chosen)
         call xtb_set_threads(preferred_scc_threads)
         if (verbose .and. unit > 0) then
            write(unit,'(a,i0,a,f6.2,a)') &
               "Thread tuner: using ", preferred_scc_threads, &
               " threads per SCC (", gain*100.0_wp, "% of baseline)."
         end if
      end if
   end subroutine xtb_record_thread_preference


   subroutine xtb_set_threads(n_threads)
      integer, intent(in) :: n_threads
      integer :: n_clamped

      if (baseline_threads == 0) baseline_threads = omp_get_max_threads()
      n_clamped = max(1, n_threads)
      call omp_set_num_threads(n_clamped)
      call xtb_setenv_int("OMP_NUM_THREADS", n_clamped)
      call xtb_setenv_int("MKL_NUM_THREADS", n_clamped)
      call xtb_setenv_int("OPENBLAS_NUM_THREADS", n_clamped)
   end subroutine xtb_set_threads


   subroutine xtb_setenv_int(name, value)
      character(len=*), intent(in) :: name
      integer, intent(in) :: value
      character(len=32) :: buff
      integer(c_int) :: ierr

      if (value <= 0) return
      write(buff,'(i0)') value
      ierr = c_setenv(trim(name)//c_null_char, trim(buff)//c_null_char, 1_c_int)
      if (ierr /= 0_c_int) then
         ! silently ignore failures, environment handling can be platform-specific
      end if
   end subroutine xtb_setenv_int

end module xtb_parallel_threads
