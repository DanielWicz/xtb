! This file is part of xtb.
!
! Helper utilities to sample the current resident set size (RSS) from
! /proc/self/status. Designed for debugging memory growth during
! geometry optimizations; guarded by the XTB_MEMLOG environment
! variable to avoid overhead in normal production runs.
module xtb_mctc_meminfo
   use iso_fortran_env, only : int64
#ifdef WITH_MKL
   ! mkl_service.mod is not always available; declare C bindings explicitly
   use iso_c_binding, only : c_int
   interface
      integer(c_int) function mkl_free_buffers() bind(C, name="mkl_free_buffers")
      end function mkl_free_buffers
      integer(c_int) function mkl_thread_free_buffers() bind(C, name="mkl_thread_free_buffers")
      end function mkl_thread_free_buffers
      integer(c_int) function mkl_disable_fast_mm() bind(C, name="mkl_disable_fast_mm")
      end function mkl_disable_fast_mm
   end interface
#endif

#ifdef WITH_OPENBLAS
   use iso_c_binding, only : c_int
   interface
      subroutine openblas_thread_cleanup() bind(C, name="openblas_thread_cleanup")
      end subroutine openblas_thread_cleanup
      subroutine openblas_set_num_threads(n) bind(C, name="openblas_set_num_threads")
         import :: c_int
         integer(c_int), value :: n
      end subroutine openblas_set_num_threads
   end interface
#endif
   implicit none
   private

   public :: rss_kb, memlog_enabled, log_memory_usage, log_memory_usage_delta
   public :: trim_memory, init_memory_management

#ifdef WITH_MKL
   logical, save :: mkl_fastmm_disabled = .false.
#endif

contains

!> One-time memory backend setup (currently disables MKL fast memory manager
!> which otherwise hoards large thread-local buffers across SCF cycles).
subroutine init_memory_management()
#ifdef WITH_MKL
   if (.not. mkl_fastmm_disabled) then
      call mkl_disable_fast_mm()
      mkl_fastmm_disabled = .true.
   end if
#endif
end subroutine init_memory_management

   !> Returns .true. once XTB_MEMLOG is set in the environment.
   logical function memlog_enabled()
      logical, save :: initialized = .false.
      logical, save :: enabled = .false.
      integer :: len, stat

      if (.not.initialized) then
         call get_environment_variable('XTB_MEMLOG', length=len, status=stat)
         enabled = (stat == 0 .and. len > 0)
         initialized = .true.
      end if

      memlog_enabled = enabled
   end function memlog_enabled

   !> Read VmRSS from /proc/self/status. Returns -1 on failure.
   integer(int64) function rss_kb()
      integer :: unit, ios, pos
      character(len=256) :: line
      character(len=32) :: suffix

      rss_kb = -1_int64
      open(newunit=unit, file='/proc/self/status', status='old', action='read', &
         & iostat=ios)
      if (ios /= 0) return

      do
         read(unit,'(A)', iostat=ios) line
         if (ios /= 0) exit
         if (line(1:5) == 'VmRSS') then
            pos = index(line, ':')
            if (pos > 0) then
               read(line(pos+1:), *, iostat=ios) rss_kb, suffix
            end if
            exit
         end if
      end do

      close(unit)
      if (ios /= 0) rss_kb = -1_int64
   end function rss_kb

   !> Conditional RSS logger (enabled via XTB_MEMLOG).
   subroutine log_memory_usage(unit, label)
      integer, intent(in) :: unit
      character(len=*), intent(in) :: label
      integer(int64) :: rss

      if (.not.memlog_enabled()) return
      rss = rss_kb()
      if (rss >= 0_int64) then
         write(unit,'(1x,"[mem]",1x,a,1x,i0," kB")') trim(label), rss
      end if
   end subroutine log_memory_usage

   !> Same as log_memory_usage, but also reports delta to previous sample.
   subroutine log_memory_usage_delta(unit, label, last_rss)
      integer, intent(in) :: unit
      character(len=*), intent(in) :: label
      integer(int64), intent(inout) :: last_rss
      integer(int64) :: rss, delta

      if (.not.memlog_enabled()) return
      rss = rss_kb()
      if (rss < 0_int64) return
      if (last_rss < 0_int64) then
         write(unit,'(1x,"[mem]",1x,a,1x,i0," kB")') trim(label), rss
      else
         delta = rss - last_rss
         write(unit,'(1x,"[mem]",1x,a,1x,i0," kB",1x,"(Δ",i0," kB)")') &
            trim(label), rss, delta
      end if
      last_rss = rss
end subroutine log_memory_usage_delta

!> Try to release unused heap pages back to the OS (requires glibc).
subroutine trim_memory()
   use iso_c_binding, only : c_int, c_size_t
   use omp_lib,        only : omp_get_max_threads
   implicit none
   interface
      function c_malloc_trim(pad) bind(C, name="malloc_trim")
         import :: c_int, c_size_t
         integer(c_size_t), value :: pad
         integer(c_int) :: c_malloc_trim
      end function c_malloc_trim
   end interface
   integer(c_int) :: ierr
#ifdef WITH_OPENBLAS
   integer :: nthreads
   ! OpenBLAS keeps per-thread caches; force cleanup and collapse to one thread
   call openblas_thread_cleanup()
   nthreads = omp_get_max_threads()
   call openblas_set_num_threads(1_c_int)
#endif
#ifdef WITH_MKL
   ! release internal MKL thread buffers that otherwise accumulate across
   ! repeated SCF/geometry steps when using the Intel/oneMKL backend
   call mkl_thread_free_buffers()
   call mkl_free_buffers()
#endif
   ierr = c_malloc_trim(0_c_size_t)
#ifdef WITH_OPENBLAS
   if (nthreads > 1) call openblas_set_num_threads(int(nthreads, c_int))
#endif
end subroutine trim_memory

end module xtb_mctc_meminfo
