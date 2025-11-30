! This file is part of xtb.
!
! Threading policy helpers to balance SCC OpenMP regions and BLAS backends.
module xtb_threading_policy
   use iso_c_binding,  only : c_int
   use xtb_mctc_accuracy, only : wp
#ifdef _OPENMP
   use omp_lib,        only : omp_get_max_threads, omp_set_num_threads, &
      &                      omp_get_max_active_levels, omp_set_max_active_levels
#endif
   implicit none
   private

   integer, parameter :: dim_threshold_default = 800
   real(wp), parameter :: work_threshold_default = 1.0e8_wp

   type :: ThreadingPolicy
      integer :: omp_threads_before         = 1
      integer :: blas_threads_before        = 1
      integer :: max_active_levels_before   = 1
      logical :: omp_changed                = .false.
      logical :: blas_changed               = .false.
      logical :: max_active_levels_changed  = .false.
   end type ThreadingPolicy

   public :: ThreadingPolicy
   public :: setup_scc_thread_policy
   public :: restore_thread_policy
   public :: should_use_scc_parallel

#if defined XTB_LAPACK_MKL
   interface
      integer(c_int) function mkl_get_max_threads() bind(C, name='mkl_get_max_threads')
         import :: c_int
      end function mkl_get_max_threads
      subroutine mkl_set_num_threads(nthreads) bind(C, name='mkl_set_num_threads')
         import :: c_int
         integer(c_int), value :: nthreads
      end subroutine mkl_set_num_threads
   end interface
#elif defined XTB_LAPACK_OPENBLAS
   interface
      subroutine openblas_set_num_threads(nthreads) bind(C, name='openblas_set_num_threads')
         import :: c_int
         integer(c_int), value :: nthreads
      end subroutine openblas_set_num_threads
      integer(c_int) function openblas_get_num_threads() bind(C, name='openblas_get_num_threads')
         import :: c_int
      end function openblas_get_num_threads
   end interface
#endif

contains

function getenv_int(name, default) result(value)
   character(len=*), intent(in) :: name
   integer,          intent(in) :: default
   integer                      :: value
   integer :: stat, lenv
   character(len=64) :: buffer

   value = default
   call get_environment_variable(name, buffer, length=lenv, status=stat)
   if (stat == 0 .and. lenv > 0) then
      read(buffer(1:lenv), *, iostat=stat) value
      if (stat /= 0 .or. value < 1) value = default
   end if
end function getenv_int

function getenv_real(name, default) result(value)
   character(len=*), intent(in) :: name
   real(wp),         intent(in) :: default
   real(wp)                      :: value
   integer :: stat, lenv
   character(len=64) :: buffer

   value = default
   call get_environment_variable(name, buffer, length=lenv, status=stat)
   if (stat == 0 .and. lenv > 0) then
      read(buffer(1:lenv), *, iostat=stat) value
      if (stat /= 0 .or. value <= 0.0_wp) value = default
   end if
end function getenv_real

logical function should_use_scc_parallel(ndim) result(use_parallel)
   integer, intent(in) :: ndim
   integer  :: dim_threshold
   real(wp) :: work_threshold
   real(wp) :: work_est

   dim_threshold   = getenv_int('XTB_SCC_DIM_THRESHOLD', dim_threshold_default)
   work_threshold  = getenv_real('XTB_SCC_WORK_THRESHOLD', work_threshold_default)
   work_est        = real(ndim, wp) * real(ndim, wp) * real(ndim, wp)

   use_parallel = .not. (ndim >= dim_threshold .or. work_est >= work_threshold)
end function should_use_scc_parallel

subroutine setup_scc_thread_policy(ndim, policy, use_parallel, blas_threads_target)
   integer, intent(in)              :: ndim
   type(ThreadingPolicy), intent(inout) :: policy
   logical, intent(out)             :: use_parallel
   integer, intent(out)             :: blas_threads_target

   integer :: outer_threads

   outer_threads = current_omp_threads()
   if (outer_threads < 1) outer_threads = 1

   use_parallel = should_use_scc_parallel(ndim)
#if defined XTB_LAPACK_MKL
   if (use_parallel) then
      blas_threads_target = 1
   else
      blas_threads_target = outer_threads
   end if
#elif defined XTB_LAPACK_OPENBLAS
   if (use_parallel) then
      blas_threads_target = 1
   else
      blas_threads_target = outer_threads
   end if
#elif defined XTB_LAPACK_CUSTOM
   blas_threads_target = outer_threads
#else
   if (use_parallel) then
      blas_threads_target = 1
   else
      blas_threads_target = outer_threads
   end if
#endif

   policy%omp_threads_before = outer_threads
   policy%max_active_levels_before = current_max_active_levels()
   call set_max_active_levels(1, policy)
   call set_blas_threads(blas_threads_target, policy)
end subroutine setup_scc_thread_policy

subroutine restore_thread_policy(policy)
   type(ThreadingPolicy), intent(inout) :: policy

   if (policy%blas_changed) then
      call reset_blas_threads(policy)
      policy%blas_changed = .false.
   end if

   if (policy%max_active_levels_changed) then
      call reset_max_active_levels(policy)
      policy%max_active_levels_changed = .false.
   end if
end subroutine restore_thread_policy

subroutine set_max_active_levels(levels, policy)
   integer, intent(in)              :: levels
   type(ThreadingPolicy), intent(inout) :: policy
#ifdef _OPENMP
#if _OPENMP >= 200805
   if (levels /= policy%max_active_levels_before) then
      call omp_set_max_active_levels(levels)
      policy%max_active_levels_changed = .true.
   end if
#else
   policy%max_active_levels_changed = .false.
#endif
#else
   policy%max_active_levels_changed = .false.
#endif
end subroutine set_max_active_levels

subroutine reset_max_active_levels(policy)
   type(ThreadingPolicy), intent(in) :: policy
#ifdef _OPENMP
#if _OPENMP >= 200805
   call omp_set_max_active_levels(policy%max_active_levels_before)
#endif
#endif
end subroutine reset_max_active_levels

subroutine set_blas_threads(nthreads, policy)
   integer, intent(in)              :: nthreads
   type(ThreadingPolicy), intent(inout) :: policy

   policy%blas_threads_before = 1

#if defined XTB_LAPACK_MKL
   policy%blas_threads_before = mkl_get_max_threads()
   call mkl_set_num_threads(int(nthreads, c_int))
   policy%blas_changed = .true.
#elif defined XTB_LAPACK_OPENBLAS
   policy%blas_threads_before = openblas_get_num_threads()
   call openblas_set_num_threads(int(nthreads, c_int))
   policy%blas_changed = .true.
#elif defined _OPENMP
   policy%blas_threads_before = current_omp_threads()
   call omp_set_num_threads(nthreads)
   policy%blas_changed = .true.
   policy%omp_changed  = .true.
#else
   policy%blas_changed = .false.
#endif
end subroutine set_blas_threads

subroutine reset_blas_threads(policy)
   type(ThreadingPolicy), intent(in) :: policy

#if defined XTB_LAPACK_MKL
   call mkl_set_num_threads(int(policy%blas_threads_before, c_int))
#elif defined XTB_LAPACK_OPENBLAS
   call openblas_set_num_threads(int(policy%blas_threads_before, c_int))
#elif defined _OPENMP
   if (policy%omp_changed) then
      call omp_set_num_threads(policy%omp_threads_before)
   end if
#endif
end subroutine reset_blas_threads

integer function current_omp_threads() result(nthreads)
#ifdef _OPENMP
   nthreads = omp_get_max_threads()
#else
   nthreads = 1
#endif
end function current_omp_threads

integer function current_max_active_levels() result(levels)
#ifdef _OPENMP
#if _OPENMP >= 200805
   levels = omp_get_max_active_levels()
#else
   levels = 1
#endif
#else
   levels = 1
#endif
end function current_max_active_levels

end module xtb_threading_policy
