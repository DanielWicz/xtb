! This file is part of xtb.
!
! Threading policy helpers to balance SCC OpenMP regions and BLAS backends.
module xtb_threading_policy
   use iso_c_binding,  only : c_int
   use xtb_mctc_accuracy, only : wp
#ifdef _OPENMP
   use omp_lib,        only : omp_get_max_threads, omp_set_num_threads, &
      &                      omp_get_max_active_levels, omp_set_max_active_levels, &
      &                      omp_get_num_procs
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
   public :: getenv_int, getenv_real

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
   logical  :: force_parallel
   logical  :: force_serial
   integer  :: outer_threads

   dim_threshold   = getenv_int('XTB_SCC_DIM_THRESHOLD', dim_threshold_default)
   work_threshold  = getenv_real('XTB_SCC_WORK_THRESHOLD', work_threshold_default)
   force_parallel  = getenv_int('XTB_SCC_FORCE_PARALLEL', 0) /= 0
   force_serial    = getenv_int('XTB_SCC_FORCE_SERIAL', 0) /= 0
   outer_threads   = current_omp_threads()

   if (force_parallel) then
      use_parallel = outer_threads > 1
      return
   end if

   if (force_serial) then
      use_parallel = .false.
      return
   end if

   if (outer_threads < 2) then
      use_parallel = .false.
      return
   end if

   work_est        = real(ndim, wp) * real(ndim, wp) * real(ndim, wp)

   use_parallel = .not. (ndim >= dim_threshold .or. work_est >= work_threshold)
end function should_use_scc_parallel

!> Decide how many BLAS threads to use for SCC builds.
!  Honors optional caps to keep work inside a socket and avoid oversubscription.
integer function select_blas_threads(use_parallel, outer_threads, allow_nested) result(target)
   logical, intent(in) :: use_parallel
   integer, intent(in) :: outer_threads
   logical, intent(in) :: allow_nested

   integer :: force_blas
   integer :: min_blas
   integer :: max_blas
   integer :: socket_cap
   integer :: hardware_cap
   integer :: user_cap

   hardware_cap = hardware_concurrency()
   socket_cap   = getenv_int('XTB_SCC_SOCKET_THREADS', hardware_cap)
   user_cap     = max(outer_threads, 1)

   force_blas   = getenv_int('XTB_SCC_FORCE_BLAS_THREADS', 0)
   min_blas     = getenv_int('XTB_SCC_MIN_BLAS_THREADS', 1)
   max_blas     = getenv_int('XTB_SCC_MAX_BLAS_THREADS', min(socket_cap, user_cap))

   if (force_blas > 0) then
      target = force_blas
   else
      if (use_parallel .and. .not. allow_nested) then
         target = 1
      else
         target = user_cap
      end if
   end if

   target = max(target, min_blas)
   target = min(target, max_blas)
   target = min(target, min(hardware_cap, socket_cap))
   target = max(target, 1)
end function select_blas_threads

subroutine setup_scc_thread_policy(ndim, policy, use_parallel, blas_threads_target)
   integer, intent(in)              :: ndim
   type(ThreadingPolicy), intent(inout) :: policy
   logical, intent(out)             :: use_parallel
   integer, intent(out)             :: blas_threads_target

   integer :: outer_threads
   logical :: allow_nested_blas

   outer_threads = current_omp_threads()
   if (outer_threads < 1) outer_threads = 1

   use_parallel = should_use_scc_parallel(ndim)
   allow_nested_blas = getenv_int('XTB_SCC_ALLOW_NESTED_BLAS', 0) /= 0
   blas_threads_target = select_blas_threads(use_parallel, outer_threads, allow_nested_blas)

   policy%omp_threads_before = outer_threads
   policy%max_active_levels_before = current_max_active_levels()
   if (.not. allow_nested_blas) then
      call set_max_active_levels(1, policy)
   else
      policy%max_active_levels_changed = .false.
   end if
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

integer function hardware_concurrency() result(nproc)
#ifdef _OPENMP
   nproc = omp_get_num_procs()
#else
   nproc = 1
#endif
   if (nproc < 1) nproc = 1
end function hardware_concurrency

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
