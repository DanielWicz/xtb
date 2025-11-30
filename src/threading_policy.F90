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
   integer, parameter :: dim_min_default       = 96
   real(wp), parameter :: work_threshold_default = 1.0e8_wp
   real(wp), parameter :: adapt_alpha          = 0.30_wp
   real(wp), parameter :: adapt_gap            = 0.15_wp
   integer, parameter  :: adapt_cooldown_iters = 2
   integer, parameter :: kernel_count          = 4
   integer, parameter :: kernel_h1             = 1
   integer, parameter :: kernel_dmat           = 2
   integer, parameter :: kernel_diag           = 3
   integer, parameter :: kernel_wiberg         = 4

   type :: ThreadingPolicy
      integer :: omp_threads_before         = 1
      integer :: blas_threads_before        = 1
      integer :: max_active_levels_before   = 1
      logical :: omp_changed                = .false.
      logical :: blas_changed               = .false.
      logical :: max_active_levels_changed  = .false.
   end type ThreadingPolicy

   type :: AdaptiveStats
      real(wp) :: avg_omp   = -1.0_wp
      real(wp) :: avg_blas  = -1.0_wp
      integer  :: n_omp     = 0
      integer  :: n_blas    = 0
      integer  :: last_mode = 0  ! 1=OMP, 2=BLAS
      integer  :: cooldown  = 0
   end type AdaptiveStats

   type(AdaptiveStats), save :: adapt_stats
   type(AdaptiveStats), save :: adapt_kernel(kernel_count)
   character(len=256), save :: adapt_state_file = ''
   logical, save :: adapt_state_loaded = .false.
   integer, save :: adapt_save_interval = 25
   integer, save :: adapt_save_counter  = 0
   logical, save :: adapt_persist       = .false.

   public :: ThreadingPolicy
   public :: setup_scc_thread_policy
   public :: restore_thread_policy
   public :: should_use_scc_parallel
    ! kernel IDs
   public :: kernel_h1, kernel_dmat, kernel_diag, kernel_wiberg
   public :: getenv_int, getenv_real, getenv_string
   public :: log_scc_iteration
   public :: log_scc_kernel

#if defined XTB_LAPACK_MKL
   interface
      integer(c_int) function mkl_get_max_threads() bind(C, name='mkl_get_max_threads')
         import :: c_int
      end function mkl_get_max_threads
      subroutine mkl_set_num_threads(nthreads) bind(C, name='mkl_set_num_threads')
         import :: c_int
         integer(c_int), value :: nthreads
      end subroutine mkl_set_num_threads
      subroutine mkl_set_dynamic(flag) bind(C, name='mkl_set_dynamic')
         import :: c_int
         integer(c_int), value :: flag
      end subroutine mkl_set_dynamic
   end interface
#elif defined XTB_LAPACK_ARMPL
   ! Arm Performance Libraries: rely on OpenMP controls (no OpenBLAS shim needed)
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

function getenv_string(name, default) result(value)
   character(len=*), intent(in) :: name
   character(len=*), intent(in) :: default
   character(len=len(default))  :: value
   integer :: stat, lenv
   character(len=len(default)) :: buffer

   value = default
   call get_environment_variable(name, buffer, length=lenv, status=stat)
   if (stat == 0 .and. lenv > 0) then
      if (lenv > len(buffer)) lenv = len(buffer)
      value = buffer(1:lenv)
   end if
end function getenv_string

logical function should_use_scc_parallel(ndim) result(use_parallel)
   integer, intent(in) :: ndim
   integer  :: dim_threshold
   integer  :: dim_min
   real(wp) :: work_threshold
   real(wp) :: work_est
   logical  :: force_parallel
   logical  :: force_serial
   integer  :: force_mode
   integer  :: outer_threads

   dim_threshold   = getenv_int('XTB_SCC_DIM_THRESHOLD', dim_threshold_default)
   dim_min         = getenv_int('XTB_SCC_MIN_DIM',       dim_min_default)
   work_threshold  = getenv_real('XTB_SCC_WORK_THRESHOLD', work_threshold_default)
   force_parallel  = getenv_int('XTB_SCC_FORCE_PARALLEL', 0) /= 0
   force_serial    = getenv_int('XTB_SCC_FORCE_SERIAL', 0) /= 0
   force_mode      = getenv_int('XTB_SCC_MODE', 0)  ! 0=auto, 1=force OMP, 2=force BLAS
   outer_threads   = current_omp_threads()

   if (force_parallel) then
      use_parallel = outer_threads > 1
      return
   end if

   if (force_serial) then
      use_parallel = .false.
      return
   end if

   select case (force_mode)
   case (1)
      use_parallel = outer_threads > 1
      return
   case (2)
      use_parallel = .false.
      return
   case default
      continue
   end select

   if (outer_threads < 2) then
      use_parallel = .false.
      return
   end if

   work_est        = real(ndim, wp) * real(ndim, wp) * real(ndim, wp)

   use_parallel = (ndim >= dim_min) .and. (ndim < dim_threshold) .and. &
      &           (work_est < work_threshold)
end function should_use_scc_parallel

integer function pick_adaptive_mode(default_parallel) result(mode)
   logical, intent(in) :: default_parallel
   real(wp) :: faster, slower, delta
   integer  :: env_adapt

   call maybe_init_adapt_state()
   env_adapt = getenv_int('XTB_SCC_ADAPT', 1)
   if (env_adapt == 0) then
      mode = merge(1, 2, default_parallel)
      return
   end if

   if (adapt_stats%last_mode == 0) adapt_stats%last_mode = merge(1, 2, default_parallel)

   if (adapt_stats%avg_omp < 0.0_wp .or. adapt_stats%avg_blas < 0.0_wp) then
      mode = adapt_stats%last_mode
      return
   end if

   faster = min(adapt_stats%avg_omp, adapt_stats%avg_blas)
   slower = max(adapt_stats%avg_omp, adapt_stats%avg_blas)
   delta  = slower - faster

   if (adapt_stats%cooldown > 0) then
      mode = adapt_stats%last_mode
      return
   end if

   if (delta > adapt_gap * faster) then
      if (adapt_stats%avg_omp < adapt_stats%avg_blas) then
         mode = 1
      else
         mode = 2
      end if
   else
      mode = adapt_stats%last_mode
   end if

   if (mode /= adapt_stats%last_mode) then
      adapt_stats%cooldown = adapt_cooldown_iters
      adapt_stats%last_mode = mode
   end if
end function pick_adaptive_mode

integer function pick_adaptive_mode_kernel(kernel_id, default_parallel) result(mode)
   integer, intent(in) :: kernel_id
   logical, intent(in) :: default_parallel
   real(wp) :: faster, slower, delta
   integer  :: env_adapt

   call maybe_init_adapt_state()
   env_adapt = getenv_int('XTB_SCC_ADAPT', 1)
   if (env_adapt == 0) then
      mode = merge(1, 2, default_parallel)
      return
   end if

   if (kernel_id < 1 .or. kernel_id > kernel_count) then
      mode = merge(1, 2, default_parallel)
      return
   end if

   if (adapt_kernel(kernel_id)%last_mode == 0) &
      adapt_kernel(kernel_id)%last_mode = merge(1, 2, default_parallel)

   if (adapt_kernel(kernel_id)%avg_omp < 0.0_wp .or. &
       adapt_kernel(kernel_id)%avg_blas < 0.0_wp) then
      mode = adapt_kernel(kernel_id)%last_mode
      return
   end if

   faster = min(adapt_kernel(kernel_id)%avg_omp, adapt_kernel(kernel_id)%avg_blas)
   slower = max(adapt_kernel(kernel_id)%avg_omp, adapt_kernel(kernel_id)%avg_blas)
   delta  = slower - faster

   if (adapt_kernel(kernel_id)%cooldown > 0) then
      mode = adapt_kernel(kernel_id)%last_mode
      return
   end if

   if (delta > adapt_gap * faster) then
      if (adapt_kernel(kernel_id)%avg_omp < adapt_kernel(kernel_id)%avg_blas) then
         mode = 1
      else
         mode = 2
      end if
   else
      mode = adapt_kernel(kernel_id)%last_mode
   end if

   if (mode /= adapt_kernel(kernel_id)%last_mode) then
      adapt_kernel(kernel_id)%cooldown = adapt_cooldown_iters
      adapt_kernel(kernel_id)%last_mode = mode
   end if
end function pick_adaptive_mode_kernel

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
   integer :: socket_count
   integer :: hardware_cap
   integer :: user_cap

   hardware_cap = hardware_concurrency()
   socket_cap   = getenv_int('XTB_SCC_SOCKET_THREADS', hardware_cap)
   socket_count = getenv_int('XTB_SCC_SOCKET_COUNT', 0)
   if (socket_count > 0) then
      ! Optional hint: provide number of sockets to cap BLAS threads per socket
      ! without hard-coding an explicit thread count. Falls back to
      ! XTB_SCC_SOCKET_THREADS when unset.
      socket_cap = max(1, hardware_cap / socket_count)
   end if
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

subroutine setup_scc_thread_policy(ndim, policy, use_parallel, blas_threads_target, kernel_id)
   integer, intent(in)              :: ndim
   type(ThreadingPolicy), intent(inout) :: policy
   logical, intent(out)             :: use_parallel
   integer, intent(out)             :: blas_threads_target
   integer, intent(in), optional    :: kernel_id

   integer :: outer_threads
   logical :: allow_nested_blas
   integer :: mode

   outer_threads = current_omp_threads()
   if (outer_threads < 1) outer_threads = 1

   use_parallel = should_use_scc_parallel(ndim)
   allow_nested_blas = getenv_int('XTB_SCC_ALLOW_NESTED_BLAS', 0) /= 0
   mode = getenv_int('XTB_SCC_MODE', 0)
   if (mode == 0) then
      if (present(kernel_id)) then
         mode = pick_adaptive_mode_kernel(kernel_id, use_parallel)
      else
         mode = pick_adaptive_mode(use_parallel)
      end if
   end if

   select case (mode)
   case (1)
      use_parallel = outer_threads > 1
      allow_nested_blas = .false.
   case (2)
      use_parallel = .false.
      allow_nested_blas = .true.
   case default
      continue
   end select

   if (present(kernel_id)) then
      adapt_kernel(kernel_id)%last_mode = merge(1, 2, use_parallel)
   else
      adapt_stats%last_mode = merge(1, 2, use_parallel)
   end if

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
   call mkl_set_dynamic(int(getenv_int('XTB_SCC_MKL_DYNAMIC', 0), c_int))
   call mkl_set_num_threads(int(nthreads, c_int))
   policy%blas_changed = .true.
#elif defined XTB_LAPACK_ARMPL
   policy%blas_threads_before = current_omp_threads()
#ifdef _OPENMP
   call omp_set_num_threads(nthreads)
   policy%omp_changed  = .true.
   policy%blas_changed = .true.
#else
   policy%blas_changed = .false.
#endif
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
#elif defined XTB_LAPACK_ARMPL
#ifdef _OPENMP
   if (policy%omp_changed) then
      call omp_set_num_threads(policy%omp_threads_before)
   end if
#endif
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

subroutine log_scc_iteration(iter_time, used_parallel)
   real(wp), intent(in) :: iter_time
   logical, intent(in)  :: used_parallel

   if (getenv_int('XTB_SCC_ADAPT', 1) == 0) return
   if (iter_time <= 0.0_wp) return

   if (used_parallel) then
      if (adapt_stats%avg_omp < 0.0_wp) adapt_stats%avg_omp = iter_time
      adapt_stats%avg_omp = (1.0_wp - adapt_alpha) * adapt_stats%avg_omp + adapt_alpha * iter_time
      adapt_stats%n_omp = adapt_stats%n_omp + 1
      adapt_stats%last_mode = 1
   else
      if (adapt_stats%avg_blas < 0.0_wp) adapt_stats%avg_blas = iter_time
      adapt_stats%avg_blas = (1.0_wp - adapt_alpha) * adapt_stats%avg_blas + adapt_alpha * iter_time
      adapt_stats%n_blas = adapt_stats%n_blas + 1
      adapt_stats%last_mode = 2
   end if

   if (adapt_stats%cooldown > 0) adapt_stats%cooldown = adapt_stats%cooldown - 1
   call maybe_save_adapt_state()
end subroutine log_scc_iteration

subroutine log_scc_kernel(kernel_id, iter_time, used_parallel)
   integer, intent(in) :: kernel_id
   real(wp), intent(in) :: iter_time
   logical, intent(in)  :: used_parallel

   if (getenv_int('XTB_SCC_ADAPT', 1) == 0) return
   if (iter_time <= 0.0_wp) return
   if (kernel_id < 1 .or. kernel_id > kernel_count) return

   if (used_parallel) then
      if (adapt_kernel(kernel_id)%avg_omp < 0.0_wp) adapt_kernel(kernel_id)%avg_omp = iter_time
      adapt_kernel(kernel_id)%avg_omp = (1.0_wp - adapt_alpha) * adapt_kernel(kernel_id)%avg_omp + adapt_alpha * iter_time
      adapt_kernel(kernel_id)%n_omp = adapt_kernel(kernel_id)%n_omp + 1
      adapt_kernel(kernel_id)%last_mode = 1
   else
      if (adapt_kernel(kernel_id)%avg_blas < 0.0_wp) adapt_kernel(kernel_id)%avg_blas = iter_time
      adapt_kernel(kernel_id)%avg_blas = (1.0_wp - adapt_alpha) * adapt_kernel(kernel_id)%avg_blas + adapt_alpha * iter_time
      adapt_kernel(kernel_id)%n_blas = adapt_kernel(kernel_id)%n_blas + 1
      adapt_kernel(kernel_id)%last_mode = 2
   end if

  if (adapt_kernel(kernel_id)%cooldown > 0) adapt_kernel(kernel_id)%cooldown = adapt_kernel(kernel_id)%cooldown - 1
   call maybe_save_adapt_state()
end subroutine log_scc_kernel

subroutine maybe_init_adapt_state()
   character(len=256) :: path
   integer :: interval

   if (adapt_state_loaded) return

   path = getenv_string('XTB_SCC_ADAPT_STATE', '')
   adapt_state_file = trim(path)
   adapt_persist = getenv_int('XTB_SCC_ADAPT_PERSIST', 0) /= 0
   interval = getenv_int('XTB_SCC_ADAPT_SAVE_INTERVAL', 25)
   if (interval < 1) interval = 25
   adapt_save_interval = interval
   adapt_save_counter  = 0

   if (adapt_persist .and. len_trim(adapt_state_file) > 0) call load_adapt_state()

   adapt_state_loaded = .true.
end subroutine maybe_init_adapt_state

subroutine load_adapt_state()
   integer :: ios, unit, k

   if (len_trim(adapt_state_file) == 0) return

   inquire(iolength=unit) unit
   open(newunit=unit, file=trim(adapt_state_file), status='old', action='read', iostat=ios)
   if (ios /= 0) return

   read(unit, *, iostat=ios) adapt_stats%avg_omp, adapt_stats%avg_blas, adapt_stats%n_omp, &
      & adapt_stats%n_blas, adapt_stats%last_mode, adapt_stats%cooldown
   if (ios /= 0) then
      close(unit)
      return
   end if

   do k = 1, kernel_count
      read(unit, *, iostat=ios) adapt_kernel(k)%avg_omp, adapt_kernel(k)%avg_blas, adapt_kernel(k)%n_omp, &
         & adapt_kernel(k)%n_blas, adapt_kernel(k)%last_mode, adapt_kernel(k)%cooldown
      if (ios /= 0) exit
   end do

   close(unit)
end subroutine load_adapt_state

subroutine save_adapt_state()
   integer :: ios, unit, k

   if (len_trim(adapt_state_file) == 0) return
   if (.not. adapt_persist) return

   inquire(iolength=unit) unit
   open(newunit=unit, file=trim(adapt_state_file), status='replace', action='write', iostat=ios)
   if (ios /= 0) return

   write(unit, '(2f18.10,2i8,2i4)') adapt_stats%avg_omp, adapt_stats%avg_blas, adapt_stats%n_omp, &
      & adapt_stats%n_blas, adapt_stats%last_mode, adapt_stats%cooldown
   do k = 1, kernel_count
      write(unit, '(2f18.10,2i8,2i4)') adapt_kernel(k)%avg_omp, adapt_kernel(k)%avg_blas, &
         & adapt_kernel(k)%n_omp, adapt_kernel(k)%n_blas, adapt_kernel(k)%last_mode, adapt_kernel(k)%cooldown
   end do
   close(unit)
end subroutine save_adapt_state

subroutine maybe_save_adapt_state()
   if (.not. adapt_persist) return
   if (len_trim(adapt_state_file) == 0) return
   adapt_save_counter = adapt_save_counter + 1
   if (adapt_save_counter >= adapt_save_interval) then
      call save_adapt_state()
      adapt_save_counter = 0
   end if
end subroutine maybe_save_adapt_state

end module xtb_threading_policy
