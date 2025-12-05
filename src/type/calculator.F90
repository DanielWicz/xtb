! This file is part of xtb.
!
! Copyright (C) 2017-2020 Stefan Grimme
!
! xtb is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! xtb is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with xtb.  If not, see <https://www.gnu.org/licenses/>.

!> abstract calculator that hides implementation details from calling codes
module xtb_type_calculator
   use xtb_mctc_accuracy, only : wp
   use xtb_solv_model, only : TSolvModel
   use xtb_type_data, only : scc_results
   use xtb_type_environment, only : TEnvironment
   use xtb_type_molecule, only : TMolecule
   use xtb_type_restart, only : TRestart
   use xtb_setparam, only : set
   use xtb_mpi_hessian, only : mpi_hessian_available, mpi_hessian_execute
   use iso_c_binding, only : c_char, c_int, c_null_char
#ifdef _OPENMP
   use omp_lib, only : omp_get_max_threads, omp_get_max_active_levels, omp_set_dynamic, &
      & omp_set_max_active_levels, omp_set_num_threads, omp_get_thread_limit, omp_set_nested
#ifdef WITH_MKL
   use mkl_service, only : mkl_set_num_threads
#endif
#endif
   implicit none

   public :: TCalculator
   private

   !> Small helper to remember environment variables we override
   type :: omp_env_state
      character(len=:), allocatable :: omp_num_threads
      character(len=:), allocatable :: omp_max_active_levels
      character(len=:), allocatable :: omp_nested
      character(len=:), allocatable :: omp_dynamic
      character(len=:), allocatable :: openblas_num_threads
      character(len=:), allocatable :: mkl_num_threads
   end type omp_env_state


   !> Base calculator
   type, abstract :: TCalculator

      real(wp) :: accuracy
      logical :: lSolv = .false.
      type(TSolvModel), allocatable :: solvation
      logical :: threadsafe = .true.

   contains

      !> Perform single point calculation
      procedure(singlepoint), deferred :: singlepoint

      !> Perform hessian calculation
      procedure :: hessian

      !> Write informative printout
      procedure(writeInfo), deferred :: writeInfo

   end type TCalculator

   interface
      integer(c_int) function c_setenv(name, value, overwrite) bind(C, name="setenv")
         import :: c_char, c_int
         character(c_char), dimension(*), intent(in) :: name
         character(c_char), dimension(*), intent(in) :: value
         integer(c_int), value :: overwrite
      end function c_setenv

      integer(c_int) function c_unsetenv(name) bind(C, name="unsetenv")
         import :: c_char, c_int
         character(c_char), dimension(*), intent(in) :: name
      end function c_unsetenv
   end interface


   abstract interface
      subroutine singlepoint(self, env, mol, chk, printlevel, restart, &
            & energy, gradient, sigma, hlgap, results)
         import :: TCalculator, TEnvironment, TMolecule, TRestart, wp
         import :: scc_results

         !> Calculator instance
         class(TCalculator), intent(inout) :: self

         !> Computational environment
         type(TEnvironment), intent(inout) :: env

         !> Molecular structure data
         type(TMolecule), intent(inout) :: mol

         !> Wavefunction data
         type(TRestart), intent(inout) :: chk

         !> Print level for IO
         integer, intent(in) :: printlevel

         !> Restart from previous results
         logical, intent(in) :: restart

         !> Total energy
         real(wp), intent(out) :: energy

         !> Molecular gradient
         real(wp), intent(out) :: gradient(:, :)

         !> Strain derivatives
         real(wp), intent(out) :: sigma(:, :)

         !> HOMO-LUMO gap
         real(wp), intent(out) :: hlgap

         !> Detailed results
         type(scc_results), intent(out) :: results

      end subroutine singlepoint


      subroutine writeInfo(self, unit, mol)
         import :: TCalculator, TMolecule

         !> Calculator instance
         class(TCalculator), intent(in) :: self

         !> Unit for I/O
         integer, intent(in) :: unit

         !> Molecular structure data
         type(TMolecule), intent(in) :: mol

      end subroutine writeInfo
   end interface


contains


!> Evaluate hessian by finite difference for all atoms
subroutine hessian(self, env, mol0, chk0, list, step, hess, dipgrad, polgrad)
   character(len=*), parameter :: source = "hessian_numdiff_numdiff2"
   !> Single point calculator
   class(TCalculator), intent(inout) :: self
   !> Computation environment
   type(TEnvironment), intent(inout) :: env
   !> Molecular structure data
   type(TMolecule), intent(in) :: mol0
   !> Restart data
   type(TRestart), intent(in) :: chk0
   !> List of atoms to displace
   integer, intent(in) :: list(:)
   !> Step size for numerical differentiation
   real(wp), intent(in) :: step
   !> Array to add Hessian to
   real(wp), intent(inout) :: hess(:, :)
   !> Array to add dipole gradient to
   real(wp), intent(inout) :: dipgrad(:, :)
   !> Array to add polarizability gradient to
   real(wp), intent(inout), optional :: polgrad(:, :)

   integer :: iat, jat, kat, ic, jc, ii, jj
   integer :: outer_threads, inner_threads, ntasks, max_threads
#ifdef _OPENMP
   integer :: thread_limit, active_levels, desired_levels
   character(len=128) :: thread_msg
   logical :: nested_ok
#endif
   real(wp) :: er, el, dr(3), dl(3), sr(3, 3), sl(3, 3), egap, step2
   real(wp) :: alphal(3, 3), alphar(3, 3)
   real(wp) :: t0, t1, w0, w1
   real(wp), allocatable :: gr(:, :), gl(:, :)
   logical :: do_parallel
   logical :: mpi_handled
   logical :: env_overridden
   type(omp_env_state) :: omp_state
#ifdef _OPENMP
   logical, parameter :: has_openmp = .true.
#else
   logical, parameter :: has_openmp = .false.
#endif

   call timing(t0, w0)

   if (step <= 0.0_wp) then
      call env%error("Step for numerical Hessian must be positive", source)
      return
   end if

   step2 = 0.5_wp / step
   ntasks = max(1, 3 * size(list))
#ifdef _OPENMP
   max_threads = omp_get_max_threads()
   thread_limit = omp_get_thread_limit()
#else
   max_threads = 1
#endif

   env_overridden = .false.

   mpi_handled = .false.
   if (mpi_hessian_available()) then
      call mpi_hessian_execute(env, mol0, chk0, list, step, hess, dipgrad, polgrad, mpi_handled)
      if (mpi_handled) return
   end if
   
   ! Logic for thread distribution:
   ! 1. If user requested --palnhess (outer), we prioritize that count.
   !    If --parallel (inner) was NOT set, we auto-scale inner threads to fill the machine.
   ! 2. If user requested --parallel (inner) but NOT --palnhess, we auto-scale outer threads.
   ! 3. If neither, we default to serial displacement (outer=1) and full parallel SCF (inner=max).
   ! 4. If both, we respect user input (allowing oversubscription if requested).

   if (set%palnhess > 0) then
      outer_threads = set%palnhess
      if (set%omp_threads > 0) then
         inner_threads = set%omp_threads
      else
         ! Auto-scale inner threads
         inner_threads = max(1, max_threads / outer_threads)
      end if
   else
      ! palnhess not set (0)
      if (set%omp_threads > 0) then
         inner_threads = set%omp_threads
         outer_threads = max(1, max_threads / inner_threads)
      else
         ! Neither set, default to standard behavior (1 task, all threads)
         inner_threads = max_threads
         outer_threads = 1
      end if
   end if

   inner_threads = max(1, inner_threads)
   outer_threads = max(1, min(outer_threads, ntasks))

   ! If calculator not threadsafe, force serial execution of displacements
   if (.not. self%threadsafe) then
      if (outer_threads > 1) then
         call env%warning("Calculator not thread-safe; running numerical Hessian serially", source)
      end if
      outer_threads = 1
   end if
   
   if (.not. has_openmp .and. outer_threads > 1) then
      call env%warning("Binary compiled without OpenMP support, running numerical Hessian serially", source)
   end if

   if (.not. has_openmp) outer_threads = 1

   do_parallel = has_openmp .and. (outer_threads > 1)

#ifdef _OPENMP
   if (do_parallel .and. thread_limit > 0 .and. outer_threads * inner_threads > thread_limit) then
      inner_threads = max(1, thread_limit / outer_threads)
      write(thread_msg,'(a,i0,a,i0,a,i0)') 'Requested palnhess*parallel (', outer_threads*inner_threads, &
         ') exceeds OpenMP thread limit (', thread_limit, '); capping inner threads to ', inner_threads
      call env%warning(trim(thread_msg), source)
   end if

   ! Stabilise thread control and enable nested OpenMP so each displacement can spawn its own SCC team.
   call omp_set_dynamic(.false.)
   if (do_parallel) then
      desired_levels = 1 + merge(1, 0, inner_threads > 1)
      call omp_set_max_active_levels(max(desired_levels, omp_get_max_active_levels()))
      call omp_set_nested(desired_levels > 1)
      active_levels = omp_get_max_active_levels()
      nested_ok = active_levels >= desired_levels
      if (.not. nested_ok .and. inner_threads > 1) then
         call env%warning("OpenMP runtime disallows nested teams; running Hessian displacements serially with inner threading", source)
         outer_threads = 1
         do_parallel = .false.
      end if
   end if

   if (do_parallel) then
      call apply_nested_env(omp_state, outer_threads, inner_threads)
      env_overridden = .true.
   end if

   if (.not. do_parallel) then
      call omp_set_num_threads(inner_threads)
#ifdef WITH_MKL
      call mkl_set_dynamic(.false.)
      call mkl_set_num_threads(inner_threads)
#endif
   end if
#endif

!$omp parallel if(do_parallel) num_threads(outer_threads) default(none) &
   !$omp shared(self, env, mol0, chk0, list, step, hess, dipgrad, polgrad, step2, t0, w0, do_parallel, outer_threads) &
   !$omp private(kat, iat, jat, jc, jj, ii, er, el, egap, gr, gl, sr, sl, dr, dl, alphar, alphal, t1, w1, ic) &
   !$omp firstprivate(inner_threads)

      ! Enforce thread limit for nested regions and libraries (BLAS/MKL)
   !$ if (do_parallel) call omp_set_num_threads(inner_threads)
#ifdef WITH_MKL
   !$ if (do_parallel) then
   !$    call mkl_set_dynamic(.false.)
   !$    call mkl_set_num_threads(inner_threads)
   !$ end if
#endif
      
      allocate(gr(3, mol0%n), gl(3, mol0%n))

      !$omp do collapse(2) schedule(runtime)
      do kat = 1, size(list)
         do ic = 1, 3

            iat = list(kat)
            ii = 3*(iat - 1) + ic
            er = 0.0_wp
            el = 0.0_wp
            gr = 0.0_wp
            gl = 0.0_wp

            call hessian_point(self, env, mol0, chk0, iat, ic, +step, er, gr, sr, egap, dr, alphar)
            call hessian_point(self, env, mol0, chk0, iat, ic, -step, el, gl, sl, egap, dl, alphal)

            if (present(polgrad)) then
               polgrad(1, ii) = (alphar(1, 1) - alphal(1, 1)) * step2
               polgrad(2, ii) = (alphar(1, 2) - alphal(1, 2)) * step2
               polgrad(3, ii) = (alphar(2, 2) - alphal(2, 2)) * step2
               polgrad(4, ii) = (alphar(1, 3) - alphal(1, 3)) * step2
               polgrad(5, ii) = (alphar(2, 3) - alphal(2, 3)) * step2
               polgrad(6, ii) = (alphar(3, 3) - alphal(3, 3)) * step2
            endif

            dipgrad(:, ii) = (dr - dl) * step2

            do jat = 1, mol0%n
               do jc = 1, 3
                  jj = 3*(jat - 1) + jc
                  hess(jj, ii) = hess(jj, ii) &
                     & + (gr(jc, jat) - gl(jc, jat)) * step2
               end do
            end do

            if (kat == 3 .and. ic == 3) then
               !$omp critical(xtb_numdiff2)
               call timing(t1, w1)
               write(*,'("estimated CPU  time",F10.2," min")') &
                  & 0.3333333_wp*size(list)*(t1-t0)/60.0_wp
               write(*,'("estimated wall time",F10.2," min")') &
                  & 0.3333333_wp*size(list)*(w1-w0)/60.0_wp
               !$omp end critical(xtb_numdiff2)
            endif

         end do
      end do
      !$omp end do

      deallocate(gr, gl)

   !$omp end parallel

   if (env_overridden) call restore_nested_env(omp_state)
end subroutine hessian

subroutine hessian_point(self, env, mol0, chk0, iat, ic, step, energy, gradient, sigma, egap, dipole, alpha)
   class(TCalculator), intent(inout) :: self
   type(TEnvironment), intent(inout) :: env
   type(TMolecule), intent(in) :: mol0
   type(TRestart), intent(in) :: chk0
   integer, intent(in) :: ic, iat
   real(wp), intent(in) :: step
   real(wp), intent(out) :: energy
   real(wp), intent(out) :: gradient(:, :)
   real(wp), intent(out) :: sigma(3, 3)
   real(wp), intent(out) :: egap
   real(wp), intent(out) :: dipole(3)
   real(wp), intent(out) :: alpha(3, 3)

   ! internal variables
   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(scc_results) :: res

   call mol%copy(mol0)
   mol%xyz(ic, iat) = mol0%xyz(ic, iat) + step
   call chk%copy(chk0)
   call self%singlepoint(env, mol, chk, -1, .true., energy, gradient, sigma, egap, res)

   dipole = res%dipole
   alpha(:, :) = res%alpha

end subroutine hessian_point

!> Convert a character string to a C NULL terminated array
pure function to_c_string(str) result(cstr)
   character(len=*), intent(in) :: str
   character(kind=c_char), allocatable :: cstr(:)
   integer :: n, i

   n = len_trim(str)
   allocate(cstr(n+1))
   do i = 1, n
      cstr(i) = str(i:i)
   end do
   cstr(n+1) = c_null_char
end function to_c_string

!> Safely capture the current value of an environment variable (if present)
subroutine capture_env(name, value)
   character(len=*), intent(in) :: name
   character(len=:), allocatable, intent(out) :: value
   integer :: stat, lenv

   call get_environment_variable(name, length=lenv, status=stat)
   if (stat == 0 .and. lenv > 0) then
      allocate(character(len=lenv) :: value)
      call get_environment_variable(name, value, status=stat)
      if (stat /= 0) then
         deallocate(value)
      end if
   end if
end subroutine capture_env

!> Set environment variable, overwriting any existing value
subroutine set_env(name, value)
   character(len=*), intent(in) :: name, value
   character(kind=c_char), allocatable :: cname(:), cvalue(:)
   integer(c_int) :: rc

   cname = to_c_string(name)
   cvalue = to_c_string(value)
   rc = c_setenv(cname, cvalue, 1_c_int)
end subroutine set_env

!> Unset environment variable if present
subroutine unset_env(name)
   character(len=*), intent(in) :: name
   character(kind=c_char), allocatable :: cname(:)
   integer(c_int) :: rc

   cname = to_c_string(name)
   rc = c_unsetenv(cname)
end subroutine unset_env

!> Apply nested OpenMP/BLAS environment for palnhess and remember previous state
subroutine apply_nested_env(state, outer_threads, inner_threads)
   type(omp_env_state), intent(out) :: state
   integer, intent(in) :: outer_threads, inner_threads
   character(len=32) :: outer_str, inner_str
   character(len=64) :: nested_list

   call capture_env('OMP_NUM_THREADS', state%omp_num_threads)
   call capture_env('OMP_MAX_ACTIVE_LEVELS', state%omp_max_active_levels)
   call capture_env('OMP_NESTED', state%omp_nested)
   call capture_env('OMP_DYNAMIC', state%omp_dynamic)
   call capture_env('OPENBLAS_NUM_THREADS', state%openblas_num_threads)
   call capture_env('MKL_NUM_THREADS', state%mkl_num_threads)

   write(outer_str, '(I0)') max(1, outer_threads)
   write(inner_str, '(I0)') max(1, inner_threads)
   nested_list = trim(outer_str)//','//trim(inner_str)

   call set_env('OMP_NUM_THREADS', nested_list)
   call set_env('OMP_MAX_ACTIVE_LEVELS', '2')
   call set_env('OMP_NESTED', 'TRUE')
   call set_env('OMP_DYNAMIC', 'FALSE')
   call set_env('OPENBLAS_NUM_THREADS', trim(inner_str))
   call set_env('MKL_NUM_THREADS', trim(inner_str))
end subroutine apply_nested_env

!> Restore environment variables saved by apply_nested_env
subroutine restore_nested_env(state)
   type(omp_env_state), intent(in) :: state

   call restore_env_var('OMP_NUM_THREADS', state%omp_num_threads)
   call restore_env_var('OMP_MAX_ACTIVE_LEVELS', state%omp_max_active_levels)
   call restore_env_var('OMP_NESTED', state%omp_nested)
   call restore_env_var('OMP_DYNAMIC', state%omp_dynamic)
   call restore_env_var('OPENBLAS_NUM_THREADS', state%openblas_num_threads)
   call restore_env_var('MKL_NUM_THREADS', state%mkl_num_threads)
end subroutine restore_nested_env

!> Helper to restore or unset a single environment variable
subroutine restore_env_var(name, value)
   character(len=*), intent(in) :: name
   character(len=:), allocatable, intent(in) :: value

   if (allocated(value)) then
      call set_env(name, value)
   else
      call unset_env(name)
   end if
end subroutine restore_env_var

end module xtb_type_calculator
