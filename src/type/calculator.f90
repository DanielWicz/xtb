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
!$   use omp_lib, only : omp_get_max_active_levels, omp_get_max_threads, &
!$      & omp_get_num_procs, omp_set_max_active_levels, omp_set_num_threads
   implicit none

   public :: TCalculator
   public :: plan_hessian_threads
   private


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

!> Compute outer/inner thread plan for Hessian without oversubscription
pure subroutine plan_hessian_threads(ndispl, env_outer, env_inner, max_procs, max_threads, &
      & outer_threads, inner_threads)
   integer, intent(in) :: ndispl, env_outer, env_inner, max_procs, max_threads
   integer, intent(out) :: outer_threads, inner_threads

   outer_threads = 1
   inner_threads = 1

   if (max_procs < 1 .or. max_threads < 1) return

   ! clamp requested threads to available pools
   ! outer caps with max_threads to respect OMP max team size
   if (env_outer > 0) then
      outer_threads = min(ndispl, min(max_threads, max_procs))
   end if

   if (env_inner > 0) then
      inner_threads = min(env_inner, max_procs)
      outer_threads = min(ndispl, min(max_threads, max(1, max_procs / inner_threads)))
   else
      outer_threads = min(ndispl, min(max_threads, max_procs))
      inner_threads = max(1, max_procs / outer_threads)
   end if

   if (outer_threads < 1) outer_threads = 1
   if (inner_threads < 1) inner_threads = 1
end subroutine plan_hessian_threads


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
   real(wp) :: er, el, dr(3), dl(3), sr(3, 3), sl(3, 3), egap, step2
   real(wp) :: alphal(3, 3), alphar(3, 3)
   real(wp) :: t0, t1, w0, w1
   real(wp), allocatable :: gr(:, :), gl(:, :)
   integer :: ndispl, outer_threads, inner_threads
   integer :: max_threads, max_procs
   integer(kind=8) :: scratch_bytes, total_scratch_bytes
   integer :: env_outer, env_inner, sep, ios
   character(len=128) :: omp_env

   call timing(t0, w0)
   step2 = 0.5_wp / step
   ndispl = 3 * size(list)

   env_outer = -1
   env_inner = -1
   omp_env = ''
   ios = 1
   ! user override: XTB_HESS_OUTER / XTB_HESS_INNER take precedence
   call get_environment_variable('XTB_HESS_OUTER', omp_env, status=ios)
   if (ios == 0 .and. len_trim(omp_env) > 0) read(omp_env, *, iostat=ios) env_outer
   omp_env = ''
   call get_environment_variable('XTB_HESS_INNER', omp_env, status=ios)
   if (ios == 0 .and. len_trim(omp_env) > 0) read(omp_env, *, iostat=ios) env_inner
   omp_env = ''
   ios = 1
   call get_environment_variable('OMP_NUM_THREADS', omp_env, status=ios)
   if (ios == 0 .and. len_trim(omp_env) > 0) then
      sep = index(omp_env, ',')
      if (sep > 1) read(omp_env(1:sep-1), *, iostat=ios) env_outer
      if (sep == 0) read(omp_env, *, iostat=ios) env_outer
      if (sep > 0 .and. sep < len_trim(omp_env)) &
         & read(omp_env(sep+1:), *, iostat=ios) env_inner
   end if

   outer_threads = 1
   inner_threads = 1
   max_procs = 1
   max_threads = 1
   ! choose nested team sizes without oversubscribing the machine
!$ max_procs = max(1, omp_get_num_procs())
!$ max_threads = max(1, omp_get_max_threads())
!$ if (env_outer > 0) max_threads = min(max_threads, env_outer)
   call plan_hessian_threads(ndispl, env_outer, env_inner, max_procs, max_threads, &
      & outer_threads, inner_threads)
!$ if (omp_get_max_active_levels() < 2) call omp_set_max_active_levels(2)

   ! small diagnostic to help detect oversubscription / memory pressure
   scratch_bytes = int(6_8 * int(mol0%n, kind=8) * (storage_size(1.0_wp)/8), kind=8)
   total_scratch_bytes = scratch_bytes * int(outer_threads, kind=8)
   if (env%unit > 0) then
      write(env%unit,'("info: Hessian threads outer=",i0,", inner=",i0, ", displacements=",i0, &
         & ", omp_procs=",i0, ", scratch_per_outer=",f8.3," MiB, total_scratch=",f9.3," MiB")') &
         outer_threads, inner_threads, ndispl, max_procs, &
         real(scratch_bytes)/(1024.0_wp**2), real(total_scratch_bytes)/(1024.0_wp**2)
   end if

   ! emit one-line note if user request was capped
   if (env%unit > 0) then
      if ( (env_outer > 0 .and. outer_threads /= env_outer) .or. &
           (env_inner > 0 .and. inner_threads /= env_inner) ) then
         write(env%unit,'(a,i0,a,i0,a,i0)') 'info: Hessian threads rescaled to outer=', &
            outer_threads, ', inner=', inner_threads, ', cores=', max_procs
      end if
   end if

   !$omp parallel if(self%threadsafe) num_threads(outer_threads) default(none) &
   !$omp shared(self, env, mol0, chk0, list, step, hess, dipgrad, polgrad, step2, t0, w0, &
   !$omp& outer_threads, inner_threads) &
   !$omp private(kat, iat, jat, jc, jj, ii, er, el, egap, gr, gl, sr, sl, dr, dl, alphar, alphal, &
   !$omp& t1, w1)

   !$ call omp_set_num_threads(inner_threads)

   ! first-touch initialize thread-local scratch to bind pages to this thread's NUMA node
   allocate(gr(3, mol0%n), source = 0.0_wp)
   allocate(gl(3, mol0%n), source = 0.0_wp)

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
   !$omp end parallel
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

end module xtb_type_calculator
