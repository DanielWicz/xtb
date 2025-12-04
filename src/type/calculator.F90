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
   integer :: thread_limit
   character(len=128) :: thread_msg
#endif
   real(wp) :: er, el, dr(3), dl(3), sr(3, 3), sl(3, 3), egap, step2
   real(wp) :: alphal(3, 3), alphar(3, 3)
   real(wp) :: t0, t1, w0, w1
   real(wp), allocatable :: gr(:, :), gl(:, :)
   logical :: do_parallel
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
   if (do_parallel .and. thread_limit > 0) then
      if (outer_threads * inner_threads > thread_limit) then
         write(thread_msg,'(a,i0,a,i0,a)') 'Requested palnhess*parallel (', outer_threads*inner_threads, &
            ') exceeds OpenMP thread limit (', thread_limit, '); runtime may cap threads'
         call env%warning(trim(thread_msg), source)
      end if
   end if

   ! Stabilise thread control and enable nested OpenMP so each displacement can spawn its own SCC team.
   call omp_set_dynamic(.false.)
   if (do_parallel) then
      call omp_set_max_active_levels(max(2, omp_get_max_active_levels()))
      call omp_set_nested(.true.)
   else
      if (set%omp_threads > 0) then
         call omp_set_num_threads(inner_threads)
#ifdef WITH_MKL
         call mkl_set_num_threads(inner_threads)
#endif
      end if
   end if
#endif

   !$omp parallel if(do_parallel) num_threads(outer_threads) default(none) &
   !$omp shared(self, env, mol0, chk0, list, step, hess, dipgrad, polgrad, step2, t0, w0, do_parallel, outer_threads) &
   !$omp private(kat, iat, jat, jc, jj, ii, er, el, egap, gr, gl, sr, sl, dr, dl, alphar, alphal, t1, w1, ic) &
   !$omp firstprivate(inner_threads)

      ! Enforce thread limit for nested regions and libraries (BLAS/MKL)
      !$ if (do_parallel) call omp_set_num_threads(inner_threads)
#ifdef WITH_MKL
      !$ if (do_parallel) call mkl_set_num_threads(inner_threads)
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
