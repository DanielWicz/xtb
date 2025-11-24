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
   use omp_lib
   implicit none

   public :: TCalculator
   private


   !> Base calculator
   type, abstract :: TCalculator

      real(wp) :: accuracy
      logical :: lSolv = .false.
      type(TSolvModel), allocatable :: solvation
      logical :: threadsafe = .true.
      logical :: high_memory_usage = .false.

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
   real(wp) :: er, el, dr(3), dl(3), sr(3, 3), sl(3, 3), egap, step2
   real(wp) :: alphal(3, 3), alphar(3, 3)
   real(wp) :: t0, t1, w0, w1, t_batch_start, t_batch_end
   real(wp), allocatable :: gr(:, :), gl(:, :)
   
   ! Adaptive load balancing variables
   integer :: n_tasks, n_threads_per_task, max_threads
   integer :: kat_start, kat_end, batch_size
   real(wp) :: throughput, best_throughput
   logical :: optimal_found, force_serial
   integer :: n_work_items

   call timing(t0, w0)
   step2 = 0.5_wp / step

   ! Environment settings
   call omp_set_max_active_levels(2)
   call omp_set_dynamic(.false.)
   max_threads = omp_get_max_threads()

   force_serial = .false.
   if (.not. self%threadsafe) force_serial = .true.

   ! Initialization for adaptive loop
   kat_start = 1
   n_tasks = 1 
   n_threads_per_task = max_threads
   best_throughput = 0.0_wp
   optimal_found = .false.

   ! Main loop over atoms (batched)
   do while (kat_start <= size(list))
      
      ! Determine batch size
      ! Ensure we have enough work items (atoms * 3 coords) to saturate n_tasks
      ! Each atom provides 3 parallelizable items (x, y, z displacements)
      ! We want at least 3 items per task to average out stragglers
      batch_size = max(1, n_tasks)
      
      ! If we are tuning, keep batch small to react quickly
      ! If optimal found, we can increase batch size for efficiency to reduce parallel overhead
      if (optimal_found) then
         batch_size = max(8, n_tasks * 4)
      endif

      kat_end = min(size(list), kat_start + batch_size - 1)
      n_work_items = (kat_end - kat_start + 1) * 3

      if (force_serial) then
         n_tasks = 1
         n_threads_per_task = max_threads
         optimal_found = .true. ! Skip tuning
      endif

      ! Set outer threads
      ! Note: omp_set_num_threads affects the NEXT parallel region
      call omp_set_num_threads(n_tasks)

      call timing(t_batch_start, w0)

      ! Parallel region for the batch
      !$omp parallel num_threads(n_tasks) default(none) &
      !$omp shared(self, env, mol0, chk0, list, step, hess, dipgrad, polgrad, step2, t0, w0, n_threads_per_task, kat_start, kat_end, t1, w1) &
      !$omp private(kat, iat, jat, jc, jj, ii, er, el, egap, gr, gl, sr, sl, dr, dl, alphar, alphal, ic)
      
      ! Set threads for inner SCF (affects subsequent parallel regions inside)
      if (n_threads_per_task > 0) call omp_set_num_threads(n_threads_per_task)

      allocate(gr(3, mol0%n), gl(3, mol0%n))

      !$omp do collapse(2) schedule(dynamic)
      do kat = kat_start, kat_end
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
                  ! Atomic update to shared array
                  !$omp atomic update
                  hess(jj, ii) = hess(jj, ii) &
                     & + (gr(jc, jat) - gl(jc, jat)) * step2
               end do
            end do

            if (kat == 3 .and. ic == 3) then
               !$omp critical(xtb_numdiff2)
               call timing(t1, w1)
               !$omp end critical(xtb_numdiff2)
            endif

         end do
      end do
      deallocate(gr, gl)
      !$omp end parallel

      call timing(t_batch_end, w0)
      
      ! Adaptive Logic
      if (.not. optimal_found .and. .not. force_serial) then
         throughput = real(n_work_items, wp) / (t_batch_end - t_batch_start + 1.0e-10_wp)
         
         if (self%threadsafe .and. set%verbose) then
            write(env%unit, '(A,I0,A,I0,A,F10.4,A)') "Hessian Tuning: Tasks=", n_tasks, &
               " Threads/Task=", n_threads_per_task, " Throughput=", throughput, " items/s"
         end if

         if (n_tasks == 1) then
            best_throughput = throughput
            ! Try doubling tasks if we have threads
            if (max_threads >= 2) then
               n_tasks = 2
               n_threads_per_task = max_threads / 2
            else
               optimal_found = .true.
            endif
         else
            ! Compare with previous
            if (throughput > best_throughput * 1.05_wp) then
               ! Improvement > 5%, keep increasing
               best_throughput = throughput
               if (n_tasks * 2 <= max_threads) then
                  n_tasks = n_tasks * 2
                  n_threads_per_task = max_threads / n_tasks
               else
                  optimal_found = .true.
               endif
            else
               ! No significant improvement or degradation
               ! Revert to previous config (which was safer/better)
               n_tasks = n_tasks / 2
               n_threads_per_task = max_threads / n_tasks
               optimal_found = .true.
               if (self%threadsafe .and. set%verbose) write(env%unit, '(A,I0,A)') "Hessian Tuning: Optimal found at ", n_tasks, " tasks."
            endif
         endif
      endif

      kat_start = kat_end + 1
   end do

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