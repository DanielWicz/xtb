! This file is part of xtb.
! SPDX-License-Identifier: LGPL-3.0-or-later
!
! MPI scaffold for numerical Hessian parallelisation.  At the moment the
! implementation is a stub that can be compiled without MPI; it defines the
! hooks we will call from the calculator and documents the intended workflow.

module xtb_mpi_hessian
   use xtb_mctc_accuracy, only : wp
   use xtb_type_environment, only : TEnvironment
   use xtb_type_molecule, only : TMolecule
   use xtb_type_restart, only : TRestart
#ifdef WITH_MPI
   use mpi
#endif
   implicit none
   private

   public :: mpi_hessian_available, mpi_hessian_execute
   public :: mpi_hessian_was_used, mpi_hessian_last_world_size, mpi_hessian_last_local_tasks

#ifdef WITH_MPI
   logical :: mpi_last_used = .false.
   integer :: mpi_last_world = 1
   integer :: mpi_last_local = 0
#endif

   abstract interface
      subroutine hess_task_iface(iat, ic, step, energy, gradient, sigma, egap, dipole, alpha)
         import :: wp
         integer, intent(in) :: iat, ic
         real(wp), intent(in) :: step
         real(wp), intent(out) :: energy
         real(wp), intent(out) :: gradient(:, :)
         real(wp), intent(out) :: sigma(3, 3)
         real(wp), intent(out) :: egap
         real(wp), intent(out) :: dipole(3)
         real(wp), intent(out) :: alpha(3, 3)
      end subroutine hess_task_iface
   end interface

contains

!> Query whether MPI Hessian path is usable (MPI library available and size>1)
logical function mpi_hessian_available()
#ifdef WITH_MPI
   integer :: flag, ierr, size

   mpi_hessian_available = .false.
   call MPI_Initialized(flag, ierr)
   if (ierr /= MPI_SUCCESS) return
   if (flag == 0) then
      call MPI_Init(ierr)
      if (ierr /= MPI_SUCCESS) return
   end if
   call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)
   if (ierr == MPI_SUCCESS .and. size > 1) mpi_hessian_available = .true.
#else
   mpi_hessian_available = .false.
#endif
end function mpi_hessian_available

!> Execute numerical Hessian using MPI task parallelism. Each rank handles a
!> round-robin subset of displacements; results are summed across all ranks.
subroutine mpi_hessian_execute(task, mol0, chk0, list, step, hess, dipgrad, polgrad, handled)
   procedure(hess_task_iface) :: task
   type(TMolecule), intent(in) :: mol0
   type(TRestart), intent(in) :: chk0
   integer, intent(in) :: list(:)
   real(wp), intent(in) :: step
   real(wp), intent(inout) :: hess(:, :)
   real(wp), intent(inout) :: dipgrad(:, :)
   real(wp), intent(inout), optional :: polgrad(:, :)
   logical, intent(out) :: handled

#ifdef WITH_MPI
   integer :: ierr, rank, size, flag
   integer :: ntasks, task_id, kat, ic, iat, ii, jat, jc, jj
   real(wp), allocatable :: hloc(:, :), dloc(:, :), ploc(:, :)
   real(wp), allocatable :: gr(:, :), gl(:, :)
   real(wp) :: er, el, egap, dr(3), dl(3), sr(3, 3), sl(3, 3)
   real(wp) :: alphal(3, 3), alphar(3, 3)
   real(wp) :: step2
   logical :: have_pol

   handled = .false.
   call MPI_Initialized(flag, ierr)
   if (flag == 0) call MPI_Init(ierr)
   if (ierr /= MPI_SUCCESS) return

   call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)
   if (ierr /= MPI_SUCCESS) return
   if (size < 2) return

   ntasks = max(1, 3 * size(list))
   step2 = 0.5_wp / step
   have_pol = present(polgrad)

   allocate(hloc(size(hess,1), size(hess,2)))
   allocate(dloc(size(dipgrad,1), size(dipgrad,2)))
   hloc = 0.0_wp
   dloc = 0.0_wp
   if (have_pol) then
      allocate(ploc(size(polgrad,1), size(polgrad,2)))
      ploc = 0.0_wp
   end if

   allocate(gr(3, mol0%n), gl(3, mol0%n))

   do task_id = 1, ntasks
      if (mod(task_id-1, size) /= rank) cycle
      kat = (task_id - 1) / 3 + 1
      ic  = mod(task_id - 1, 3) + 1
      iat = list(kat)
      ii  = 3*(iat - 1) + ic

      er = 0.0_wp; el = 0.0_wp
      gr = 0.0_wp; gl = 0.0_wp

      call task(iat, ic, +step, er, gr, sr, egap, dr, alphar)
      call task(iat, ic, -step, el, gl, sl, egap, dl, alphal)

      if (have_pol) then
         ploc(1, ii) = (alphar(1, 1) - alphal(1, 1)) * step2
         ploc(2, ii) = (alphar(1, 2) - alphal(1, 2)) * step2
         ploc(3, ii) = (alphar(2, 2) - alphal(2, 2)) * step2
         ploc(4, ii) = (alphar(1, 3) - alphal(1, 3)) * step2
         ploc(5, ii) = (alphar(2, 3) - alphal(2, 3)) * step2
         ploc(6, ii) = (alphar(3, 3) - alphal(3, 3)) * step2
      end if

      dloc(:, ii) = (dr - dl) * step2

      do jat = 1, mol0%n
         do jc = 1, 3
            jj = 3*(jat - 1) + jc
            hloc(jj, ii) = hloc(jj, ii) + (gr(jc, jat) - gl(jc, jat)) * step2
         end do
      end do
   end do

   call MPI_Allreduce(hloc, hess, size(hess), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
   call MPI_Allreduce(dloc, dipgrad, size(dipgrad), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
   if (have_pol) then
      call MPI_Allreduce(ploc, polgrad, size(polgrad), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      deallocate(ploc)
   end if

   deallocate(hloc, dloc, gr, gl)

   mpi_last_used  = .true.
   mpi_last_world = size
   mpi_last_local = ntasks/size + merge(1,0, mod(ntasks,size) > rank)
   handled = .true.

   call MPI_Barrier(MPI_COMM_WORLD, ierr)
#else
   handled = .false.
#endif
end subroutine mpi_hessian_execute

!> Query flag set when MPI Hessian path was taken in the last call
logical function mpi_hessian_was_used()
#ifdef WITH_MPI
   mpi_hessian_was_used = mpi_last_used
#else
   mpi_hessian_was_used = .false.
#endif
end function mpi_hessian_was_used

integer function mpi_hessian_last_world_size()
#ifdef WITH_MPI
   mpi_hessian_last_world_size = mpi_last_world
#else
   mpi_hessian_last_world_size = 1
#endif
end function mpi_hessian_last_world_size

integer function mpi_hessian_last_local_tasks()
#ifdef WITH_MPI
   mpi_hessian_last_local_tasks = mpi_last_local
#else
   mpi_hessian_last_local_tasks = 0
#endif
end function mpi_hessian_last_local_tasks

end module xtb_mpi_hessian
