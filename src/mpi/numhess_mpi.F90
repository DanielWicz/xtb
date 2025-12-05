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

contains

!> Check whether an MPI-backed Hessian path is available at runtime
logical function mpi_hessian_available()
#ifdef WITH_MPI
   integer :: flag, ierr
   mpi_hessian_available = .false.
   call MPI_Initialized(flag, ierr)
   if (ierr == MPI_SUCCESS) then
      if (flag == 0) then
         call MPI_Init(ierr)
      end if
      call MPI_Initialized(flag, ierr)
      if (flag /= 0) then
         mpi_hessian_available = .true.
      end if
   end if
#else
   mpi_hessian_available = .false.
#endif
end function mpi_hessian_available

!> Placeholder MPI execution path.  The intent is to dispatch displacements
!> to worker ranks or spawned processes and accumulate the Hessian by
!> reduction.  Currently this routine simply reports "not handled" so the
!> OpenMP fallback is used.
subroutine mpi_hessian_execute(env, mol0, chk0, list, step, hess, dipgrad, polgrad, handled)
   type(TEnvironment), intent(inout) :: env
   type(TMolecule), intent(in) :: mol0
   type(TRestart), intent(in) :: chk0
   integer, intent(in) :: list(:)
   real(wp), intent(in) :: step
   real(wp), intent(inout) :: hess(:, :)
   real(wp), intent(inout) :: dipgrad(:, :)
   real(wp), intent(inout), optional :: polgrad(:, :)
   logical, intent(out) :: handled

#ifdef WITH_MPI
   ! Future implementation sketch (not active yet):
   !  * Split MPI_COMM_WORLD so rank 0 orchestrates the job.
   !  * Broadcast molecular/restart data once to all ranks.
   !  * Distribute displacement indices round-robin across ranks.
   !  * Each rank calls a provided single-point callback for its share and
   !    accumulates local Hessian slices.
   !  * Use MPI_Allreduce to sum hess/dipgrad/polgrad contributions.
   !  * Rank 0 keeps wall-clock estimates; all ranks respect set%omp_threads
   !    for intra-SCC threading.
   handled = .false.
#else
   handled = .false.
#endif
end subroutine mpi_hessian_execute

end module xtb_mpi_hessian
