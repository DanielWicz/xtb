! This file is part of xtb.
! SPDX-Identifier: LGPL-3.0-or-later
!
! MPI regression test for numerical Hessian path. The test is only active
! when compiled with WITH_MPI and will skip otherwise.

module test_mpi_hessian
   use testdrive, only : new_unittest, unittest_type, error_type, check
   use xtb_mctc_accuracy, only : wp
   use xtb_type_environment
   use xtb_type_molecule
   use xtb_type_restart
   use xtb_type_calculator
   use xtb_type_data, only : scc_results
   use xtb_setparam, only : set
   use xtb_mpi_hessian, only : mpi_hessian_was_used, mpi_hessian_last_world_size, mpi_hessian_last_local_tasks
#ifdef WITH_MPI
   use mpi
#endif
   implicit none
   private

   public :: collect_mpi_hessian

   type, extends(TCalculator) :: TMpiDummyCalc
   contains
      procedure :: singlepoint => dummy_sp
      procedure :: writeInfo   => dummy_info
   end type TMpiDummyCalc

contains

subroutine collect_mpi_hessian(testsuite)
   type(unittest_type), allocatable, intent(out) :: testsuite(:)
#ifdef WITH_MPI
   testsuite = [ new_unittest("mpi_hessian_parallel", test_mpi_hessian_parallel) ]
#else
   allocate(testsuite(0))
#endif
end subroutine collect_mpi_hessian

#ifdef WITH_MPI
subroutine dummy_sp(self, env, mol, chk, printlevel, restart, energy, gradient, sigma, hlgap, results)
   class(TMpiDummyCalc), intent(inout) :: self
   type(TEnvironment), intent(inout) :: env
   type(TMolecule), intent(inout) :: mol
   type(TRestart), intent(inout) :: chk
   integer, intent(in) :: printlevel
   logical, intent(in) :: restart
   real(wp), intent(out) :: energy
   real(wp), intent(out) :: gradient(:, :)
   real(wp), intent(out) :: sigma(:, :)
   real(wp), intent(out) :: hlgap
   type(scc_results), intent(out) :: results

   energy = 0.0_wp
   gradient = 0.0_wp
   sigma = 0.0_wp
   hlgap = 0.0_wp
   results = scc_results()
end subroutine dummy_sp

subroutine dummy_info(self, unit, mol)
   class(TMpiDummyCalc), intent(in) :: self
   integer, intent(in) :: unit
   type(TMolecule), intent(in) :: mol
end subroutine dummy_info

subroutine test_mpi_hessian_parallel(error)
   type(error_type), allocatable, intent(out) :: error
   type(TMpiDummyCalc) :: calc
   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TRestart) :: chk
   real(wp), allocatable :: hess(:, :), dipgrad(:, :)
   integer, allocatable :: list(:)
   integer :: rank, size, ierr
   real(wp) :: step

   call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)
   if (size < 2) then
      return
   end if

   call init(env)
   call init(mol, ["H", "H"], reshape([0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.74_wp], [3, 2]))

   set%palnhess = 2
   set%omp_threads = 1

   step = 1.0e-3_wp
   allocate(hess(3*mol%n, 3*mol%n), dipgrad(3, 3*mol%n))
   hess = 0.0_wp
   dipgrad = 0.0_wp
   list = [(i, i=1, mol%n)]

   call calc%hessian(env, mol, chk, list, step, hess, dipgrad)

   if (rank == 0) then
      call check(error, mpi_hessian_was_used())
      call check(error, mpi_hessian_last_world_size() == size)
   end if
end subroutine test_mpi_hessian_parallel
#else
subroutine dummy_sp(self, env, mol, chk, printlevel, restart, energy, gradient, sigma, hlgap, results)
   class(TMpiDummyCalc), intent(inout) :: self
   type(TEnvironment), intent(inout) :: env
   type(TMolecule), intent(inout) :: mol
   type(TRestart), intent(inout) :: chk
   integer, intent(in) :: printlevel
   logical, intent(in) :: restart
   real(wp), intent(out) :: energy
   real(wp), intent(out) :: gradient(:, :)
   real(wp), intent(out) :: sigma(:, :)
   real(wp), intent(out) :: hlgap
   type(scc_results), intent(out) :: results
   energy = 0.0_wp; gradient = 0.0_wp; sigma = 0.0_wp; hlgap = 0.0_wp; results = scc_results()
end subroutine dummy_sp
subroutine dummy_info(self, unit, mol)
   class(TMpiDummyCalc), intent(in) :: self
   integer, intent(in) :: unit
   type(TMolecule), intent(in) :: mol
end subroutine dummy_info
subroutine test_mpi_hessian_parallel(error)
   type(error_type), allocatable, intent(out) :: error
   call skip(error, "MPI not enabled in build")
end subroutine test_mpi_hessian_parallel
#endif

end module test_mpi_hessian
