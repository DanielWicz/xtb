! This file is part of xtb.
! SPDX-Identifier: LGPL-3.0-or-later
!
! Regression test to ensure palnhess keeps requested inner OpenMP team size
! when nested parallelism is used inside numerical Hessian.

module test_palnhess_threads
   use testdrive, only : new_unittest, unittest_type, error_type, check
   use xtb_mctc_accuracy, only : wp
   use xtb_type_environment
   use xtb_type_molecule
   use xtb_type_restart
   use xtb_type_calculator
   use xtb_type_data, only : scc_results
   use xtb_setparam, only : set
   use omp_lib, only : omp_get_num_threads, omp_get_max_threads, omp_get_thread_limit
   implicit none
   private

   public :: collect_palnhess_threads

   integer :: observed_threads = 0

   type, extends(TCalculator) :: TDummyCalc
   contains
      procedure :: singlepoint => dummy_singlepoint
      procedure :: writeInfo   => dummy_writeinfo
   end type TDummyCalc

contains

!> Collect exported unit tests
subroutine collect_palnhess_threads(testsuite)
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ new_unittest("nested_inner_threads", test_nested_inner_threads) ]
end subroutine collect_palnhess_threads

subroutine reset_observed()
   observed_threads = 0
end subroutine reset_observed

subroutine dummy_singlepoint(self, env, mol, chk, printlevel, restart, &
      & energy, gradient, sigma, hlgap, results)
   class(TDummyCalc), intent(inout) :: self
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

   ! Determine the team size the runtime actually spawns for inner regions
   ! (hessian sets omp_set_num_threads(inner_threads) before calling us).
   block
      integer :: team
      team = 0
      !$omp parallel default(shared) reduction(max:team)
      team = omp_get_num_threads()
      !$omp end parallel

      !$omp critical
      observed_threads = max(observed_threads, team)
      !$omp end critical
   end block
end subroutine dummy_singlepoint

subroutine dummy_writeinfo(self, unit, mol)
   class(TDummyCalc), intent(in) :: self
   integer, intent(in) :: unit
   type(TMolecule), intent(in) :: mol
end subroutine dummy_writeinfo

subroutine test_nested_inner_threads(error)
   type(error_type), allocatable, intent(out) :: error
   type(TDummyCalc) :: calc
   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(TEnvironment) :: env
   real(wp), allocatable :: hess(:, :), dipgrad(:, :)
   integer, allocatable :: list(:)
   integer :: saved_palnhess, saved_omp
   integer :: requested_outer, requested_inner, expected_inner
   integer :: thread_limit, max_threads, i
   real(wp) :: step

   call init(env)
   call init(mol, ["H"], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))

   saved_palnhess = set%palnhess
   saved_omp = set%omp_threads

   requested_outer = 2
   requested_inner = 3
   max_threads = omp_get_max_threads()
   thread_limit = omp_get_thread_limit()

   if (requested_outer > max_threads) requested_outer = max_threads

   set%palnhess = requested_outer
   set%omp_threads = requested_inner

   expected_inner = requested_inner
   if (thread_limit > 0 .and. requested_outer * expected_inner > thread_limit) then
      expected_inner = max(1, thread_limit / requested_outer)
   end if

   allocate(hess(3*mol%n, 3*mol%n), dipgrad(3, 3*mol%n))
   hess = 0.0_wp
   dipgrad = 0.0_wp
   list = [(i, i = 1, mol%n)]
   step = 1.0e-3_wp

   call reset_observed()
   call calc%hessian(env, mol, chk, list, step, hess, dipgrad)

   call check(error, observed_threads >= max(1, expected_inner))
   call check(error, observed_threads <= expected_inner)

   set%palnhess = saved_palnhess
   set%omp_threads = saved_omp
end subroutine test_nested_inner_threads

end module test_palnhess_threads
