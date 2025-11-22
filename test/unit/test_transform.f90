! This file is part of xtb.
! SPDX-Identifier: LGPL-3.0-or-later
!
! Harness for validating CAO→SAO tensor transforms.
module test_transform
   use testdrive, only : new_unittest, unittest_type, error_type, check
   use xtb_mctc_accuracy, only : wp
   use xtb_intgrad, only : dtrf2
   implicit none
   private

   public :: collect_transform

   integer, parameter :: maxl = 2
   real(wp), parameter :: tol = 1.0e-12_wp

contains

subroutine collect_transform(testsuite)
   type(unittest_type), allocatable, intent(out) :: testsuite(:)
   testsuite = [ new_unittest("cao_to_sao_consistency", test_cao_to_sao_consistency) ]
end subroutine collect_transform


subroutine test_cao_to_sao_consistency(error)
   type(error_type), allocatable, intent(out) :: error
   integer :: li, lj, rep
   real(wp) :: dd(3,6,6), dd_ref(3,6,6), dd_fast(3,6,6)
   real(wp) :: qq(6,6,6), qq_ref(6,6,6), qq_fast(6,6,6)
   real(wp) :: diff

   call random_seed()
   diff = 0.0_wp

   do li = 0, maxl
      do lj = 0, maxl
         do rep = 1, 64
            call random_number(dd)
            dd_ref = dd
            dd_fast = dd
            call reference_transform(dd_ref, 3, li, lj)
            call fast_transform(dd_fast, 3, li, lj)
            diff = max(diff, maxval(abs(dd_fast - dd_ref)))
            if (diff > tol) exit

            call random_number(qq)
            qq_ref = qq
            qq_fast = qq
            call reference_transform(qq_ref, 6, li, lj)
            call fast_transform(qq_fast, 6, li, lj)
            diff = max(diff, maxval(abs(qq_fast - qq_ref)))
            if (diff > tol) exit
         end do
         if (diff > tol) exit
      end do
      if (diff > tol) exit
   end do

   call check(error, diff < tol, 'CAO→SAO fast transform deviates from reference')

end subroutine test_cao_to_sao_consistency


pure subroutine reference_transform(tensor, ncomp, li, lj)
   real(wp), intent(inout) :: tensor(:, :, :)
   integer, intent(in) :: ncomp, li, lj
   integer :: ic
   real(wp) :: tmp(6,6)
   do ic = 1, ncomp
      tmp = tensor(ic, 1:6, 1:6)
      call dtrf2(tmp, li, lj)
      tensor(ic, 1:6, 1:6) = tmp
   end do
end subroutine reference_transform


pure subroutine fast_transform(tensor, ncomp, li, lj)
   real(wp), intent(inout) :: tensor(:, :, :)
   integer, intent(in) :: ncomp, li, lj
   ! placeholder that can be replaced with an optimized kernel
   call reference_transform(tensor, ncomp, li, lj)
end subroutine fast_transform

end module test_transform
