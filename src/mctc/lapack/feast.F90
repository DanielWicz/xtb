! This file is part of xtb.
!
! Copyright (C) 2024
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

!> Lightweight wrapper around MKL FEAST to replace DSYEVD when available.
module xtb_mctc_lapack_feast
   use xtb_mctc_accuracy, only : dp
#ifdef WITH_MKL
   use mkl_solvers_ee
#endif
   implicit none
   private

   public :: feast_available, feast_syevd_dp, feast_sygvd_dp

   logical, parameter :: feast_available = &
#ifdef WITH_MKL
      .true.
#else
      .false.
#endif

contains

   subroutine feast_syevd_dp(jobz, uplo, n, a, lda, w, info)
      !! Attempt to solve the symmetric eigenproblem with FEAST. On success the
      !! semantics mirror LAPACK's DSYEVD: eigenvalues in `w`, eigenvectors in `a`.
      character(len=1), intent(in)    :: jobz
      character(len=1), intent(in)    :: uplo
      integer,          intent(in)    :: n
      integer,          intent(in)    :: lda
      real(dp),         intent(inout) :: a(lda, *)
      real(dp),         intent(out)   :: w(*)
      integer,          intent(out)   :: info

#ifdef WITH_MKL
      integer :: fpm(128)
      integer :: loop, m0, m_found
      real(dp) :: epsout
      real(dp) :: emin, emax
      real(dp), allocatable :: x(:, :)
      real(dp), allocatable :: res(:)
      logical :: need_vec

      info = 0

      call feastinit(fpm)
      fpm(1) = 0          ! silence
      fpm(2) = 12         ! contour points
      fpm(3) = 12         ! 10^-12 tolerance target
      fpm(4) = 8          ! max refinement loops
      fpm(5) = 1          ! recycle input subspace

      call compute_bounds(n, a, lda, emin, emax)

      m0 = n
      allocate(x(n, m0))
      x = a
      allocate(res(m0))

      call dfeast_syev(uplo, n, a, lda, fpm, epsout, loop, emin, emax, m0, w, x, m_found, res, info)

      if (info == 0 .and. m_found == n) then
         need_vec = (jobz == 'V' .or. jobz == 'v')
         if (need_vec) then
            a(:, 1:m0) = x(:, 1:m0)
         end if
      else
         if (info == 0) info = -7776 ! FEAST converged but did not return all eigenpairs
      end if

      deallocate(res)
      deallocate(x)
#else
      call touch_unused(jobz, uplo, n, lda, a, w)
      info = -1
#endif
   end subroutine feast_syevd_dp


   subroutine feast_sygvd_dp(itype, jobz, uplo, n, a, lda, b, ldb, w, info)
      !! Generalised symmetric-definite eigenproblem using FEAST. Mirrors DSYGVD
      !! semantics: on success, eigenvalues in w and eigenvectors in a.
      integer,          intent(in)    :: itype
      character(len=1), intent(in)    :: jobz
      character(len=1), intent(in)    :: uplo
      integer,          intent(in)    :: n
      integer,          intent(in)    :: lda, ldb
      real(dp),         intent(inout) :: a(lda, *)
      real(dp),         intent(inout) :: b(ldb, *)
      real(dp),         intent(out)   :: w(*)
      integer,          intent(out)   :: info

#ifdef WITH_MKL
      integer :: fpm(128)
      integer :: loop, m0, m_found
      real(dp) :: epsout
      real(dp) :: emin, emax
      real(dp), allocatable :: x(:, :)
      real(dp), allocatable :: res(:)
      logical :: need_vec

      info = 0

      call feastinit(fpm)
      fpm(1) = 0
      fpm(2) = 12
      fpm(3) = 12
      fpm(4) = 8
      fpm(5) = 1

      call compute_bounds(n, a, lda, emin, emax)

      m0 = n
      allocate(x(n, m0))
      x = a
      allocate(res(m0))

      call dfeast_sygv(uplo, n, a, lda, b, ldb, fpm, epsout, loop, emin, emax, m0, w, x, m_found, res, info, itype)

      if (info == 0 .and. m_found == n) then
         need_vec = (jobz == 'V' .or. jobz == 'v')
         if (need_vec) then
            a(:, 1:m0) = x(:, 1:m0)
         end if
      else
         if (info == 0) info = -7777
      end if

      deallocate(res)
      deallocate(x)
#else
      call touch_unused(jobz, uplo, n, lda, a, w)
      info = -1
#endif

   end subroutine feast_sygvd_dp

#ifdef WITH_MKL
   subroutine compute_bounds(n, a, lda, emin, emax)
      !! Gershgorin bounds used to span the full spectrum for FEAST.
      integer,  intent(in)  :: n
      integer,  intent(in)  :: lda
      real(dp), intent(in)  :: a(lda, *)
      real(dp), intent(out) :: emin, emax

      integer :: i
      real(dp) :: radius, diag, pad, maxabs

      emin = huge(0.0_dp)
      emax = -huge(0.0_dp)
      maxabs = 0.0_dp

      do i = 1, n
         radius = sum(abs(a(i, 1:n))) - abs(a(i, i))
         diag = a(i, i)
         emin = min(emin, diag - radius)
         emax = max(emax, diag + radius)
         maxabs = max(maxabs, abs(diag) + radius)
      end do

      pad = 1.0e-6_dp * max(1.0_dp, maxabs)
      emin = emin - pad
      emax = emax + pad

   end subroutine compute_bounds

#else

   subroutine touch_unused(jobz, uplo, n, lda, a, w)
      !! Helper to silence unused dummy warnings when MKL is not available.
      character(len=1), intent(in)    :: jobz, uplo
      integer,          intent(in)    :: n, lda
      real(dp),         intent(inout) :: a(lda, *)
      real(dp),         intent(inout) :: w(*)

      if (.false.) then
         a(1, 1) = a(1, 1)
         w(1) = w(1)
         if (jobz == uplo .and. n == lda) continue
      end if
   end subroutine touch_unused

#endif

end module xtb_mctc_lapack_feast
