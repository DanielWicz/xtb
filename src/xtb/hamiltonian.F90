! This file is part of xtb.
!
! Copyright (C) 2019-2020 Sebastian Ehlert
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

!> Implementation of the xTB core Hamiltonian
module xtb_xtb_hamiltonian
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_constants, only : pi
   use xtb_mctc_convert, only : evtoau
   use xtb_xtb_data, only : THamiltonianData
   use xtb_intgrad
   use xtb_lin
   use xtb_scc_core, only : shellPoly, h0scal
   use xtb_grad_core, only : dshellPoly
   use omp_lib, only : omp_get_wtime
   implicit none
   private

   public :: getSelfEnergy, build_SDQH0, build_dSDQH0, build_dSDQH0_noreset
   public :: count_dpint, count_qpint


   interface getSelfEnergy
      module procedure :: getSelfEnergyFlat
      module procedure :: getSelfEnergy2D
   end interface getSelfEnergy


   integer,private, parameter :: llao (0:3) = (/ 1, 3, 6,10/)
   integer,private, parameter :: llao2(0:3) = (/ 1, 3, 5, 7/)
   logical, private, save :: timing_init = .false.
   logical, private, save :: timing_enabled = .false.

contains

subroutine init_sdqh0_timing()
   character(len=8) :: env
   integer :: stat
   if (timing_init) return
   timing_init = .true.
   call get_environment_variable('XTB_SDQH0_TIMING', env, status=stat)
   if (stat == 0 .and. len_trim(env) > 0) timing_enabled = .true.
end subroutine init_sdqh0_timing


subroutine getSelfEnergyFlat(hData, nShell, at, cn, qat, selfEnergy, dSEdcn, dSEdq)
   type(THamiltonianData), intent(in) :: hData
   integer, intent(in) :: nShell(:)
   integer, intent(in) :: at(:)
   real(wp), intent(in), optional :: cn(:)
   real(wp), intent(in), optional :: qat(:)
   real(wp), intent(out) :: selfEnergy(:)
   real(wp), intent(out), optional :: dSEdcn(:)
   real(wp), intent(out), optional :: dSEdq(:)

   integer :: ind, iAt, iZp, iSh, nat

   nat = size(at)

   selfEnergy(:) = 0.0_wp
   if (present(dSEdcn)) dSEdcn(:) = 0.0_wp
   if (present(dSEdq)) dSEdq(:) = 0.0_wp
   ind = 0
   do iAt = 1, nat
      iZp = at(iAt)
      do iSh = 1, nShell(iZp)
         selfEnergy(ind+iSh) = hData%selfEnergy(iSh, iZp)
      end do
      ind = ind + nShell(iZp)
   end do
   if (present(dSEdcn) .and. present(cn)) then
      ind = 0
      do iAt = 1, nat
         iZp = at(iAt)
         do iSh = 1, nShell(iZp)
            selfEnergy(ind+iSh) = selfEnergy(ind+iSh) &
               & - hData%kCN(iSh, iZp) * cn(iAt)
            dSEdcn(ind+iSh) = -hData%kCN(iSh, iZp)
         end do
         ind = ind + nShell(iZp)
      end do
   end if
   if (present(dSEdq) .and. present(qat)) then
      ind = 0
      do iAt = 1, nat
         iZp = at(iAt)
         do iSh = 1, nShell(iZp)
            selfEnergy(ind+iSh) = selfEnergy(ind+iSh) &
               & - hData%kQShell(iSh,iZp)*qat(iAt) - hData%kQAtom(iZp)*qat(iAt)**2
            dSEdq(ind+iSh) = -hData%kQShell(iSh,iZp) - hData%kQAtom(iZp)*2*qat(iAt)
         end do
         ind = ind + nShell(iZp)
      end do
   end if

end subroutine getSelfEnergyFlat


subroutine getSelfEnergy2D(hData, nShell, at, cn, qat, selfEnergy, dSEdcn, dSEdq)
   type(THamiltonianData), intent(in) :: hData
   integer, intent(in) :: nShell(:)
   integer, intent(in) :: at(:)
   real(wp), intent(in), optional :: cn(:)
   real(wp), intent(in), optional :: qat(:)
   real(wp), intent(out) :: selfEnergy(:, :)
   real(wp), intent(out), optional :: dSEdcn(:, :)
   real(wp), intent(out), optional :: dSEdq(:, :)

   integer :: iAt, iZp, iSh, nat

   nat = size(at)

   selfEnergy(:, :) = 0.0_wp
   if (present(dSEdcn)) dSEdcn(:, :) = 0.0_wp
   if (present(dSEdq)) dSEdq(:, :) = 0.0_wp
   do iAt = 1, nat
      iZp = at(iAt)
      do iSh = 1, nShell(iZp)
         selfEnergy(iSh, iAt) = hData%selfEnergy(iSh, iZp)
      end do
   end do
   if (present(dSEdcn) .and. present(cn)) then
      do iAt = 1, nat
         iZp = at(iAt)
         do iSh = 1, nShell(iZp)
            selfEnergy(iSh, iAt) = selfEnergy(iSh, iAt) &
               & - hData%kCN(iSh, iZp) * cn(iAt)
            dSEdcn(iSh, iAt) = -hData%kCN(iSh, iZp)
         end do
      end do
   end if
   if (present(dSEdq) .and. present(qat)) then
      do iAt = 1, nat
         iZp = at(iAt)
         do iSh = 1, nShell(iZp)
            selfEnergy(iSh, iAt) = selfEnergy(iSh, iAt) &
               & - hData%kQShell(iSh,iZp)*qat(iAt) - hData%kQAtom(iZp)*qat(iAt)**2
            dSEdq(iSh, iAt) = -hData%kQShell(iSh,iZp) &
               & - hData%kQAtom(iZp)*2*qat(iAt)
         end do
      end do
   end if

end subroutine getSelfEnergy2D


!> Computes the dipole and quadrupole integrals and performs screening to
!  determine, which contribute to potential
subroutine build_SDQH0(nShell, hData, nat, at, nbf, nao, xyz, trans, selfEnergy, &
      & intcut, caoshell, saoshell, nprim, primcount, alp, cont, &
      & sint, dpint, qpint, H0, H0_noovlp)
   implicit none
   integer, intent(in) :: nShell(:)
   type(THamiltonianData), intent(in) :: hData
   !> # of atoms
   integer, intent(in)  :: nat
   !> # of spherical AOs (SAOs)
   integer, intent(in)  :: nao
   !> # of Cartesian AOs (CAOs)
   integer, intent(in)  :: nbf
   integer, intent(in)  :: at(nat)
   !> Cartesian coordinates
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(in)    :: trans(:, :)
   real(wp), intent(in) :: selfEnergy(:, :)
   !> Integral cutoff according to prefactor from Gaussian product theorem
   real(wp),intent(in)  :: intcut
   !> Map shell of atom to index in CAO space (lowest Cart. component is taken)
   integer, intent(in)  :: caoshell(:,:)
   !> Map shell of atom to index in SAO space (lowest m_l component is taken)
   integer, intent(in)  :: saoshell(:,:)
   integer, intent(in)  :: nprim(:)
   !> Index of first primitive (over entire system) of given CAO
   integer, intent(in)  :: primcount(:)
   real(wp),intent(in)  :: alp(:)
   real(wp),intent(in)  :: cont(:)
   !> Overlap integral matrix
   real(wp),intent(out) :: sint(nao,nao)
   !> Dipole integral matrix
   real(wp),intent(out) :: dpint(3,nao,nao)
   !> Quadrupole integral matrix
   real(wp),intent(out) :: qpint(6,nao,nao)
   !> Core Hamiltonian
   real(wp),intent(out) :: H0(:)
   !> Core Hamiltonian without overlap contribution
   real(wp),intent(out) :: H0_noovlp(:)


   integer i,j,k,l,m,ii,jj,ll,mm,kk,ki,kj,kl,mi,mj,ij
   real(wp) tmp1,tmp2,tmp3,tmp4,step,step2
   real(wp) dx,dy,dz,s00r,s00l,s00,alpj
   real(wp) skj,r1,r2,tt,t1,t2,t3,t4,thr2,f,ci,cc,cj,alpi,rab2,ab,est

   real(wp)  ra(3),rb(3),f1,f2,point(3)
   real(wp) dtmp(3),qtmp(6),ss(6,6),dd(3,6,6),qq(6,6,6),tmp(6,6)
   real(wp) :: sblk(7,7), hblk(7,7), dblk(3,7,7), qblk(6,7,7)
   integer ip,jp,iat,jat,izp,jzp,ish,jsh,icao,jcao,iao,jao,jshmax
   integer ishtyp,jshtyp,iptyp,jptyp,naoi,naoj,mli,mlj,iprim,jprim
   integer :: il, jl, itr
   real(wp) :: zi, zj, zetaij, km, hii, hjj, hav, shpoly, shpoly_acc, hav_shpoly
   integer :: ntrans
   real(wp) :: trans_time_sum, transform_time_sum, screen_time_sum, diag_time_sum
   real(wp) :: t_block, t_local, t_screen
   integer itt(0:3)
   parameter(itt  =(/0,1,4,10/))
   real(wp), parameter :: shpoly_tol = 1.0e-12_wp
   real(wp) :: saw(10)
   logical :: skip_pair
   integer, allocatable :: rowStart(:)


   call init_sdqh0_timing()
   trans_time_sum = 0.0_wp
   transform_time_sum = 0.0_wp
   screen_time_sum = 0.0_wp
   diag_time_sum = 0.0_wp

   ! integrals
   H0(:) = 0.0_wp
   H0_noovlp(:) = 0.0_wp
   sint = 0.0_wp
   dpint = 0.0_wp
   qpint = 0.0_wp
   ! --- Aufpunkt for moment operator
   point = 0.0_wp
   ntrans = size(trans, dim=2)
   allocate(rowStart(nao))
   rowStart(1) = 0
   do i = 2, nao
      rowStart(i) = rowStart(i-1) + (i-1)
   end do

   !$omp parallel do default(none) &
   !$omp shared(nat, xyz, at, nShell, hData, selfEnergy, caoshell, saoshell, &
   !$omp& nprim, primcount, alp, cont, intcut, trans, point, ntrans, rowStart) &
   !$omp private (iat,jat,izp,ci,ra,rb,saw,dx,dy,dz, &
   !$omp& rab2,jzp,ish,ishtyp,icao,naoi,iptyp, &
   !$omp& jsh,jshmax,jshtyp,jcao,naoj,jptyp,ss,dd,qq,shpoly, shpoly_acc, hav_shpoly,&
   !$omp& est,alpi,alpj,ab,iprim,jprim,ip,jp,il,jl,hii,hjj,km,zi,zj,zetaij,hav, &
   !$omp& mli,mlj,tmp,tmp1,tmp2,iao,jao,ii,jj,k,ij,itr,sblk,hblk,dblk,qblk, &
   !$omp& skip_pair,t_block, t_local, t_screen) &
   !$omp reduction(+:trans_time_sum, transform_time_sum, screen_time_sum) &
   !$omp shared(sint,dpint,qpint,H0,H0_noovlp,timing_enabled) &
   !$omp collapse(2) schedule(dynamic,32)
   do iat = 1, nat
      do jat = 1, nat
         if (jat >= iat) cycle
         ra(1:3) = xyz(1:3,iat)
         izp = at(iat)
         jzp = at(jat)
         do ish = 1, nShell(izp)
            ishtyp = hData%angShell(ish,izp)
            icao = caoshell(ish,iat)
            naoi = llao(ishtyp)
            iptyp = itt(ishtyp)
            do jsh = 1, nShell(jzp)
               jshtyp = hData%angShell(jsh,jzp)
               jcao = caoshell(jsh,jat)
               naoj = llao(jshtyp)
               jptyp = itt(jshtyp)
               mli = llao2(ishtyp)
               mlj = llao2(jshtyp)

               il = ishtyp+1
               jl = jshtyp+1
               ! diagonals are the same for all H0 elements
               hii = selfEnergy(ish, iat)
               hjj = selfEnergy(jsh, jat)

               ! we scale the two shells depending on their exponent
               zi = hData%slaterExponent(ish, izp)
               zj = hData%slaterExponent(jsh, jzp)
               zetaij = (2 * sqrt(zi*zj)/(zi+zj))**hData%wExp
               call h0scal(hData,il,jl,izp,jzp,hData%valenceShell(ish, izp).ne.0, &
                  & hData%valenceShell(jsh, jzp).ne.0,km)

               hav = 0.5_wp * km * (hii + hjj) * zetaij

               skip_pair = .false.
               shpoly_acc = 0.0_wp
               sblk(1:mlj,1:mli) = 0.0_wp
               hblk(1:mlj,1:mli) = 0.0_wp
               dblk(:,1:mlj,1:mli) = 0.0_wp
               qblk(:,1:mlj,1:mli) = 0.0_wp
               block
                  real(wp) :: shpoly_cache(ntrans)
                  real(wp) :: rb_cache(3,ntrans)
                  integer :: n_active, itr_loc
                  n_active = 0
                  if (timing_enabled) t_block = omp_get_wtime()
                  do itr = 1, ntrans
                     rb(1:3) = xyz(1:3,jat) + trans(:, itr)
                     dx = rb(1) - ra(1)
                     dy = rb(2) - ra(2)
                     dz = rb(3) - ra(3)
                     rab2 = dx*dx + dy*dy + dz*dz

                     ! distance dependent polynomial
                     if (timing_enabled) t_screen = omp_get_wtime()
                     shpoly=shellPoly(hData%shellPoly(il,izp),hData%shellPoly(jl,jzp),&
                        &             hData%atomicRad(izp),hData%atomicRad(jzp),ra,rb)
                     if (timing_enabled) screen_time_sum = screen_time_sum + (omp_get_wtime()-t_screen)

                     if (abs(shpoly) <= shpoly_tol) cycle
                     n_active = n_active + 1
                     rb_cache(:, n_active) = rb
                     shpoly_cache(n_active) = shpoly
                  end do

                  skip_pair = (n_active == 0)
                  if (.not. skip_pair) then
                     do itr_loc = 1, n_active
                        rb = rb_cache(:, itr_loc)
                        shpoly = shpoly_cache(itr_loc)
                        ss = 0.0_wp
                        dd = 0.0_wp
                        qq = 0.0_wp
                        call get_multiints(icao,jcao,naoi,naoj,ishtyp,jshtyp,ra,rb,point, &
                           &               intcut,nprim,primcount,alp,cont,ss,dd,qq)
                        !transform from CAO to SAO
                        if (timing_enabled) t_local = omp_get_wtime()
                        call dtrf2(ss,ishtyp,jshtyp)
                        if (timing_enabled) transform_time_sum = transform_time_sum + (omp_get_wtime()-t_local)
                        do k = 1,3
                           if (timing_enabled) t_local = omp_get_wtime()
                           tmp(1:6,1:6) = dd(k,1:6,1:6)
                           call dtrf2(tmp,ishtyp,jshtyp)
                           dd(k,1:6,1:6) = tmp(1:6,1:6)
                           if (timing_enabled) transform_time_sum = transform_time_sum + (omp_get_wtime()-t_local)
                        enddo
                        do k = 1,6
                           if (timing_enabled) t_local = omp_get_wtime()
                           tmp(1:6,1:6) = qq(k,1:6,1:6)
                           call dtrf2(tmp,ishtyp,jshtyp)
                           qq(k,1:6,1:6) = tmp(1:6,1:6)
                           if (timing_enabled) transform_time_sum = transform_time_sum + (omp_get_wtime()-t_local)
                        enddo
                        hblk(1:mlj,1:mli) = hblk(1:mlj,1:mli) + shpoly * ss(1:mlj,1:mli)
                        sblk(1:mlj,1:mli) = sblk(1:mlj,1:mli) + ss(1:mlj,1:mli)
                        dblk(:,1:mlj,1:mli) = dblk(:,1:mlj,1:mli) + dd(:,1:mlj,1:mli)
                        qblk(:,1:mlj,1:mli) = qblk(:,1:mlj,1:mli) + qq(:,1:mlj,1:mli)
                        shpoly_acc = shpoly_acc + shpoly
                     enddo
                  end if
                  if (timing_enabled) trans_time_sum = trans_time_sum + (omp_get_wtime()-t_block)
               end block
               if (skip_pair) cycle
               hav_shpoly = hav * shpoly_acc
               do ii = 1,mli
                  iao = ii+saoshell(ish,iat)
                  do jj = 1,mlj
                     jao = jj+saoshell(jsh,jat)
                     if (iao >= jao) then
                        ij = rowStart(iao) + jao
                     else
                        ij = rowStart(jao) + iao
                     end if
                     H0(ij) = H0(ij) + hav * hblk(jj, ii)
                     H0_noovlp(ij) = H0_noovlp(ij) + hav_shpoly
                     sint(jao, iao) = sint(jao, iao) + sblk(jj, ii)
                     dpint(:, jao, iao) = dpint(:, jao, iao) + dblk(:, jj, ii)
                     qpint(:, jao, iao) = qpint(:, jao, iao) + qblk(:, jj, ii)
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo
   !$omp parallel do default(none) shared(nao, sint, dpint, qpint) private(iao, jao)
   do iao = 1, nao
      !$omp simd
      do jao = 1, iao - 1
         sint(iao, jao) = sint(jao, iao)
         dpint(:, iao, jao) = dpint(:, jao, iao)
         qpint(:, iao, jao) = qpint(:, jao, iao)
      end do
   end do

   ! diagonal elements
   !$omp parallel do default(none) schedule(dynamic) &
   !$omp shared(H0, H0_noovlp, sint, dpint, qpint, timing_enabled) &
   !$omp shared(nat, xyz, at, nShell, hData, saoshell, selfEnergy, caoshell, &
   !$omp& point, intcut, nprim, primcount, alp, cont, rowStart) &
   !$omp private(iat, ra, izp, ish, ishtyp, iao, i, ii, icao, naoi, iptyp, &
   !$omp& jsh, jshtyp, jcao, ss, dd, qq, jao, jj, naoj, jptyp, mli, mlj, dblk, qblk, tmp, &
   !$omp& t_block, t_local) &
   !$omp reduction(+:diag_time_sum, transform_time_sum)
   do iat = 1, nat
      if (timing_enabled) t_block = omp_get_wtime()
      ra = xyz(:, iat)
      izp = at(iat)
      do ish = 1, nShell(izp)
         ishtyp = hData%angShell(ish,izp)
         do iao = 1, llao2(ishtyp)
            i = iao+saoshell(ish,iat)
            ii = rowStart(i) + i
            sint(i,i) = 1.0_wp + sint(i,i)
            H0(ii) = H0(ii) + selfEnergy(ish, iat)
            H0_noovlp(ii) = H0_noovlp(ii) + selfEnergy(ish, iat)
         end do

         icao = caoshell(ish,iat)
         naoi = llao(ishtyp)
         iptyp = itt(ishtyp)
         do jsh = 1, ish
            jshtyp = hData%angShell(jsh,izp)
            jcao = caoshell(jsh,iat)
            naoj = llao(jshtyp)
            jptyp = itt(jshtyp)
            ss = 0.0_wp
            dd = 0.0_wp
            qq = 0.0_wp
            call get_multiints(icao,jcao,naoi,naoj,ishtyp,jshtyp,ra,ra,point, &
               &               intcut,nprim,primcount,alp,cont,ss,dd,qq)
            do k = 1,3
               if (timing_enabled) t_local = omp_get_wtime()
               tmp(1:6, 1:6) = dd(k,1:6, 1:6)
               call dtrf2(tmp, ishtyp, jshtyp)
               dd(k, 1:6, 1:6) = tmp(1:6, 1:6)
               if (timing_enabled) transform_time_sum = transform_time_sum + (omp_get_wtime()-t_local)
            enddo
            do k = 1,6
               if (timing_enabled) t_local = omp_get_wtime()
               tmp(1:6, 1:6) = qq(k, 1:6, 1:6)
               call dtrf2(tmp, ishtyp, jshtyp)
               qq(k, 1:6, 1:6) = tmp(1:6, 1:6)
               if (timing_enabled) transform_time_sum = transform_time_sum + (omp_get_wtime()-t_local)
            enddo
            mli = llao2(ishtyp)
            mlj = llao2(jshtyp)
            dblk(:,1:mlj,1:mli) = dd(:,1:mlj,1:mli)
            qblk(:,1:mlj,1:mli) = qq(:,1:mlj,1:mli)
            do ii = 1, mli
               iao = ii + saoshell(ish,iat)
               do jj = 1, mlj
                  jao = jj + saoshell(jsh,iat)
                  if (jao > iao .and. ish ==  jsh) cycle
                  dpint(1:3, iao, jao) = dpint(1:3, iao, jao) + dblk(1:3, jj, ii)
                  if (iao /= jao) then
                     dpint(1:3, jao, iao) = dpint(1:3, jao, iao) + dblk(1:3, jj, ii)
                  end if
                  qpint(1:6, iao, jao) = qpint(1:6, iao, jao) + qblk(1:6, jj, ii)
                  if (jao /= iao) then
                     qpint(1:6, jao, iao) = qpint(1:6, jao, iao) + qblk(1:6, jj, ii)
                  end if
               end do
            end do
         end do
         if (timing_enabled) diag_time_sum = diag_time_sum + (omp_get_wtime()-t_block)
      end do
   end do

   if (timing_enabled) then
      write(*,'(a,3f12.6)') 'build_SDQH0 timing [translation, transform, screening] (s):', &
         trans_time_sum, transform_time_sum, screen_time_sum
      write(*,'(a,f12.6)') 'build_SDQH0 diagonal (s):', diag_time_sum
   end if
   deallocate(rowStart)

end subroutine build_SDQH0


!> Computes the gradient of the dipole/qpole integral contribution
subroutine build_dSDQH0(nShell, hData, selfEnergy, dSEdcn, intcut, nat, nao, nbf, &
      & at, xyz, trans, caoshell, saoshell, nprim, primcount, alp, cont, &
      & p, Pew, ves, vs, vd, vq, dhdcn, g, sigma)
   integer, intent(in) :: nShell(:)
   type(THamiltonianData), intent(in) :: hData
   real(wp), intent(in) :: selfEnergy(:, :)
   real(wp), intent(in) :: dSEdcn(:, :)
   !> # of atoms
   integer, intent(in)    :: nat
   !> # of spherical AOs (SAOs)
   integer, intent(in)    :: nao
   !> # of Cartesian AOs (CAOs)
   integer, intent(in)    :: nbf
   !> Atomic numbers of atoms
   integer, intent(in)    :: at(nat)
   !> Integral cutoff according to prefactor from Gaussian product theorem
   real(wp),intent(in)    :: intcut
   real(wp),intent(in)    :: ves(:, :)
   real(wp),intent(in)    :: vs(nat)
   real(wp),intent(in)    :: vd(3,nat)
   real(wp),intent(in)    :: vq(6,nat)
   !> Cartesian coordinates
   real(wp),intent(in)    :: xyz(3,nat)
   real(wp),intent(in)    :: trans(:, :)
   !> Map shell of atom to index in CAO space (lowest Cart. component is taken)
   integer, intent(in)    :: caoshell(:,:)
   !> Map shell of atom to index in SAO space (lowest m_l component is taken)
   integer, intent(in)    :: saoshell(:,:)
   integer, intent(in)    :: nprim(:)
   !> Index of first primitive (over entire system) of given CAO, dimension: nbf
   integer, intent(in)    :: primcount(:)
   real(wp),intent(in)    :: alp(:)
   real(wp),intent(in)    :: cont(:)
   !> Density matrix
   real(wp),intent(in) :: p(:, :)
   !> Energy weighted density matrix
   real(wp),intent(in) :: Pew(:, :)
   real(wp),intent(inout) :: g(:, :)
   real(wp),intent(inout) :: sigma(:, :)
   real(wp),intent(inout) :: dhdcn(:)

   integer itt(0:3)
   parameter (itt=(/0,1,4,10/))
   real(wp) tmp1,tmp2,tmp3,step,step2,step3,s00r,s00l,s00,alpj
   real(wp) skj,r1,r2,tt,t1,t2,t3,t4,thr2,f,ci,cj,alpi,rij2,ab,est
   real(wp) f1,f2,point(3),tmp(6,6),rij(3),ri(3),rj(3)
   real(wp) stmp,ral(3,3),rar(3,3),rbl(3,3),pre
   real(wp) dtmp,qtmp,rbr(3,3),r2l(3),r2r(3),qqa(6,6,6,3)
   real(wp)  ss(6,6,3),ddc(3,6,6,3),qqc(6,6,6,3),dda(3,6,6,3)
   integer i,j,k,l,m,ii,jj,ll,mm,kk,ki,kj,kl,mi,mj,ij,jshmax
   integer ip,jp,iat,jat,izp,jzp,ish,jsh,icao,jcao,iao,jao,ixyz
   integer ishtyp,jshtyp,iptyp,jptyp,naoi,naoj,mli,mlj,iprim,jprim
   real(wp) :: dumdum(3,10),dum(10),sdq(10,6,6),sdqg(3,19,6,6)
   integer :: il, jl, itr
   real(wp) :: zi, zj, zetaij, km, hii, hjj, hav, shpoly, dshpoly(3)
   real(wp) :: Pij, Hij, HPij, g_xyz(3)
   real(wp), parameter :: rthr = 1600.0_wp

   thr2 = intcut
   point = 0.0_wp
   ! call timing(t1,t3)
   !$omp parallel do default(none) &
   !$omp shared(nat, at, xyz, trans, nShell, hData, selfEnergy, dSEdcn, P, Pew, &
   !$omp& ves, vs, vd, vq, intcut, nprim, primcount, caoshell, saoshell, alp, cont) &
   !$omp private(iat,jat,ixyz,izp,ci,rij2,jzp,ish,ishtyp, &
   !$omp& icao,naoi,iptyp,jsh,jshmax,jshtyp,jcao,naoj,jptyp, &
   !$omp& sdq,sdqg,est,alpi,alpj,ab,iprim,jprim,ip,jp,ri,rj,rij,km,shpoly,dshpoly, &
   !$omp& mli,mlj,dum,dumdum,tmp,stmp,dtmp,qtmp,il,jl,zi,zj,zetaij,hii,hjj,hav, &
   !$omp& iao,jao,ii,jj,k,pij,hij,hpij,g_xyz,itr) &
   !$omp reduction(+:g,sigma,dhdcn) &
   !$omp collapse(2) schedule(dynamic,32)
   do iat = 1,nat
      do jat = 1,nat
         if (jat >= iat) cycle
         ri = xyz(:,iat)
         jzp = at(jat)
         izp = at(iat)

         do ish = 1,nShell(izp)
            ishtyp = hData%angShell(ish,izp)
            icao = caoshell(ish,iat)
            naoi = llao(ishtyp)
            iptyp = itt(ishtyp)
            jshmax = nShell(jzp)
            !              if(iat.eq.jat) jshmax = ish
            do jsh = 1,jshmax ! jshells
               jshtyp = hData%angShell(jsh,jzp)
               jcao = caoshell(jsh,jat)
               naoj = llao(jshtyp)
               jptyp = itt(jshtyp)

               il = ishtyp+1
               jl = jshtyp+1
               ! diagonals are the same for all H0 elements
               hii = selfEnergy(ish, iat)
               hjj = selfEnergy(jsh, jat)

               ! we scale the two shells depending on their exponent
               zi = hData%slaterExponent(ish, izp)
               zj = hData%slaterExponent(jsh, jzp)
               zetaij = (2 * sqrt(zi*zj)/(zi+zj))**hData%wExp
               call h0scal(hData,il,jl,izp,jzp,hData%valenceShell(ish, izp).ne.0, &
                  & hData%valenceShell(jsh, jzp).ne.0,km)

               ! averaged H0 element (without overlap contribution!)
               hav = 0.5_wp * km * (hii + hjj) * zetaij * evtoau

               do itr = 1, size(trans, dim=2)
                  rj = xyz(:,jat) + trans(:, itr)
                  rij = ri - rj
                  rij2 =  sum( rij**2 )

                  if (rij2 > rthr) cycle

                  ! distance dependent polynomial
                  call dshellPoly(hData%shellPoly(il,izp),hData%shellPoly(jl,jzp),&
                     & hData%atomicRad(izp),hData%atomicRad(jzp),rij2,ri,rj,&
                     & shpoly,dshpoly)

                  sdqg = 0;sdq = 0
                  call get_grad_multiint(icao,jcao,naoi,naoj,ishtyp,jshtyp,ri,rj, &
                     &                   intcut,nprim,primcount,alp,cont,sdq,sdqg)
                  tmp(1:6,1:6) = sdq(1,1:6,1:6)
                  call dtrf2(tmp,ishtyp,jshtyp)
                  sdq(1,1:6,1:6) = tmp(1:6,1:6)
                  do k = 1,19 ! 1 S, 2-4 D, 5-10 Q, 11-13 D, 14-19 Q
                     do ixyz = 1,3
                        ! transform from CAO to SAO
                        !                        call dtrf2(sdqg(ixyz,k,1:6,1:6),ishtyp,jshtyp)
                        tmp(1:6,1:6) = sdqg(ixyz,k,1:6,1:6)
                        call dtrf2(tmp,ishtyp,jshtyp)
                        sdqg(ixyz,k,1:6,1:6) = tmp(1:6,1:6)
                     enddo
                  enddo
                  g_xyz(:) = 0.0_wp
                  do ii = 1,llao2(ishtyp)
                     iao = ii+saoshell(ish,iat)
                     do jj = 1,llao2(jshtyp)
                        jao = jj+saoshell(jsh,jat)
                        Pij = p(jao,iao)

                        ! Hamiltonian element without overlap
                        Hij  = hav * shpoly
                        HPij = Hij * Pij

                        g_xyz(:) = g_xyz + 2*HPij*sdq(1,jj,ii)*dshpoly/shpoly

                        do ixyz = 1,3
                           stmp = sdqg(ixyz,1,jj,ii)*(2*HPij - 2*Pew(jao, iao) &
                              & -Pij*(ves(ish,iat)+ves(jsh,jat)) &
                              & +Pij*(vs(iat)+vs(jat)))
                           dtmp = Pij*sum(sdqg(ixyz,11:13,jj,ii)*vd(1:3,iat) &
                              & +sdqg(ixyz, 2:4, jj,ii)*vd(1:3,jat) )
                           qtmp = Pij*sum( sdqg(ixyz,14:19,jj,ii)*vq(1:6,iat) &
                              & +sdqg(ixyz, 5:10,jj,ii)*vq(1:6,jat) )
                           g_xyz(ixyz) = g_xyz(ixyz)+stmp+dtmp+qtmp

                        enddo ! ixyz

                        ! Hamiltonian without Hav
                        HPij = km * zetaij * shpoly * Pij * sdq(1,jj,ii) * evtoau
                        ! save dE/dCN for CNi
                        dhdcn(iat) = dhdcn(iat) + HPij*dSEdcn(ish, iat)
                        ! save dE/dCN for CNj
                        dhdcn(jat) = dhdcn(jat) + HPij*dSEdcn(jsh, jat)
                     enddo
                  enddo
                  g(:,iat) = g(:,iat)+g_xyz
                  g(:,jat) = g(:,jat)-g_xyz
                  sigma(:, 1) = sigma(:, 1) + g_xyz(1) * rij
                  sigma(:, 2) = sigma(:, 2) + g_xyz(2) * rij
                  sigma(:, 3) = sigma(:, 3) + g_xyz(3) * rij
               enddo ! lattice translations
            enddo ! jsh : loop over shells on jat
         enddo  ! ish : loop over shells on iat
      enddo ! jat
   enddo  ! iat

   ! diagonal contributions
   !$omp parallel do default(none) schedule(dynamic) reduction(+:dhdcn) &
   !$omp shared(nat, at, nshell, hData, saoshell, P, dSEdcn) &
   !$omp private(iat, izp, ish, ishtyp, iao, i, Pij)
   do iat = 1, nat
      izp = at(iat)
      do ish = 1, nShell(izp)
         ishtyp = hData%angShell(ish,izp)
         do iao = 1, llao2(ishtyp)
            i = iao+saoshell(ish,iat)

            Pij = P(i,i)
            ! save dE/dCN for CNi
            dhdcn(iat) = dhdcn(iat) + Pij*dSEdcn(ish, iat)*evtoau
         end do
      end do
   end do

end subroutine build_dSDQH0


!> Computes the gradient of the dipole/qpole integral contribution
subroutine build_dSDQH0_noreset(nShell, hData, selfEnergy, dSEdcn, intcut, &
      & nat, nao, nbf, at, xyz, caoshell, saoshell, nprim, primcount, &
      & alp, cont, H0, S, p, Pew, ves, vs, vd, vq, dhdcn, g, sigma)
   integer, intent(in) :: nShell(:)
   type(THamiltonianData), intent(in) :: hData
   real(wp), intent(in) :: selfEnergy(:, :)
   real(wp), intent(in) :: dSEdcn(:, :)
   !> # of atoms
   integer, intent(in)    :: nat
   !> # of spherical AOs (SAOs)
   integer, intent(in)    :: nao
   !> # of Cartesian AOs (CAOs)
   integer, intent(in)    :: nbf
   !> Atomic numbers of atoms
   integer, intent(in)    :: at(nat)
   !> Integral cutoff according to prefactor from Gaussian product theorem
   real(wp),intent(in)    :: intcut
   real(wp),intent(in)    :: ves(:, :)
   real(wp),intent(in)    :: vs(nat)
   real(wp),intent(in)    :: vd(3,nat)
   real(wp),intent(in)    :: vq(6,nat)
   !> Cartesian coordinates
   real(wp),intent(in)    :: xyz(3,nat)
   !> Map shell of atom to index in CAO space (lowest Cart. component is taken)
   integer, intent(in)    :: caoshell(:,:)
   !> Map shell of atom to index in SAO space (lowest m_l component is taken)
   integer, intent(in)    :: saoshell(:,:)
   integer, intent(in)    :: nprim(:)
   !> Index of first primitive (over entire system) of given CAO, dimension: nbf
   integer, intent(in)    :: primcount(:)
   real(wp),intent(in)    :: alp(:)
   real(wp),intent(in)    :: cont(:)
   !> Molecular Hamiltonian (packed lower-triangular storage, in eV)
   real(wp),intent(in) :: H0(:)
   !> Molecular Overlap
   real(wp),intent(in) :: S(:, :)
   !> Density matrix
   real(wp),intent(in) :: p(:, :)
   !> Energy weighted density matrix
   real(wp),intent(in) :: Pew(:, :)
   real(wp),intent(inout) :: g(:, :)
   real(wp),intent(inout) :: sigma(:, :)
   real(wp),intent(inout) :: dhdcn(:)

   integer itt(0:3)
   parameter (itt=(/0,1,4,10/))
   real(wp) tmp1,tmp2,tmp3,step,step2,step3,s00r,s00l,s00,alpj
   real(wp) skj,r1,r2,tt,t1,t2,t3,t4,thr2,f,ci,cj,alpi,rij2,ab,est
   real(wp) f1,f2,point(3),tmp(6,6),rij(3),ri(3),rj(3)
   real(wp) ral(3,3),rar(3,3),rbl(3,3),pre
   real(wp) dtmp,qtmp,rbr(3,3),r2l(3),r2r(3),qqa(6,6,6,3)
   real(wp)  ss(6,6,3),ddc(3,6,6,3),qqc(6,6,6,3),dda(3,6,6,3)
   integer i,j,k,l,m,ii,jj,ll,mm,kk,ki,kj,kl,mi,mj,ij,jshmax
   integer ip,jp,iat,jat,izp,jzp,ish,jsh,icao,jcao,iao,jao,ixyz
   integer ishtyp,jshtyp,iptyp,jptyp,naoi,naoj,mli,mlj,iprim,jprim
   real(wp) :: dumdum(3,10),dum(10),sdq(10,6,6),sdqg(3,19,6,6)
   integer :: il, jl, itr
   real(wp) :: zi, zj, zetaij, km, hii, hjj, hav, shpoly, dshpoly(3), dCN
   real(wp) :: Pij, Hij, HPij, g_xyz(3)
   real(wp), parameter :: rthr = 1600.0_wp
   integer, allocatable :: rowStart(:)

   ! local OpenMP variables
!$ real(wp), allocatable :: g_omp(:, :), sigma_omp(:, :), dhdcn_omp(:)

   thr2 = intcut
   point = 0.0_wp
   allocate(rowStart(nao))
   rowStart(1) = 0
   do i = 2, nao
      rowStart(i) = rowStart(i-1) + (i-1)
   end do
   ! call timing(t1,t3)
   !$omp parallel default(none) &
   !$omp shared(nat, at, xyz, nShell, hData, selfEnergy, dSEdcn, P, Pew, &
   !$omp& H0, S, ves, vs, vd, vq, intcut, nprim, primcount, caoshell, saoshell, &
   !$omp& alp, cont, g, sigma, dhdcn, rowStart) &
   !$omp private(iat,jat,ixyz,izp,ci,rij2,jzp,ish,ishtyp,ij,i, &
   !$omp& icao,naoi,iptyp,jsh,jshmax,jshtyp,jcao,naoj,jptyp,dCN, &
   !$omp& sdq,sdqg,est,alpi,alpj,ab,iprim,jprim,ip,jp,ri,rj,rij,km,shpoly,dshpoly, &
   !$omp& mli,mlj,dum,dumdum,tmp,dtmp,qtmp,il,jl,zi,zj,zetaij,hii,hjj,hav, &
   !$omp& iao,jao,ii,jj,k,pij,hij,hpij,g_xyz,itr, g_omp, sigma_omp, dhdcn_omp)

!$ allocate(g_omp(size(g, dim=1), size(g, dim=2)), source = 0.0_wp)
!$ allocate(sigma_omp(size(sigma, dim=1), size(sigma, dim=2)), source = 0.0_wp)
!$ allocate(dhdcn_omp(size(dhdcn, dim=1)), source = 0.0_wp)

#ifndef _OPENMP
   associate(g_omp => g, sigma_omp => sigma, dhdcn_omp => dhdcn)
#endif

   !$omp do collapse(2) schedule(dynamic,32)
   do iat = 1,nat
      do jat = 1,nat
         if (jat >= iat) cycle
         izp = at(iat)
         jzp = at(jat)

         ri = xyz(:,iat)
         rj = xyz(:,jat)
         rij = ri - rj
         rij2 =  sum( rij**2 )

         if (rij2 > rthr) cycle
         do ish = 1,nShell(izp)
            ishtyp = hData%angShell(ish,izp)
            icao = caoshell(ish,iat)
            naoi = llao(ishtyp)
            iptyp = itt(ishtyp)
            jshmax = nShell(jzp)
            !              if(iat.eq.jat) jshmax = ish
            do jsh = 1,jshmax ! jshells
               jshtyp = hData%angShell(jsh,jzp)
               jcao = caoshell(jsh,jat)
               naoj = llao(jshtyp)
               jptyp = itt(jshtyp)

               il = ishtyp+1
               jl = jshtyp+1
               ! diagonals are the same for all H0 elements
               hii = selfEnergy(ish, iat)
               hjj = selfEnergy(jsh, jat)

               ! distance dependent polynomial
               call dshellPoly(hData%shellPoly(il,izp),hData%shellPoly(jl,jzp),&
                  & hData%atomicRad(izp),hData%atomicRad(jzp),rij2,ri,rj,&
                  & shpoly,dshpoly)

               ! averaged H0 element (without overlap contribution!)
               hav = 0.5_wp * (hii + hjj)

               sdqg = 0;sdq = 0
               call get_grad_multiint(icao,jcao,naoi,naoj,ishtyp,jshtyp,ri,rj, &
                  &                   intcut,nprim,primcount,alp,cont,sdq,sdqg)
               do k = 1,19 ! 1 S, 2-4 D, 5-10 Q, 11-13 D, 14-19 Q
                  do ixyz = 1,3
                     ! transform from CAO to SAO
                     !                        call dtrf2(sdqg(ixyz,k,1:6,1:6),ishtyp,jshtyp)
                     tmp(1:6,1:6) = sdqg(ixyz,k,1:6,1:6)
                     call dtrf2(tmp,ishtyp,jshtyp)
                     sdqg(ixyz,k,1:6,1:6) = tmp(1:6,1:6)
                  enddo
               enddo
               g_xyz(:) = 0.0_wp
               dCN = 0.0_wp
               do ii = 1,llao2(ishtyp)
                  iao = ii+saoshell(ish,iat)
                  do jj = 1,llao2(jshtyp)
                     jao = jj+saoshell(jsh,jat)
                     Pij = p(jao,iao)

                     ! Hamiltonian element without overlap
                     if (jao >= iao) then
                        ij = rowStart(jao) + iao
                     else
                        ij = rowStart(iao) + jao
                     end if
                     Hij  = H0(ij) * evtoau
                     HPij = Hij * Pij

                     g_xyz(:) = g_xyz + 2*HPij*S(jao,iao)*dshpoly/shpoly &
                        & + sdqg(:,1,jj,ii)*(2*HPij - 2*Pew(jao, iao) &
                        & - Pij*(ves(ish,iat)+ves(jsh,jat)) &
                        & + Pij*(vs(iat)+vs(jat)))

                     do ixyz = 1,3
                        dtmp = Pij*sum(sdqg(ixyz,11:13,jj,ii)*vd(1:3,iat) &
                           & +sdqg(ixyz, 2:4, jj,ii)*vd(1:3,jat) )
                        qtmp = Pij*sum( sdqg(ixyz,14:19,jj,ii)*vq(1:6,iat) &
                           & +sdqg(ixyz, 5:10,jj,ii)*vq(1:6,jat) )
                        g_xyz(ixyz) = g_xyz(ixyz)+dtmp+qtmp

                     enddo ! ixyz

                     ! Hamiltonian without Hav
                     dCN = dCN + HPij / hav * S(jao, iao)
                  enddo
               enddo
               ! save dE/dCN for CNi
               dhdcn_omp(iat) = dhdcn_omp(iat) + dCN*dSEdcn(ish, iat)
               ! save dE/dCN for CNj
               dhdcn_omp(jat) = dhdcn_omp(jat) + dCN*dSEdcn(jsh, jat)
               g_omp(:,iat) = g_omp(:,iat)+g_xyz
               g_omp(:,jat) = g_omp(:,jat)-g_xyz
               sigma_omp(:, 1) = sigma_omp(:, 1) + g_xyz(1) * rij
               sigma_omp(:, 2) = sigma_omp(:, 2) + g_xyz(2) * rij
               sigma_omp(:, 3) = sigma_omp(:, 3) + g_xyz(3) * rij
            enddo ! jsh : loop over shells on jat
         enddo  ! ish : loop over shells on iat
      enddo ! jat
   enddo  ! iat
   !$omp end do nowait

   ! diagonal contributions
   !$omp do schedule(dynamic)
   do iat = 1, nat
      izp = at(iat)
      do ish = 1, nShell(izp)
         ishtyp = hData%angShell(ish,izp)
         do iao = 1, llao2(ishtyp)
            i = iao+saoshell(ish,iat)

            Pij = P(i,i)
            ! save dE/dCN for CNi
            dhdcn_omp(iat) = dhdcn_omp(iat) + Pij*dSEdcn(ish, iat)*evtoau
         end do
      end do
   end do
   !$omp end do nowait

#ifndef _OPENMP
   end associate
#endif

   !$omp critical (g_crt)
!$ g(:,:) = g + g_omp
   !$omp end critical (g_crt)
   !$omp critical (sigma_crt)
!$ sigma(:,:) = sigma + sigma_omp
   !$omp end critical (sigma_crt)
   !$omp critical (dhdcn_crt)
!$ dhdcn(:) = dhdcn + dhdcn_omp
   !$omp end critical (dhdcn_crt)

   !$omp end parallel

   deallocate(rowStart)

end subroutine build_dSDQH0_noreset


!> Count number of significant dipole integrals
subroutine count_dpint(ndp, dpint, thr)

   !> Number of significant dipole integrals
   integer, intent(out) :: ndp

   !> Dipole integrals
   real(wp), intent(in) :: dpint(:, :, :)

   !> Neglect threshold for dipole integrals
   real(wp), intent(in) :: thr

   integer :: i, j
   real(wp) :: tmp1, thr2

   ndp = 0
   thr2 = (thr*1.0e-2_wp)-thr*1.0e-12_wp

!$omp parallel do default(shared) private(j,tmp1) reduction(+:ndp) schedule(static)
   do i = 1, size(dpint, dim=3)
      do j = 1, i
         tmp1 = dot_product(dpint(1:3, j, i), dpint(1:3, j, i))
         if (tmp1 > thr2) ndp = ndp + 1
      enddo
   enddo
!$omp end parallel do

end subroutine count_dpint


!> Count number of significant quadrupole integrals
subroutine count_qpint(nqp, qpint, thr)

   !> Number of significant quadrupole integrals
   integer, intent(out) :: nqp

   !> Quadrupole integrals
   real(wp), intent(in) :: qpint(:, :, :)

   !> Neglect threshold for quadrupole integrals
   real(wp), intent(in) :: thr

   integer :: i, j
   real(wp) :: tmp2, thr2

   nqp = 0
   thr2 = (thr*1.0e-2_wp)-thr*1.0e-12_wp

!$omp parallel do default(shared) private(j,tmp2) reduction(+:nqp) schedule(static)
   do i = 1, size(qpint, dim=3)
      do j = 1, i
         tmp2 = dot_product(qpint(1:3, j, i), qpint(1:3, j, i)) &
            & + 2.0_wp*dot_product(qpint(4:6, j, i), qpint(4:6, j, i))
         if (tmp2 > thr2) nqp = nqp + 1
      enddo
   enddo
!$omp end parallel do

end subroutine count_qpint


end module xtb_xtb_hamiltonian
