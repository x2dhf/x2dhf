! ***************************************************************************
! *                                                                         *
! *   Copyright (C) 1996       Leif Laaksonen, Dage Sundholm                *
! *   Copyright (C) 1996-2010  Jacek Kobus <jkob@fizyka.umk.pl>             *
! *   Copyright (C) 2018-      Susi Lehtola                                 *
! *                                                                         *
! *   This program is free software; you can redistribute it and/or modify  *
! *   it under the terms of the GNU General Public License version 2 as     *
! *   published by the Free Software Foundation.                            *
! *                                                                         *
! ***************************************************************************
! ### dointerp ###
!
!	Performs interpolations of orbitals and potentials between grids.
!       Orbitals and potentials are treated separately.

subroutine dointerp (ic,nmuall_p,nmuall,fbefore,fafter)
  use params
  use discret
  use commons8
  use commons16

  implicit none
  integer :: i,ic,im,imu,imu_bext,in,nmuall,nmuall_p

  real (PREC), dimension(nni_p,*)           ::  fbefore
  real (PREC), dimension(nni,*)             ::  fafter
  real (PREC), dimension(:,:), allocatable  ::  fmiddle
  real (PREC) :: gtol

  logical :: muchange, nuchange, gridchange, rinfchange
  logical :: muinterp, nuinterp, usemiddle

  ! Relative tolerance for change in bond length or rinf
  gtol = 1e2*precis

  ! Figure out what has changed
  muchange = nmuall .ne. nmuall_p
  nuchange = nni .ne. nni_p
  gridchange = ngrids .ne. ngrids_p
  rinfchange = abs(rinf-rinf_p).gt.(gtol*rinf_p)

  ! Do we need to interpolate in mu?
  muinterp=(muchange .or. gridchange .or. rinfchange)

  ! For nu, we only do interpolation if the size of the nu grid changes.
  nuinterp=(nuchange)

  ! Allocate helper if necessary
  if(muinterp .and. nuinterp) then
     usemiddle = .true.
     allocate(fmiddle(nni_p,nmuall))
  else
     usemiddle = .false.
  end if

  ! The previous grid contains nni_p points in \nu (\eta) direction and
  ! nmuall_pt points in \mu (\xi) direction
  if (muinterp) then
     ! Prepare data for routines evaluating coefficients of the
     ! Lagrange polynomial

     if (ic.eq.1) then
        iord=iord_mu_orb
     elseif (ic.eq.2) then
        iord=iord_mu_coul
     elseif (ic.eq.3) then
        iord=iord_mu_exch
     else
        stop "invalid argument ic"
     end if

     ! Interpolate in mu variable
     kbeg=1
     kend=iord
     iord2=iord/2
     do imu=1,nmuall_p
        vmuq(imu)=vmu_p(imu)
     end do

     if (usemiddle) then
        call dointerp_mu (nni_p,nmuall_p,nmuall,fbefore,fmiddle)
     else
        call dointerp_mu (nni_p,nmuall_p,nmuall,fbefore,fafter)
     end if

     ! Extrapolation (needed when R_infy is being increased seems to
     ! cause degradation of accuracy. That is why these values are
     ! taken equal (for each nu) to those corresponding to nmuall_p,
     ! i.e to the values of the last column of fbefore

     if(vmu(nmuall).gt.vmu_p(nmuall_p)) then
        imu_bext=nmuall
        do imu=1,nmuall
           if (vmu(imu).ge.vmu_p(nmuall_p)) then
              imu_bext=imu
              exit
           end if
        end do

        if(usemiddle) then
           do in=1,nni_p
              do im=imu_bext,nmuall
                 fmiddle(in,im)=fbefore(in,nmuall_p)
              end do
           end do
        else
           do in=1,nni_p
              do im=imu_bext,nmuall
                 fafter(in,im)=fbefore(in,nmuall_p)
              end do
           end do
        end if
     end if
  end if

  if(nuinterp) then
     !        interpolate in nu variable

     if (nni_p.gt.maxmu) then
        print *,"dointerp: error! array vmuq is too short"
        stop 'dointerp'
     end if

     if (ic.eq.1) then
        iord=iord_nu_orb
     elseif (ic.eq.2) then
        iord=iord_nu_coul
     elseif (ic.eq.3) then
        iord=iord_nu_exch
     else
        stop "invalid argument ic"
     end if

     ! Interpolate in nu variable
     kbeg=1
     kend=iord
     iord2=iord/2
     do i=1,nni_p
        vmuq(i)=vni_p(i)
     end do

     if(usemiddle) then
        call dointerp_nu (nmuall_p,fmiddle,fafter)
     else
        call dointerp_nu (nmuall_p,fbefore,fafter)
     end if
  end if

  if(usemiddle) deallocate(fmiddle)

  ! If grids are the same, just use the same values
  if(.not. muinterp .and. .not. nuinterp) then
     do in=1,nni
        do im=1,nmuall
           fafter(in,im)=fbefore(in,im)
        end do
     end do
  end if

end subroutine dointerp
