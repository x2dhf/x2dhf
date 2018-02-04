! ***************************************************************************
! *                                                                         *
! *   Copyright (C) 1996 Leif Laaksonen, Dage Sundholm                      *
! *   Copyright (C) 1996-2010 Jacek Kobus <jkob@fizyka.umk.pl>              *
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
  integer :: i,ic,im,imu,imu_bext,in,nall,nall_p,nstart,nstart_p,nmuall,nmuall_p

  real (PREC), dimension(nni_p,*) ::  fbefore
  real (PREC), dimension(nni,*)   ::  fafter
  real (PREC) :: rinftol
  logical :: muchange, nuchange, gridchange, rinfchange

  ! Tolerance for change in rinf
  rinftol = 1e2*precis

  ! Figure out what has changed
  muchange = nmuall .ne. nmuall_p
  nuchange = nni .ne. nni_p
  gridchange = ngrids .ne. ngrids_p
  rinfchange = abs(rinf-rinf_p).gt.(rinftol*rinf_p)

  ! The previous grid contains nni_p points in \nu (\eta) direction and
  ! nmuall_pt points in \mu (\xi) direction

  !      print *,nmuall,ngrids,rinf,nni
  !      print *,nmuall_p,ngrids_p,rinf_p,nni_p
  !      write (*,'(2e30.16)') rinf,rinf_p
  if ((muchange .or. gridchange .or. rinfchange) .and. .not. nuchange) then
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

     nstart_p=1
     nstart  =1
     nall_p=nmuall_p
     nall  =nmuall
     call dointerp_mu (nni_p,nstart_p,nall_p,nstart,nall,fbefore,fafter)

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

        do in=1,nni_p
           do im=imu_bext,nmuall
              fafter(in,im)=fbefore(in,nmuall_p)
           end do
        end do
     end if

  elseif (nuchange .and. (.not. muchange .and. .not. gridchange .and. .not. rinfchange)) then
     !        interpolate in nu variable

     if (nni_p.gt.maxmu) then
        print *,"dointerp: error! array vmuq is too short"
        stop 'dointerp'
     end if

     do i=1,nni_p
        vmuq(i)=vni_p(i)
     end do

     if (ic.eq.1) then
        iord=iord_nu_orb
        kbeg=1
        kend=iord
        iord2=iord/2
        !            stop "dointerp: error check array dimension"
     elseif (ic.eq.2) then
        iord=iord_nu_coul
        kbeg=1
        kend=iord
        iord2=iord/2
     elseif (ic.eq.3) then
        iord=iord_nu_exch
        kbeg=1
        kend=iord
        iord2=iord/2
     else
        stop "invalid argument ic"
     end if
     call dointerp_nu (nmuall_p,fbefore,fafter)

  else
     if (nuchange .and. rinfchange) then
        write(*,'(/,1x,"dointerp: cannot perform interpolation since both the number of nu grid points " &
             & /11x "and the value of practical infinity change")')
        stop "dointerp"
     end if
     write(*,'(/1x,"dointerp: no interpolation performed; check input data")')
     stop "dointerp"
  end if
end subroutine dointerp
