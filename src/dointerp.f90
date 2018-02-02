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
  integer :: i,ic,im,imu,imu_bext,in,nall,nall_p,ninner,ninner_p,nmuall,nmuall_p,nstart,nstart_p

  real (PREC), dimension(nni_p,*) ::  fbefore,fafter
  real (PREC) :: rinftol

  ! Tolerance for change in rinf
  rinftol = 1e4*precis

  !  The previous grid contains nni_p points in \nu (\eta) direction and
  !  nmuall_pt points in \mu (\xi) direction

  !      print *,nmuall,ngrids,rinf,nni
  !      print *,nmuall_p,ngrids_p,rinf_p,nni_p
  !      write (*,'(2e30.16)') rinf,rinf_p

  if ((nmuall.ne.nmuall_p.or.ngrids.ne.ngrids_p.or. &
       abs(rinf-rinf_p).gt.rinftol).and.nni.eq.nni_p) then

     ! rinf and rinf_p might differ a bit because of round-off errors

     if (nmuall.ne.nmuall_p.and.abs(rinf/rinf_p-one).gt.rinftol) then
        write(*,'(/,1x,"dointerp: cannot perform interpolation" &
             & " since both the number of mu grid points "/11x    &
             & "and the value of practical infinity change")')
        stop "dointerp"
     endif

     !  prepare data for routines evaluating coefficients of the Lagrange
     !  polynomial

     if (ic.eq.1) then

        !  Interpolate in mu variable

        !  Interpolation of orbitals is done separately in two regions in mu
        !  variable. The inner region extends from mu=1 to the last maximum
        !  and the outer (tail) region extends from that maximum to the practical
        !  infinity. The values in the latter are logarithmed before
        !  the interpolation is carried out.

        !  Find the tail region: the maximum is sought along the line which is
        !  almost perpendicular to the intermolecular axis (to avoid picking up
        !  zeros for ungerade orbitals)

        !  The above scheme turned out to be unnecessary

        iord=iord_mu_orb
        kbeg=1
        kend=iord
        iord2=iord/2
        do imu=1,nmuall_p
           vmuq(imu)=vmu_p(imu)
        enddo

        nall_p=nmuall_p
        nall  =nmuall
        nstart_p  =1
        nstart    =1
        call dointerp_mu (nni_p,nstart_p,nall_p,nstart,nall,fbefore,fafter)

     elseif (ic.eq.2) then
        iord=iord_mu_coul
        kbeg=1
        kend=iord
        iord2=iord/2
        do imu=1,nmuall_p
           vmuq(imu)=vmu_p(imu)
        enddo

        ninner_p=1
        ninner  =1
        nall_p=nmuall_p
        nall  =nmuall
        call dointerp_mu (nni_p,ninner_p,nall_p,ninner,nall,fbefore,fafter)

     elseif (ic.eq.3) then
        iord=iord_mu_exch
        kbeg=1
        kend=iord
        iord2=iord/2
        do imu=1,nmuall_p
           vmuq(imu)=vmu_p(imu)
        enddo

        ninner_p=1
        ninner  =1
        nall_p=nmuall_p
        nall  =nmuall
        call dointerp_mu (nni_p,ninner_p,nall_p,ninner,nall,fbefore,fafter)
     endif

     !        Extrapolation (needed when R_infy is being increased seems to
     !        cause degradation of accuracy. That is why these values are
     !        taken equal (for each nu) to those corresponding to nmuall_p,
     !        i.e to the values of the last column of fbefore

     imu_bext=nmuall
     do imu=1,nmuall
        if (vmu(imu).ge.vmu_p(nmuall_p)) then
           imu_bext=imu
           goto 1234
        endif
     enddo
01234 continue
     do in=1,nni_p
        do im=imu_bext,nmuall
           fafter(in,im)=fbefore(in,nmuall_p)
        enddo
     enddo

  elseif (nni.ne.nni_p.and.(nmuall.eq.nmuall_p.and.ngrids.eq.ngrids_p.and. &
       abs(rinf-rinf_p).lt.rinftol)) then

     !        interpolate in nu variable

     if (nni_p.gt.maxmu) then
        print *,"dointerp: error! array vmuq is too short"
        stop 'dointerp'
     endif

     do i=1,nni_p
        vmuq(i)=vni_p(i)
     enddo

     if (ic.eq.1) then
        iord=iord_nu_orb
        kbeg=1
        kend=iord
        iord2=iord/2
        !            stop "dointerp: error check array dimension"
        call dointerp_nu (nmuall_p,fbefore,fafter)
     elseif (ic.eq.2) then
        iord=iord_nu_coul
        kbeg=1
        kend=iord
        iord2=iord/2
        call dointerp_nu (nmuall_p,fbefore,fafter)
     elseif (ic.eq.3) then
        iord=iord_nu_exch
        kbeg=1
        kend=iord
        iord2=iord/2
        call dointerp_nu (nmuall_p,fbefore,fafter)
     endif
  else
     if (nni.ne.nni_p.and.abs(rinf-rinf_p).gt.rinftol) then
        write(*,'(/,1x,"dointerp: cannot perform interpolation since both the number of nu grid points " &
             & /11x "and the value of practical infinity change")')
        stop "dointerp"
     endif
     write(*,'(/1x,"dointerp: no interpolation performed; check input data")')
     stop "dointerp"
  endif
end subroutine dointerp
