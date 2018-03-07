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
! ### propet3 ####

!     Calculates dipole, quadrupole, octopole and hexadecapole moments
!     relative to the geometrical centre and the centre of mass.

subroutine propet3 (cw_orb,cw_coul,cw_exch,f4,wgt2,wk0,wk1,wk2,wk3,wk4,wk5)
  use params
  use discret
  use scf
  use commons8

  use blas_m
  use plegendg_m
  
  implicit none

  integer :: i1beg,igp,imu,in,inioff,iorb,izz1,izz2,k,ngorb
  real (PREC) :: atw1,atw2,cm2zz,costh,qq,qtot,qktot,qkz1,qkz2,rr,&
       sum1,sum2,sum3,xxplusyy,z,zcm

  real (PREC), dimension(*) :: cw_orb,cw_coul,cw_exch,f4,wgt2,wk0,wk1,wk2,wk3,wk4,wk5

  !     1 bohr electron = 2.541765 Debye -- Gaussian94 User's Reference
  !     data au2Debye /2.5417650_PREC/
  !     data bohr2ang /0.5291772490_PREC/


  do k=1,maxmpole
     total(k)=0.0_PREC
  enddo

  zcm=0.0_PREC
  qtot=0.0_PREC


  izz1=nint(z1)
  izz2=nint(z2)
  atw1=atweight(izz1)
  atw2=atweight(izz2)
  zcm=(-atw1+atw2)/(atw1+atw2)*r/2.0_PREC

  ! multipole moments relative to the centre of mass


  izz1=nint(z1)
  izz2=nint(z2)
  atw1=atweight(izz1)
  atw2=atweight(izz2)
  zcm=(-atw1+atw2)/(atw1+atw2)*r/2.0_PREC

  ! multipole moments relative to the centre of mass

  do k=1,mpole

     !  loop over grid points

     do imu=1,mxnmu
        inioff=(imu-1)*nni
        do in=1,nni
           igp=inioff+in

           !              rr=(r/2.0_PREC)*sqrt(vxisq(imu)+vetasq(in)-1.0_PREC)
           !              x=(r/2.0_PREC)*vxi1(imu)*veta1(in)

           z=(r/2.0_PREC)*vxi(imu)*veta(in)-zcm
           xxplusyy=((r/2.0_PREC)*vxi1(imu)*veta1(in))**2
           rr=sqrt(xxplusyy+z*z)

           if (abs(rr).lt.precis) then
              costh=0.0_PREC
           else
              costh=z/rr
           endif
           wk0(igp)=rr**k*plegendg(k,izero,costh)
        enddo
     enddo

     qtot=0.0_PREC
     cm2zz=0.0_PREC
     do iorb=1,norb
        i1beg=i1b(iorb)
        ngorb=i1si(iorb)
        qq=-occ(iorb)
        qtot=qtot+qq
        call prod2 (ngorb,cw_orb(i1beg),cw_orb(i1beg),wk3)

        nni2=nni/2
        sum1=0.0_PREC
        sum2=0.0_PREC
        do imu=mxnmu,1,-1
           inioff=(imu-1)*nni
           do in=1,nni2
              igp=inioff+in
              sum1=sum1+wk0(igp)*wk3(igp)*f4(igp)*wgt2(igp)
           enddo

           do in=nni2+1,nni
              igp=inioff+in
              sum2=sum2+wk0(igp)*wk3(igp)*f4(igp)*wgt2(igp)
           enddo
        enddo

        sum3=0.0_PREC
        do imu=1,mxnmu
           inioff=(imu-1)*nni
           do in=1,nni
              igp=inioff+in
              wk5(igp)=wk0(igp)*wk3(igp)*f4(igp)*wgt2(igp)
              sum3=sum3+wk5(igp)
           enddo
        enddo

        if (iprint(214).ne.0) then
           write(iout6,'(i5,i4,1x,a8,a1,e28.16,e10.2)') k, iorn(iorb),bond(iorb),gut(iorb),sum3,sum1+sum2-sum3
        endif

        cm2zz=cm2zz+qq*(sum1+sum2)
     enddo

     qkz1=z1*abs(r/2.0_PREC+zcm)**k*plegendg(k,izero,-1.0_PREC)
     qkz2=z2*abs(r/2.0_PREC-zcm)**k*plegendg(k,izero,1.0_PREC)
     qktot=cm2zz+qkz1+qkz2
     total(k)=qktot
  enddo

end subroutine propet3
