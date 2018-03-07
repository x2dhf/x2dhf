! ***************************************************************************
! *                                                                         *
! *   Copyright (C) 2005-2010 Jacek Kobus <jkob@fizyka.umk.pl>              *
! *                                                                         *
! *   This program is free software; you can redistribute it and/or modify  *
! *   it under the terms of the GNU General Public License version 2 as     *
! *   published by the Free Software Foundation.                            *
! *                                                                         *
! ***************************************************************************
! ### nuclder ####

!     This subroutine calculates the derivative of the orbital density
!     over z (along the internuclear distance) at (0,0,-R/2).  At this
!     point the derivative cannot be calculated in (\nu,\mu) variables
!     since the denominators of the derivative expresion are zero though
!     the derivative over z is finite at the point.  The derivative has
!     to be calculated for non -1 values of \nu and then extrapolated to
!     \nu=-1.

subroutine nuclder(psi,e,f0,wgt1,wgt2,dermu,dernu,dznu,wk2,wk3,wk4,wk5)
  use params
  use discret
  use commons8

  use blas_m
  use diffmu_m
  use diffnu_m
  use factor_m
  
  implicit none
  integer :: i,i1beg,ib,imu,inu,iorb,iprt,isym,l1,m1,n1,ngorb,ngpot,nmut
  real (PREC) :: ez1,fn,fn1,shn1,tderz
  real (PREC), dimension(*) :: e,f0,wgt1,wgt2,dznu,wk2,wk3,wk4,wk5
  real (PREC), dimension(nni,*) :: psi,dermu,dernu
  real (PREC), dimension(6) :: dzmu
  real (PREC), dimension(maxorb) :: tderzorb1,tderzorb2,tderzorb3,tderzorb4
  iprt=0

  !     prepare differentialtion arrays used by diffnu and difmu

  if (nni.gt.5000) then
     write(*,*) 'nuclder: nni > 5000 '
     stop 'nuclder'
  endif

  call prepdiff1

  !     orbital density functions are of even symmetry

  write(*,*)
  write(*,*) 'derivative and derivative/(N*N) of density over z at (0,0,-R/2) for orbital '
  write(*,*)

  !     calculate the derivative of density as d(rho)/dz = 2 psi d(psi)/dz

  isym=1
  tderz=.0_PREC
  nmut=mxnmu

  do iorb=1,norb

     !        for hydrogen-like systems the accuracy of the derivatives can
     !        be checked since the analytical formulae are known. For example,
     !        the 1s orbital density is N^2*exp(-2Zr) and the derivative is
     !        -2Z up to the normalization constant. This procedure calculates
     !        the derivative and prints out both its value and the value
     !        divided by N squared.

     n1=mgx(1,iorb)
     l1=mgx(2,iorb)
     m1=mgx(3,iorb)
     ez1=eza1(iorb)

     ! normalization factor for Laguere polynomials

     fn1=(2.0_PREC*ez1/dble(n1))**(3.0_PREC/2.0_PREC+dble(l1)) &
          *sqrt(factor(n1+l1)/(2.0_PREC*dble(n1)*factor(n1-l1-1)))/factor(2*l1+1)

     ! normalization factor for spherical harmonics

     shn1=(-1.0_PREC)**dble(m1)/sqrt(4.0_PREC*pii)*sqrt((2*l1+1)*factor(l1-m1)/factor(l1+m1))
     if (m1.eq.0) shn1=1.0_PREC/sqrt(4.0_PREC*pii)*sqrt((2*l1+1)*factor(l1-m1)/factor(l1+m1))

     fn=fn1*shn1

     i1beg=i1b(iorb)
     nmut=i1mu(iorb)
     ib=(iorb-1)*nmut+1

     call nuclder1(iorb,psi(1,ib),dznu,dzmu,dernu,dermu,wk2,wk3,wk4,wk5)

     tderz=tderz+dznu(nni)+dzmu(1)
     tderzorb1(iorb)=dznu(nni)
     tderzorb2(iorb)=dzmu(1)
     tderzorb3(iorb)=dznu(nni)+dzmu(1)
     tderzorb4(iorb)=dznu(nni)/(fn*fn)
  enddo

  write(*,*)
  write(*,*) 'derivative of the total density',tderz
  write(*,*) 'and its components '
  do iorb=1,norb
     write(*,1000) iorn(iorb),bond(iorb),tderzorb1(iorb),tderzorb2(iorb),tderzorb3(iorb),tderzorb4(iorb)
  enddo

  write(*,*)
  write(*,*)

  !     calculate the dervivative of density directly

  tderz=0.0_PREC
  do iorb=1,norb

     !        for hydrogen-like systems the accuracy of the derivatives can
     !        be checked since the analytical formulae are known For example,
     !        the 1s orbital density is N^2*exp(-2Zr) and the derivative is
     !        -2Z up to the normalization constant. This procedure calculates
     !        the derivative and prints outs both its value and the value
     !        divided by N squared.

     n1=mgx(1,iorb)
     l1=mgx(2,iorb)
     m1=mgx(3,iorb)
     ez1=eza1(iorb)

     !        normalization factor for Laguere polynomials

     fn1=(2.0_PREC*ez1/dble(n1))**(3.0_PREC/2.0_PREC+dble(l1)) &
          *sqrt(factor(n1+l1)/(2.0_PREC*dble(n1)*factor(n1-l1-1)))/factor(2*l1+1)

     !        normalization factor for spherical harmonics

     shn1=(-1.0_PREC)**dble(m1)/sqrt(4.0_PREC*pii)*sqrt((2*l1+1)*factor(l1-m1)/factor(l1+m1))
     if (m1.eq.0) shn1=1.0_PREC/sqrt(4.0_PREC*pii)*sqrt((2*l1+1)*factor(l1-m1)/factor(l1+m1))

     fn=fn1*shn1

     i1beg=i1b(iorb)
     nmut=i1mu(iorb)
     ngorb=i1si(iorb)
     ngpot=i2si(iorb)

     call tden(iorb,ngorb,psi,wk2)
     call copy(ngorb,wk2,ione,psi(ione,ione),ione)

     call putin (nni,nmut,isym,psi(ione,ione),wk3)
     call diffnu (nmut,wk3,dermu,dernu,wk2)
     call putout (nni,nmut,dernu,dermu)

     call diffmu (nmut,wk3,wk2)
     call putout (nni,nmut,dermu,wk2)

     !  derivatives at the right nuclei (0,0,R/2)

     !       dernu(1,1)= exeven(1)*dernu(2,1)+exeven(2)*dernu(3,1)
     !     &          +exeven(3)*dernu(4,1)
     !     &          +exeven(4)*dernu(5,1)+exeven(5)*dernu(6,1)

     !        dermu(1,1)= exeven(1)*dermu(1,2)+exeven(2)*dermu(1,3)
     !     &          +exeven(3)*dermu(1,4)
     !     &          +exeven(4)*dermu(1,5)+exeven(5)*dermu(1,6)

     !        write(*,*) '(dernu(i,1),i=1,6) after extrapolation'
     !        write(*,*) (dernu(i,1),i=1,6)

     !        write(*,*) 'difmu'
     !        write(*,*) (dermu(1,i),i=1,6)

     dernu(nni,1)= exeven(1)*dernu(nni-1,1)+exeven(2)*dernu(nni-2,1)+exeven(3)*dernu(nni-3,1) &
          +exeven(4)*dernu(nni-4,1)+exeven(5)*dernu(nni-5,1)

     dermu(nni,1)= exeven(1)*dermu(nni,2)+exeven(2)*dermu(nni,3)+exeven(3)*dermu(nni,4) &
          +exeven(4)*dermu(nni,5)+exeven(5)*dermu(nni,6)

     !  vxi1(i)=sqrt(xi^2-1)
     !  veta1(i)=sqrt(1-eta^2)

     !  Having derivative over nu calculate the derivative over z near
     !  (0,0,-R/2)

     do inu=nni-5,nni
        dznu(inu)=(-1.0_PREC)*2.0_PREC/(r*veta1(inu))*dernu(inu,1)
     enddo

     dzmu(1)=.0_PREC
     do imu=2,6
        dzmu(imu)=(-1.0_PREC)*2.0_PREC/(r*vxi1(imu))*dermu(nni,imu)
     enddo

     !       calculate the derivative over z at (0,0,-R/2+) by extrapolation

     dznu(nni)= exeven(1)*dznu(nni-1)+exeven(2)*dznu(nni-2)+exeven(3)*dznu(nni-3) &
          +exeven(4)*dznu(nni-4)+exeven(5)*dznu(nni-5)

     !   calculate the derivative over z at (0,0,-R/2-) by extrapolation

     dzmu(1)= exeven(1)*dzmu(2)+exeven(2)*dzmu(3)+exeven(3)*dzmu(4) &
          +exeven(4)*dzmu(5)+exeven(5)*dzmu(6)

     iprt=1
     if (iprt.eq.0) then
        write(*,*) 'dznu(nni) ',(dznu(i),i=nni-5,nni)
        write(*,*) 'dzmu(1)   ',(dzmu(i),i=1,6)

        write(*,*) iorn(iorb),'  ',bond(iorb),dznu(nni),dznu(nni)/(fn*fn),dznu(nni)*(fn*fn)
     endif

     tderz=tderz+dznu(nni)+dzmu(1)
     tderzorb1(iorb)=dznu(nni)
     tderzorb2(iorb)=dzmu(1)
     tderzorb3(iorb)=dznu(nni)+dzmu(1)
     tderzorb4(iorb)=dznu(nni)/(fn*fn)
  enddo

  write(*,*)
  write(*,*) 'derivative of the total density',tderz
  write(*,*) 'and its components '
  do iorb=1,norb
     write(*,1000) iorn(iorb),bond(iorb),tderzorb1(iorb),tderzorb2(iorb),tderzorb3(iorb),tderzorb4(iorb)
  enddo

01000 format(1x,i2,a6,2x,4e20.10)
  write(*,*)
end subroutine nuclder



