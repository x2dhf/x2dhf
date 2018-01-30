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
! ### vexch ###

!     Calculates the value of exchange potential from the multipole
!     expansion.

subroutine vexch(i,j,mt,ipc,pe)
  use params
  use discret
  use scf
  use commons8

  implicit none
  integer :: i,ipc,j,mt,n
  real (PREC) :: costh,costh2,costh4,costh6,costh8,rr,sini,sini2,sini3,xrn
  real (PREC), dimension(maxmpole) :: dome,pe

  pe(1)=0.0_PREC
  
  if (mt.gt.4) return
  
  rr=sqrt(vxisq(j)+vetasq(i)-1.0_PREC)		
  costh=0.0_PREC
  if (abs(rr).gt.precis) then
     costh=veta(i)*vxi(j)/rr	 
  endif
  
  rr=rr*r2
  
  if (mt.eq.0) then
     !        m=0
     
     dome(1)=costh
     dome(2)=(3.0_PREC*costh*costh-1.0_PREC)*0.50_PREC
     do n=2,mpole-1
        dome(n+1)=(dble(2*n+1)*costh*dome(n)-dble(n)*dome(n-1))/dble(n+1)
     enddo
     
!      for m not equal zero the associated Legendre polynomials 
!      are multiplied by sqrt((k-|q|)!/(k+|q|)!) factor 

  elseif (mt.eq.1) then
     !      m=1

     costh2=costh*costh
     sini2=abs(1.0_PREC-costh2)
     sini =sqrt(sini2)
     
     dome(1)=-sini/sqrt(2.0_PREC)
     dome(2)=-3.0_PREC*sini*costh/sqrt(6.0_PREC)
     
     if (mpole.lt.3) goto 100
     dome(3)=(3.0_PREC/2.0_PREC)*sini*(1.0_PREC-5.0_PREC*costh2)/sqrt(12.0_PREC)
     
     if (mpole.lt.4) goto 100
     dome(4)=(5.0_PREC/2.0_PREC)*sini*costh*(3.0_PREC-7.0_PREC*costh2)/sqrt(20.0_PREC)
     
     if (mpole.lt.5) goto 100
     costh4=costh2*costh2
     dome(5)=(15.0_PREC/8.0_PREC)*sini*(-1.0_PREC+14.0_PREC*costh2-21.0_PREC*costh4)/sqrt(30.0_PREC)

     if (mpole.lt.6) goto 100
     dome(6)=(21.0_PREC/8.0_PREC)*sini*costh*(-5.0_PREC+30.0_PREC*costh2-33.0_PREC*costh4)/sqrt(42.0_PREC)

     if (mpole.lt.7) goto 100
     costh6=costh4*costh2
     dome(7)=(7.0_PREC/16.0_PREC)*sini*(5.0_PREC-135.0_PREC*costh2+495.0_PREC*costh4-429.0_PREC*costh6)/sqrt(56.0_PREC)

     if (mpole.lt.8) goto 100
     dome(8)=(9.0_PREC/16.0_PREC)*sini*costh*(35.0_PREC-385.0_PREC*costh2+1001.0_PREC*costh4-715.0_PREC*costh6)/sqrt(72.0_PREC)
     
     !      minus sign for odd values of mt
     !       pe=-pe
     
  elseif (mt.eq.2) then
     !      m=2
     costh2=costh*costh
     sini2=abs(1.0_PREC-costh2)
     
     dome(1)=0.0_PREC
     dome(2)=3.0_PREC*sini2/sqrt(24.0_PREC)
     
     if (mpole.lt.3) goto 100
     dome(3)=15.0_PREC*costh*sini2/sqrt(120.0_PREC)
     
     if (mpole.lt.4) goto 100
     costh4=costh2*costh2
     dome(4)=(15.0_PREC/2.0_PREC)*(-1.0_PREC+8.0_PREC*costh2-7.0_PREC*costh4)/sqrt(360.0_PREC)

     if (mpole.lt.5) goto 100
     dome(5)=(105.0_PREC/2.0_PREC)*costh*(-1.0_PREC+4.0_PREC*costh2-3.0_PREC*costh4)/sqrt(840.0_PREC)
     
     if (mpole.lt.6) goto 100
     costh6=costh4*costh2
     dome(6)=(105.0_PREC/8.0_PREC)*(1.0_PREC-19.0_PREC*costh2+51.0_PREC*costh4-33.0_PREC*costh6)/sqrt(1680.0_PREC)
     
     if (mpole.lt.7) goto 100
     dome(7)=(63.0_PREC/8.0_PREC)*costh*(15.0_PREC-125.0_PREC*costh2+253.0_PREC*costh4-143.0_PREC*costh6)/sqrt(3024.0_PREC)
     
     if (mpole.lt.8) goto 100
     costh8=costh6*costh2
     dome(8)=(315.0_PREC/8.0_PREC)*(-1.0_PREC+34.0_PREC*costh2-176.0_PREC*costh4+&
          286.0_PREC*costh6-143.0_PREC*costh8)/sqrt(5040.0_PREC)
     
  elseif (mt.eq.3) then
     !      m=3
     costh2=costh*costh
     sini2=abs(1.0_PREC-costh2)
     ! FIXME sini =(sini2) or sini =sqrt(sini2)
     sini =(sini2)
     sini3=sini2*sini
     
     dome(1)=0.0_PREC
     dome(2)=0.0_PREC
     
     if (mpole.lt.3) goto 100
     dome(3)=-15.0_PREC*sini3/sqrt(720.0_PREC)
     
     if (mpole.lt.4) goto 100
     dome(4)=-105.0_PREC*sini3*costh/sqrt(5040.0_PREC)
     
     if (mpole.lt.5) goto 100
     dome(5)=-(105.0_PREC/2.0_PREC)*sini3*(-1.0_PREC+9.0_PREC*costh2)/sqrt(20160.0_PREC)
     
     if (mpole.lt.6) goto 100
     dome(6)=-(315.0_PREC/2.0_PREC)*sini3*costh*(-3.0_PREC+11.0_PREC*costh2)/sqrt(60480.0_PREC)
     
     if (mpole.lt.7) goto 100
     costh4=costh2*costh2
     dome(7)=-(315.0_PREC/8.0_PREC)*sini3*(3.0_PREC-66.0_PREC*costh2+143.0_PREC*costh4)/sqrt(151200.0_PREC)
     
     if (mpole.lt.8) goto 100
     dome(8)=-(3465.0_PREC/8.0_PREC)*sini3*costh*(3.0_PREC-26.0_PREC*costh2+39.0_PREC*costh4)/sqrt(332640.0_PREC)
     
     !      minus sign for odd values of mt
     !       pe=-pe
     
  elseif (mt.eq.4) then
     !      m=4
     
     dome(1)=0.0_PREC
     dome(2)=0.0_PREC
     dome(3)=0.0_PREC
     
     if (mpole.lt.4) goto 100
     costh2=costh*costh
     costh4=costh2*costh2
     dome(4)=105.0_PREC*(1.0_PREC-2.0_PREC*costh2+costh4)/sqrt(40320.0_PREC)
     
     if (mpole.lt.5) goto 100
     dome(5)=945.0_PREC*costh*(1.0_PREC-2.0_PREC*costh2+costh4)/sqrt(362880.0_PREC)
     
     if (mpole.lt.6) goto 100
     costh6=costh4*costh2
     dome(6)=(945.0_PREC/2.0_PREC)*(-1.0_PREC+13*costh2-23.0_PREC*costh4+11.0_PREC*costh6)/sqrt(1814400.0_PREC)
     
     if (mpole.lt.7) goto 100
     dome(7)=(3465.0_PREC/2.0_PREC)*costh*(-3.0_PREC+19*costh2-29.0_PREC*costh4+13.0_PREC*costh6)/sqrt(6652800.0_PREC)
     
     if (mpole.lt.8) goto 100
     costh8=costh6*costh2
     dome(8)=(10395.0_PREC/8.0_PREC)*(1.0_PREC-28.0_PREC*costh2+118.0_PREC*costh4 &
          -156.0_PREC*costh6+65.0_PREC*costh8)/sqrt(19958400.0_PREC)
  endif
  
100 continue
  
  xrn=rr*rr
  pe(1)=        dome(1)*excdi(ipc)/xrn
  xrn=xrn*rr
  pe(2)=pe(1) + dome(2)*excqu(ipc)/xrn
  
  if (mpole.ge.3) then
     xrn=xrn*rr
     pe(3)=pe(2) + dome(3)*excoc(ipc)/xrn
  endif
  
  if (mpole.ge.4) then
     xrn=xrn*rr
     pe(4)=pe(3) + dome(4)*exche(ipc)/xrn
  endif
  
  if (mpole.ge.5) then
     xrn=xrn*rr
     pe(5)=pe(4) + dome(5)*exc5(ipc)/xrn
  endif
  
  if (mpole.ge.6) then
     xrn=xrn*rr
     pe(6)=pe(5) + dome(6)*exc6(ipc)/xrn
  endif
  
  if (mpole.ge.7) then
     xrn=xrn*rr
     pe(7)=pe(6) + dome(7)*exc7(ipc)/xrn
  endif
  
  if (mpole.ge.8) then
     xrn=xrn*rr
     pe(8)=pe(7) + dome(8)*exc8(ipc)/xrn
  endif
end subroutine vexch
