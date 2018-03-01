! ***************************************************************************
! *                                                                                                                     *
! *   Copyright (C) 1996 Leif Laaksonen, Dage Sundholm                                   *
! *   Copyright (C) 1996-2010 Jacek Kobus <jkob@fizyka.umk.pl>                     *
! *                                                                                                                     *     
! *   This program is free software; you can redistribute it and/or modify            *
! *   it under the terms of the GNU General Public License version 2 as              *
! *   published by the Free Software Foundation.                                               *
! *                                                                                                                     *
! ***************************************************************************
! ### mulex ###

!     This routine generates values of the associate Legendre functions
!     P(k,q) multiplied by [(k-|q|)!/(k+|q|)!]^{1/2} for k=1,..,8 and
!     q=mt.  See also routine which employs the same definition of the
!     functions multi.

!     Note that the original routine (mulex) generates values scaled by
!     (-1)**mt

subroutine mulex(i,j,mt,dome)
  use params
  use discret
  use commons8

  implicit none
  integer :: i,j,mt,n
  real (PREC) :: costh,costh2,costh4,costh6,costh8,rr,sini,sini2,sini3
  real (PREC), dimension(10) :: dome

  dome(1)=0.0_PREC

  if (mt.gt.4) return
  
  rr=sqrt(vxisq(j)+vetasq(i)-1.0_PREC)
  costh=0.0_PREC
  if (abs(rr).gt.precis) then
     costh=veta(i)*vxi(j)/rr
  endif
  
  
  if (mt.eq.0) then
     !        m=0
     dome(1)=costh
     dome(2)=(3.0_PREC*costh*costh-1.0_PREC)*0.50_PREC
     do n=2,mpole-1
        dome(n+1)=(dble(2*n+1)*costh*dome(n)-dble(n)*dome(n-1))/dble(n+1)
     enddo
     return
     
  elseif (mt.eq.1) then
     !        m=1
     
     costh2=costh*costh
     sini2=abs(1.0_PREC-costh2)
     sini =sqrt(sini2)
     !        P_1^1
     dome(1)=-sini/sqrt(2.0_PREC)
     !        P_2^1
     dome(2)=-3.0_PREC*sini*costh/sqrt(6.0_PREC)
     
     if (mpole.lt.3) return
     !        P_3^1
     dome(3)=(3.0_PREC/2.0_PREC)*sini*(1.0_PREC-5.0_PREC*costh2)/sqrt(12.0_PREC)
     
     if (mpole.lt.4) return
     !        P_4^1
     dome(4)=(5.0_PREC/2.0_PREC)*sini*costh*(3.0_PREC-7.0_PREC*costh2)/sqrt(20.0_PREC)
     
     if (mpole.lt.5) return
     !        P_5^1
     costh4=costh2*costh2
     dome(5)=(15.0_PREC/8.0_PREC)*sini*(-1.0_PREC+14.0_PREC*costh2-21.0_PREC*costh4)/sqrt(30.0_PREC)
     if (mpole.lt.6) return
     !        P_6^1
     dome(6)=(21.0_PREC/8.0_PREC)*sini*costh*(-5.0_PREC+30.0_PREC*costh2-33.0_PREC*costh4)/sqrt(42.0_PREC)
     
     if (mpole.lt.7) return
     !        P_7^1
     costh6=costh4*costh2
     dome(7)=(7.0_PREC/16.0_PREC)*sini*(5.0_PREC-135.0_PREC*costh2+495.0_PREC*costh4-429.0_PREC*costh6)/sqrt(56.0_PREC)
     
     if (mpole.lt.8) return
     !        P_8^1
     dome(8)=(9.0_PREC/16.0_PREC)*sini*costh*(35.0_PREC-385.0_PREC*costh2+1001.0_PREC*costh4-715.0_PREC*costh6)/sqrt(72.0_PREC)
     return
     
  elseif (mt.eq.2) then
     !        m=2
     costh2=costh*costh
     sini2=abs(1.0_PREC-costh2)
     
     !        P_1^2
     dome(1)=0.0_PREC
     !        P_2^2
     dome(2)=3.0_PREC*sini2/sqrt(24.0_PREC)
     
     if (mpole.lt.3) return
     !        P_3^2
     dome(3)=15.0_PREC*costh*sini2/sqrt(120.0_PREC)
     
     if (mpole.lt.4) return
     !        P_4^2
     costh4=costh2*costh2
     dome(4)=(15.0_PREC/2.0_PREC)*(-1.0_PREC+8.0_PREC*costh2-7.0_PREC*costh4)/sqrt(360.0_PREC)

     if (mpole.lt.5) return
     !        P_5^2
     dome(5)=(105.0_PREC/2.0_PREC)*costh*(-1.0_PREC+4.0_PREC*costh2-3.0_PREC*costh4)/sqrt(840.0_PREC)
     
     if (mpole.lt.6) return
     !        P_6^2
     costh6=costh4*costh2
     dome(6)=(105.0_PREC/8.0_PREC)*(1.0_PREC-19.0_PREC*costh2+51.0_PREC*costh4-33.0_PREC*costh6)/sqrt(1680.0_PREC)
     
     if (mpole.lt.7) return
     !        P_7^2
     dome(7)=(63.0_PREC/8.0_PREC)*costh*(15.0_PREC-125.0_PREC*costh2+253.0_PREC*costh4-143.0_PREC*costh6)/sqrt(3024.0_PREC)
     
     if (mpole.lt.8) return
     !        P_8^2
     costh8=costh6*costh2
     dome(8)=(315.0_PREC/8.0_PREC)*(-1.0_PREC+34.0_PREC*costh2-176.0_PREC*costh4+ &
          286.0_PREC*costh6-143.0_PREC*costh8)/sqrt(5040.0_PREC)
     return
     
  elseif (mt.eq.3) then 
     !        m=3
     
     costh2=costh*costh
     sini2=abs(1.0_PREC-costh2)
     sini =sqrt(sini2)
     sini3=sini2*sini
     
     dome(1)=0.0_PREC
     dome(2)=0.0_PREC
     
     if (mpole.lt.3) return
     !        P_3^3
     dome(3)=-15.0_PREC*sini3/sqrt(720.0_PREC)
     
     if (mpole.lt.4) return
     !        P_4^3
     dome(4)=-105.0_PREC*sini3*costh/sqrt(5040.0_PREC)
     
     if (mpole.lt.5) return
     !        P_5^3
     dome(5)=-(105.0_PREC/2.0_PREC)*sini3*(-1.0_PREC+9.0_PREC*costh2)/sqrt(20160.0_PREC)

     if (mpole.lt.6) return
     !        P_6^3
     dome(6)=-(315.0_PREC/2.0_PREC)*sini3*costh*(-3.0_PREC+11.0_PREC*costh2)/sqrt(60480.0_PREC)
     
     if (mpole.lt.7) return
     !        P_7^3
     costh4=costh2*costh2
     dome(7)=-(315.0_PREC/8.0_PREC)*sini3*(3.0_PREC-66.0_PREC*costh2+143.0_PREC*costh4)/sqrt(151200.0_PREC)
     
     if (mpole.lt.8) return
     dome(8)=-(3465.0_PREC/8.0_PREC)*sini3*costh*(3.0_PREC-26.0_PREC*costh2+39.0_PREC*costh4)/sqrt(332640.0_PREC)
     return
     
  elseif (mt.eq.4) then
     !        m=4
     
     dome(1)=0.0_PREC
     dome(2)=0.0_PREC
     dome(3)=0.0_PREC
     
     if (mpole.lt.4) return
     !        P_4^4
     costh2=costh*costh
     costh4=costh2*costh2
     dome(4)=105.0_PREC*(1.0_PREC-2.0_PREC*costh2+costh4)/sqrt(40320.0_PREC)
     
     if (mpole.lt.5) return
     !        P_5^4
     dome(5)=945.0_PREC*costh*(1.0_PREC-2.0_PREC*costh2+costh4)/sqrt(362880.0_PREC)
     
     if (mpole.lt.6) return
     !        P_6^4
     costh6=costh4*costh2
     dome(6)=(945.0_PREC/2.0_PREC)*(-1.0_PREC+13*costh2-23.0_PREC*costh4+11.0_PREC*costh6)/sqrt(1814400.0_PREC)
     
     if (mpole.lt.7) return
     !        P_7^4
     dome(7)=(3465.0_PREC/2.0_PREC)*costh*(-3.0_PREC+19*costh2-29.0_PREC*costh4+13.0_PREC*costh6)/sqrt(6652800.0_PREC)
     
     if (mpole.lt.8) return
     !        P_8^4
     costh8=costh6*costh2
     dome(8)=(10395.0_PREC/8.0_PREC)*(1.0_PREC-28.0_PREC*costh2+118.0_PREC*costh4- &
          156.0_PREC*costh6+65.0_PREC*costh8)/sqrt(19958400.0_PREC)
     return
  endif
  
end subroutine mulex

