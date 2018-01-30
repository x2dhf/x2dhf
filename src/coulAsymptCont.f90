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
! ### coulAsymptCont ###
!     evaluates multipole moment contributions to the boundary values of
!     Coulomb potential for a given orbital at imu=mxnmu, inu=(nn1-1)/2 (Pi/2)

subroutine coulAsymptCont(iorb,pot)
  use params
  use discret
  use scf
  use commons8

  implicit none
  integer :: i,iorb,itt,j,kk,kxk,m,n
  real (PREC) :: costh,pe,rr,xr,xrr 
  real (PREC), dimension(10) :: dome,pottmp
  real (PREC), dimension(*) :: pot

!     potentials are calculated for ni=Pi/2

  j=mxnmu
  itt=(j-1)*nni
  i=(nni-1)/2
  kk=i+itt
  rr=sqrt(vxisq(j)+vetasq(i)-1.0_PREC)
  costh=veta(i)*vxi(j)/rr
  xr=1.0_PREC/(rr*r2)
  
  dome(1)=costh
  dome(2)=(3.0_PREC*costh*costh-1.0_PREC)*0.50_PREC
  do n=2,mpole-1
     dome(n+1)=(dble(2*n+1)*costh*dome(n)-dble(n)*dome(n-1))/dble(n+1)
  enddo
  
  pe=0.0_PREC
  xrr=xr
  do m=1,mpole
     xrr=xrr*xr
     kxk=iorb+(m-1)*norb
     pe=pe+cmulti(kxk)*dome(m)*xrr
     pottmp(m)=r2*vxi(j)*(pe+xr)	  
  enddo
  
  write(*,1000) iorn(iorb),bond(iorb),gut(iorb),pottmp(1),(pottmp(m)-pottmp(m-1),m=2,mpole),pottmp(mpole)
1000 format(/i4,1x,a8,a1,3x,/4e13.5/5e13.5)
end subroutine coulAsymptCont






