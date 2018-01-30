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
! ### coulMom ###
!
!     Calculates multipole moments up to order k_max=7.

subroutine coulMom (psi,f4,wgt2,d1,d2,d3,d4,d5,d6,d7,d8,wk1,wk2)
  use params
  use discret
  use scf
  use commons8

  implicit none
  integer :: i,ibeg,iorb,mu,n,ni,ngrid
  real (PREC) :: costh,rr,xr,xw
  real (PREC), dimension(*) ::  psi,f4,wgt2,wk1,wk2
  real (PREC), dimension(10) ::  dome
  real (PREC), dimension(nni,mxnmu) :: d1(nni,mxnmu),d2(nni,mxnmu),d3(nni,mxnmu),d4(nni,mxnmu),&
       d5(nni,mxnmu),d6(nni,mxnmu),d7(nni,mxnmu),d8(nni,mxnmu)

  real (PREC), external :: dot


  do iorb=1,norb
     if (itouch(iorb).eq.0) goto 500	    	
     ibeg = i1b (iorb)
     ngrid= i1si(iorb)	  
     
     do mu=1,mxnmu
        do ni=1,nni
           rr=sqrt(vxisq(mu)+vetasq(ni)-1.0_PREC)
           xr=r2*rr
           
           !              xr=r2*rr  r2=R/2  costh == xi*eta/rr
           !              r==xr=(R/2)rr see eq.13 (CPC 98 (1996) 346)
           
           if (abs(rr).lt.precis) then
              costh=0.0_PREC
           else
              costh=vxi(mu)*veta(ni)/rr
           endif
           !              domei=P_k (Legendre polynomial of order k)
           
           dome(1)=costh
           dome(2)=(3.0_PREC*costh*costh-1.0_PREC)*0.50_PREC
           do n=2,(mpole-1)
              dome(n+1)=(dble(2*n+1)*costh*dome(n)-dble(n)*dome(n-1))/dble(n+1)
           enddo
           
           xw=xr
           d1(ni,mu)=dome(1)*xw
           xw=xw*xr
           d2(ni,mu)=dome(2)*xw
           
           if (mpole.ge.3) then
              xw=xw*xr
              d3(ni,mu)=dome(3)*xw
           endif
           
           if (mpole.ge.4) then
              xw=xw*xr
              d4(ni,mu)=dome(4)*xw
           endif
           
           if (mpole.ge.5) then
              xw=xw*xr
              d5(ni,mu)=dome(5)*xw
           endif
           if (mpole.ge.6) then
              xw=xw*xr
              d6(ni,mu)=dome(6)*xw
           endif
           
           if (mpole.ge.7) then
              xw=xw*xr
              d7(ni,mu)=dome(7)*xw
           endif
           
           if (mpole.ge.8) then
              xw=xw*xr
              d8(ni,mu)=dome(8)*xw
           endif
        enddo
     enddo
     
     if (iprint(160).ne.0) then
        write(*,*) 'momen0: d1 '
        call pmtx(nni,mxnmu,d1,ione,ione,incrni,incrmu)
        write(*,*) 'momen0: d3 '
        call pmtx(nni,mxnmu,d3,ione,ione,incrni,incrmu)
        write(*,*) 'momen0: d5 '
        call pmtx(nni,mxnmu,d5,ione,ione,incrni,incrmu)
     endif
     
     !        Now d1, d2, and d3 are equal to P1*r, P2*r**2 and P3*r**2,
     !        respectively (r==xr), where P1, P2, ...  are Legendre polynomials of a
     !        specified order.
     
     !        Prepare integrands and calculate moments. Multiplication by f4 is 
     !        necessary because the factor 2/(r*cosh(mu)) is incorporated in wgt2.
     
     call prod2 (ngrid,psi(ibeg),psi(ibeg),wk1)
     call prod  (ngrid,f4,wk1)
     
     if (ihomon.ne.0) then	
        cmulti (iorb	  ) = 0.0_PREC
     else
        call prod2 (ngrid,d1,wk1,wk2)
        cmulti (iorb	  ) =dot (ngrid,wgt2,ione,wk2,ione)
     endif
     
     call prod2 (ngrid,d2,wk1,wk2)
     cmulti (iorb+norb) =dot (ngrid,wgt2,ione,wk2,ione)
     
     if (mpole.ge.3) then
        if (ihomon.ne.0) then	
           cmulti (iorb+2*norb) = 0.0_PREC
        else
           call prod2 (ngrid,d3,wk1,wk2)
           cmulti (iorb+2*norb) =dot (ngrid,wgt2,ione,wk2,ione)
        endif
     endif
     
     if (mpole.ge.4) then
        call prod2 (ngrid,d4,wk1,wk2)
        cmulti (iorb+3*norb) =dot (ngrid,wgt2,ione,wk2,ione)
     endif
     
     if (mpole.ge.5) then
        if (ihomon.ne.0) then	
           cmulti (iorb+4*norb) = 0.0_PREC
        else
           call prod2 (ngrid,d5,wk1,wk2)
           cmulti (iorb+4*norb) =dot (ngrid,wgt2,ione,wk2,ione)
        endif
     endif
     
     if (mpole.ge.6) then
        call prod2 (ngrid,d6,wk1,wk2)
        cmulti (iorb+5*norb) =dot (ngrid,wgt2,ione,wk2,ione)
     endif
     
     if (mpole.ge.7) then
        if (ihomon.ne.0) then	
           cmulti (iorb+6*norb) = 0.0_PREC
        else
           call prod2 (ngrid,d7,wk1,wk2)
           cmulti (iorb+6*norb) =dot (ngrid,wgt2,ione,wk2,ione)
        endif
     endif
     
     if (mpole.ge.8) then
        call prod2 (ngrid,d8,wk1,wk2)
        cmulti (iorb+7*norb) =dot (ngrid,wgt2,ione,wk2,ione)
     endif
     
     if (iprint(161).ne.0) then
        write(*,*)
        write(*,*) 'momen0 - multipole moments for orbital ',iorb
        write(*,1111) (cmulti (iorb+i*norb),i=0,mpole-1)
1111    format(4e25.14)
     endif
500  continue
  enddo

end subroutine coulMom
