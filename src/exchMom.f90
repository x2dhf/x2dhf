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
! c ### exchMom ###	
! c
! c     Calculates multipole moments to k=8 and m=4 (iorb1<=iorb2)

subroutine exchMom (iorb1,iorb2,psi,f4,wgt2,d1,d2,d3,d4,d5,d6,d7,d8,wk1,wk2)
  use params
  use discret
  use scf
  use commons8

  implicit none
  integer :: i,ibeg1,ibeg2,iorb1,iorb2,idel,ido,ipc,mu,ni,ngrid
  
  real (PREC) ::  xrr,xw
  real (PREC), dimension(*) ::  psi,f4,wgt2,wk1,wk2
  real (PREC), dimension(10) ::  dome
  real (PREC), dimension(nni,mxnmu) :: d1(nni,mxnmu),d2(nni,mxnmu),d3(nni,mxnmu),d4(nni,mxnmu),&
       d5(nni,mxnmu),d6(nni,mxnmu),d7(nni,mxnmu),d8(nni,mxnmu)

  real (PREC), external :: dot

  !     iorb2>=iorb1
  
  if (itouch(iorb1).eq.0.and.itouch(iorb2).eq.0) return
  if (iorb1.eq.iorb2.and.mgx(6,iorb2).eq.0) return
  
  ido=0
1234  ido=ido+1
  
  !     calculate m
  
  idel=abs(mgx(6,iorb2)-mgx(6,iorb1))
  if (iorb1.eq.iorb2) idel=2*mgx(6,iorb2)
  !      ipc=iorb1+iorb2*(iorb2-1)/2
  ipc=i3xk(iorb1,iorb2)
  
  !      if ((iorb1.eq.iorb2).and.(ilc(ipc).lt.1)) return
  
  if (ipc.eq.0) return
  
  !     second values of expansion coefficients for the same pair of
  !     orbitals are stored in the second half of the arrays.
  
  if (ido.eq.2) then
     idel=mgx(6,iorb2)+mgx(6,iorb1)
     ipc=ipc+norb*(norb+1)/2
  endif
  
  !   homonuclear case is treated as a heteronuclear one
  
  ibeg1=i1b (iorb1)
  ibeg2=i1b (iorb2)
  
  do mu=1,mxnmu
     do ni=1,nni
        xrr=r2*sqrt (vxisq(mu)+vetasq(ni)-1.0_PREC)
        call mulex(ni,mu,idel,dome)
        
        xw=xrr
        d1(ni,mu)=dome(1)*xw
        xw=xw*xrr
        d2(ni,mu)=dome(2)*xw
        
        if (mpole.ge.3) then
           xw=xw*xrr
           d3(ni,mu)=dome(3)*xw
        endif
        
        if (mpole.ge.4) then
           xw=xw*xrr
           d4(ni,mu)=dome(4)*xw
        endif
        
        if (mpole.ge.5) then
           xw=xw*xrr
           d5(ni,mu)=dome(5)*xw
        endif
        
        if (mpole.ge.6) then
           xw=xw*xrr
           d6(ni,mu)=dome(6)*xw
        endif
        
        if (mpole.ge.7) then
           xw=xw*xrr
           d7(ni,mu)=dome(7)*xw
        endif
        
        if (mpole.ge.8) then
           xw=xw*xrr
           d8(ni,mu)=dome(8)*xw
        endif
        
        if (iprint(165).ne.0) then
           if ((iorb1.eq.1).and.(iorb2.eq.1).and.(ni.eq.nni/2).and.(mu.eq.mxnmu/2)) then
              write(*,'(10e16.7)') xrr,xw
              write(*,'(10e16.7)') (dome(i),i=1,8)
              print *,'iorb1,iorb2,ni,mu,nnu,mxnmu',iorb1,iorb2,ni,mu,nni,mxnmu
           endif
        endif
     enddo
  enddo
  
  !     Now d1,d2,... are equal to N*P(1,q)*r, N*P(2,q)*r**2, ...,
  !     respectively. P(1,q), P(2,q),...  are associate Legendre polynomials
  !     of a specified order (k) and degree q (q==idel). N is equal to
  !     [(k-|q|)!/(k+|q|)!]^{1/2}.
  !     Note that N*N is present in eq. (19). One N is taken care of when 
  !     evaluating Q_kq (Q_{k\Delta m}^{ab}) and the other when generating 
  !     Pkm (P_{k|\Delta m|).   
  
  !     Prepare integrands and calculate moments. Multiplication by f4 is 
  !     necessary because of the factor 2/(r*cosh(mu)) incorporated in wgt2.
  
  ngrid = min (i1si(iorb1),i1si(iorb2))
  call prod2 (ngrid,psi(ibeg1),psi(ibeg2),wk1)
  call prod  (ngrid,f4,wk1)
  
  call prod2 (ngrid,d1,wk1,wk2)
  excdi(ipc) =dot (ngrid,wgt2,ione,wk2,ione)
  
  call prod2 (ngrid,d2,wk1,wk2)
  excqu(ipc) =dot (ngrid,wgt2,ione,wk2,ione)
  
  if (mpole.ge.3) then
     call prod2 (ngrid,d3,wk1,wk2)
     excoc(ipc) =dot (ngrid,wgt2,ione,wk2,ione)
  endif
  
  if (mpole.ge.4) then
     call prod2 (ngrid,d4,wk1,wk2)
     exche(ipc) =dot (ngrid,wgt2,ione,wk2,ione)
  endif
  
  if (mpole.ge.5) then
     call prod2 (ngrid,d5,wk1,wk2)
     exc5(ipc) =dot (ngrid,wgt2,ione,wk2,ione)
  endif
  
  if (mpole.ge.6) then
     call prod2 (ngrid,d6,wk1,wk2)
     exc6(ipc) =dot (ngrid,wgt2,ione,wk2,ione)
  endif
  
  if (mpole.ge.7) then
     call prod2 (ngrid,d7,wk1,wk2)
     exc7(ipc) =dot (ngrid,wgt2,ione,wk2,ione)
  endif
  
  if (mpole.ge.8) then
     call prod2 (ngrid,d8,wk1,wk2)
     exc8(ipc) =dot (ngrid,wgt2,ione,wk2,ione)
  endif
  
  if (iprint(166).ne.0) then
     write(*,1000) iorn(iorb1),bond(iorb1),gut(iorb1),iorn(iorb2),bond(iorb2),gut(iorb2),idel
1000 format(/,i4,1x,a8,a1,3x,i4,1x,a8,a1,3x,i5)
     
     write(*,1010) excdi(ipc),excqu(ipc),excoc(ipc),exche(ipc),exc5(ipc),exc6(ipc),exc7(ipc),exc8(ipc)
1010 format(4d25.14)
  endif
  
  !     do not replace iorb1+iorb2*(iorb2-1)/2 by ipc!
  if (ilc(iorb1+iorb2*(iorb2-1)/2).eq.2.and.ido.eq.1) go to 1234
  
end subroutine exchMom
