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
! ### Eab1HF ###

!     Evaluates off-diagonal Lagrange multipliers in cases
!     occ(iorb1)=occ(iorb2)

!      ipc12=iorb1+(iorb2-1)*norb
!      ipc21=iorb2+(iorb1-1)*norb
!      e(iorb2,iorb1) = e(2,1) = engo(ipc21)=<1|T_k +V_n+V_C-V_x|2>
!      e(iorb1,iorb2) = e(1,2) = engo(ipc12)=<2|T_k +V_n+V_C-V_x|1>

subroutine Eab1HF(iorb1,iorb2,psi,pot,excp,e,f0,wgt1,wgt2,wk0,wk1,wk2,wk3)
  use params
  use discret
  use scf
  use commons8

  implicit none
  integer :: i,i1beg,i1beg1,i1beg2,i2beg,i3beg,ihc,iorb,iorb1,iorb2,ipc12,isym,kex,ngex,&
       ngorb,ngorb1,ngorb2,ngpot,ngpot1,ngpot2,ngrid,ngrid2,norb2,nmut

  real (PREC) :: coo,w,woneel,wtwoel
  real (PREC), dimension(*) :: psi,pot,excp,e,f0,wgt1,wgt2,wk0,wk1,wk2,wk3
  real (PREC), external :: dot

  ipc12=iorb1+(iorb2-1)*norb
  engo(ipc12)=zero
  
  if (mm(iorb1).ne.mm(iorb2)) return
  if (ige(iorb1).ne.ige(iorb2)) return
  
  i1beg1=i1b(iorb1)
  ngorb1=i1si(iorb1)
  ngpot1=i2si(iorb1)
  
  i1beg2=i1b(iorb2)
  ngorb2=i1si(iorb2)
  ngpot2=i2si(iorb2)
  
  norb2=norb*norb
  
  nmut=i1mu(iorb1)
  isym=isymOrb(iorb1)
  
  !    calculate derivatives over mu and ni 

  call putin (nni,nmut,isym,psi(i1beg1),wk3)
  call diffnu (nmut,wk3,wk0,wk1,wk2)
  call putout (nni,nmut,wk1,wk0)
  
  call diffmu (nmut,wk3,wk2)
  call putout (nni,nmut,wk0,wk2)
  
  !     add contributions from derivatives over mu and ni
  ngorb=min(ngorb1,ngorb2)
  call add (ngorb,wk0,wk1)
  
  if (mm(iorb1).eq.0) then
     call copy (ngorb,f0,ione,wk0,ione)
  else
     
     !        nuclear energy for non-sigma orbitals contains contribution
     !        from e term (in totalHF this term is correctly added to the
     !        kinetic energy contribution); e enters the expression with
     !        minus sign which is already incorporated in e
     
     w=dble(mm(iorb1)*mm(iorb1))
     
     call copy (ngorb,f0,ione,wk0,ione)
     call axpy (ngorb,w,e,ione,wk0,ione)
  endif
  
  call proda (ngorb,psi(i1beg1),wk0,wk1)
  call prod (ngorb,psi(i1beg2),wk1)
  
  woneel=dot(ngorb,wgt1,ione,wk1,ione)
  
  !     add contributions from Coulomb and exchange potentials
  
  do iorb=1,norb
     i1beg=i1b(iorb)
     ngorb=i1si(iorb)
     ngpot=i2si(iorb)

     kex=iorb2+norb*(iorb-1)
     if (iorb2.le.iorb) then
        ihc=iorb2+iorb*(iorb-1)/2
     else
        ihc=iorb+iorb2*(iorb2-1)/2
     endif
     
     i2beg=i2b(iorb)
     i3beg=i3b(ihc)
     ngex =i3si(ihc)
     
     ngrid=min(ngorb2,ngpot)
     do i=1,mxsize
        wk0(i)=0.0_PREC
     enddo
     
     coo=occ(iorb)
     if (iorb2.eq.iorb) coo=coo-one
     
     call copy (ngrid,pot(i2beg),ione,wk0,ione)
     call prod  (ngrid,psi(i1beg2),wk0)
     
     if (coo.ne.one) then
        call scal (ngrid,coo,wk0,ione)
     endif
     
     ngrid2=min(ngorb2,ngex,ngorb)
     if (iorb.ne.iorb2)  then
        coo=gec(kex)
        call prodas (ngrid2,-coo,psi(i1beg),excp(i3beg),wk0)
        if (ilc(ihc).gt.1) then
           coo=gec(kex+norb*norb)
           call prodas (ngrid2,-coo,psi(i1beg),excp(i3beg+i3si(ihc)),wk0)
        endif
     else
        if ((mm(iorb1).gt.0).and.(ilc(ihc).gt.0)) then
           coo=gec(kex)
           call prodas (ngrid2,-coo,psi(i1beg),excp(i3beg),wk0)
        endif
     endif
     
     if (iorb.eq.1) then
        call copy (ngrid,wk0,ione,wk1,ione)
     else
        call add   (ngrid,wk0,wk1)
     endif
  enddo
  
  !     to complete the integrand wk1 has to be multiplied by psi(i1beg1)
  call prod (ngrid,psi(i1beg1),wk1)
  wtwoel=dot(ngrid,wgt2,ione,wk1,ione)
  
  engo(ipc12)=woneel+wtwoel
  
end subroutine Eab1HF




