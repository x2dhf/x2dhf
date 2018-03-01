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
! ### fock ###
!
!     Calculates the Fock potential for a given orbital.

subroutine fock (iorb,psi,pot,excp,e,f0,f1,f2,f4,fock1,fock2,wk1,wk2)
  use params
  use discret
  use scf
  use commons8

  implicit none
  integer :: i,ibexp,iborb,iborb1,ibpot,ibpot1,idexp,iorb,iorb1,ipc,kex,ngexp,ngorb,ngorb1,ngpot,ngpot1,norb2,ngrid
  real (PREC) :: coo,coo0,coo1,coo2,w
  real (PREC), dimension(*) :: psi,pot,excp,e,f0,f1,f2,f4,fock1,fock2,wk1,wk2
  real (PREC), external :: checkev

  iborb=i1b(iorb)
  ngorb=i1si(iorb)
  ibpot=i2b(iorb)
  ngpot=i2si(iorb)
  norb2=norb*norb

  call copy (mxsize,f0,ione,fock1,ione)
  call axpy (mxsize,eng(iorb),f1,ione,fock1,ione)

  if (mm(iorb).ne.0) then
     !      e enter the expression with minus sign which is already incorporated in e
     w=dble(mm(iorb)*mm(iorb))
     call axpy (mxsize,w,e,ione,fock1,ione)
  endif

  !   in the asymptotic region E-V should be negative for orbitals to
  !   converge; it may be not the case when external electric field is
  !   too strong or the practical infinity is too small

  if (checkev(nni,mxnmu,fock1).lt.zero) then
     write(6,'("Warning! E-V > 0 for orbital",i3)') iorb
  endif

  if (iprint(95).ne.0) then
     write(6,*) "Checking (E-V) values for orbital",iorb
     call checkasymp (nni,mxnmu,fock1,ione,ione,incrni,incrmu)
  endif

  do i=1,mxsize
     fock2(i)=0.0_PREC
  enddo

  if (iprint(91).ne.0) then
     write(*,*) 'fock: one-electron contribution.'
     call pmtx (nni,nmu(1),fock1,ione,ione,incrni,incrmu)
  endif

  if (nel.eq.1) return

  do i=1,mxsize
     wk1 (i)=0.0_PREC
  enddo

  if (norb.eq.1) then

     !      FIXME what about one non-sigma orbital
     !      contribution from coulomb potential of one sigma orbital
     call copy (ngpot,pot(ibpot),ione,wk1,ione)
  else

     !      add contributions from Coulomb and exchange potentials

     do iorb1=1,norb
        iborb1=i1b(iorb1)
        ngorb1=i1si(iorb1)
        kex=iorb+norb*(iorb1-1)
        ipc=iorb1+norb*(iorb-1)

        if (iorb.le.iorb1) then
           idexp=iorb+iorb1*(iorb1-1)/2
        else
           idexp=iorb1+iorb*(iorb-1)/2
        endif

        ibpot1=i2b(iorb1)
        ngpot1=i2si(iorb1)
        ibexp=i3b(idexp)
        ngexp=i3si(idexp)

        coo=occ(iorb1)

        !         FIXME
        if (iorb.eq.iorb1) coo=coo-1.0_PREC
        coo0=coo

        call axpy (ngpot1,coo,pot(ibpot1),ione,wk1,ione)

        coo1=0.0_PREC
        coo2=0.0_PREC

        if (iorb1.ne.iorb)  then
           coo=gec(kex)
           coo1=coo
           call copy (ngexp,excp(ibexp),ione,wk2,ione)
           call scal (ngexp,coo,wk2,ione)
           if (ilc(idexp).gt.1) then
              coo=gec(kex+norb2)
              coo2=coo
              call axpy (ngexp,coo,excp(ibexp+ngexp),ione,wk2,ione)
           endif

           if (engo(ipc).ne.zero) then
              call axpy (ngexp,engo(ipc),f4,ione,wk2,ione)
           endif

           call prod (ngorb1,psi(ibpot1),wk2)
           call add  (ngorb1,wk2,fock2)
        else
           if ((mm(iorb).gt.0).and.(ilc(idexp).gt.0)) then
              coo=gec(kex)
              coo2=coo
              call copy (ngexp,excp(ibexp),ione,wk2,ione)
              ngrid=max(ngexp,ngorb1)
              call prod  (ngorb1,psi(iborb1),wk2)
              call scal (ngrid,coo,wk2,ione)
              call add (ngrid,wk2,fock2)
           endif
        endif

        if (iprint(92).ne.0) then
           write(*,1000) iorn(iorb),bond(iorb),gut(iorb),iorn(iorb1),bond(iorb1),gut(iorb1),coo0,coo1,coo2,engo(ipc)
1000       format(4x,i2,1x,a5,1x,a1,1x,i2,1x,a5,1x,a1,1x,4f8.3)
        endif
     enddo
  endif

  if (iprint(93).ne.0) then
     write(*,*) 'fock: Coulomb potential contribution'
     call pmtx (nni,nmu(1),wk1,ione,ione,incrni,incrmu)
  endif

  if (iprint(94).ne.0) then
     write(*,*) 'fock: Exchange potential contribution'
     call pmtx (nni,nmu(1),fock2,ione,ione,incrni,incrmu)
  endif

  !   multiply Coulomb and exchange potentials by f2

  call prod (mxsize,f2,wk1)
  call prod (mxsize,f2,fock2)

  !  add one electron and coulomb pot. contributions
  !  are coulomb and exchange pot zero at infinity?

  call add   (mxsize,wk1,fock1)
  !    call copy (mxsize,wk1,ione,fock1,ione)

  return
end subroutine fock





subroutine checkasymp (n1,n2,a,n1st,n2st,incrn1,incrn2)
  use params

  integer :: i,j,n1,n2,n1st,n2st,incrn1,incrn2
  real (PREC), dimension(n1,n2) :: a

  write(6,*) '(E-V) values along eta=1 axis (R/2,+infty)'
  write(6,1000) (j, j=n2st,n2,incrn2)
  i=1
  write(6,1010) i, (a(i,j),j=n2st,n2,incrn2)

  write(6,*) '(E-V) values along eta=-1 axis (-infty,-R/2)'
  write(6,1000) (j, j=n2st,n2,incrn2)
  i=n1
  write(6,1010) i, (a(i,j),j=n2st,n2,incrn2)

01000 format(4i25)
01010 format(/,1x,i4,4e25.16,/(5x,4e25.16))

end subroutine checkasymp

function checkev(n1,n2,a)
  use params


  integer :: n1,n2
  real (PREC) :: checkev
  real (PREC), dimension(n1,n2) :: a

!  parameter (ione=1)

  checkev=0.0_PREC

  if (a(ione,n2)*a(n1,n2).lt.0.0_PREC) then
     checkev=-1.0_PREC
  endif

end function checkev



