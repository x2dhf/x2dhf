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
! ### EaDFT ###

!     Evaluates an eigenvalue of the Fock equation for a given orbital
!     using a DFT exchange potential
!
subroutine EaDFT(iorb,psi,pot,excp,e,f0,wgt1,wgt2,wk0,wk1,wk2,wk3)
  use params
  use discret
  use memory
  use scf
  use commons8
  use diffmu_m
  use diffnu_m
  use blas_m

  implicit none
  integer :: iorb,iorb1,i1beg,i1beg1,i2beg,i2beg1,ipc,isym,nmut,ngorb,ngpot,ngorb1,ngpot1
  real (PREC) :: w,woneel,wtwoel
  real (PREC), dimension(*) :: psi,pot,excp,e,f0,wgt1,wgt2,wk0,wk1,wk2,wk3

  if (ifix(iorb).eq.1) return

  call zeroArray(mxsize,wk0)
  call zeroArray(mxsize,wk1)

  wtwoel=zero

  i1beg=i1b(iorb)
  i2beg=i2b(iorb)
  nmut=i1mu(iorb)
  ngorb=i1si(iorb)
  ngpot=i2si(iorb)

  if (ngorb.ne.ngpot) then
     write(*,*) 'Orbitals and corresponding Coulomb potentials have to be defined on the same number of grid points'
     stop 'EaDFT'
  endif

  ipc=iorb+(iorb-1)*norb
  engo(ipc)=zero
  isym=isymOrb(iorb)

  !     calculate derivatives over mu and ni

  call putin (nni,nmut,isym,psi(i1beg),wk3)
  call diffnu (nmut,wk3,wk0,wk1,wk2)
  call putout (nni,nmut,wk1,wk0)

  call diffmu (nmut,wk3,wk2)
  call putout (nni,nmut,wk0,wk2)

  !     add contribution from derivatives over mu and ni

  call add (ngorb,wk0,wk1)

  call copy (ngorb,f0,ione,wk0,ione)

  if (mm(iorb).ne.0) then
     !        e enters the expression with minus sign which is already
     !        incorporated in e

     w=dble(mm(iorb)*mm(iorb))
     call axpy (ngorb,w,e,ione,wk0,ione)
  endif

  call proda (ngorb,psi(i1beg),wk0,wk1)
  call prod  (ngorb,psi(i1beg),wk1)

  woneel=dot(ngorb,wgt1,ione,wk1,ione)

  if(iprint(67).ne.0) then
     write(*,*) 'woneel', woneel
  endif
  ! (jk)
  call zeroArray(mxsize,wk2)

  if (nel.gt.1) then

     !        Add contributions from coulomb and local exchange potential.
     !        In the local exchange approximation the coulomb potential
     !        includes also the contribution from the orbital

     do iorb1=1,norb
        if (inhyd(iorb1).eq.1) goto 10
        i1beg1=i1b(iorb1)
        i2beg1=i2b(iorb1)
        ngorb1=i1si(iorb1)
        ngpot1=i2si(iorb1)
        call axpy (ngpot1,occ(iorb1),pot(i2beg1),ione,wk2,ione)
10      continue
     enddo
     call prod  (ngorb,psi(i1beg),wk2)

     !        multiply the local exchange potential by psi(i1beg) and
     !        add the result to the coulomb potential (multiplied by
     !        psi(i1beg))

     call prod2 (ngorb,psi(i1beg),excp(length3-mxsize),wk1)
     call add (ngorb,wk1,wk2)
  endif

  !     to complete the integrand wk2 has to be multiplied by psi(i1beg)

  call prod (ngorb,psi(i1beg),wk2)
  wtwoel=dot(ngorb,wgt2,ione,wk2,ione)

  engo(ipc) = woneel+wtwoel
  eng(iorb) = engo(ipc)

  if (iprint(68).ne.0) then
     write(*,*) 'EaDFT: woneel, wtwoel, eng ',iorb,woneel,wtwoel,eng(iorb)
  endif

  if (iprint(69).ne.0) then
     write(*,*) 'EaDFT: ',iorn(iorb),' ',bond(iorb)
     call pmtx (nni,nmu(1),psi(i1beg),ione,ione,incrni,incrmu)
  endif

  if (iprint(70).ne.0) then
     write(*,*) 'EaDFT: wgt1'
     call pmtx (nni,nmu(1),wgt1,ione,ione,incrni,incrmu)
     write(*,*) 'EaDFT: wgt2'
     call pmtx (nni,nmu(1),wgt2,ione,ione,incrni,incrmu)
  endif

  if (iprint(71).ne.0) then
     write(*,*) 'EaDFT: f0'
     call pmtx (nni,nmu(1),f0,ione,ione,incrni,incrmu)
  endif


  if (iprint(72).ne.0) then
     write(*,*) 'EaDFT: wk1'
     call pmtx (nni,nmu(1),wk1,ione,ione,incrni,incrmu)
     write(*,*) 'EaDFT: wk2'
     call pmtx (nni,nmu(1),wk2,ione,ione,incrni,incrmu)
  endif

end subroutine EaDFT
