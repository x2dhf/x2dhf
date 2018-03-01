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
! ### Eab2DFT ###

!     Evaluates the off-diagonal Lagrange multipliers in case of a local
!     exchange approximation

subroutine Eab2DFT(iorb1,iorb2,psi,pot,excp,wgt2,wk0,wk1)
  use params
  use discret
  use scf
  use commons8

  implicit none

  integer :: length 
  integer :: i,iborb,ibpot,iorb1,iorb2,iborb1,isiorb1,iborb2,iorb,ipc1,ipc2,&
       isiorb2,isipot,ngrid,nmut
  real (PREC) :: engo1,engo2,engoprv1,engoprv2,oc,wtwoel
  real (PREC), dimension(*) :: psi,pot,excp,wgt2,wk0,wk1
  real (PREC), external :: dot

  if (mgx(6,iorb1).ne.mgx(6,iorb2)) return
  if (ige(iorb1).ne.ige(iorb2)) return

! FIXME 
  ipc1=iorb1+(iorb2-1)*norb
  ipc2=iorb2+(iorb1-1)*norb
  engo(ipc1)=zero
  
  iborb1=i1b(iorb1)
  isiorb1=i1si(iorb1)
  iborb2=i1b(iorb2)
  isiorb2=i1si(iorb2)
  
  !     Add contributions from the Coulomb and local exchange potentials.
  !     In the local exchange approximtion the Coulomb potential includes
  !     also the contribution from a given orbital
  
  !     Calculate the Coulomb potential contribution from all the orbitals
  
  do i=1,mxsize
     wk0(i)=0.0_PREC
  enddo
  
  do iorb=1,norb
     iborb=i1b(iorb)
     ibpot=i2b(iorb)
     nmut=i1mu(iorb)
     isipot=i2si(iorb)
     oc=occ(iorb)
     call axpy (isipot,oc,pot(ibpot),ione,wk0,ione)
  enddo
  
  call prod  (isiorb2,psi(iborb2),wk0)
  
  !     Multiply the local exchange potential by psi(iborb2) 
  !     and add the result to the Coulomb potential

  ! FIXME length must be initialized
  call prod2 (isiorb2,psi(iborb2),excp(length-mxsize),wk1)
  call add (isiorb2,wk0,wk1)
  
  !     To complete the integrand wk1 has to be multiplied by psi(iborb1)
  
  ngrid=min(isiorb1,isiorb2)
  call prod (ngrid,psi(iborb1),wk1)
  wtwoel=dot(ngrid,wgt2,ione,wk1,ione)

  ! FIXME 
  engo(ipc1)= wtwoel
  engo(ipc2)= wtwoel
  
  ! FIXME
  !        Damping factors are set to 0 by default (obsolete).
  ! engo1/2 must be initialized
  engo(ipc1)=sflagra*((1.0_PREC-dflagra)*engo1+dflagra*engoprv1)
  engo(ipc2)=sflagra*((1.0_PREC-dflagra)*engo2+dflagra*engoprv2)
  
  if (iprint(48).ne.0) then
     write(*,'(a8,i4,e16.6,i4,a8,i4,a8,2e16.8)') 'Eab2DFT: ', &
          lmtype,sflagra,iorn(iorb),bond(iorb1),iorn(iorb1),bond(iorb1),engo(ipc1),engo(ipc2)
  endif
  
end subroutine Eab2DFT
