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
! ### Eab1DFT ###

!     Evaluates off-diagonal Lagrange multipliers in cases of a local exchange
!     approximation

!      ipc12=iorb1+(iorb2-1)*norb
!      ipc21=iorb2+(iorb1-1)*norb
!      e(iorb2,iorb1) = e(2,1) = engo(ipc21)=<1|T_k +V_n+V_C-V_x|2>
!      e(iorb1,iorb2) = e(1,2) = engo(ipc12)=<2|T_k +V_n+V_C-V_x|1>

subroutine Eab1DFT(iorb1,iorb2,psi,pot,excp,e,f0,wgt1,wgt2,wk0,wk1,wk2,wk3)
  use params
  use discret
  use memory
  use scf
  use commons8

  implicit none
  integer :: i1beg1,i1beg2,iborb,iborb1,iborb2,ibpot,iorb,iorb1,iorb2,ipc1,ipc2,&
       ipc12,isiorb1,isiorb2,isipot,isym,ngorb,ngorb1,ngorb2,ngpot1,&
       ngpot2,ngrid,norb2,nmut

  real (PREC) :: oc,w,woneel,wtwoel
  real (PREC), dimension(*) :: psi,pot,excp,e,f0,wgt1,wgt2,wk0,wk1,wk2,wk3
  real (PREC), external :: dot


! FIXME!!!
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
  
  !     calculate derivatives over mu and ni 
  
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
     !        from e term (in toten this term is correctly added to the
     !        kinetic energy contribution); e enters the expression with
     !        minus sign which is already incorporated in e
     
     w=dble(mm(iorb1)*mm(iorb1))
     
     call copy (ngorb,f0,ione,wk0,ione)
     call axpy (ngorb,w,e,ione,wk0,ione)
  endif
  
  call proda (ngorb,psi(i1beg1),wk0,wk1)
  call prod (ngorb,psi(i1beg2),wk1)
  
  woneel=dot(ngorb,wgt1,ione,wk1,ione)
  
  ipc1=iorb1+(iorb2-1)*norb
  ipc2=iorb2+(iorb1-1)*norb
  
  iborb1=i1b(iorb1)
  isiorb1=i1si(iorb1)
  iborb2=i1b(iorb2)
  isiorb2=i1si(iorb2)

  !     Add contributions from the Coulomb and local exchange potentials.
  !     In the local exchange approximtion the Coulomb potential includes
  !     also the contribution from a given orbital
  
  !     Calculate the Coulomb potential contribution from all the orbitals
  
  call zeroArray(mxsize,wk0)
  
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
  
  call prod2 (isiorb2,psi(iborb2),excp(length3-mxsize),wk1)
  call add (isiorb2,wk0,wk1)
  
  !     To complete the integrand wk1 has to be multiplied by psi(iborb1)
  
  ngrid=min(isiorb1,isiorb2)
  call prod (ngrid,psi(iborb1),wk1)
  wtwoel=dot(ngrid,wgt2,ione,wk1,ione)
  
  engo(ipc12)=woneel+wtwoel
  
end subroutine Eab1DFT




