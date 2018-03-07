! ***************************************************************************
! *                                                                         *
! *   Copyright (C) 2010 Jacek Kobus <jkob@fizyka.umk.pl>                   *
! *                                                                         *
! *   This program is free software; you can redistribute it and/or modify  *
! *   it under the terms of the GNU General Public License version 2 as     *
! *   published by the Free Software Foundation.                            *
! *                                                                         *
! ***************************************************************************
! ### flypcs ###

!     Calculates (and returns in wk2) the correlation potential of Lee,
!     Yang and Parr (PRB 37 (1988) 785, eq.23, i.e. using the
!     closed-shell formula.

subroutine flypcs (psi,f4,rhot,rhotup,rhotdown,grhot,grhotup,grhotdown, &
     wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7)
  use params
  use discret
  use commons8
  use blas_m

  implicit none
  integer :: iborb,iorb,isiorb,nmut
  real (PREC) :: ocdown,ocup
  real (PREC), dimension(*) :: psi,f4,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7, &
       rhot,rhotup,rhotdown,grhot,grhotup,grhotdown

  call zeroArray(mxsize,rhotup)
  call zeroArray(mxsize,rhotdown)

  !     calculate total densities due to up and down spins
  do iorb=1,norb
     if (inhyd(iorb).eq.1) goto 10
     iborb=i1b(iorb)
     isiorb=i1si(iorb)
     nmut=i1mu(iorb)

     call exocc (iorb,ocup,ocdown)

     if (ocup.ne.ocdown) then
        write(*,*) "Warning: This implementation of LYP potential is only valid for closed shell systems."
        stop 'flypcs'
     endif

     call prod2 (isiorb,psi(iborb),psi(iborb),wk1)
     call scal (isiorb,ocup,wk1,ione)

     call prod2 (isiorb,psi(iborb),psi(iborb),wk2)
     call scal (isiorb,ocdown,wk2,ione)

     !        store total densities
     call add(isiorb,wk1,rhotup)
     call add(isiorb,wk2,rhotdown)
10   continue
  enddo

  call copy(mxsize,rhotup,ione,rhot,ione)
  call add(mxsize,rhotdown,rhot)

  call flypsupcs (rhot,rhotup,rhotdown,grhot,grhotup,grhotdown,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7)

  call prod (mxsize,f4,wk7)
  call copy(mxsize,wk7,ione,wk2,ione)

end subroutine flypcs

