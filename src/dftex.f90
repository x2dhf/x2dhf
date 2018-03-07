! ***************************************************************************
! *                                                                         *
! *   Copyright (C) 2010 Jacek Kobus <jkob@fizyka.umk.pl>                   *
! *                                                                         *
! *   This program is free software; you can redistribute it and/or modify  *
! *   it under the terms of the GNU General Public License version 2 as     *
! *   published by the Free Software Foundation.                            *
! *                                                                         *
! ***************************************************************************
! ### dftex ###

!     Calculates exchange contribution due to various DFT functionals

module dftex_m
  implicit none
contains
  subroutine dftex (psi,pot,wgt2,wk0,wk1,wk2,wk3,rhot,rhotup,rhotdown, &
       grhot,grhotup,grhotdown,wk10,wk11,wk12,wk13)
    use params
    use discret
    use commons8
    use util

    use blas_m
    use exxalpha_m
    use exbe88_m
    use expw86_m

    implicit none
    integer :: i,iborb,ibpot,iorb,isiorb,isipot
    real (PREC) :: oc,w,wex,wndc
    real (PREC), dimension(*) :: psi,pot,wgt2,wk0,wk1,wk2,wk3,rhot,rhotup,rhotdown, &
         grhot,grhotup,grhotdown,wk10,wk11,wk12,wk13

    wndc =0.0_PREC
    wex=0.0_PREC

    !     calculate the coulomb potential contribution from all orbitals
    !     (include 1/2 factor )

    do i=1,mxsize
       wk2(i)=0.0_PREC
    enddo

    do iorb=1,norb
       ibpot=i2b(iorb)
       isipot=i2si(iorb)
       oc=occ(iorb)/two
       call axpy (isipot,oc,pot(ibpot),ione,wk2,ione)
    enddo

    !     contribution from the Coulomb interaction

    do iorb=1,norb
       iborb=i1b(iorb)
       isiorb=i1si(iorb)
       call prod2 (isiorb,psi(iborb),psi(iborb),wk0)
       call prod (isiorb,wk2,wk0)
       call scal (isiorb,occ(iorb),wk0,ione)
       w=dot(isiorb,wgt2,ione,wk0,ione)
       wndc=wndc+w
    enddo

    !     DFT exchange energy corrections

    edftex=0.0_PREC
    if     (idftex.eq.1) then
       edftex=exxalpha(psi,wgt2,rhot,rhotup,rhotdown,grhot,grhotup,grhotdown, &
            wk0,wk1,wk2,wk3,wk10,wk11,wk12,wk13)
    elseif (idftex.eq.2) then
       edftex=exbe88(psi,wgt2,rhot,rhotup,rhotdown,grhot,grhotup,grhotdown,   &
            wk0,wk1,wk2,wk3,wk10,wk11,wk12,wk13)
    elseif (idftex.eq.3) then
       edftex=expw86(psi,wgt2,rhot,rhotup,rhotdown,grhot,grhotup,grhotdown,   &
            wk0,wk1,wk2,wk3,wk10,wk11,wk12,wk13)
    endif

  end subroutine dftex
end module dftex_m

