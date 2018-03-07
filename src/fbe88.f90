! ***************************************************************************
! *                                                                         *
! *   Copyright (C) 2010 Jacek Kobus <jkob@fizyka.umk.pl>                   *
! *                                                                         *
! *   This program is free software; you can redistribute it and/or modify  *
! *   it under the terms of the GNU General Public License version 2 as     *
! *   published by the Free Software Foundation.                            *
! *                                                                         *
! ***************************************************************************
! ### fbe88 ###

!     Calculates the Becke gradient-corrected exchange potential (PRA 88
!     (1988) 3098-3100)

module fbe88_m
  implicit none
contains
  subroutine fbe88 (psi,f4,rhot,rhotup,rhotdown,grhot,grhotup,grhotdown, &
       wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7)
    use params
    use discret
    use commons8
    use util

    use blas_m
    use exocc_m
    use fbe88sup_m

    implicit none
    integer :: i,iborb,iorb,isiorb,nmut
    real (PREC) :: ocdown,ocup
    real (PREC), dimension(*) :: psi,f4,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7,rhot,rhotup,rhotdown, &
         grhot,grhotup,grhotdown

    do i=1,mxsize
       rhotup(i)   =0.0_PREC
       rhotdown(i) =0.0_PREC
    enddo

    !     calculate total densities due to up and down spins
    do iorb=1,norb
       if (inhyd(iorb).eq.1) goto 10
       iborb=i1b(iorb)
       isiorb=i1si(iorb)
       nmut=i1mu(iorb)

       call exocc (iorb,ocup,ocdown)

       call prod2 (isiorb,psi(iborb),psi(iborb),wk1)
       call scal (isiorb,ocup,wk1,ione)

       call prod2 (isiorb,psi(iborb),psi(iborb),wk2)
       call scal (isiorb,ocdown,wk2,ione)

       !        store total densities
       call add(isiorb,wk1,rhotup)
       call add(isiorb,wk2,rhotdown)
10     continue
    enddo

    call fbe88sup (rhotup,grhotup,rhot,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7)

    call fbe88sup (rhotdown,grhotdown,rhot,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7)

    call copy(mxsize,grhotup,ione,wk2,ione)
    call add (mxsize,grhotdown,wk2)

    call prod (mxsize,f4,wk2)

  end subroutine fbe88
end module fbe88_m
