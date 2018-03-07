! ***************************************************************************
! *                                                                         *
! *   Copyright (C) 2010 Jacek Kobus <jkob@fizyka.umk.pl>                   *
! *                                                                         *
! *   This program is free software; you can redistribute it and/or modify  *
! *   it under the terms of the GNU General Public License version 2 as     *
! *   published by the Free Software Foundation.                            *
! *                                                                         *
! ***************************************************************************
! ### exbe88 ###

!     Calculates exchange energy according to Becke eq.8 (PRA 38
!     (1988) 3098-3100)

module exbe88_m
  implicit none
contains
  function exbe88 (psi,wgt2,rhot,rhotup,rhotdown,grhot,grhotup,grhotdown, &
       wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7)
    use params
    use discret
    use commons8
    use util

    use blas_m
    use exocc_m
    use fdften_m
    use exbe88sup_m
    use multf4_m
    use nfng_m

    implicit none
    real (PREC) :: exbe88
    integer :: i,iborb,iorb,isiorb,nmut
    real (PREC) ::  const43,ocdown,ocup,wxalpha
    real (PREC), dimension(*) :: psi,wgt2,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7,rhot,rhotup,rhotdown,grhot,grhotup,grhotdown

    parameter (const43=4.0_PREC/3.0_PREC)

    do i=1,mxsize
       rhotup(i)  =0.0_PREC
       rhotdown(i)=0.0_PREC
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

       !        store total spin densities
       call add(isiorb,wk1,rhotup)
       call add(isiorb,wk2,rhotdown)
10     continue
    enddo

    do i=1,mxsize
       wk1(i)=rhotup(i)**const43
       wk2(i)=rhotdown(i)**const43
    enddo

    call add (mxsize,wk1,wk2)

    !    take care of f4 factor
    call multf4(wk2)

    wxalpha=fdften(alphaf)*dot(mxsize,wgt2,ione,wk2,ione)

    !    calculate (nabla rho nabla rho)

    ! FIXME
    call nfng (rhotup,rhotup,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7)
    call copy(mxsize,wk7,ione,grhotup,ione)

    call nfng (rhotdown,rhotdown,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7)
    call copy(mxsize,wk7,ione,grhotdown,ione)

    exbe88=wxalpha+exbe88sup(wgt2,rhotup,grhotup,wk0,wk1)+exbe88sup(wgt2,rhotdown,grhotdown,wk0,wk1)

  end function exbe88
end module exbe88_m
