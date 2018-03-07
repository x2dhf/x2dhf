! ***************************************************************************
! *                                                                         *
! *   Copyright (C) 2010 Jacek Kobus <jkob@fizyka.umk.pl>                   *
! *                                                                         *
! *   This program is free software; you can redistribute it and/or modify  *
! *   it under the terms of the GNU General Public License version 2 as     *
! *   published by the Free Software Foundation.                            *
! *                                                                         *
! ***************************************************************************
! ### exxalpha ###

!     Calculates exchange energy according to Xalpha scheme

module exxalpha_m
  implicit none
contains
  function exxalpha (psi,wgt2,rhot,rhotup,rhotdown,grhot,grhotup,grhotdown,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7)
    use params
    use discret
    use commons8
    use util

    use blas_m
    use exocc_m
    use fdften_m
    use multf4_m
    
    implicit none
    integer :: i,iborb,iorb,isiorb,nmut
    real (PREC) :: exxalpha
    real (PREC) :: const43,ocdown,ocup,w
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
       rhotup(i)=rhotup(i)**const43
       rhotdown(i)=rhotdown(i)**const43
    enddo

    call add (mxsize,rhotdown,rhotup)

    !    take care of F4 factor
    call multf4(rhotup)

    w=dot(mxsize,wgt2,ione,rhotup,ione)

    exxalpha=fdften(alphaf)*w

  end function exxalpha
end module exxalpha_m
