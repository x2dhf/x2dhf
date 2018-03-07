! ***************************************************************************
! *                                                                         *
! *   Copyright (C) 2010 Jacek Kobus <jkob@fizyka.umk.pl>                   *
! *                                                                         *
! *   This program is free software; you can redistribute it and/or modify  *
! *   it under the terms of the GNU General Public License version 2 as     *
! *   published by the Free Software Foundation.                            *
! *                                                                         *
! ***************************************************************************
! c ### eclyptot ###

!     Calculates exchange correlation energy according to Vosko, Wilk, Nusair
!     Can. J. Phys. 58 (1980) 1200

module ecvwntot_m
  implicit none
contains
  function ecvwntot (psi,wgt2,rhot,rhotup,rhotdown, &
       grhot,grhotup,grhotdown,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7)
    use params
    use discret
    use commons8
    use util
    
    use blas_m
    use exocc_m
    use dftauxfun_m
    use multf4_m

    implicit none
    real (PREC) :: ecvwntot
    integer :: i,iborb,iorb,isiorb,nmut
    real (PREC) :: ck1,cl1,cm1,cn1,const16,constx,ocdown,ocup,x
    real (PREC), dimension(*) :: psi,wgt2,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7, &
         rhot,rhotup,rhotdown,grhot,grhotup,grhotdown

    parameter(ck1=0.03109070_PREC,cl1=-0.104980_PREC,cm1=3.727440_PREC,cn1=12.93520_PREC, &
         const16=1.0_PREC/6.0_PREC)

    constx=(three/(four*pii))**const16

    do i=1,mxsize
       rhotup(i)  =0.0_PREC
       rhotdown(i)=0.0_PREC
    enddo

    !     calculate total densities for up and down spins
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
       rhot(i)=rhotup(i)+rhotdown(i)
    enddo

    do i=1,mxsize
       if (rhot(i).lt.precis) then
          wk0(i)=0.0_PREC
       else
          x=constx*rhot(i)**(-const16)
          wk0(i)=rhot(i)*qvwn(x,ck1,cl1,cm1,cn1)
       endif

    enddo

    !     take care of f4 factor
    call multf4(wk0)

    ecvwntot=dot(mxsize,wgt2,ione,wk0,ione)

  end function ecvwntot
end module ecvwntot_m
