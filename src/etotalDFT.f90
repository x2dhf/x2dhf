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
! ### etotalDFT ###

!     Calculates total energy using several DFT functionals

module etotalDFT_m
  implicit none
contains
  subroutine etotalDFT (psi,pot,excp,e,f0,wgt1,wgt2,wk0,wk1,wk2,wk3,rhot,rhotup,rhotdown,&
       grhot,grhotup,grhotdown,wk10,wk11,wk12,wk13)
    use params
    use discret
    use commons8
    use util

    use blas_m
    use diffmu_m
    use diffnu_m
    use eclyptot_m
    use ecvwntot_m
    use exxalpha_m
    use exbe88_m
    use expw86_m
    use expw91_m
    use n2f_m
    use putin_m
    use putout_m
    use testn2f_m
    use testnfng_m

    implicit none
    integer :: i,iborb,ibpot,iorb,isiorb,isipot,isym,nmut,norb2

    real (PREC) :: oc,w,wcorr,wex,wndc,woneel
    real (PREC), dimension(*) :: psi,pot,excp,e,f0,wgt1,wgt2,wk0,wk1,wk2,wk3,rhot,rhotup,rhotdown,&
         grhot,grhotup,grhotdown,wk10,wk11,wk12,wk13

    if (iprint(75).ne.0) then
       !        test n2f and nfnf routines
       print *,'Testing n2f'
       call testn2f(rhot,wk0,wk1,wk2,wk3)

       print *,' '
       print *,'Testing nfng'
       call testnfng(rhot,rhotup,rhotdown,grhotdown,grhotup,wk0,wk1,wk2,wk3,wk10)

       call n2f(psi,wk0,wk1,wk2,wk3)
       call prod(mxsize,psi,wk3)
    endif

    !     calculate first contributions from one particle operators and
    !     Coulomb potential contributions within the same shell

    norb2=norb*norb

    vnt=zero
    vkt=zero
    woneel=zero
    wndc =zero
    wex=zero
    wcorr=zero

    do iorb=1,norb
       if (inhyd(iorb).eq.1) goto 10
       iborb=i1b(iorb)
       isiorb=i1si(iorb)
       nmut=i1mu(iorb)
       ibpot=i2b(iorb)
       isipot=i2si(iorb)
       oc=occ(iorb)

       isym=isymOrb(iorb)

       !        calculate derivatives over mu and ni variables by means of matrix
       !        multiplication

       call putin  (nni,nmut,isym,psi(iborb),wk3)
       call diffnu (nmut,wk3,wk0,wk1,wk2)
       call putout (nni,nmut,wk1,wk0)
       call diffmu (nmut,wk3,wk2)
       call putout (nni,nmut,wk0,wk2)

       !        add derivatives over mu and ni

       call add (isiorb,wk0,wk1)

       !        add contribution from phi part of laplasian
       !        e enters the expression with minus sign which is already
       !        incorporated in e

       if (mm(iorb).ne.0) then
          w=dble(mm(iorb)*mm(iorb))
          call copy (isiorb,e,ione,wk0,ione)
          if (mm(iorb).ne.1) then
             call scal (isiorb,w,wk0,ione)
          endif
          call prod  (isiorb,psi(iborb),wk0)
          call add   (isiorb,wk0,wk1)
       endif

       call copy (isiorb,wk1,ione,wk2,ione)
       call prod  (isiorb,psi(iborb),wk2)
       w=dot(isiorb,wgt1,ione,wk2,ione)
       vk(iorb)=w

       call copy (isiorb,f0,ione,wk0,ione)
       call prod (isiorb,psi(iborb),wk0)

       call prod2 (isiorb,psi(iborb),wk0,wk2)
       w=dot(isiorb,wgt1,ione,wk2,ione)
       vn(iorb)=w

       call add (isiorb,wk0,wk1)
       call prod (isiorb,psi(iborb),wk1)

       w =dot(isiorb,wgt1,ione,wk1,ione)
       woneel=woneel+oc*w

       vnt=vnt+oc*vn(iorb)
       vkt=vkt+oc*vk(iorb)
10     continue
    enddo

    evt=woneel
    etot=evt+z1*z2/r
    virrat=evt/vkt
    enkin=vkt
    ennucel=vnt
    encoul=wndc
    enexch=wex

    encouldft=wndc
    enexchdft=wex

    edftcorr=wcorr

    !     contribution from coulomb interaction within the same shell
    !     calculate the coulomb potential contribution from all orbitals
    !     (include 1/2 factor )

    do i=1,mxsize
       wk2(i)=zero
    enddo

    do iorb=1,norb
       if (inhyd(iorb).eq.1) goto 20
       ibpot=i2b(iorb)
       isipot=i2si(iorb)
       oc=occ(iorb)/two
       call axpy (isipot,oc,pot(ibpot),ione,wk2,ione)
20     continue
    enddo

    !     contribution from the Coulomb interaction

    do iorb=1,norb
       if (inhyd(iorb).eq.1) goto 30
       iborb=i1b(iorb)
       isiorb=i1si(iorb)
       call prod2 (isiorb,psi(iborb),psi(iborb),wk0)
       call prod (isiorb,wk2,wk0)
       call scal (isiorb,occ(iorb),wk0,ione)
       w=dot(isiorb,wgt2,ione,wk0,ione)
       wndc=wndc+w
30     continue
    enddo
    encoul=wndc
    encouldft=wndc    

    !    if (nel.eq.1) return
    
    !     DFT exchange energy corrections

    if     (idftex.eq.1) then
       wex=exxalpha(psi,wgt2,rhot,rhotup,rhotdown,grhot,grhotup,grhotdown,wk0,wk1,wk2,wk3,wk10,wk11,wk12,wk13)
    elseif (idftex.eq.2) then
       wex=exbe88(psi,wgt2,rhot,rhotup,rhotdown,grhot,grhotup,grhotdown,wk0,wk1,wk2,wk3,wk10,wk11,wk12,wk13)
    elseif (idftex.eq.3) then
       wex=expw86(psi,wgt2,rhot,rhotup,rhotdown,grhot,grhotup,grhotdown,wk0,wk1,wk2,wk3,wk10,wk11,wk12,wk13)
    elseif (idftex.eq.4) then
       wex=expw91(psi,wgt2,rhot,rhotup,rhotdown,grhot,grhotup,grhotdown,wk0,wk1,wk2,wk3,wk10,wk11,wk12,wk13)
    endif

    !     DFT correlation energy corrections

    if (idftcorr.eq.1) then
       wcorr=eclyptot(psi,wgt2,rhot,rhotup,rhotdown,grhot,grhotup,grhotdown,wk0,wk1,wk2,wk3,wk10,wk11,wk12,wk13)
    elseif (idftcorr.eq.2) then
       wcorr=ecvwntot(psi,wgt2,rhot,rhotup,rhotdown,grhot,grhotup,grhotdown,wk0,wk1,wk2,wk3,wk10,wk11,wk12,wk13)
    endif
    evt=woneel+wndc+wex+wcorr

    etot=evt+z1*z2/r
    virrat=evt/vkt

    enkin=vkt
    ennucel=vnt
    encoul=wndc
    enexch=wex

    encouldft=wndc
    enexchdft=wex+wcorr

    edftcorr=wcorr

    entot=enkin+ennucel+encoul+enexch+edftcorr+z1*z2/r

    if (iprint(77).ne.0) then
       write(*,*)'etotalDFT - woneel,wndc,wex,etotal ', woneel,wndc,wex,evt
    endif

  end subroutine etotalDFT
end module etotalDFT_m
