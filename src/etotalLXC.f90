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
! ### etotalLXC ###

! Calculates total energy using several LXC functionals

module etotalLXC_m
  implicit none
contains
  subroutine etotalLXC (psi,pot,excp,e,f0,wgt1,wgt2,wk0,wk1,wk2,wk3,rhot,rhotup,rhotdown,&
       grhot,grhotup,grhotdown,wk10,wk11,wk12,wk13)
    use params
    use discret
    use commons8
    use util

    use blas_m
    use diffmu_m
    use diffnu_m
    use multf4_m
    use n2f_m
    use nfng_m
    use putin_m
    use putout_m
    use testn2f_m
    use testnfng_m
    use zeroArray_m
    
    use xc_f90_types_m
    use xc_f90_lib_m

    implicit none
    integer :: i,iborb,ibpot,iorb,isiorb,isipot,isym,nmut,norb2

    real (PREC) :: etsum,oc,w,wcorr,wex,wndc,woneel
    real (PREC), dimension(*) :: psi,pot,excp,e,f0,wgt1,wgt2,wk0,wk1,wk2,wk3,rhot,rhotup,rhotdown,&
         grhot,grhotup,grhotdown,wk10,wk11,wk12,wk13

    TYPE(xc_f90_pointer_t) :: xc_func
    TYPE(xc_f90_pointer_t) :: xc_info
    integer :: func_id = 1

    
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

    ! calculate first contributions from one particle operators and Coulomb potential
    ! contributions within the same shell

    norb2=norb*norb

    vnt=zero
    vkt=zero
    woneel=zero
    
    do iorb=1,norb
       if (inhyd(iorb).eq.1) goto 10
       iborb=i1b(iorb)
       isiorb=i1si(iorb)
       nmut=i1mu(iorb)
       ibpot=i2b(iorb)
       isipot=i2si(iorb)
       oc=occ(iorb)

       isym=isymOrb(iorb)

       ! calculate derivatives over mu and ni variables by means of matrix multiplication

       call putin  (nni,nmut,isym,psi(iborb),wk3)
       call diffnu (nmut,wk3,wk0,wk1,wk2)
       call putout (nni,nmut,wk1,wk0)
       call diffmu (nmut,wk3,wk2)
       call putout (nni,nmut,wk0,wk2)

       ! add derivatives over mu and ni

       call add (isiorb,wk0,wk1)

       ! add contribution from phi part of laplasian e enters the expression with minus
       ! sign which is already incorporated in e

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

       if (iprint(78).ne.0) then
          etsum=etsum+oc*vk(iorb)
          write(*,7028) iorn(iorb),bond(iorb),gut(iorb),iorn(iorb),bond(iorb),gut(iorb),vk(iorb),oc,etsum
7028      format('<',i4,1x,a5,a1,'| T |',i4,1x,a5,a1,' >',26x,d25.16, f8.2, d25.16)

          etsum=etsum+oc*vn(iorb)
          write(*,7030) iorn(iorb),bond(iorb),gut(iorb),iorn(iorb),bond(iorb),gut(iorb),vn(iorb),oc,etsum
7030      format('<',i4,1x,a5,a1,'| V |',i4,1x,a5,a1,' >',26x,d25.16, f8.2, d25.16)
       endif


10     continue
    enddo

    evt=woneel
    etot=evt+z1*z2/r
    
    if (nel.eq.1) return

    wndc =zero
    wex=zero
    wcorr=zero
    
    ! contribution from coulomb interaction within the same shell

    ! calculate the coulomb potential contribution from all orbitals (include 1/2 factor )

    do i=1,mxsize
       wk2(i)=zero
       wk12(i)=zero
    enddo

    do iorb=1,norb
       if (inhyd(iorb).eq.1) goto 20
       ibpot=i2b(iorb)
       isipot=i2si(iorb)
       oc=occ(iorb)/two
       call axpy (isipot,oc,pot(ibpot),ione,wk2,ione)
20     continue
    enddo

    ! contribution from the Coulomb interaction

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

    ! DFT exchange energy corrections

    call zeroarray(mxsize,rhot)    
    call zeroarray(mxsize,wk12)
    call zeroarray(mxsize,wk13)    
    
    do iorb=1,norb
       if (inhyd(iorb).eq.1) cycle
       iborb =i1b (iorb)
       call prodas (mxsize,occ(iorb),psi(iborb),psi(iborb),rhot)
    enddo

    do func_id=1,lxcFuncs
       
       call xc_f90_func_init(xc_func, xc_info, lxcFuncs2use(func_id), XC_UNPOLARIZED)
       select case (xc_f90_info_family(xc_info))
       case(XC_FAMILY_LDA)
          call xc_f90_lda_exc(xc_func, mxsize, rhot(1), wk12(1))
       case(XC_FAMILY_GGA)
          ! calculate nabla rho nabla rho 
          call nfng (rhot,rhot,wk0,wk1,wk2,wk3,rhotup,rhotdown,grhotdown,grhotup)
          ! do i=1,mxsize
          !    print *,i,grhotup(i)
          ! enddo
          call xc_f90_gga_exc(xc_func, mxsize, rhot(1), grhotup(1), wk12(1))
       case(XC_FAMILY_HYB_GGA)
          ! calculate nabla rho nabla rho 
          call nfng (rhot,rhot,wk0,wk1,wk2,wk3,rhotup,rhotdown,grhotdown,grhotup)
          call xc_f90_gga_exc(xc_func, mxsize, rhot(1), grhotup(1), wk12(1))

       case default
          write(*,'("Error! Undefined libxc functional.")')
          stop 'etotalLXC'
       end select
       call xc_f90_func_end(xc_func)
       
       call add (mxsize,wk12,wk13)
    enddo

    call prod (mxsize,rhot,wk13)
    
    ! take care of F4 factor
    call multf4(wk13)

    wex=dot(mxsize,wgt2,ione,wk13,ione)
    
    ! DFT correlation energy corrections
    evt=woneel+wndc+wex+wcorr
    engt(1)=evt
    etot=evt+z1*z2/r
    virrat=evt/vkt
    enkin=vkt
    ennucel=vnt
    encoul=wndc
    enexch=wex

    encouldft=wndc
    enexchdft=wex

    edftcorr=wcorr

    entot=enkin+ennucel+encoul+enexch+edftcorr+z1*z2/r

    if (iprint(77).ne.0) then
       write(*,*)'etotalLXC - woneel,wndc,wex,etotal ', woneel,wndc,wex,evt
    endif

  end subroutine etotalLXC
end module etotalLXC_m
