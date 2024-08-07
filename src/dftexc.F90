! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2010-2023  Jacek Kobus 

module dftexc
  implicit none
contains
  ! ### eclyptot ###
  !
  !     Calculates correlation energy according to Lee, Yang, Prr (PRB 37
  !     (1988) 785-789)
  !
  function eclyptot (psi,wgt2,rhot,rhotup,rhotdown, &
       grhot,grhotup,grhotdown,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7)
    use params
    use discrete
    use commons
    use blas
    use elocc
    use nabla
    use utils

    implicit none
    real (PREC) :: eclyptot
    integer (KIND=IPREC) :: i,iborb,iorb,isiorb
    real (PREC) :: a,b,c,cf,d,const13,const83,const23,const53,const19,const118,ocdown,ocup,t1,t2

    real (PREC), dimension(*) :: psi,wgt2,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7,  &
         rhot,rhotup,rhotdown,grhot,grhotup,grhotdown

    parameter (const13=1.0_PREC/3.0_PREC,const83=8.0_PREC/3.0_PREC,const23=2.0_PREC/3.0_PREC, &
         const53=5.0_PREC/3.0_PREC,const19=1.0_PREC/9.0_PREC,const118=1.0_PREC/18.0_PREC)

    ! coefficients of the Colle-Salvetti formula
    parameter(a=0.049180_PREC,b=0.1320_PREC,c=0.25330_PREC,d=0.3490_PREC)

#ifdef BLAS    
    real (PREC) ddot
    external ddot
#endif

    cf=0.30_PREC*(three*pii*pii)**const23

    call zeroArray(mxsize,rhotup)
    call zeroArray(mxsize,rhotdown)

    !     calculate total densities for up and down spins
    do iorb=1,norb
       !if (inhyd(iorb).eq.1) goto 10
       if (inDFT(iorb).eq.0) cycle
       iborb=i1b(iorb)
       isiorb=i1si(iorb)

       call exocc (iorb,ocup,ocdown)

       call prod2 (isiorb,psi(iborb),psi(iborb),wk1)
       call dscal (isiorb,ocup,wk1,ione)

       call prod2 (isiorb,psi(iborb),psi(iborb),wk2)
       call dscal (isiorb,ocdown,wk2,ione)

       !        store total spin densities
       call add(isiorb,wk1,rhotup)
       call add(isiorb,wk2,rhotdown)
    enddo

    ! calculate nabla rho nabla rho (up)
    call nfng (rhotup,rhotup,wk0,wk1,wk2,wk3,wk4,wk5,wk6,grhotup)

    ! calculate nabla rho nabla rho (down)
    call nfng (rhotdown,rhotdown,wk0,wk1,wk2,wk3,wk4,wk5,wk6,grhotdown)


    ! calculate nabla^2 rho (up)
    call n2f (rhotup,wk0,wk1,wk2,rhot)

    ! calculate nabla^2 rho (down)
    call n2f (rhotdown,wk0,wk1,wk2,grhot)

    ! tw (up)
    do i=1,mxsize
       if (abs(rhotup(i)).lt.precis) then
          wk1(i)=0.0_PREC
       else
          wk1(i)=(grhotup(i)/rhotup(i)-rhot(i))/8.0_PREC
       endif
    enddo

    ! tw (down)
    do i=1,mxsize
       if (abs(rhotdown(i)).lt.precis) then
          wk2(i)=0.0_PREC
       else
          wk2(i)=(grhotdown(i)/rhotdown(i)-grhot(i))/8.0_PREC
       endif
    enddo

    ! wk7 - total density, wk6 - tw(total)
    do i=1,mxsize
       wk7(i)=rhotup(i)+rhotdown(i)
       if (abs(wk7(i)).lt.precis) then
          ! FIXME
          wk6(i)=two
       else
          wk6(i)=two*(one-(rhotup(i)*rhotup(i)+rhotdown(i)*rhotdown(i))/(wk7(i)*wk7(i)) )
       endif
    enddo

    ! wk4 - Eq.22 terms from square brackets except for rhot tw
    do i=1,mxsize
       wk4(i)=two**const23*cf*(rhotup(i)**const83+rhotdown(i)**const83)  &
            +const19*(rhotup(i)*wk1(i)+rhotdown(i)*wk2(i))               &
            +const118*(rhotup(i)*rhot(i)+rhotdown(i)*grhot(i))
    enddo

    ! wk7->rhot
    call dcopy(mxsize,wk7,ione,rhot,ione)

    ! wk4->grhotup
    call dcopy(mxsize,wk4,ione,grhot,ione)

    ! calculate nabla rho nabla rho (up+down)
    call nfng (rhot,rhot,wk0,wk1,wk2,wk3,wk4,wk5,wk6,grhotup)

    ! calculate nabla^2 rho (up+down)
    call n2f (rhot,wk0,wk1,wk2,grhotdown)

    ! wk5 - tw
    ! wk6 - gamma (total)
    do i=1,mxsize
       if (rhot(i).lt.precis) then
          ! FIXME
          wk5(i)=0.0_PREC
          wk6(i)=two
       else
          wk5(i)=(grhotup(i)/rhot(i)-grhotdown(i))/8.0_PREC
          wk6(i)=two*(one-(rhotup(i)*rhotup(i)+rhotdown(i)*rhotdown(i))/(rhot(i)*rhot(i)))
       endif
    enddo

    do i=1,mxsize
       if (rhot(i).lt.precis) then
          wk0(i)=0.0_PREC
       else
          ! wk4/grhot +  rhot tw
          t1=grhot(i)-rhot(i)*wk5(i)
          t2=two*b*rhot(i)**(-const53)*t1*exp(-c*rhot(i)**(-const13))
          wk0(i)=-a*(rhot(i)+t2)*wk6(i)/(one+d*rhot(i)**(-const13))
       endif
    enddo

    ! take care of f4 factor
    call multf4(wk0)

    eclyptot=ddot(mxsize,wgt2,ione,wk0,ione)

  end function eclyptot

  ! ### eclyptot ###
  !
  !     Calculates exchange correlation energy according to Vosko, Wilk, Nusair
  !     Can. J. Phys. 58 (1980) 1200
  !
  function ecvwntot (psi,wgt2,rhot,rhotup,rhotdown, &
       grhot,grhotup,grhotdown,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7)
    use params
    use discrete
    use commons
    use utils

    use blas
    use elocc
    use dftvxc

    implicit none
    real (PREC) :: ecvwntot
    integer (KIND=IPREC) :: i,iborb,iorb,isiorb,nmut
    real (PREC) :: ck1,cl1,cm1,cn1
    !real (PREC) :: ck2,cl2,cm2,cn2
    real (PREC) :: const16,constx,ocdown,ocup,x
    real (PREC), dimension(*) :: psi,wgt2,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7, &
         rhot,rhotup,rhotdown,grhot,grhotup,grhotdown

    ! see math.nist.gov/DFTdata/atomdata/node5.html#SECTION00021200000000000000
    ! paramagnetic parametrization:
    parameter(ck1=0.0310907_PREC,cl1=-0.10498_PREC,cm1=3.72744_PREC,cn1=12.9352_PREC)

    ! ferromagnetic parametrization (no needed for 0-spin states): 
    !parameter(ck2=0.01554535_PREC,cl2=-0.325000_PREC,cm2=7.060420_PREC,cn2=18.05780_PREC

    parameter(const16=1.0_PREC/6.0_PREC)
#ifdef BLAS    
    real (PREC) ddot
    external ddot
#endif
    
    constx=(three/(four*pii))**const16

    do i=1,mxsize
       rhotup(i)  =0.0_PREC
       rhotdown(i)=0.0_PREC
    enddo

    !     calculate total densities for up and down spins
    do iorb=1,norb
       !if (inhyd(iorb).eq.1) goto 10
       if (inDFT(iorb).eq.0) cycle       
       iborb=i1b(iorb)
       isiorb=i1si(iorb)
       nmut=i1mu(iorb)

       call exocc (iorb,ocup,ocdown)

       call prod2 (isiorb,psi(iborb),psi(iborb),wk1)
       call dscal (isiorb,ocup,wk1,ione)

       call prod2 (isiorb,psi(iborb),psi(iborb),wk2)
       call dscal (isiorb,ocdown,wk2,ione)

       !        store total spin densities
       call add(isiorb,wk1,rhotup)
       call add(isiorb,wk2,rhotdown)
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

    ! take care of f4 factor
    call multf4(wk0)

    ecvwntot=ddot(mxsize,wgt2,ione,wk0,ione)

  end function ecvwntot

  ! ### exxalpha ###
  !
  !     Calculates exchange energy according to Xalpha scheme
  !
  function exxalpha (psi,wgt2,rhot,rhotup,rhotdown,grhot,grhotup,grhotdown,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7)
    use params
    use discrete
    use commons
    use blas
    use elocc
    use dftvxc
    use utils

    implicit none
    integer (KIND=IPREC) :: i,iborb,iorb,isiorb,nmut
    real (PREC) :: exxalpha
    real (PREC) :: const43,ocdown,ocup,w
    real (PREC), dimension(*) :: psi,wgt2,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7,rhot,rhotup,rhotdown,grhot,grhotup,grhotdown

    parameter (const43=4.0_PREC/3.0_PREC)
#ifdef BLAS    
    real (PREC) ddot
    external ddot
#endif

    call zeroArray(mxsize,rhotup)
    call zeroArray(mxsize,rhotdown)

    ! calculate total densities due to up and down spins
    do iorb=1,norb
       !if (inhyd(iorb).eq.1) goto 10
       if (inDFT(iorb).eq.0) cycle
       iborb=i1b(iorb)
       isiorb=i1si(iorb)
       nmut=i1mu(iorb)

       call exocc (iorb,ocup,ocdown)
       
       call prod2 (isiorb,psi(iborb),psi(iborb),wk1)
       call dscal (isiorb,ocup,wk1,ione)

       call prod2 (isiorb,psi(iborb),psi(iborb),wk2)
       call dscal (isiorb,ocdown,wk2,ione)

       ! store total spin densities
       call add(isiorb,wk1,rhotup)
       call add(isiorb,wk2,rhotdown)
    enddo

    do i=1,mxsize
       rhotup(i)=rhotup(i)**const43
       rhotdown(i)=rhotdown(i)**const43
    enddo

    call add (mxsize,rhotdown,rhotup)

    ! take care of F4 factor
    call multf4(rhotup)

    w=ddot(mxsize,wgt2,ione,rhotup,ione)
    exxalpha=fdften(alphaf)*w

  end function exxalpha
  
  ! ### exbe88 ###
  !
  !     Calculates exchange energy according to Becke eq.8 (PRA 38
  !     (1988) 3098-3100)
  !
  function exbe88 (psi,wgt2,rhot,rhotup,rhotdown,grhot,grhotup,grhotdown, &
       wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7)
    use params
    use discrete
    use commons
    use blas
    use elocc
    use dftvxc
    use nabla
    use utils
    implicit none
    real (PREC) :: exbe88
    integer (KIND=IPREC) :: i,iborb,iorb,isiorb,nmut
    real (PREC) ::  const43,ocdown,ocup,wxalpha
    real (PREC), dimension(*) :: psi,wgt2,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7,rhot,rhotup,rhotdown,grhot,grhotup,grhotdown

    parameter (const43=4.0_PREC/3.0_PREC)
#ifdef BLAS    
    real (PREC) ddot
    external ddot
#endif

    do i=1,mxsize
       rhotup(i)  =0.0_PREC
       rhotdown(i)=0.0_PREC
    enddo

    ! calculate total densities due to up and down spins
    do iorb=1,norb
       !if (inhyd(iorb).eq.1) goto 10
       if (inDFT(iorb).eq.0) cycle
       iborb=i1b(iorb)
       isiorb=i1si(iorb)
       nmut=i1mu(iorb)

       call exocc (iorb,ocup,ocdown)

       call prod2 (isiorb,psi(iborb),psi(iborb),wk1)
       call dscal (isiorb,ocup,wk1,ione)

       call prod2 (isiorb,psi(iborb),psi(iborb),wk2)
       call dscal (isiorb,ocdown,wk2,ione)

       ! store total spin densities
       call add(isiorb,wk1,rhotup)
       call add(isiorb,wk2,rhotdown)

    enddo

    do i=1,mxsize
       wk1(i)=rhotup(i)**const43
       wk2(i)=rhotdown(i)**const43
    enddo

    call add (mxsize,wk1,wk2)

    ! take care of f4 factor
    call multf4(wk2)

    wxalpha=fdften(alphaf)*ddot(mxsize,wgt2,ione,wk2,ione)

    ! calculate (nabla rho nabla rho)

    ! FIXME
    call nfng (rhotup,rhotup,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7)
    call dcopy(mxsize,wk7,ione,grhotup,ione)

    call nfng (rhotdown,rhotdown,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7)
    call dcopy(mxsize,wk7,ione,grhotdown,ione)

    exbe88=wxalpha+exbe88sup(wgt2,rhotup,grhotup,wk0,wk1)+exbe88sup(wgt2,rhotdown,grhotdown,wk0,wk1)
    
  end function exbe88

  ! ### exbe88sup ###
  !
  !     Calculates exchange energy according to Becke's gradient formula
  !     (PRA 88 (1988) 3098-3100)
  !
  function exbe88sup (wgt2,rho,grho,wk0,wk1)
    use params
    use discrete
    use commons
    use blas
    use utils
    implicit none
    real (PREC) :: exbe88sup
    integer (KIND=IPREC) :: i
    real (PREC) :: ash,bbeta,const43,rho43,s,s2

    real (PREC), dimension(*) :: wgt2,rho,grho,wk0,wk1

    parameter (const43=4.0_PREC/3.0_PREC,bbeta=0.0042_PREC)
#ifdef BLAS    
    real (PREC) ddot
    external ddot
#endif
    ! sinh^{-1}(s)
    ash(s)=log(s+sqrt(1.0_PREC+s*s))

    ! grho = nabla rho  nabla rho
    ! |nabla rho| = sqrt(grho)

    do i=1,mxsize
       if (abs(rho(i)).lt.precis) then
          wk0(i)=0.0_PREC
       else
          rho43=rho(i)**const43
          s2=grho(i)/(rho43*rho43)
          s=sqrt(s2)
          ! fgga=-bbeta*s2/(one+6.0_PREC*bbeta*s*asinh(s))
          wk0(i)=-bbeta*s2/(one+6.0_PREC*bbeta*s*ash(s))*rho43
       endif
    enddo

    !     take care of f4 factor
    call multf4(wk0)

    exbe88sup=ddot(mxsize,wgt2,ione,wk0,ione)

  end function exbe88sup

  ! ### expw86 ###
  !
  !     Calculates exchange energy according to a formula of Parr and Wang
  !     Yue (PRB 33 (1986) 8800)
  !
  function expw86 (psi,wgt2,rhot,rhotup,rhotdown,grhot,grhotup,grhotdown, &
       wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7)
    use params
    use discrete
    use commons
    use utils
    use blas
    use elocc
    use dftvxc
    use nabla

    implicit none
    integer (KIND=IPREC) :: i,iborb,iorb,isiorb,nmut
    real (PREC) :: expw86
    real (PREC) :: ocdown,ocup
    real (PREC), dimension(*) :: psi,wgt2,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7,  &
         rhot,rhotup,rhotdown,grhot,grhotup,grhotdown
#ifdef BLAS    
    real (PREC) ddot
    external ddot
#endif

    do i=1,mxsize
       rhotup(i)  =0.0_PREC
       rhotdown(i)=0.0_PREC
    enddo

    ! calculate total densities due to up and down spins
    do iorb=1,norb
       !if (inhyd(iorb).eq.1) goto 10
       if (inDFT(iorb).eq.0) cycle
       iborb=i1b(iorb)
       isiorb=i1si(iorb)
       nmut=i1mu(iorb)

       call exocc (iorb,ocup,ocdown)

       call prod2 (isiorb,psi(iborb),psi(iborb),wk1)
       call dscal (isiorb,ocup,wk1,ione)

       call prod2 (isiorb,psi(iborb),psi(iborb),wk2)
       call dscal (isiorb,ocdown,wk2,ione)

       ! store total spin densities
       call add(isiorb,wk1,rhotup)
       call add(isiorb,wk2,rhotdown)
    enddo

    ! calculate (nabla rho nabla rho)
    ! FIXME
    call nfng (rhotup,rhotup,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7)
    call dcopy(mxsize,wk7,ione,grhotup,ione)

    call nfng (rhotdown,rhotdown,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7)
    call dcopy(mxsize,wk7,ione,grhotdown,ione)

    do i=1,mxsize
       rhotup(i)   =two*rhotup(i)
       rhotdown(i) =two*rhotdown(i)
       grhotup(i)  =four*grhotup(i)
       grhotdown(i)=four*grhotdown(i)
    enddo

    ! total exchange energy is calculated as
    ! Ex(rhoup,rhodown)=(1/2)Ex(2*rhoup)+(1/2)Ex(2*rhodown)

    expw86=fdften(alphaf)*(  half*expw86sup(wgt2,rhotup,grhotup,wk0,wk1) &
         + half*expw86sup(wgt2,rhotdown,grhotdown,wk0,wk1))

  end function expw86

  ! ### expw86sup ###
  !
  !     Calculates exchange energy according to a formula of Parr and Wang
  !     Yue (PRB 33 (1986) 8800)
  !
  function expw86sup (wgt2,rho,grho,wk0,wk1)
    use params
    use discrete
    use commons
    use blas
    use utils

    implicit none
    integer (KIND=IPREC) :: i
    real (PREC) :: expw86sup
    real (PREC) :: a,b,c,const,const23,const43,const83,const115,fgga,s,s2
    real (PREC), dimension(*) :: wgt2,rho,grho,wk0,wk1

    parameter (a=1.2960_PREC,b=14.0_PREC,c=0.20_PREC,const23=2.0_PREC/3.0_PREC, &
         const43=4.0_PREC/3.0_PREC,const83=8.0_PREC/3.0_PREC,const115=1.0_PREC/15.0_PREC)
#ifdef BLAS    
    real (PREC) ddot
    external ddot
#endif

    ! grho = nabla rho  nabla rho
    ! |nabla rho| = sqrt(grho)

    const=(24.0_PREC*pii*pii)**(-const23)
    do i=1,mxsize
       if (abs(rho(i)).lt.precis) then
          wk0(i)=0.0_PREC
       else
          s2=const*grho(i)/rho(i)**const83
          s=sqrt(s2)
          fgga=(one+a*s2+b*s2*s2+c*s2*s2*s2)**const115

          ! Parr & Yue 1986
          wk0(i)=rho(i)**const43*fgga

          ! Langreth-Mehl
          ! wk0(i)=rho(i)**const43*(one+1.521*0.0864*s2)

          ! GEA
          ! wk0(i)=rho(i)**const43*(one+0.0864*s2)

          ! LDA
          ! wk0(i)=rho(i)**const43
       endif
    enddo

    ! take care of F4 factor
    call multf4(wk0)
    expw86sup=ddot(mxsize,wgt2,ione,wk0,ione)
  end function expw86sup

  ! ### expw91 ###
  !
  !     Calculates exchange energy according to a formula of Parr and Wang
  !     Yue Phys. Rev. B 54, 16 533 (1996), Phys. Rev. B 45, 13 244
  !     (1992).
  !
  function expw91 (psi,wgt2,rhot,rhotup,rhotdown, &
       grhot,grhotup,grhotdown,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7)
    use params
    use discrete
    use commons
    use blas
    use elocc
    use dftvxc
    use nabla
    use utils
    
    implicit none
    integer (KIND=IPREC) :: i,iborb,iorb,isiorb,nmut
    real (PREC) :: expw91
    real (PREC) :: ocdown,ocup

    real (PREC), dimension(*) :: psi,wgt2,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7,  &
         rhot,rhotup,rhotdown,grhot,grhotup,grhotdown
#ifdef BLAS    
    real (PREC) ddot
    external ddot
#endif


    do i=1,mxsize
       rhotup(i)  =0.0_PREC
       rhotdown(i)=0.0_PREC
    enddo

    !     calculate total densities due to up and down spins
    do iorb=1,norb
       !if (inhyd(iorb).eq.1) goto 10
       if (inDFT(iorb).eq.0) cycle
       iborb=i1b(iorb)
       isiorb=i1si(iorb)
       nmut=i1mu(iorb)

       call exocc (iorb,ocup,ocdown)

       call prod2 (isiorb,psi(iborb),psi(iborb),wk1)
       call dscal (isiorb,ocup,wk1,ione)

       call prod2 (isiorb,psi(iborb),psi(iborb),wk2)
       call dscal (isiorb,ocdown,wk2,ione)

       !        store total spin densities
       call add(isiorb,wk1,rhotup)
       call add(isiorb,wk2,rhotdown)
    enddo

    ! calculate (nabla rho nabla rho)
    ! FIXME
    call nfng (rhotup,rhotup,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7)
    call dcopy(mxsize,wk7,ione,grhotup,ione)

    call nfng (rhotdown,rhotdown,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7)
    call dcopy(mxsize,wk7,ione,grhotdown,ione)

    do i=1,mxsize
       rhotup(i)   =two*rhotup(i)
       rhotdown(i) =two*rhotdown(i)
       grhotup(i)  =four*grhotup(i)
       grhotdown(i)=four*grhotdown(i)
    enddo

    ! total exchange energy is calculated as
    ! Ex(rhoup,rhodown)=(1/2)Ex(2*rhoup)+(1/2)Ex(2*rhodown)

    expw91=fdften(alphaf)*(  half*expw91sup(wgt2,rhotup,grhotup,wk0,wk1) &
         + half*expw91sup(wgt2,rhotdown,grhotdown,wk0,wk1))

  end function expw91

  !  ### expw86sup ###
  !
  !     Calculates exchange energy according to a formula of Parr and Wang Yue
  !     Phys. Rev. B 54, 16 533 (1996), Phys. Rev. B 45, 13 244 (1992).
  !
  !     J. P. Perdew, J. A. Chevary, S. H. Vosko, K. A. Jackson, M. R. Pederson
  !     and C. Fiolhais, Phys. Rev. B 46, 6671 (1992)
  !
  function expw91sup (wgt2,rho,grho,wk0,wk1)
    use params
    use discrete
    use commons
    use blas
    use utils

    implicit none
    integer (KIND=IPREC) :: i
    real (PREC) :: expw91sup
    real (PREC) :: a,a1,a2,a3,a4,ash,b1,const,const23,const43,rho43,s,s2

    real (PREC), dimension(*) :: wgt2,rho,grho,wk0,wk1

    parameter(a=7.795600_PREC,a1=0.1964500_PREC,a2=0.274300_PREC,a3=0.1508400_PREC, &
         a4=100.00_PREC,b1=0.00400_PREC,const23=2.0_PREC/3.0_PREC,const43=4.0_PREC/3.0_PREC)
#ifdef BLAS    
    real (PREC) ddot
    external ddot
#endif

    ! arcsinh
    ash(s)=log(s+sqrt(one+s*s))

    const=(24.00_PREC*pii*pii)**(-const23)

    ! grho = nabla rho  nabla rho
    ! |nabla rho| = sqrt(grho)
    do i=1,mxsize
       if (abs(rho(i)).lt.precis) then
          wk0(i)=0.0_PREC
       else
          rho43=rho(i)**const43
          s2=const*grho(i)/(rho43*rho43)
          s=sqrt(s2)
          wk0(i)=rho43*(one+a1*s*ash(a*s)+s2*(a2-a3*exp(-a4*s2)))/(one+a1*s*ash(a*s)+b1*s2*s2)
       endif
    enddo

    ! take care of F4 factor
    call multf4(wk0)

    expw91sup=ddot(mxsize,wgt2,ione,wk0,ione)

  end function expw91sup
  
end module dftexc
