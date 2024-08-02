! SPDX-License-Identifier: GPL-2.0-or-later

! Copyright (C) 2010-2024  Jacek Kobus 

module dftvxc
  implicit none
contains

  ! ### fdften ###
  !
  !     Calculates coefficients appearing in LDA energy formulea.
  !
  function fdften (alpha)
    use params
    use commons

    implicit none

    real (PREC) :: fdften
    real (PREC) :: alpha,const13,const34,const98
    parameter (const13=1.0_PREC/3.0_PREC,const34=3.0_PREC/4.0_PREC, &
         const98=9.0_PREC/8.0_PREC)

    if (idftex.eq.1.or.idftex.eq.2) then
       ! L.Laaksonen, P.Pyykko, D.Sundholm (Int. J.Quantum Chem. 27 (1985) 601)
       ! additional factor 2**(1/3) is due to the way the total density is calculated
       ! (see fldapot)

       ! Becke (idftex=2)
       fdften = -alpha*const98*(three/pii)**const13*two**(const13)
    elseif (idftex.eq.3) then
       ! generalized gradient approximation (GGA) (Perdew & Wang 86 )
       fdften = -const34*(three/pii)**const13

    elseif (idftex.eq.4) then
       ! generalized gradient approximation (GGA) (Perdew & Wang 91 )
       fdften = -const34*(three/pii)**const13
    else
       stop "Invalid argument idftex to fdften"
    endif
  end function fdften

  ! ### fdftpot ###
  !
  !     Calculates the coefficients appearing in the LDA potential formulea.
  !
  function fdftpot (alpha)
    use params
    use commons

    implicit none
    real (PREC) :: fdftpot
    real (PREC) :: alpha,const13,const23,const34

    parameter (const13=1.0_PREC/3.0_PREC,const23=2.0_PREC/3.0_PREC,const34=3.0_PREC/4.0_PREC)

    if (idftex.eq.1) then
       ! fldapot=-alphaf*three/two*(three/pii)**const13

       ! In case of the spin local density approximation where
       !   \rho_{total}^{1/3}=\sum_{spinorbital i}
       !                [rho_{\alpha}^{1/3}(i) + \rho_{\beta}^{1/3}(i)]
       ! an additional factor 2**(-2/3) shows up
       !fdftpot=-alpha*three/two*(three/pii)**const13/two**(const23)
       fdftpot=-alpha*three/two*(two*three/pii)**const13

       !fdftpot=-two**(-const23)*nine/eight*alpha*three/two*(three/pii)**const13
       ! according to ACES
       !fdftpot=-alpha*three/two*(three/pii)**const13*two**(const13)

    elseif (idftex.eq.2) then
       ! D.Becke Phys. Rev. A 38 (1988) 3098
       fdftpot=-alpha*three/two*(three/pii)**const13/two**(const23)

    elseif (idftex.eq.3) then
       ! generalized gradient approximation (GGA)
       ! Perdew & Wang Phys. Rev. B 33 (1986) 8800
       fdftpot=-const34*(three/pii)**const13
    else
       stop "Invalid argument idftex to fdftpot."
    endif
  end function fdftpot

  ! ### fbe88 ###
  !
  !     Calculates the Becke gradient-corrected exchange potential (PRA 88
  !     (1988) 3098-3100)
  !
  subroutine fbe88 (psi,f4,rhot,rhotup,rhotdown,grhot,grhotup,grhotdown, &
       wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7)
    use params
    use discrete
    use commons
    use blas
    use elocc
    use utils
    
    implicit none
    integer (KIND=IPREC) :: i,iborb,iorb,isiorb,nmut
    real (PREC) :: ocdown,ocup
    real (PREC), dimension(*) :: psi,f4,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7,rhot,rhotup,rhotdown, &
         grhot,grhotup,grhotdown

    do i=1,mxsize
       rhotup(i)   =0.0_PREC
       rhotdown(i) =0.0_PREC
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

       !        store total densities
       call add(isiorb,wk1,rhotup)
       call add(isiorb,wk2,rhotdown)
    enddo

    call fbe88sup (rhotup,grhotup,rhot,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7)

    call fbe88sup (rhotdown,grhotdown,rhot,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7)

    call dcopy(mxsize,grhotup,ione,wk2,ione)
    call add (mxsize,grhotdown,wk2)

    call prod (mxsize,f4,wk2)

  end subroutine fbe88

  ! ### fbe88sup ###
  !
  !     Calculates the Becke exchange potential for a given density and
  !     returns it in grhot array (see Johnson, Gill, Pople, JCP 98 (1993)
  !     p.5623) and L.Fan and T.Ziegler, JCP 94 (1991) 6057)
  !
  subroutine fbe88sup (rhot,grhot,wk,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7)
    use params
    use discrete
    use commons

    use blas
    use nabla

    implicit none
    integer (KIND=IPREC) :: i
    real (PREC) :: ash,bbeta,const13,const43,const53, &
         g1,g3,f,fm1,s,s2,t1,t2,t3

    real (PREC), dimension(*) :: wk,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7,rhot,grhot
    parameter (const13=1.0_PREC/3.0_PREC,const43=4.0_PREC/3.0_PREC, &
         const53=5.0_PREC/3.0_PREC,bbeta=0.0042_PREC)

    ! arcsinh
    ash(s)=log(s+sqrt(one+s*s))

    ! Calculate |(nabla rhot nabla rhot)| (grhot)
    call nfng (rhot,rhot,wk0,wk1,wk2,wk3,wk4,wk5,wk6,grhot)

    ! calculate gamma (see JGP, p.5623)
    do i=1,mxsize
       if (abs(rhot(i)).lt.precis) then
          wk(i)=0.0_PREC
       else
          wk(i)=sqrt(grhot(i))/rhot(i)**const43
       endif
    enddo

    ! nabla^2 rho
    call n2f (rhot,wk0,wk1,wk2,grhot)

    ! nabla gamma nabla rho (wk7)
    call nfng (wk,rhot,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7)

    ! calculate g and g prime
    do i=1,mxsize
       g1=fdftpot(alphaf)*rhot(i)**const13

       s=wk(i)
       s2=s*s

       if (abs(rhot(i)).lt.precis) then
          wk1(i)=g1
       else
          fm1=(one+6.0_PREC*bbeta*s*ash(s))
          if (fm1.lt.precis) then
             f=0.0_PREC
          else
             f=one/fm1
          endif

          g3=s/sqrt(one+s2)

          t1=const43*s2*rhot(i)**const53
          !           nabla^2 rho
          t2=one+f*(one-6.0_PREC*bbeta*s*g3)

          !           nabla rho nabla gamma
          t3=6.0_PREC*bbeta*f*((one+two*f)*ash(s)+g3*(one/(one+s2)+two*f*(two-6.0_PREC*bbeta*s*g3)))

          wk1(i)=g1-bbeta*f/rhot(i)**const43*(t1-grhot(i)*t2+wk7(i)*t3)

       endif
    enddo

    call dcopy(mxsize,wk1,ione,grhot,ione)

  end subroutine fbe88sup

  ! ### flypcs ###
  !
  !     Calculates (and returns in wk2) the correlation potential of Lee,
  !     Yang and Parr (PRB 37 (1988) 785, eq.23, i.e. using the
  !     closed-shell formula.
  !
  subroutine flypcs (psi,f4,rhot,rhotup,rhotdown,grhot,grhotup,grhotdown, &
       wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7)
    use params
    use discrete
    use commons
    use blas
    use elocc
    use nabla
    use utils
    
    implicit none
    integer (KIND=IPREC) :: iborb,iorb,isiorb,nmut
    real (PREC) :: ocdown,ocup
    real (PREC), dimension(*) :: psi,f4,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7, &
         rhot,rhotup,rhotdown,grhot,grhotup,grhotdown

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

       if (abs(ocup-ocdown)>epsilon(zero)) then
          write(*,*) "Warning: This implementation of LYP potential is only valid for closed shell systems."
          !stop 'flypcs'
       endif

       call prod2 (isiorb,psi(iborb),psi(iborb),wk1)
       call dscal (isiorb,ocup,wk1,ione)

       call prod2 (isiorb,psi(iborb),psi(iborb),wk2)
       call dscal (isiorb,ocdown,wk2,ione)

       ! store total densities
       call add(isiorb,wk1,rhotup)
       call add(isiorb,wk2,rhotdown)
    enddo

    call dcopy(mxsize,rhotup,ione,rhot,ione)
    call add(mxsize,rhotdown,rhot)

    call flypsupcs (rhot,rhotup,rhotdown,grhot,grhotup,grhotdown,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7)

    call prod (mxsize,f4,wk7)
    call dcopy(mxsize,wk7,ione,wk2,ione)

  end subroutine flypcs

  ! ### flypsupcs ###
  !
  !     Calculates the correlation potential of Lee, Yang and Parr (PRB 37
  !     (1988) 785, eq.23, i.e. using the closed-shell formula); it is
  !     returned in wk10
  !
  subroutine flypsupcs (rhot,rhotup,rhotdown,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7,wk8,wk9,wk10)
    use params
    use discrete
    use commons

    use blas
    use nabla

    implicit none
    integer (KIND=IPREC) :: i
    real (PREC) :: a,b,c,cf,d,const13,const23,const53,const83, &
         f1,g0,g1,g2,g3,rh
    real (PREC), dimension(*) :: wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7,wk8,wk9,wk10, &
         rhot,rhotup,rhotdown

    parameter (const13=1.0_PREC/3.0_PREC,const23=2.0_PREC/3.0_PREC,const53=5.0_PREC/3.0_PREC, &
         const83=8.0_PREC/3.0_PREC)

    ! coefficients of the Colle-Salveti formula
    parameter(a=0.049180_PREC,b=0.1320_PREC,c=0.25330_PREC,d=0.3490_PREC)

    !call dcopy (mxsize,wk6,ione,rhotup,ione)

    ! total density in rhot
    ! nabla^2 rho (rhotup)
    call n2f(rhot,wk0,wk1,wk2,rhotup)

    ! nabla rho nabla rho (rhotdown)
    call nfng (rhot,rhot,rhotdown,wk3,wk4,wk5,wk6,wk7,wk8,wk9)
    call dcopy (mxsize,wk9,ione,rhotdown,ione)

    cf=0.30_PREC*(three*pii*pii)**const23

    do i=1,mxsize
       if (abs(rhot(i)).lt.precis) then
          wk10(i)=0.0_PREC
       else
          rh=rhot(i)**(-const13)
          f1=one/(one+d*rh)
          g0=rhot(i)**(-const53)*exp(-c*rh)
          g1=-const53+c/three*rhot(i)**(-const13)
          g2=g1+d/three*f1*rhot(i)**(-const13)
          g3=d*d*f1*f1*rhot(i)**(-const13) -d*f1-c

          wk10(i)=                                                     &

                                ! See documentation for the derivation of these terms
                                !        term (1)
               -a*d/three*rh*f1*f1-a*f1                                &

                                !        term (2)
               -a*b*cf*f1*g0*rhot(i)**(const53)*(g2+const83)           &

                                !        term (3)
               -a*b*19.0_PREC/18.0_PREC*f1*g0*rhotup(i)                        &

                                !        term (4)
               -a*b*1.0_PREC/72.0_PREC*f1*g0*g2*(42.0_PREC*rhotup(i)+59.0_PREC*rhotdown(i)/rhot(i)) &

                                !        term (5)
               -a*b*7.0_PREC/24.0_PREC*f1*g0/rhot(i)*rhotdown(i)*( g2*(d/three*f1*rh+g1-1.0_PREC)+g3*rh/9.0_PREC )

       endif
    enddo

    return
  end subroutine flypsupcs

  ! ### fvwncs ###
  !
  !     Calculates (and returns in wk2) the correlation potential of VWN
  !     using the closed-shell formula.
  !
  subroutine fvwncs (psi,f4,rhot,rhotup,rhotdown,grhot,grhotup,grhotdown, &
       wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7)
    use params
    use discrete
    use commons
    use utils
    use elocc
    use blas

    implicit none
    integer (KIND=IPREC) :: i,iborb,iorb,isiorb,nmut
    real (PREC) :: ocdown,ocup
    real (PREC), dimension(*) :: psi,f4,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7, &
         rhot,rhotup,rhotdown,grhot,grhotup,grhotdown

    do i=1,mxsize
       rhotup(i)   =0.0_PREC
       rhotdown(i) =0.0_PREC
    enddo

    !     calculate total densities due to up and down spin
    do iorb=1,norb
       !if (inhyd(iorb).eq.1) cycle
       if (inDFT(iorb).eq.0) cycle
       iborb=i1b(iorb)
       isiorb=i1si(iorb)
       nmut=i1mu(iorb)

       call exocc (iorb,ocup,ocdown)

       if (abs(ocup-ocdown)>epsilon(zero)) then
          write(*,*) "Warning: This implementation of VWN potential is only valid for closed shell systems."
          stop 'fvwncs'
       endif

       call prod2 (isiorb,psi(iborb),psi(iborb),wk1)
       call dscal (isiorb,ocup,wk1,ione)

       call prod2 (isiorb,psi(iborb),psi(iborb),wk2)
       call dscal (isiorb,ocdown,wk2,ione)

       !        store total densities
       call add(isiorb,wk1,rhotup)
       call add(isiorb,wk2,rhotdown)
    enddo

    call dcopy(mxsize,rhotup,ione,rhot,ione)
    call add(mxsize,rhotdown,rhot)

    call fvwnsupcs (rhot,rhotup,rhotdown,grhot,grhotup,grhotdown,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7)

    call prod (mxsize,f4,wk7)
    call dcopy(mxsize,wk7,ione,wk2,ione)

  end subroutine fvwncs

  ! ### fvwnsupcs ###
  !
  !     Calculates the correlation potential of VWN using the closed-shell
  !     formula; it is returned in wk10 (see also fvwncs)
  !
  subroutine fvwnsupcs (rhot,rhotup,rhotdown,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7,wk8,wk9,wk10)
    use params
    use discrete
    use commons
    use nabla

    implicit none
    integer (KIND=IPREC) :: i
    real (PREC) :: ck1,cl1,cm1,cn1,const16,const76,const56,constx,constxp,g1,g2,x,xderrhot
    real (PREC), dimension(*) :: wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7,wk8,wk9,wk10,rhot,rhotup,rhotdown

    parameter(ck1=0.03109070_PREC,cl1=-0.104980_PREC,cm1=3.727440_PREC,cn1=12.93520_PREC, &
         const16=1.0_PREC/6.0_PREC,const76=7.0_PREC/6.0_PREC,const56=5.0_PREC/6.0_PREC)

    !     total density in rhot

    constx=(three/(four*pii))**const16
    constxp=-(four*pii/three)**const56/(eight*pii)

    ! dx/d(rho)=xderrhot
    do i=1,mxsize
       if (abs(rhot(i)).lt.precis) then
          wk10(i)=0.0_PREC
       else
          x=constx*rhot(i)**(-const16)
          xderrhot=constxp*rhot(i)**(-const76)
          g1=qvwn(x,ck1,cl1,cm1,cn1)
          g2=qvwnderx(x,ck1,cl1,cm1,cn1)
          !            print *,'2',i,constx,constxp,x,xderrhot,g1,g2
          wk10(i)=g1+rhot(i)*xderrhot*g2
       endif
    enddo
  end subroutine fvwnsupcs

  ! ### qvwn ###
  !
  !     Calculates q function of VWN functional (see
  !     http://wild.life.nctu.edu.tw/~jsyu/molpro2002.1/doc/manual/node184.html)
  !
  function qvwn(x,a,p,c,d)
    use params
    use commons
    
    implicit none
    real (PREC) :: qvwn
    real (PREC) :: a,arctan,c,d,p,qcap,t1,t2,t3,xcap,x

    xcap=x*x+c*x+d
    qcap=sqrt(four*d-c*c)
    arctan=atan(qcap/(two*x+c))

    t1=log(x*x/xcap)
    t2=two*c/qcap*arctan
    t3=-c*p/xcap*(log((x-p)*(x-p)/xcap) +two*(c+two*p)/qcap*arctan)
    
    ! modified formula give better agreement with Libxc lda_c_vwn functional
    ! if used the derivative formula should be modified accordingly
    ! see libxc/maple/vwn.mpl 
    ! t3=-c*p/xcap*(log((x-p)*(x-p)/xcap))

    qvwn=a*(t1+t2+t3)
    return
  end function qvwn

  ! ### qvwnderx ###
  !
  !     Calculates dq/dx of the VWN functional
  !
  function qvwnderx(x,a,p,c,d)
    use params
    use commons

    implicit none
    real (PREC) qvwnderx
    real (PREC) :: a,c,d,p,qcap,qcap2,t1,t2,t3,t4,x,xcap,xcapder,xcapder2,xcapp

    ! X=xcap, dX/dx=xcapder
    xcap=x*x+c*x+d
    xcapder=two*x+c
    xcapder2=xcapder*xcapder

    xcapp=p*p+c*p+d

    qcap=sqrt(four*d-c*c)
    qcap2=(four*d-c*c)

    t1=two/x-xcapder/xcap
    t2=-four*c/(qcap2+xcapder2)
    t3=four*c*p*(c+two*p)/xcapp/(qcap2+xcapder2)
    t4=-c*p/xcapp*(two/(x-p)-xcapder/xcap)

    qvwnderx=a*(t1+t2+t3+t4)

  end function qvwnderx

  ! ### ffbar ###
  !
  !     Calculates Fbar(s) (see fpw86sup)
  !
  function ffbar (s)
    use params
    implicit none
    real (PREC) :: ffbar
    real (PREC) ::  a,b,c,s,s2
    parameter (a=1.2960_PREC,b=14.0_PREC,c=0.20_PREC)

    s2=s*s
    ffbar=(1.0_PREC+a*s2+b*s2*s2+c*s2*s2*s2)

  end function ffbar

  ! ### ffbarp ###
  !
  !     Calculates Fbar(s) prime (see fpw86sup)
  !
  function ffbarp (s)
    use params
    implicit none
    real (PREC) :: ffbarp
    real (PREC) :: a,b,c,s,s2
    parameter (a=1.2960_PREC,b=14.0_PREC,c=0.20_PREC)

    s2=s*s
    ffbarp=(2.0_PREC*a*s+4.0_PREC*b*s*s2+6.0_PREC*c*s*s2*s2)

  end function ffbarp

  ! ### ff ###
  !
  !     Calculates F(s) (see fpw86sup)
  !
  function ff (s)
    use params

    implicit none
    real (PREC) :: ff
    real (PREC) ::  const115,s

    parameter (const115=1.0_PREC/15.0_PREC)

    ff=ffbar(s)**const115

  end function ff

  ! ### ffp ###
  !
  !     Calculates dF(s)/ds (see fpw86sup)
  !
  function ffp (s)
    use params
    implicit none
    real (PREC) :: ffp
    real (PREC) :: const115,s,s2
    parameter (const115=1.0_PREC/15.0_PREC)

    s2=s*s
    ffp=const115*ffbar(s)**(const115-1.0_PREC)*ffbarp(s)
  end function ffp


  ! ### ffdp ###
  !
  !     Calculates d(1/s dF(s)/ds)/ds (see fpw86sup)
  !
  function ffdp (s)
    use params
    implicit none
    real (PREC) :: ffdp
    real (PREC) :: b,c,const115,s,s2

    parameter (b=14.0_PREC,c=0.20_PREC,const115=1.0_PREC/15.0_PREC)

    s2=s*s
    ffdp= const115*(8.0_PREC*b*s+24.0_PREC*c*s*s2)*ffbar(s)**(const115-1.0_PREC) &
         +const115*(const115-1.0_PREC)*ffbarp(s)**2.0_PREC/s*ffbar(s)**(const115-2.0_PREC)

  end function ffdp
end module dftvxc
