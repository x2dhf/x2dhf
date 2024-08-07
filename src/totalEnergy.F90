! SPDX-License-Identifier: GPL-2.0-or-later

! Copyright (C) 1996       Leif Laaksonen, Dage Sundholm               
! Copyright (C) 1996-2023  Jacek Kobus 

module totalEnergy
  implicit none
contains
  ! ### etotal ###
  !
  !    Calculates total HF energy
  !
  subroutine etotalHF 
    use params
    use discrete
    use commons
    use utils
    use blas
    use exchContribs
    use diskInterface
    use inout
    use integrals
    use sharedMemory
    use utils
    implicit none
    integer (KIND=IPREC) :: ibex,iborb,iborb1,iborb2,ibpot,ibpot1,ibpot2,iex,iex1,iorb,iorb1,iorb2,isiex,isiex1,isiorb,&
         isiorb1,isiorb2,isipot,isipot1,isipot2,isym,ngrid,nmut,nmut1,nmut2

    real (PREC) :: epscharge,etsum,oc,ocx2,ocx1,w,wcouldft,wdcoul,wex1,wex2,wndc,woneel
    !real (PREC), dimension(*) :: psi,pot,excp,e,f0,wgt1,wgt2,wk0,wk1,wk2,wk3,wk4,wk5,wk6, &
    !wk7,wk8,wk9,wk10,wk11,wk12,wk13

    real (PREC), dimension(:), pointer :: psi,excp,e,f0,wgt1,wgt2,&
              wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7,wk8,wk9,wk10,wk11,wk12,wk13
    e=>supplptr(i4b(4):)
    excp=>exchptr
    f0=>supplptr(i4b(5):)
    psi=>orbptr
    wgt1=>supplptr(i4b(13):)
    wgt2=>supplptr(i4b(14):)

    wk0 =>scratchptr(          1:   mxsize8)
    wk1 =>scratchptr(   mxsize8+1: 2*mxsize8)
    wk2 =>scratchptr( 2*mxsize8+1: 3*mxsize8)
    wk3 =>scratchptr( 3*mxsize8+1: 4*mxsize8)
    wk4 =>scratchptr( 4*mxsize8+1: 5*mxsize8)            
    wk5 =>scratchptr( 5*mxsize8+1: 6*mxsize8)
    wk6 =>scratchptr( 6*mxsize8+1: 7*mxsize8)
    wk7 =>scratchptr( 7*mxsize8+1: 8*mxsize8)            
    wk8 =>scratchptr( 8*mxsize8+1: 9*mxsize8)
    wk9 =>scratchptr( 9*mxsize8+1:10*mxsize8)
    wk10=>scratchptr(10*mxsize8+1:11*mxsize8)
    wk11=>scratchptr(11*mxsize8+1:12*mxsize8)
    wk12=>scratchptr(12*mxsize8+1:13*mxsize8)
    wk13=>scratchptr(13*mxsize8+1:14*mxsize8)


    data epscharge /1.e-7_PREC/

    engt(1)=ee(ione,ione)

    ! FIXME
    engt(1)=zero

    ! calculate first contributions from one particle operators and
    ! coulomb potential contributions within the same shell

    vkt=zero
    vnt=zero
    woneel=zero
    wdcoul=zero
    wcouldft=zero
    etsum=zero

#ifdef PRINT
! print= 66: etotalHF: HF total energy individual contributions and their sum
    if (iprint(66).ne.0) then
       write(*,*)
       write(*,'("-------------------------------------------------------------------------"$)')
       write(*,'("-----------------------------------------")')
       write(*,'("HF total energy: individual contributions and their sum: ")')
    endif
#endif
    
    do iorb=1,norb
       woneel=woneel+occ(iorb)*oneelii(iorb)
       vkt=vkt+occ(iorb)*vk(iorb)
       vnt=vnt+occ(iorb)*vn(iorb)
       oc=occ(iorb)

#ifdef PRINT
       if (iprint(66).ne.0) then
          etsum=etsum+oc*vk(iorb)
          write(*,7028) iorn(iorb),bond(iorb),gusym(iorb),iorn(iorb),bond(iorb),gusym(iorb),vk(iorb),oc,etsum
7028      format('<',i4,1x,a5,a1,'| T |',i4,1x,a5,a1,' >',26x,1Pd25.16, 0Pf8.2, 1Pd25.16)

          etsum=etsum+oc*vn(iorb)
          write(*,7030) iorn(iorb),bond(iorb),gusym(iorb),iorn(iorb),bond(iorb),gusym(iorb),vn(iorb),oc,etsum
7030      format('<',i4,1x,a5,a1,'| V |',i4,1x,a5,a1,' >',26x,1Pd25.16, 0Pf8.2, 1Pd25.16)
       endif
#endif
       
       if (occ(iorb).gt.one) then
          w=coulij (iorb,iorb,psi,excp,wgt2,wk0)
          wdcoul=wdcoul+occ(iorb)*(occ(iorb)-one)/two*w
          ! To calculate DFT exchange energy the Coulomb energy must include
          ! J_{ii}=K_{ii} terms

          wcouldft=wcouldft+occ(iorb)*occ(iorb)/two*w

#ifdef PRINT          
          if (iprint(66).ne.0) then
             etsum=etsum+oc*(oc-one)/two*w
             write(*,7031) iorn(iorb),bond(iorb),gusym(iorb),iorn(iorb),bond(iorb),gusym(iorb), &
                  iorn(iorb),bond(iorb),gusym(iorb),iorn(iorb),bond(iorb),gusym(iorb),w,oc*(oc-one)/two,etsum
7031         format('<',i4,1x,a5,a1,'| J1(',i4,1x,a5,a1,i4,1x,a5,a1,' ) |',i4,1x,a5,a1,' >',1Pd25.16, 0Pf8.2, 1Pd25.16)
          endif
#endif          
       endif
    enddo

    ! contribution from Coulomb and exchange interaction between shells
    wndc=zero
    wex1=zero
    wex2=zero

   
    do iorb1=1,norb
       ! calculate exchange interaction within pi, delta, etc. shell
       if (mm(iorb1).gt.0.and.abs(occ(iorb1)-one).gt.epscharge) then
          call exint (iorb1,ocx1)
          w=exchij(0,iorb1,iorb1,psi,excp,wgt2,wk0)
          wex1=wex1+ocx1*w

#ifdef PRINT                    
          if (iprint(66).ne.0) then
             etsum=etsum-ocx1*w
             write(*,7032) iorn(iorb1),bond(iorb1),gusym(iorb1),iorn(iorb1),bond(iorb1),gusym(iorb1), &
                  iorn(iorb1),bond(iorb1),gusym(iorb1),iorn(iorb1),bond(iorb1),gusym(iorb1),w,-ocx1,etsum
7032         format('<',i4,1x,a5,a1,'| K (',i4,1x,a5,a1,i4,1x,a5,a1,' ) |',i4,1x,a5,a1,' >',1Pd25.16, 0Pf8.2, 1Pd25.16)
          endif
#endif
       endif

       do iorb2=iorb1+1,norb
          iex=iorb1+iorb2*(iorb2-1)/2
          
          ! Coulomb interaction between shells
          w=coulij(iorb2,iorb1,psi,excp,wgt2,wk0)
          wndc=wndc+occ(iorb1)*occ(iorb2)*w
          wcouldft=wcouldft+occ(iorb1)*occ(iorb2)*w

#ifdef PRINT                    
          if (iprint(66).ne.0) then
             etsum=etsum+occ(iorb1)*occ(iorb2)*w
             write(*,7034) iorn(iorb1),bond(iorb1),gusym(iorb1),iorn(iorb2),bond(iorb2),gusym(iorb2), &
                  iorn(iorb2),bond(iorb2),gusym(iorb2),iorn(iorb1),bond(iorb1),gusym(iorb1),&
                  w,occ(iorb1)*occ(iorb2),etsum
7034         format('<',i4,1x,a5,a1,'| J (',i4,1x,a5,a1,i4,1x,a5,a1,' ) |',i4,1x,a5,a1,' >',1Pd25.16, 0Pf8.2, 1Pd25.16)
          endif
#endif
          ! exchange interaction between shells (same lambda)
          call excont (iorb1,iorb2,ocx1,ocx2)
          w=exchij(0,iorb1,iorb2,psi,excp,wgt2,wk0)
          wex1=wex1+ocx1*w

          if (iprint(66).ne.0) then
             etsum=etsum-ocx1*w
             write(*,7036) iorn(iorb1),bond(iorb1),gusym(iorb1),iorn(iorb1),bond(iorb1),gusym(iorb1), &
                  iorn(iorb2),bond(iorb2),gusym(iorb2),iorn(iorb2),bond(iorb2),gusym(iorb2),w,-ocx1,etsum
7036         format('<',i4,1x,a5,a1,'| K (',i4,1x,a5,a1,i4,1x,a5,a1,' ) |',i4,1x,a5,a1,' >',1Pd25.16, 0Pf8.2, 1Pd25.16)
          endif

          ! exchange interaction between shells (different lambda)
          if (ilc(iex).gt.1) then
             w=exchij(1,iorb1,iorb2,psi,excp,wgt2,wk0)
             wex2=wex2+ocx2*w

#ifdef PRINT                       
             if (iprint(66).ne.0) then
                etsum=etsum-ocx2*w
                write(*,7038) iorn(iorb1),bond(iorb1),gusym(iorb1),iorn(iorb1),bond(iorb1),gusym(iorb1), &
                     iorn(iorb2),bond(iorb2),gusym(iorb2),iorn(iorb2),bond(iorb2),gusym(iorb2),w,-ocx2,etsum
7038            format('<',i4,1x,a5,a1,'| K1(',i4,1x,a5,a1,i4,1x,a5,a1,' ) |',i4,1x,a5,a1,' >',1Pd25.16, 0Pf8.2, 1Pd25.16)
             endif
#endif
          endif
       enddo
    enddo

    evt=woneel+wdcoul+wndc-wex1-wex2
    epott=vnt+wdcoul+wndc-wex1-wex2
    virrat=(epott+z1*z2/r)/vkt
    etot=evt+z1*z2/r

    enkin=vkt
    ennucel=vnt
    encoul=wdcoul+wndc
    enexch=-wex1-wex2
    entot=enkin+ennucel+encoul+enexch+z1*z2/r

    ! to be able to compare HF and DFT exchange energies these contributions have to be
    ! calculated as
    encouldft=wcouldft
    enexchdft=evt-woneel-wcouldft

    engt(1)=evt

    etsum=zero

#ifdef PRINT              
    if (iprint(66).ne.0) then
       write(*,*)
       !write(*,'(" total energy contributions and their sum: ")')
       etsum=etsum+vkt
       write(*,'("  kinetic:                  ",1Pd25.16,1Pd25.16)') vkt,etsum
       etsum=etsum+vnt
       write(*,'("  nuclear attraction:       ",1Pd25.16,1Pd25.16)') vnt,etsum
       write(*,'("  one-electron:             ",1Pd25.16,1Pd25.16)') woneel,etsum
       etsum=etsum+wdcoul+wndc
       write(*,'("  Coulomb:                  ",1Pd25.16,1Pd25.16)') wdcoul+wndc,etsum
       etsum=etsum-wex1-wex2
       write(*,'("  exchange:                 ",1Pd25.16,1Pd25.16)') -wdcoul-wex1-wex2,etsum
       write(*,'("  two-electron:             ",1Pd25.16,1Pd25.16)') wdcoul+wndc-wex1-wex2
       write(*,'("  nuclear repulsion:        ",1Pd25.16,1Pd25.16)') z1*z2/r,etsum+z1*z2/r
       write(*,'("  total energy:             ",25x,1Pd25.16)') entot
       write(*,'("-------------------------------------------------------------------------"$)')
       write(*,'("-----------------------------------------")')
       write(*,*)
    endif
#endif
  end subroutine etotalHF

  ! ### etotalDFT ###
  !
  !     Calculates total energy using several DFT functionals
  !
  subroutine etotalDFT 
    use params
    use discrete
    use commons
    use utils
    use blas
    use dftexc
    use integrals
    use nabla
    use inout
    use sharedMemory
    use utils
    implicit none
    integer (KIND=IPREC) :: i,iborb,ibpot,iorb,isiorb,isipot,isym,nmut

    real (PREC) :: oc,w,wcorr,wex,wcoul,woneel,etsum

    real (PREC), dimension(:), pointer :: psi,psi1,e,excp,excp1,f0,wgt1,wgt2,&
         rhot,rhotup,rhotdown,grhot,grhotup,grhotdown,&
         wk0,wk1,wk2,wk3,wk10,wk11,wk12,wk13

#ifdef BLAS    
    real (PREC) ddot
    external ddot
#endif
    
    e=>supplptr(i4b(4):)
    excp=>exchptr
    f0=>supplptr(i4b(5):)
    psi=>orbptr
    wgt1=>supplptr(i4b(13):)
    wgt2=>supplptr(i4b(14):)

    wk0      =>scratchptr(          1:   mxsize8)
    wk1      =>scratchptr(   mxsize8+1: 2*mxsize8)
    wk2      =>scratchptr( 2*mxsize8+1: 3*mxsize8)
    wk3      =>scratchptr( 3*mxsize8+1: 4*mxsize8)
    rhot     =>scratchptr( 4*mxsize8+1: 5*mxsize8)            
    rhotup   =>scratchptr( 5*mxsize8+1: 6*mxsize8)
    rhotdown =>scratchptr( 6*mxsize8+1: 7*mxsize8)
    grhot    =>scratchptr( 7*mxsize8+1: 8*mxsize8)            
    grhotup  =>scratchptr( 8*mxsize8+1: 9*mxsize8)
    grhotdown=>scratchptr( 9*mxsize8+1:10*mxsize8)
    wk10     =>scratchptr(10*mxsize8+1:11*mxsize8)
    wk11     =>scratchptr(11*mxsize8+1:12*mxsize8)
    wk12     =>scratchptr(12*mxsize8+1:13*mxsize8)
    wk13     =>scratchptr(13*mxsize8+1:14*mxsize8)

#ifdef DEBUG
! debug= 75: etotalDFT: test of n2f and nfnf routines   
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
#endif
    ! calculate first contributions from one particle operators and
    ! Coulomb potential contributions within the same shell

    vnt=zero
    vkt=zero
    woneel=zero

    do iorb=1,norb
       !if (inhyd(iorb).eq.1) goto 10
       if (inDFT(iorb).eq.0) cycle
       woneel=woneel+occ(iorb)*oneelii(iorb)
       vkt=vkt+occ(iorb)*vk(iorb)
       vnt=vnt+occ(iorb)*vn(iorb)
    enddo
    
    evt=woneel
    etot=evt+z1*z2/r
    enkin=vkt
    ennucel=vnt

    if (nel.eq.1) return    
    
    ! contribution from coulomb interaction within the same shell
    ! calculate the coulomb potential contribution from all orbitals
    ! (include 1/2 factor )

    call zeroArray(mxsize,wk2)

    wcorr=zero
    wcoul=zero
    wex=zero

    do iorb=1,norb
       if (inDFT(iorb).eq.0) cycle
       call daxpy (mxsize,occ(iorb)/two,exchptr(i2b(iorb):),ione,wk2,ione)
    enddo

    ! contribution from the Coulomb interaction (DFT style, i.e. plus <i|K(i,i)|i> term)
    do iorb=1,norb
       if (inDFT(iorb).eq.0) cycle
       !psi1=>psi(i1b(iorb):)
       iborb=i1b(iorb)
       psi1=>orbptr(iborb:)
       call prod2 (mxsize,psi1,psi1,wk0)
       call prod (mxsize,wk2,wk0)
       call dscal (mxsize,occ(iorb),wk0,ione)
       w=ddot(mxsize,wgt2,ione,wk0,ione)
       wcoul=wcoul+w
    enddo
    encoul=wcoul
    encouldft=wcoul
    
    ! DFT exchange energy corrections
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
    if (idftcorr.eq.1) wcorr=eclyptot(psi,wgt2,rhot,rhotup,rhotdown,grhot,grhotup,grhotdown,&
         wk0,wk1,wk2,wk3,wk10,wk11,wk12,wk13)
    if (idftcorr.eq.2) wcorr=ecvwntot(psi,wgt2,rhot,rhotup,rhotdown,grhot,grhotup,grhotdown,&
         wk0,wk1,wk2,wk3,wk10,wk11,wk12,wk13)
    evt=woneel+wcoul+wex+wcorr
    etot=evt+z1*z2/r
    virrat=(evt-vkt)/vkt

    enexch=wex
    enexchdft=wex+wcorr
    edftcorr=wcorr

    entot=enkin+ennucel+encoul+enexch+edftcorr+z1*z2/r

    etsum=zero
    
#ifdef PRINT
! print= 68: etotalDFT: DFT total energy contributions and their sum
    if (iprint(68).ne.0) then
       write(*,*)
       !write(*,'(" total DFT energy contributions and their sum: ")')
       etsum=etsum+vkt
       write(*,'("  kinetic:                  ",1Pd25.16,1Pd25.16)') vkt,etsum
       etsum=etsum+vnt
       write(*,'("  nuclear attraction:       ",1Pd25.16,1Pd25.16)') vnt,etsum
       write(*,'("  one-electron:             ",1Pd25.16,1Pd25.16)') woneel,etsum
       etsum=etsum+wcoul
       write(*,'("  Coulomb (DFT)             ",1Pd25.16,1Pd25.16)') wcoul,etsum
       etsum=etsum+wex
       write(*,'("  exchange (DFT)            ",1Pd25.16,1Pd25.16)') wex,etsum
       etsum=etsum+wcorr
       write(*,'("  correlation (DFT):        ",1Pd25.16,1Pd25.16)') wcorr,etsum
       write(*,'("  two-electron (DFT):       ",1Pd25.16,1Pd25.16)') wcoul+wex+wcorr
       write(*,'("  nuclear repulsion:        ",1Pd25.16,1Pd25.16)') z1*z2/r,etsum+z1*z2/r
       write(*,'("  total energy:             ",25x,1Pd25.16)') entot
       write(*,'("-------------------------------------------------------------------------"$)')
       write(*,'("-----------------------------------------")')
       write(*,*)
    endif
#endif
  end subroutine etotalDFT
  
  ! ### etotalGauss ###
  !
  !     Evaluate the interaction energy between two nuclei with Gaussian
  !     charge distributions Z1*Z2/r lowerGamma(half,eta12*r*r)/sqrt(pi),
  !     1/eta12=1/eta1+1/eta2
  !
  subroutine etotalGauss
    use params
    use discrete
    use commons

    implicit none

    real (PREC) :: at3,atw,eta1,eta12,eta2,fmtoau,gammaLower,rrms,rrmsfm
    real (PREC), external :: dgamit,dgamma

    fmtoau=1.0e-13_PREC/ainfcm

    if (abs(z1)<epsilon(zero)) then
       eta1=1.e35_PREC
    endif

    if (abs(z2)<epsilon(zero)) then
       eta2=1.e35_PREC
    endif

    atw=z1atmass
    at3=atw**(one/three)
    rrmsfm = 0.8360_PREC*at3+0.5700_PREC
    rrms = rrmsfm*fmtoau
    eta1=3.0_PREC/(2.0_PREC*rrms*rrms)

    atw=z2atmass
    at3=atw**(one/three)
    rrmsfm = 0.8360_PREC*at3+0.5700_PREC
    rrms = rrmsfm*fmtoau
    eta2=3.0_PREC/(2.0_PREC*rrms*rrms)

    eta12=eta1*eta2/(eta1+eta2)
    ! lower incomplete gamma function
    gammaLower=(eta12*r*r)**(half)*dgamma(half)*dgamit(half,eta12*r*r)
    etotFN=etot-z1*z2/r+z1*z2*gammaLower/r/sqrt(pii)

  end subroutine etotalGauss


end module totalEnergy
