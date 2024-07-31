! SPDX-License-Identifier: GPL-2.0-or-later

! Copyright (C) 1996       Leif Laaksonen, Dage Sundholm               
! Copyright (C) 1996-2023  Jacek Kobus 
module printInitData
  implicit none
contains

  ! ### banner ###
  !
  !     Prints a banner of the program
  subroutine printBanner
    implicit none

    write(*,'(  "//////////////////////////// "$)')    
    write(*,'(" FINITE DIFFERENCE 2D HARTREE-FOCK "$)')
    write(*,'(" ////////////////////////////// ")')    
    write(*,'(  "//////////////////////////// "$)')    
    write(*,'("            version 3.0            "$)')
    write(*,'(" ////////////////////////////// ")')    

  end subroutine printBanner

  subroutine separator
    write(*,'("///////////////////////////////////////////////////////////////////////////////////////////////")')
  end subroutine separator

  ! ### fermi ###
  !
  !     Analyzes Fermi charge distribution.
  !
  subroutine fermi
    use params
    use discrete
    use commons

    implicit none

    integer (KIND=IPREC) :: iprt,mu,nc1,nc2,ni
    real (PREC) :: a,at3,atw,c,facto1,fmtoau,rr,rrms,rrmsfm,t,tfm,xr2,xrr

    iprt=0

#ifdef PRINT
! print=230: finite nuclei with the Fermi charge distribution
    if (iprint(230).ne.0) then
       write(*,*) ' '
       write(*,*) '  finite nuclei with the Fermi charge distribution'
    endif

    !     calculate parameters of the finite nucleus charge distribution
    !     (Fremi distribution)

    fmtoau=1.0e-13_PREC/ainfcm

    if (abs(z1)>epsilon(zero)) then

       ! set atomic weight for centre A

       atw=z1atmass
       at3=atw**(1.0_PREC/3.0_PREC)
       rrmsfm = 0.8360_PREC*at3+0.5700_PREC
       tfm = 2.30_PREC

       ! change units from fm into bohr

       rrms = rrmsfm*fmtoau
       t = tfm*fmtoau
       a = t/(4.00_PREC*log(3.00_PREC))
       facto1 = rrms**2-(7.00_PREC/5.00_PREC)*(pii**2)*(a**2)
       c = sqrt (5.00_PREC/3.00_PREC) * sqrt (facto1)

       ni=nni
       nc1=0
       xr2=r/2.00_PREC
       do mu=1,mxnmu
          rr=sqrt(vxisq(mu)+vetasq(ni)-1.0_PREC)
          xrr=xr2*rr
          if (abs(xrr-xr2).lt.c) then
             nc1=nc1+1
             if (iprint(231).ne.0) write(*,*) 'mu,xrr',mu,xrr
          endif
       enddo

       mu=1
       nc2=0
       do ni=nni,nni/2,-1
          rr=sqrt(vxisq(mu)+vetasq(ni)-1.0_PREC)
          xrr=xr2*rr
          if (abs(xrr-xr2).lt.c) then
             nc2=nc2+1
             if (iprint(231).ne.0) write(*,*) 'ni,xrr',ni,xrr
          endif
       enddo

       if (iprt.ne.0) then
          write(*,*) '    center A: '
          write(*,*) '    c = ',c,' bohr', '  a = ',a,' bohr'
          write(*,*) '    nu=-1:',nc1, ' points inside radius c'
          write(*,*) '    mu= 1:',nc2, ' points inside radius c'
          write(*,*) ' '
       endif
    elseif (abs(z2)>epsilon(zero)) then
       ! set atomic weight for centre B

       atw=z2atmass
       at3=atw**(1.0_PREC/3.0_PREC)
       rrmsfm = 0.8360_PREC*at3+0.5700_PREC
       tfm = 2.300_PREC

       ! change units from fm into bohr

       rrms = rrmsfm*fmtoau
       t = tfm*fmtoau
       a = t/(4.00_PREC*log (3.00_PREC))
       facto1 = rrms**2-(7.00_PREC/5.00_PREC)*(pii**2)*(a**2)
       c = sqrt (5.00_PREC/3.00_PREC) * sqrt (facto1)

       ni=1
       nc1=0
       xr2=r/2.00_PREC
       do mu=1,mxnmu
          rr=sqrt(vxisq(mu)+vetasq(ni)-1.0_PREC)
          xrr=xr2*rr
          if (abs(xrr-xr2).lt.c) then
             nc1=nc1+1
             if (iprint(231).ne.0) write(*,*) 'mu,xrr',mu,xrr
          endif
       enddo

       mu=1
       nc2=0
       do ni=nni,nni/2,-1
          rr=sqrt(vxisq(mu)+vetasq(ni)-1.0_PREC)
          xrr=xr2*rr
          if (abs(xrr-xr2).lt.c) then
             nc2=nc2+1
             if (iprint(231).ne.0) write(*,*) 'ni,xrr',ni,xrr
          endif
       enddo

       if (iprint(230).ne.0) then
          write(*,*) '     Center B: '
          write(*,*) '     c = ',c,' bohr', '  a = ',a,' bohr'
          write(*,*) '     nu=1:',nc1, ' points inside radius c'
          write(*,*) '     mu= 1:',nc2, ' points inside radius c'
       endif
    endif
#endif
  end subroutine fermi

  ! ### fockform ###
  !
  !     Print 2-electron part of the Fock equation for every orbital.
  !
  subroutine fockform
    use params
    use discrete
    use scfshr
    use commons

    implicit none
    integer (KIND=IPREC) :: iorb,iorb1,ipc,kex,orborder
    real (PREC) :: coo
    parameter (orborder=1)

    ! exchange contributions due to different pair of orbitals
    ! Ka:  sigma i - nonsigma j
    ! Kb:  nonsigma i - nonsigma i
    ! Kc:  nonsigma i - nonsigma j

    if (nel.eq.1) return

    if (norb.eq.1) then
       iorb=1
       write(*,'(/5x,"2-electron part of Fock operator for orbital",i4,1x,a8,a1)') &
            iorn(iorb),bond(iorb),gusym(iorb)
       write(*,'(15x,f6.2, "  J (",i4,1x,a8,a1,",",i4,1x,a8,a1,")" )') &
            1.0,iorn(iorb),bond(iorb),gusym(iorb),iorn(iorb),bond(iorb),gusym(iorb)
       return
    endif

    if (orborder.eq.1) goto 1000

    do iorb=1,norb
       write(*,'(/5x,"2-electron part of Fock operator for orbital",i4,1x,a8,a1)') &
            iorn(iorb),bond(iorb),gusym(iorb)

       ! add contributions from Coulomb potentials
       do iorb1=1,norb
          coo=occ(iorb1)
          if (iorb.eq.iorb1) coo=coo-1.0_PREC
          !           write(*,'(2i5," J1  " ,1Pe12.2)') iorb,iorb1,coo
          write(*,'(15x,f6.2, "  J (",i4,1x,a8,a1,",",i4,1x,a8,a1,")" )') &
               coo,iorn(iorb1),bond(iorb1),gusym(iorb1),iorn(iorb1),bond(iorb1),gusym(iorb1)
       enddo

       ! add contributions from exchange potentials
       do iorb1=1,norb
          kex=iorb+norb*(iorb1-1)
          if (iorb1.ne.iorb)  then
             coo=gec(kex)
             !              write(*,'(2i5," Ka  " ,1Pe12.2)') iorb,iorb1,-coo
             write(*,'(15x,f6.2, "  Ka(",i4,1x,a8,a1,",",i4,1x,a8,a1,")" )') &
                  -coo,iorn(iorb1),bond(iorb1),gusym(iorb1),iorn(iorb),bond(iorb),gusym(iorb)

             if (ilc(k2(iorb,iorb1))==2) then
                coo=gec(kex+norb*norb)
                !              call daxpy (ngexp,coo,excp(ibexp+ngexp),ione,wk2,ione)
                !                 write(*,'(2i5," Kc  " ,1Pe12.2)') iorb,iorb1,-coo
                write(*,'(15x,f6.2, "  Kc(",i4,1x,a8,a1,",",i4,1x,a8,a1,")" )') &
                     -coo,iorn(iorb1),bond(iorb1),gusym(iorb1),iorn(iorb),bond(iorb),gusym(iorb)

             endif
          else
             if ((mm(iorb).gt.0).and.(ilc(k2(iorb,iorb1)).gt.0)) then
                coo=gec(kex)
                write(*,'(15x,f6.2, "  Kb(",i4,1x,a8,a1,",",i4,1x,a8,a1,")" )') &
                     -coo,iorn(iorb1),bond(iorb1),gusym(iorb1),iorn(iorb),bond(iorb),gusym(iorb)
             endif
          endif
       enddo
    enddo

1000 continue

    do iorb=norb,1,-1
       write(*,'(/5x,"2-electron part of Fock operator for orbital",i4,1x,a8,a1)') &
            iorn(iorb),bond(iorb),gusym(iorb)

       ! add contributions from Coulomb potentials
       do iorb1=norb,1,-1
          kex=iorb+norb*(iorb1-1)
          coo=occ(iorb1)

          if (iorb.eq.iorb1) coo=coo-1.0_PREC
          write(*,'(15x,f6.2, "  J (",i4,1x,a8,a1,",",i4,1x,a8,a1,")" )') &
               coo,iorn(iorb1),bond(iorb1),gusym(iorb1),iorn(iorb1),bond(iorb1),gusym(iorb1)
       enddo

       ! add contributions from exchange potentials
       do iorb1=norb,1,-1
          kex=iorb+norb*(iorb1-1)
          ipc=iorb1+norb*(iorb-1)
          coo=occ(iorb1)

          if (iorb.eq.iorb1) coo=coo-1.0_PREC

          if (iorb1.ne.iorb)  then
             coo=gec(kex)
             write(*,'(15x,f6.2, "  Ka(",i4,1x,a8,a1,",",i4,1x,a8,a1,")" )') &
                  -coo,iorn(iorb1),bond(iorb1),gusym(iorb1),iorn(iorb),bond(iorb),gusym(iorb)

             if (ilc(k2(iorb,iorb1))==2) then
                coo=gec(kex+norb*norb)
                write(*,'(15x,f6.2, "  Kc(",i4,1x,a8,a1,",",i4,1x,a8,a1,")" )') &
                     -coo,iorn(iorb1),bond(iorb1),gusym(iorb1),iorn(iorb),bond(iorb),gusym(iorb)
             endif
          else
             if ((mm(iorb).gt.0).and.(ilc(k2(iorb,iorb1)).gt.0)) then
                coo=gec(kex)
                write(*,'(15x,f6.2, "  Kb(",i4,1x,a8,a1,",",i4,1x,a8,a1,")" )') &
                     -coo,iorn(iorb1),bond(iorb1),gusym(iorb1),iorn(iorb),bond(iorb),gusym(iorb)
             endif
          endif
       enddo
    enddo

    return
  end subroutine fockform

  ! ### printCase ###
  !
  subroutine printElConf
    use params
    use discrete
    use commons

    implicit none

    integer (KIND=IPREC) :: i,iorb,iorbr,iput,izz1,izz2,&
         nsigma,nsigmag,nsigmau,npi,ndelta,nphi,nspins

    izz1=nint(z1)
    izz2=nint(z2)

    write(iout6,*)
    write(iout6,*) '  Electronic configuration:'
    write(iout6,*)

    nsigma=0
    nsigmag=0
    nsigmau=0
    npi   =0
    ndelta=0
    nphi  =0

    do iorb=1,norb
       if (orbsym(iorb).eq.'sigma') nsigma=nsigma+1
       !        if (orbsym(iorb).eq.'sigma g') nsigmag=nsigmag+1
       !        if (orbsym(iorb).eq.'sigma u') nsigmau=nsigmau+1
       if (orbsym(iorb).eq.'pi'   ) npi   =npi   +1
       if (orbsym(iorb).eq.'delta') ndelta=ndelta+1
       if (orbsym(iorb).eq.'phi'  ) nphi  =nphi  +1
    enddo

    iorbr=0
    do iorb=1,norb
       if (iorbr.eq.0.and.orbsym(iorb).eq.'sigma') iorbr=nsigma
       !        if (iorbr.eq.0.and.orbsym(iorb).eq.'sigma g') iorbr=nsigmag
       !        if (iorbr.eq.0.and.orbsym(iorb).eq.'sigma u') iorbr=nsigmau
       if (iorbr.eq.0.and.orbsym(iorb).eq.'pi') iorbr=npi
       if (iorbr.eq.0.and.orbsym(iorb).eq.'delta') iorbr=ndelta
       if (iorbr.eq.0.and.orbsym(iorb).eq.'phi') iorbr=nphi
       nspins=4
       if (orbsym(iorb).eq.'sigma') nspins=2
       iput=4*(iorb-1)
       write(iout6,1020) iorbr, orbsym(iorb),gusym(iorb),(spin(i+iput),i=1,nspins)
       iorbr=iorbr-1
    enddo

    !     number of electrons and total charge

    write(iout6,*)
    write(iout6,1010) izz1+izz2-nel,nel,norb,norb,nexch
    !write(iout6,1010) izz1+izz2-nel,nel,norb,norb,nexch,nexchpots(1)

    if (linitFuncsHydrogen.or.initFuncsOED) then
       write(iout6,*)
       write(iout6,*) '  LCAO via hydrogenic functions:'
       write(iout6,*)
       write(iout6,1030)

       iorbr=0
       do iorb=1,norb
          if (iorbr.eq.0.and.orbsym(iorb).eq.'sigma') iorbr=nsigma
          if (iorbr.eq.0.and.orbsym(iorb).eq.'pi') iorbr=npi
          if (iorbr.eq.0.and.orbsym(iorb).eq.'delta') iorbr=ndelta
          if (iorbr.eq.0.and.orbsym(iorb).eq.'phi') iorbr=nphi

          iput=4*(iorb-1)
          write(iout6,1032) iorbr, orbsym(iorb),gusym(iorb),mgx(1,iorb),mgx(2,iorb),eza1(iorb),co1(iorb),&
               mgx(4,iorb),mgx(5,iorb),eza2(iorb),co2(iorb)
          iorbr=iorbr-1
       enddo
    endif

1010 format(/10x,'total charge            =',i3 &
         /10x,'number of',                    &
         /14x,'electrons           =',i3      &
         /14x,'orbitals            =',i3      &
         /14x,'Coulomb potentials  =',i3      &
         /14x,'exchange potentials =',i3)
         !/14x,'exchange potentials =',i3,' ('i2,')')

    !  1010 format(10x,'number of electrons =',i5
    !      &     /10x,'total charge        =',i5
    !      &     /10x,'number of orbitals, Coulomb and exchange potentials',
    !      &     3i4)

1020 format(10x,i2,2x,a8,a1,2x,4a2)
1030 format(10x,' orbital           n1 l1   Z1    c1','       n2 l2   Z2    c2 '/)
1032 format(10x,i2,2x,a8,a1,5x,2i3,2f6.2,5x,2i3,2f6.2)

  end subroutine printElConf
  
  ! ### printCase ###
  !
  !     This routine can be used to check if input data have been
  !     correctly transformed into values of numerous variable used by the
  !     program.
  !
  subroutine printCase
    use commons
    use discrete
    use memory
    use params
    use scfshr
    use solver

#ifdef LIBXC
    use, intrinsic :: iso_c_binding
    use xc_f90_lib_m
#endif
    
    implicit none
    integer (KIND=IPREC) :: i,ib,ie,iorb,iorb1,iorb2,ip,izz1,izz2,isizeint,isizereal
    integer (KIND=IPREC) :: maxorb2
    integer (KIND=8) :: lengtht,length0

    real (PREC) ::  heta,hxibeg,hxiend,omb,rinfig

#ifdef LIBXC
    type(xc_f90_func_info_t) :: xc_info
    type(xc_f90_func_t) :: xc_func
    integer(c_int) :: vmajor, vminor, vmicro, family_id, func_id, kind_id, number, err, refnumber
    character(len=120) :: name, kind1, family, ref
    character(len=14) :: nspin
    type(xc_f90_func_reference_t) :: reference 
#endif
    
    ! calculate the total length of working arrays (in bytes)
    ! system under consideration

    izz1=nint(z1)
    izz2=nint(z2)

    write(iout6,*)
    write(iout6,*) '  Atomic/molecular system: '
    write(iout6,*)
    write(iout6,1000) element(izz1), z1, element(izz2),z2 ,r,r*0.5291772490_PREC

    if (lpotFermi) then    
       write(iout6,1401) z1atmass,element(izz1),z2atmass,element(izz2)
    elseif (lpotGauss) then
       write(iout6,1402) z1atmass,element(izz1),z2atmass,element(izz2)
    endif

    write(iout6,*)
    if (HF)   write(iout6,'(3x,"Method: HF")')
    if (OED)  write(iout6,'(3x,"Method: OED")')
    if (TED)  write(iout6,'(3x,"Method: TED")')
    if (SCMC) write(iout6,'(3x,"Method: SCMC")')

    if (DFT.or.HFS) then
       write(iout6,'(3x,"Method: DFT     ")')
       if (.not.LXC) then
          if (idftex.eq.1) then
             write(iout6,'(20x,a4," functional (alpha = ",f7.5,")")') cdftex(idftex),alphaf
          elseif(idftex>1) then
             write(iout6,'(20x,a4," functional")') cdftex(idftex)
          endif
          if (idftcorr.ne.0) then
             write(iout6,'(20x,a4," functional")') cdftcorr(idftcorr)
          endif
       endif
    endif

#ifdef LIBXC
    If (LXC) then
       write(iout6,'(3x,"Method: DFT     ")')
       call xc_f90_version(vmajor, vminor, vmicro)
       do i=1,lxcFuncs
          if (lxcPolar) then
             nspin="XC_POLARIZED"
             call xc_f90_func_init(xc_func, lxcFuncs2use(i), XC_POLARIZED, err)
          else
             nspin="XC_UNPOLARIZED"
             call xc_f90_func_init(xc_func, lxcFuncs2use(i), XC_UNPOLARIZED, err)
          endif
          ! It turns out that for all the functionals used xc_f90_func_init returns err=0
          ! so the following piece of code is redundant
          ! if (err>0) then
          !    write(*,'("Error: xc_f90_func_init failed with err=",i5)') err
          !    stop "printCase"
          ! endif
          
          xc_info=xc_f90_func_get_info(xc_func)
          family_id=xc_f90_func_info_get_family(xc_info)
          kind_id=xc_f90_func_info_get_kind(xc_info)
          name=xc_f90_functional_get_name(lxcFuncs2use(i))
          
          select case(kind_id)
          case (XC_EXCHANGE)
             write(kind1, '(a)') 'XC_EXCHANGE'
          case (XC_CORRELATION)
             write(kind1, '(a)') 'XC_CORRELATION'
          case (XC_EXCHANGE_CORRELATION)
             write(kind1, '(a)') 'XC_EXCHANGE-CORRELATION'
          case (XC_KINETIC)
             write(kind1, '(a)') 'XC_KINETIC'
          case default
             write(kind1, '(a)') 'unknown'
          end select
          
          select case (family_id)
          case (XC_FAMILY_LDA);
             write(family,'(a)') "LDA"
          case (XC_FAMILY_HYB_LDA);
             write(family,'(a)') "Hybrid LDA"
             alphaflxc=xc_f90_hyb_exx_coef(xc_func)
             if (alphaflxc>zero) then
                write(*,"(3x,'           alpha = ',f7.5,&
                     ' (from xc_hyb_exx_coef)')") alphaflxc
                alphaf=alphaflxc
             else
                write(*,"(3x,'           alpha = ',f7.5)") alphaf
             endif
          case (XC_FAMILY_GGA);
             write(family,'(a)') "GGA"
             
          case (XC_FAMILY_HYB_GGA);
             write(family,'(a)') "Hybrid GGA"
             alphaflxc=xc_f90_hyb_exx_coef(xc_func)
             if (alphaflxc>zero) then
                write(*,"(3x,'           alpha = ',f7.5,&
                     ' (from xc_hyb_exx_coef)')") alphaflxc
                alphaf=alphaflxc
             else
                write(*,"(3x,'           alpha = ',f7.5)") alphaf
             endif
             
          case (XC_FAMILY_MGGA);
             write(family,'(a)') "MGGA"
          case (XC_FAMILY_HYB_MGGA);
             write(family,'(a)') "Hybrid MGGA"
          case default;
             write(family,'(a)') "unknown"
          end select
          
          write(*,'(/"      Libxc version = ",I1,".",I1,".",I1)') vmajor, vminor, vmicro             
          write(*,'( "               name = ", a)') trim(name)
          write(*,'( "               kind = ", a)') trim(kind1)
          write(*,'( "             family = ", a)') trim(family)
          write(*,'( "              nspin = ", a)') nspin
          write(*,'( "       reference(s) = ",$)') 
          
          number=0
          refnumber=0
          do while (number<=4)
             reference=xc_f90_func_info_get_references(xc_info, number)
             ref=xc_f90_func_reference_get_ref(reference)
             refnumber=refnumber+1
             if     (number==-1) then
                if (refnumber==1) then
                   write(*, '(a1,i1,2a)') '[', refnumber, '] ', trim(ref)
                else
                   write(*, '(22x,a,i1,2a)') '[', refnumber, '] ', trim(ref)
                endif
                exit
             else
                if (refnumber==1) then
                   write(*, '(a1,i1,2a)') '[', refnumber, '] ', trim(ref)
                else
                   write(*, '(22x,a,i1,2a)') '[', refnumber, '] ', trim(ref)
                endif
             endif
          end do
          call xc_f90_func_end(xc_func)             
          write (*,*)
       enddo
    endif
#endif

    write(iout6,*)
    if (lpotCoulomb.or.lpotFermi.or.lpotGauss) write(iout6,'(3x,"Nuclear potential: Coulomb")')
    if (lpotCoul2)   write(iout6,'(3x,"Nuclear potential: smoothed Coulomb 2D")')
    if (lpotCoul3)   write(iout6,'(3x,"Nuclear potential: smoothed Coulomb 3D")')
    if (lpotGSZ)     write(iout6,'(3x,"Nuclear potential: Green-Sellin-Zachor")')
    if (lpotGSZG)    write(iout6,'(3x,"Nuclear potential: Green-Sellin-Zachor + finite Gauss nuclei")')
    if (lpotHarm2)   write(iout6,'(3x,"Potential: harmonic 2D"$)')
    if (lpotHarm3.or.lextracule) write(iout6,'(3x,"Nuclear potential: harmonic 3D"$)')                   
    if (lpotKH)      write(iout6,'(3x,"Potential: Kramers-Henneberger")')
    if (lintracule)  write(iout6,'(3x,"Potential: harmonic 3D + 1/2r"$)') 
    if (lpotHooke)   write(iout6,'(3x,"Potential: harmonium"$)')

    if (lpotHooke.or.lpotHarm2.or.lpotHarm3.or.lextracule.or.lintracule) then
       write(iout6,'(6x,"force constant: ",1Pe9.2," au")') hooke
    endif
    if (lfefield) write(iout6,1300) ffield
    
    ! its electronic configuration

    call printElConf

#ifdef PRINT
! print=175: print 2-electron part of Fock equation for every orbital
    if (iprint(175).eq.1) call fockform
#endif
    
    write(iout6,*)
    write(iout6,*) '  Grid:'
    write(iout6,1040) nni,hni,nmu(ngrids),hmu(ngrids),rinf

    write(iout6,*)
    write(iout6,*) '  SCF: '
    write(iout6,1050) maxscf,ienterm,inoterm,recalcMMfactor,mpole

    write(iout6,*)
    if       (.not.lfixorb) then
       write (iout6,'(10x,"orbitals are relaxed")')
    else
       write (iout6,'(10x,"orbitals are kept frozen")')
    endif
    if       (.not.lfixcoul) then
       write (iout6,'(10x,"Coulomb potentials are relaxed")')
    else
       write (iout6,'(10x,"Coulomb potentials are kept frozen")')
    endif

    if (loffDiagLM) then
       write(*,'(10x,"non-zero off-diagonal Lagrange multipliers")')
       do iorb1=1,norb
          do iorb2=iorb1+1,norb
             if (offDiagLM(iorb1,iorb2)) then
                write(*,'(12x,i2,1x,a5,1x,a1,": ",(i2,1x,a5,1x,a1))') &
                     iorn(iorb1),bond(iorb1),gusym(iorb1),iorn(iorb2),bond(iorb2),gusym(iorb2)
             endif
          enddo
       enddo
       
    endif

    if (HF) then
       if (lfixexch) then
          write (iout6,'(10x,"exchange potentials are kept frozen")')
       elseif (lfastexch) then
          write (iout6,'(10x,"exchange potential for each pair of orbitals is relaxed once per single scf iteration")')
       elseif (.not.lfastexch) then
          write (iout6,'(10x,"exchange potential for each pair of orbitals is relaxed twice per single scf iteration")')
       endif
    endif

    write(iout6,1058) mpole

    if (openmp) write(iout6,'(/3x,"MT support: OpenMP")')
    if (mcsorpt) then
       if (lpthreadpoolq) then
          write(iout6,'(/3x,"MT support: pthread (a queue of threads)")')    
       else
          write(iout6,'(/3x,"MT support: pthread")')
       endif
    endif
    
    write(iout6,'(/3x,"(MC)SOR:")')

    if     (lorbmcsor) then
       write(iout6,'(10x,"MCSOR method used for relaxing orbitals (",i1," threads)")') nthreads
    else
       write(iout6,'(10x,"SOR method used for relaxing orbitals")')
    endif

    if (lpotmcsor) then
       if (openmp) then
          write(iout6,'(10x,"MCSOR method used for relaxing Coulomb/exchange potentials (",i2," threads)")') nthreads
       elseif (pthread.or.tpool) then
          write(iout6,'(10x,"SOR method used for relaxing Coulomb/exchange potentials")')
       endif
    endif
    write(iout6,'(/10x,"maximal number of Coulomb+exchange potentials per orbital =",i3)') nexchpots(1)
    
    write(iout6,'(/10x,"micro and macro SOR iterations for orbitals   = ",2i3)') maxsororb(2),maxsororb(1)
    write(iout6,'( 10x,"micro and macro SOR iterations for potentials = ",2i3)') maxsorpot(2),maxsorpot(1)    
    
    if (altSweeps) then
       if (meshOrdering=='col-wise')  &
            write(iout6,'(/10x,"ordering: column-wise + forward/backward sweeps")') 
       if (meshOrdering=='middle')    &
            write(iout6,'(/10x,"ordering: middle + forward/backward sweeps")')
       if (meshOrdering=='middlemt')    &
            write(iout6,'(/10x,"ordering: middle + forward/backward sweeps")')
       if (meshOrdering=='row-wise')  &
            write(iout6,'(/10x,"ordering: row-wise + forward/backward sweeps")')
       if (meshOrdering=='rrow-wis') &
            write(iout6,'(/10x,"ordering: row-wise + forward/backward sweeps")')
    else
       if (meshOrdering=='col-wise')  &
            write(iout6,'(/10x,"ordering: column-wise")')
       if (meshOrdering=='middle')    &
            write(iout6,'(/10x,"ordering: middle")')
       if (meshOrdering=='middlemt')    &
            write(iout6,'(/10x,"ordering: middle for MCSOR")')
       if (meshOrdering=='row-wise')  &
            write(iout6,'(/10x,"ordering: row-wise")')
       if (meshOrdering=='rrow-wis') &
            write(iout6,'(/10x,"ordering: reversed row-wise")')
    endif

    write(iout6,*)
    write(iout6,'(10x,"overrelaxation parameters:   orbitals       potentials ")')
    write(iout6,'(38x,2x,f5.3,7x,f5.3,3x,f5.3)') ovforb,ovfcoul,ovfexch

    write(iout6,*)
    write(iout6,'("   Machine accuracy      = ",1Pe11.2)') precis
    write(iout6,*)

    if (lengthfp.eq.8) then
       write(iout6,'("   Constants: "/&
            & "               pi        = ",1Pe25.16/&
            & "               bohr      = ",1Pe25.16," angstroms"/&
            & )') pii,bohr2ang
    else
       write(iout6,'("   constants: "/ &
            & "               pi        = ",1Pe45.34/ &
            & "               bohr      = ",1Pe45.34," angstroms"/ &
            & )') pii,bohr2ang
    endif


    write(iout6,*)
    write(iout6,1110)

    !     Text and data segments of the program have about 700 and 210 KB,
    !     respectively.

    !     The amount of memory required by the program depends on a given
    !     case and the memory is allocated dynamically (if alloc routine
    !     is supported by the system).

    !     The amount of 'statically' allocated memeory (i.e. for a given
    !     choice of maxnu and maxmu values) is due to several integer and
    !     real arrays depending on the current values of maxnu and maxmu.

    maxorb2=(maxorb*maxorb)/2+maxorb

    !     real arrays
    isizereal= (35+9*maxorb+4*maxorb2)*maxnu+(41)*maxmu+(10+4*maxorb+maxorb/2)*maxorb &
         +18*maxorb+8*maxorb*maxorb+(maxbasis+maxorb+5)*maxbasis+44000+maxorb*maxorb

    isizeint=6*maxorb*maxorb+11*maxorb2 +(25)*maxorb+4*maxbasis

    length0=lengthfp*isizereal+lengthint*isizeint

    lengtht=8*(length1+length2+length3+length4+length5)+lengthint*(length6)

    if (OED.or.HFS.or.DFT.or.SCMC) then
       nexch=1
    endif

    !     1MB=2^20B
    omb=2.0_PREC**20
    write(iout6,1200) dble(length0)/omb,dble(lengtht)/omb
1200 format(10x,'text+data                           ',9x,'  0.9 MB '/,&
         10x,'bss (common block contributions)    ',6x,f8.1,' MB'/,&
         10x,'dynamical allocation                ',6x,f8.1,' MB ')
    write(iout6,1211) dble((length1*lengthfp)/dble(omb))
    write(iout6,1212) dble((length2*lengthfp)/dble(omb))
    write(iout6,1213) dble((length3*lengthfp)/dble(omb))

1211 format(14x,'orbitals                              ',f8.1,' MB')
1212 format(14x,'Coulomb potentials                    ',f8.1,' MB')
1213 format(14x,'exchange potentials                   ',f8.1,' MB')

    call separator

1000 format(10x,a2,'(',f6.2,')',3x,a2,' (',f6.2,')',3x,'R =',f9.5,' bohr =',f8.5,' angstroms')
1040 format(10x,'nu (h_nu)  = ',i4,'  (',f7.5,')'/10x,'mu (h_mu)  = ',i4,'  (',f7.5,')'/,10x,'R_infty    =',f7.2)
1042 format(10x,'nu (h_nu) ',i6,'    (',f7.5,')'/10x,'mu (h_mu) ',i6,'              '/,10x,'R_infty   ',f6.2)
1050 format(10x,'thresholds',/14x,'scf iterations           =',i6/14x,'orbital energy           = 1.00E-',i2/ &
         14x,'orbital norm             = 1.00E-',i2/14x,'multipole moments recalc =',1p,e9.2&
         '  (mpole=',i1,')')

1058 format(/,10x,'multipole expansion coefficients = ',i2)

1060 format(28x,'(mc)sor iterations',/27x,' orbital  potentials')
1061 format(10x,i4,1x,a8,a1,1x,i8,i10)

1110 format(3x,'Memory usage:')

1300 format(11x,'finite electric field: ',1Pe9.2,' au')
    
1401 format(/3x,'Finite nuclei (Fermi nuclear charge distribution):',&
          /10x,'atomic masses:', 1Pe16.9,' ( ',a2,')',1Pe16.9,' ( ',a2,')')
1402 format(/3x,'Finite nuclei (Gauss nuclear charge distribution):',&
          /10x,'atomic masses:', 1Pe16.9,' ( ',a2,')',1Pe16.9,' ( ',a2,')')
    
    return
  end subroutine printCase

end module printInitdata
