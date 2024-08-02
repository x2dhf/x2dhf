! SPDX-License-Identifier: GPL-2.0-or-later

! Copyright (C) 1996       Leif Laaksonen, Dage Sundholm               
! Copyright (C) 1996-2024  Jacek Kobus 

module totalEnergyLXC
  implicit none
contains

  ! ### etotalLXC ###
  !
  !     Calculates total energy using several LXC functionals
  !
  subroutine etotalLXC 
    use params
    use detect
    use discrete
    use commons
    use utils
    use blas
    use elocc
    use integrals
    use memory
    use nabla
    use inout
    use sharedMemory
    use utils
#ifdef LIBXC
    use, intrinsic :: iso_c_binding
    use xc_f90_lib_m
    use libxc_funcs_m
#endif

    implicit none
    integer (KIND=IPREC) :: i,j,k,iborb,ibpot,iorb,isym,n,nmut,nan
    
    real (PREC) :: etsum,ocdown,ocup,w,wcorr,wcoul,wex,wndc,woneel,wtwoelLXC1,wtwoelLXC2
    ! real (PREC), dimension(*) :: psi,pot,excp,e,f0,wgt1,wgt2,wk0,wk1,wk2,wk3,rhot,rhotup,rhotdown,&
    !      grhot,grhotup,grhotdown,wk10,wk11,wk12,wk13

    !real (PREC), target :: work(length7)
    real (PREC), pointer :: w3(:),w4(:),w5(:),w6(:),w7(:),w8(:),w9(:),&
         w10(:),w11(:),w12(:),w13(:),w16(:)

#ifdef LIBXC
    type(xc_f90_func_info_t) :: xc_info
    type(xc_f90_func_t) :: xc_func
    integer(c_int) :: vmajor, vminor, vmicro, family_id, func_id, kind_id, err
    character(len=120) :: name, kind1, family, ref
    integer(c_size_t) :: size
    integer (KIND=IPREC) :: allocstat    
    real (PREC), allocatable :: w0(:),w1(:),w2(:)
    !integer (KIND=IPREC),external :: detectNaN
#endif

    real (PREC), dimension(:), pointer :: psi,excp,e,f0,wgt1,wgt2,&
         wk0,wk1,wk2,wk3,wk10,wk11,wk12,wk13,&
         rhot,rhotup,rhotdown,grhot,grhotup,grhotdown
    real (PREC), dimension(:), pointer :: psi1
 
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

    !(psi,pot,excp,e,f0,wgt1,wgt2,wk0,wk1,wk2,wk3,rhot,rhotup,rhotdown,&
    !   grhot,grhotup,grhotdown,wk10,wk11,wk12,wk13,work)

    
    wk0       =>scratchptr(          1:   mxsize8)
    wk1       =>scratchptr(   mxsize8+1: 2*mxsize8)
    wk2       =>scratchptr( 2*mxsize8+1: 3*mxsize8)
    wk3       =>scratchptr( 3*mxsize8+1: 4*mxsize8)
    rhot      =>scratchptr( 4*mxsize8+1: 5*mxsize8)            
    rhotup    =>scratchptr( 5*mxsize8+1: 6*mxsize8)
    rhotdown  =>scratchptr( 6*mxsize8+1: 7*mxsize8)
    grhot     =>scratchptr( 7*mxsize8+1: 8*mxsize8)            
    grhotup   =>scratchptr( 8*mxsize8+1: 9*mxsize8)
    grhotdown =>scratchptr( 9*mxsize8+1:10*mxsize8)
    wk10      =>scratchptr(10*mxsize8+1:11*mxsize8)
    wk11      =>scratchptr(11*mxsize8+1:12*mxsize8)
    wk12      =>scratchptr(12*mxsize8+1:13*mxsize8)
    wk13      =>scratchptr(13*mxsize8+1:14*mxsize8)

    w3  => scratch4lxcptr(mxsize8+1:2*mxsize8)
    w4  => scratch4lxcptr(2*mxsize8+1:3*mxsize8)    
    w5  => scratch4lxcptr(3*mxsize8+1:4*mxsize8)    
    w6  => scratch4lxcptr(4*mxsize8+1:5*mxsize8)
    w7  => scratch4lxcptr(5*mxsize8+1:6*mxsize8)
    w8  => scratch4lxcptr(6*mxsize8+1:7*mxsize8)
    w9  => scratch4lxcptr(7*mxsize8+1:8*mxsize8)
    w10 => scratch4lxcptr(8*mxsize8+1:9*mxsize8)    
    w11 => scratch4lxcptr(9*mxsize8+1:10*mxsize8)
    w12 => scratch4lxcptr(10*mxsize8+1:11*mxsize8)
    w13 => scratch4lxcptr(11*mxsize8+1:14*mxsize8)
    w16 => scratch4lxcptr(14*mxsize8+1:16*mxsize8)        

#ifdef DEBUG
! debug= 75: etotalLXC: testing testn2f routine     
    if (idebug(75).ne.0) then
       ! test n2f and nfnf routines
       print *,'Testing n2f'
       call testn2f(rhot,wk0,wk1,wk2,wk3)
    endif
#endif
    
#ifdef DEBUG
! debug= 76: etotalLXC: testing testnfng routine     
    if (idebug(76).ne.0) then
       print *,' '
       print *,'Testing nfng'
       call testnfng(rhot,rhotup,rhotdown,grhotdown,grhotup,wk0,wk1,wk2,wk3,wk10)
       call n2f(psi,wk0,wk1,wk2,wk3)
       call prod(mxsize,psi,wk3)
    endif
#endif
    ! calculate first contributions from one particle operators and Coulomb potential
    ! contributions within the same shell

    vnt=zero
    vkt=zero
    woneel=zero

    do iorb=1,norb
       if (inDFT(iorb)==0) cycle
       woneel=woneel+occ(iorb)*oneelii(iorb)
       vkt=vkt+occ(iorb)*vk(iorb)
       vnt=vnt+occ(iorb)*vn(iorb)
    enddo

    evt=woneel
    etot=evt+z1*z2/r
    enkin=vkt
    ennucel=vnt
    
    if (nel==1) return

    wcoul=zero
    wex=zero
    wcorr=zero
    wexhyb=zero
    
    ! contribution from coulomb interaction within the same shell
    ! calculate the coulomb potential contribution from all orbitals (include 1/2 factor )

    call zeroArray(mxsize,wk2)
    call zeroArray(mxsize,wk12)

    do iorb=1,norb
       if (inDFT(iorb).eq.0) cycle
       call daxpy (mxsize,occ(iorb)/two,exchptr(i2b(iorb):),ione,wk2,ione)
    enddo

    ! contribution from the Coulomb interaction (DFT style, i.e. plus <i|K(i,i)|i> term)
    do iorb=1,norb
       if (inDFT(iorb).eq.0) cycle
       iborb=i1b(iorb)
       psi1=>orbptr(iborb:)
       call prod2 (mxsize,psi1,psi1,wk0)
       call prod (mxsize,wk2,wk0)
       call dscal (mxsize,occ(iorb),wk0,ione)
       w=ddot(mxsize,wgt2,ione,wk0,ione)
       wcoul=wcoul+w
    enddo

    ! DFT exchange energy corrections
    call zeroarray(mxsize,rhot)    
    call zeroarray(mxsize,wk12)
    call zeroarray(mxsize,wk13)    

#ifdef LIBXC
    size=mxsize

    if (lxcPolar) then
       allocate(w0(2*mxsize),stat=allocstat)
       if (allocstat/=0) print *,'etotallxcHyb: allocate fails for w0'

       allocate(w1(3*mxsize),stat=allocstat)
       if (allocstat/=0) print *,'etotallxcHyb: allocate fails for w1'

       allocate(w2(2*mxsize),stat=allocstat)
       if (allocstat/=0) print *,'etotallxcHyb: allocate fails for w2'
       
       call zeroarray(mxsize,rhotup)
       call zeroarray(mxsize,rhotdown)
       call zeroarray(mxsize,wk13)    
    
       do iorb=1,norb
          !if (inhyd(iorb).eq.1) cycle
          if (inDFT(iorb).eq.0) cycle
          call exocc (iorb,ocup,ocdown)
          call prodas (mxsize,ocup,  psi(i1b(iorb):),psi(i1b(iorb):),rhotup)
          call prodas (mxsize,ocdown,psi(i1b(iorb):),psi(i1b(iorb):),rhotdown)
       enddo

       ! Libxc DFT functionals are implemented in C and expect 2*size rho
       ! array with spin-up and spin-down densities packed row-wise
       do i=1,mxsize
           w0(2*i-1)=rhotup(i)
           w0(2*i  )=rhotdown(i)
       enddo

       
       do n=1,lxcFuncs
          call xc_f90_func_init(xc_func, lxcFuncs2use(n), XC_POLARIZED, err)

          ! It turns out that for all the functionals used xc_f90_func_init returns err=0
          ! so the following piece of code is redundant
          ! if (err>0) then
          !    write(*,'("Error: xc_f90_func_init failed with err=",i5)') err
          !    stop "etotalLXC"
          ! endif

          xc_info=xc_f90_func_get_info(xc_func)

          select case (xc_f90_func_info_get_family(xc_info))
          case(XC_FAMILY_LDA)
             call xc_f90_lda_exc(xc_func, size, w0(1), wk12(1:))
             if (ldetectNaN) then
                nan=detectNaN(mxsize,wk12)
                if (nan>0) then
                   write(*,'("Error: NaN detected in (XC_FAMILY_LDA) xc_f90_lda_exc at",i5)') nan
                endif
             endif
             if (lfixNaN) then          
                call fixNaN(nni,mxnmu,wk12)
             endif

          case(XC_FAMILY_HYB_LDA)
             call xc_f90_lda_exc(xc_func, size, w0(1), wk12(1:))
             if (ldetectNaN) then
                nan=detectNaN(mxsize,wk12)
                if (nan>0) then
                   write(*,'("Error: NaN detected in (XC_FAMILY_LDA) xc_f90_lda_exc at",i5)') nan
                endif
             endif
             if (lfixNaN) then          
                call fixNaN(nni,mxnmu,wk12)
             endif

          case(XC_FAMILY_GGA)
             !print *,'GGA is supported for open-shell cases'
             call nfng(rhotup,rhotup,wk0,wk1,wk2,wk3,wk10,wk11,w2,grhotup)
             call nfng(rhotup,rhotdown,wk0,wk1,wk2,wk3,wk10,wk11,w2,grhot)       
             call nfng(rhotdown,rhotdown,wk0,wk1,wk2,wk3,wk10,wk11,w2,grhotdown)       
                  
             ! Libxc DFT functionals are implemented in C and expect 3*size
             ! array with \nabla rho(spin-up) \nabla rho(spin-up), 
             ! \nabla rho(spin-up) \nabla rho(spin-down)
             ! and array with \nabla rho(spin-down) \nabla rho(spin-down)
             ! contracted gradients of density packed row-wise
             do i=1,mxsize
                w1(3*i-2)=grhotup(i)
                w1(3*i-1)=grhot(i)
                w1(3*i  )=grhotdown(i)
             enddo
             
             call xc_f90_gga_exc(xc_func, size, w0(1), w1(1), wk12(1:))
             
             if (ldetectNaN) then
                nan=detectNaN(mxsize,wk12)
                if (nan>0) then
                   write(*,'("Error: NaN detected in (XC_FAMILY_LDA) xc_f90_lda_exc at",i5)') nan
                endif
             endif
             if (lfixNaN) then          
                call fixNaN(nni,mxnmu,wk12)
             endif

          case(XC_FAMILY_MGGA)
             ! calculate nabla rho nabla rho 
             call nfng(rhotup,rhotup,wk0,wk1,wk2,wk3,wk10,wk11,w2,grhotup)
             call nfng(rhotup,rhotdown,wk0,wk1,wk2,wk3,wk10,wk11,w2,grhot)       
             call nfng(rhotdown,rhotdown,wk0,wk1,wk2,wk3,wk10,wk11,w2,grhotdown)       
                  
             ! Libxc DFT functionals are implemented in C and expect 3*size
             ! array with \nabla rho(spin-up) \nabla rho(spin-up), 
             ! \nabla rho(spin-up) \nabla rho(spin-down)
             ! and array with \nabla rho(spin-down) \nabla rho(spin-down)
             ! contracted gradients of density packed row-wise
             do i=1,mxsize
                w1(3*i-2)=grhotup(i)
                w1(3*i-1)=grhot(i)
                w1(3*i  )=grhotdown(i)
             enddo

             !double *rho, double *sigma, double *lapl, double *tau, double
             !*exc)
             ! lapl[]: the laplacian of the density
             ! tau[]: the kinetic energy density

             ! w16==nabla^2 rho
             call n2f (rhotup,w3,w4,w5,w7)
             call n2f (rhotdown,w3,w4,w5,w8)
             do i=1,mxsize
                w16(2*i-1)=w7(i)
                w16(2*i  )=w8(i)
             enddo

             call zeroArray(mxsize,wk10)
             call tau(psi,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,grhotup)
             call dscal (mxsize,half,grhotup,ione)
             call dscal (mxsize,half,grhotup,ione)
             call zeroArray(mxsize,wk10)             
             call tau(psi,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,grhotdown)
             call dscal (mxsize,half,grhotdown,ione)
             call dscal (mxsize,half,grhotdown,ione)
             do i=1,mxsize
                w13(2*i-1)=grhotup(i)
                w13(2*i  )=grhotdown(i)
             enddo
             
             call xc_f90_mgga_exc(xc_func, size, w0, w1, w16, w13, wk12(1:))

          case(XC_FAMILY_HYB_GGA)
             print *,'HYB_GGA are not supported for open-shell cases'
             stop "etotalLXC" 
          case default
             write(*,'("Error! Unsupported family of libxc functionals.")')
             stop 'etotalLXC'
          end select
          call xc_f90_func_end(xc_func)
          
          call add (mxsize,wk12,wk13)
       enddo
       call add (mxsize,rhotdown,rhotup)
       call prod (mxsize,rhotup,wk13)

       ! take care of F4 factor
       call multf4(wk13)
       wtwoelLXC=ddot(mxsize,wgt2,ione,wk13,ione)
       wex=wtwoelLXC

       deallocate(w0,stat=allocstat)
       if (allocstat/=0) print *,'etotallxcHyb: deallocate fails for w0'

       deallocate(w1,stat=allocstat)
       if (allocstat/=0) print *,'etotallxcHyb: deallocate fails for w1'

       deallocate(w2,stat=allocstat)
       if (allocstat/=0) print *,'etotallxcHyb: deallocate fails for w2'

    else
       call zeroarray(mxsize,rhot)
       call zeroarray(mxsize,wk13)
       
       do iorb=1,norb
          !if (inhyd(iorb).eq.1) cycle
          if (inDFT(iorb).eq.0) cycle
          call prodas (mxsize,occ(iorb),psi(i1b (iorb):),psi(i1b (iorb):),rhot)
       enddo

       do n=1,lxcFuncs
          call xc_f90_func_init(xc_func, lxcFuncs2use(n), XC_UNPOLARIZED, err)
          
          ! xc_f90_func_init returns err=0 so the following piece of code is redundant
          ! if (err>0) then
          !    write(*,'("Error: xc_f90_func_init failed with err=",i5)') err
          !    stop "etotalLXHybC"
          ! endif

          xc_info=xc_f90_func_get_info(xc_func)
          select case (xc_f90_func_info_get_family(xc_info))
          case(XC_FAMILY_LDA)
             call xc_f90_lda_exc(xc_func, size, rhot(1:), wk12(1:))
          case(XC_FAMILY_HYB_LDA)
             call xc_f90_lda_exc(xc_func, size, rhot(1:), wk12(1:))
          case(XC_FAMILY_GGA)
             ! calculate nabla rho nabla rho 
             call nfng (rhot,rhot,wk0,wk1,wk2,wk3,rhotup,rhotdown,grhotdown,grhotup)

             ! nfng returns values without NaNs, and therefore the
             ! following piece of code is redundant
             ! if (nan>0) then
             !    write(*,'("Error: NaN detected in grhotup at",i5)') nan
             !    stop "etotalLXC"
             ! endif

             call xc_f90_gga_exc(xc_func, size, rhot(1:), grhotup(1:), wk12(1:))

             ! xc_f90_gga_exc returns values without NaNs, and therefore
             ! the following piece of code is redundant
             ! nan=detectNaN(mxsize,wk12)
             ! if (nan>0) then
             !    write(*,'("Error: NaN detected in wk12 at",i5)') nan
             !    stop "etotalLXC"
             ! endif
          case(XC_FAMILY_MGGA)
             ! calculate nabla rho nabla rho 
             call nfng (rhot,rhot,wk0,wk1,wk2,wk3,rhotup,rhotdown,grhotdown,grhotup)

             !double *rho, double *sigma, double *lapl, double *tau, double
             !*exc)
             ! lapl[]: the laplacian of the density
             ! tau[]: the kinetic energy density

             ! w6==nabla^2 rho
             call n2f (rhot,w3,w4,w5,w6)

             call zeroArray(mxsize,wk10)
             !call tau(psi,wk3,wk4,wk5,wk6,wk7,wk8,wk9,wk10,wk11,wk12,wk13)
             call tau(psi,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13)
             call dscal (mxsize,half,w13,ione)
             
             call xc_f90_mgga_exc(xc_func, size, rhot(1:), grhotup(1:), w6, w13, wk12(1:))

             ! It turns out that for all the functionals used xc_f90_gga_exc returns
             ! proper values, i.e. without NaNs, and therefore the following piece of
             ! code is redundant
             ! nan=detectNaN(mxsize,wk12)
             ! if (nan>0) then
             !    write(*,'("Error: NaN detected in wk12 at",i5)') nan
             !    stop "etotalLXC"
             ! endif
             
          case(XC_FAMILY_HYB_GGA)
             call nfng (rhot,rhot,wk0,wk1,wk2,wk3,rhotup,rhotdown,grhotdown,grhotup)
             call xc_f90_gga_exc(xc_func, size, rhot(1:), grhotup(1:), wk12(1:))
          case default
             write(*,'("Error! Unsupported family of libxc functionals.")')             
             stop 'etotalLXC'
          end select
          call xc_f90_func_end(xc_func)
          
          call add (mxsize,wk12,wk13)
       enddo
       call prod (mxsize,rhot,wk13)
       ! take care of F4 factor
       call multf4(wk13)
       wtwoelLXC=ddot(mxsize,wgt2,ione,wk13,ione)
    endif
    wex=wtwoelLXC

    if (lxcHyb) then
       call eExchangeLXC
    endif
#endif

    ! DFT correlation energy corrections
    evt=woneel+wcoul+wex+wcorr+wexhyb

    engt(1)=evt
    etot=evt+z1*z2/r
    virrat=(evt-vkt)/vkt
    enkin=vkt
    ennucel=vnt
    encoul=wcoul
    enexch=wex+wexhyb

    encouldft=wcoul
    enexchdft=wex+wexhyb
    edftcorr=wcorr
    entot=enkin+ennucel+encoul+enexch+edftcorr+z1*z2/r

#ifdef PRINT
! print= 66: etotalLXC: LXC total energy contributions and their sum
    if (iprint(66).ne.0) then
       write(*,'(2x,"etotalLXC: ")')        
       etsum=zero
       write(*,*)
       !write(*,'(" total energy contributions and their sum: ")')
       etsum=etsum+vkt
       write(*,'("  kinetic:                    ",1Pd25.16,1Pd25.16)') vkt,etsum
       etsum=etsum+vnt
       write(*,'("  nuclear attraction:         ",1Pd25.16,1Pd25.16)') vnt,etsum
       write(*,'("  one-electron:               ",1Pd25.16,1Pd25.16)') woneel,etsum
       etsum=etsum+wcoul
       write(*,'("  Coulomb (LXC):              ",1Pd25.16,1Pd25.16)') wcoul,etsum
       etsum=etsum+wex+wexhyb
       write(*,'("  exchange-correlation (LXC): ",1Pd25.16,1Pd25.16)') wex+wexhyb,etsum
       write(*,'("  two-electron:               ",1Pd25.16,1Pd25.16)') wcoul+wex+wexhyb
       write(*,'("  nuclear repulsion:          ",1Pd25.16,1Pd25.16)') z1*z2/r,etsum+z1*z2/r
       write(*,'("  total energy:               ",25x,1Pd25.16)') entot
       write(*,'("-------------------------------------------------------------------------"$)')
       write(*,'("----------------------")')
       write(*,*)
    endif
#endif
  end subroutine etotalLXC

  subroutine eExchangeLXC
    use params 
    use discrete 
    use commons 
    use utils
    use blas
    use exchContribs
    use diskInterface
    use integrals
    use inout
    use sharedMemory
    use utils
    
    implicit none
    integer (KIND=IPREC) :: ibex,iborb,iborb1,iborb2,ibpot,ibpot1,ibpot2,iorb,iorb1,iorb2,&
         isym,nmut,nmut1,nmut2

    real (PREC) :: epscharge,etsum,oc,oc1,oc2,ocx2,ocx1,w,wcouldft,wdcoul,&
         wex0,wex1,wex2,woneel
    real (PREC), dimension(:), pointer :: psi,excp,wgt2,wk0,wk1,wk2

    data epscharge /1.e-7_PREC/
#ifdef BLAS    
    real (PREC) ddot
    external ddot
#endif

    excp=>exchptr
    psi=>orbptr
    wgt2=>supplptr(i4b(14):)

    wk0       =>scratchptr(          1:   mxsize8)
    wk1       =>scratchptr(   mxsize8+1: 2*mxsize8)
    wk2       =>scratchptr( 2*mxsize8+1: 3*mxsize8)

    wex0=zero
    wex1=zero
    wex2=zero
    wdcoul=zero
    do iorb1=1,norb
      if (occ(iorb1).gt.one) then
         w=coulij (iorb1,iorb1,psi,excp,wgt2,wk0)
         wex0=wex0+occ(iorb1)*w/two
      endif
       ! calculate exchange interaction within pi, delta, etc. shells
      if (mm(iorb1).gt.0.and.abs(occ(iorb1)-one).gt.epscharge) then
         call exint (iorb1,ocx1)
         w=exchij(0,iorb1,iorb1,psi,excp,wgt2,wk0)
         wex1=wex1+ocx1*w
#ifdef PRINT
! print= 68: eExchangeLXC: exchange contributions          
         if (iprint(68).ne.0) then
            etsum=etsum-ocx1*w
            write(*,7032) iorn(iorb1),bond(iorb1),gusym(iorb1),&
                 iorn(iorb1),bond(iorb1),gusym(iorb1),&
                 iorn(iorb1),bond(iorb1),gusym(iorb1),&
                 iorn(iorb1),bond(iorb1),gusym(iorb1),w,-ocx1,etsum
7032        format('<',i4,1x,a5,a1,'| K (',i4,1x,a5,a1,i4,1x,a5,a1,&
                 ' ) |',i4,1x,a5,a1,' >',1Pd25.16, 0Pf8.2, 1Pd25.16)
         endif
#endif         
      endif
      
      do iorb2=iorb1+1,norb
          ! exchange interaction between shells (same lambda)
          call excont (iorb1,iorb2,ocx1,ocx2)
          w=exchij(0,iorb1,iorb2,psi,excp,wgt2,wk0)
          wex1=wex1+ocx1*w

#ifdef PRINT
! print= 68: eExchangeLXC: K contribution                     
          if (iprint(68).ne.0) then
             etsum=etsum-ocx1*w
             write(*,7036) iorn(iorb1),bond(iorb1),gusym(iorb1),&
                  iorn(iorb1),bond(iorb1),gusym(iorb1),&
                  iorn(iorb2),bond(iorb2),gusym(iorb2),&
                  iorn(iorb2),bond(iorb2),gusym(iorb2),w,-ocx1,etsum
7036         format('<',i4,1x,a5,a1,'| K (',i4,1x,a5,a1,i4,1x,a5,a1,&
                  ' ) |',i4,1x,a5,a1,' >',1Pd25.16, 0Pf8.2, 1Pd25.16)
          endif
#endif
          ! exchange interaction between shells (different lambda)
          !if (ilc(iex).gt.1) then
          if (ilc2(iorb1,iorb2)==2) then          
             w=exchij(1,iorb1,iorb2,psi,excp,wgt2,wk0)
             wex2=wex2+ocx2*w

#ifdef PRINT
! print= 68: eExchangeLXC: K1 contribution                     
             if (iprint(68).ne.0) then
                etsum=etsum-ocx2*w
                write(*,7038) iorn(iorb1),bond(iorb1),gusym(iorb1),&
                     iorn(iorb1),bond(iorb1),gusym(iorb1),&
                     iorn(iorb2),bond(iorb2),gusym(iorb2),&
                     iorn(iorb2),bond(iorb2),gusym(iorb2),w,-ocx2,etsum
7038            format('<',i4,1x,a5,a1,'| K1(',i4,1x,a5,a1,i4,1x,a5,a1,&
                     ' ) |',i4,1x,a5,a1,' >',1Pd25.16, 0Pf8.2, 1Pd25.16)
             endif
#endif             
          endif
       enddo
    enddo
    wexhyb=-alphaf*(wex0+wex1+wex2)

#ifdef PRINT          
! print= 69: eExchangeLXC: summary of exchange contributions          
    if (iprint(69).ne.0) then
       write(*,'(2x,"eExchangeLXC: ")') 
       write(*,*)
       write(*,'("  alpha (HYB):    ",2d25.16)') alphaf
       write(*,'("  exchange (HYB): ",2d25.16)') wexhyb
       write(*,'("  total:        : ",2d25.16)') entot
       write(*,*)
    endif
#endif
  end subroutine eExchangeLXC

  
end module totalEnergyLXC
