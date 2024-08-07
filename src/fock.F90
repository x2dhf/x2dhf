! SPDX-License-Identifier: GPL-2.0-or-later

! Copyright (C) 1996       Leif Laaksonen, Dage Sundholm               
! Copyright (C) 1996-2023  Jacek Kobus 

module fock
  implicit none
contains

  ! ### fock ###
  !
  !     Calculates the Fock potential for a given orbital.
  !
  subroutine checkasymp (n1,n2,a,n1st,n2st,incrn1,incrn2)
    use params

    integer (KIND=IPREC) :: i,j,n1,n2,n1st,n2st,incrn1,incrn2
    real (PREC), dimension(n1,n2) :: a

    write(6,*) '(E-V) values along eta=1 axis (R/2,+infty)'
    write(6,1000) (j, j=n2st,n2,incrn2)
    i=1
    write(6,1010) i, (a(i,j),j=n2st,n2,incrn2)

    write(6,*) '(E-V) values along eta=-1 axis (-infty,-R/2)'
    write(6,1000) (j, j=n2st,n2,incrn2)
    i=n1
    write(6,1010) i, (a(i,j),j=n2st,n2,incrn2)

01000 format(4i25)
01010 format(/,1x,i4,4e25.16,/(5x,4e25.16))

  end subroutine checkasymp

  function checkev(n1,n2,a)
    use params


    integer (KIND=IPREC) :: n1,n2
    real (PREC) :: checkev
    real (PREC), dimension(n1,n2) :: a

    !  parameter (ione=1)

    checkev=0.0_PREC

    if (a(ione,n2)*a(n1,n2).lt.0.0_PREC) then
       checkev=-1.0_PREC
    endif

  end function checkev

  ! ### fockHF ###
  !
  !     Calculates the Fock potential for a given orbital.
  !
  subroutine fockHF (iorb)
    use commons
    use params
    use discrete
    use scfshr
    use sharedMemory
    use utils
    use blas
    implicit none
    integer (KIND=IPREC) :: i,ibexp,iborb,iborb1,ibpot,ibpot1,idexp,&
         iorb,iorb1,ipc,kex,istart,istop
    real (PREC) :: coo,w
    real (PREC), dimension(:), pointer :: psi,excp,e,f0,f1,f2,f4,&
         fock1,fock2,wk1,wk2

    e=>supplptr(i4b(4):)
    excp=>exchptr
    f0=>supplptr(i4b(5):)
    f1=>supplptr(i4b(6):)
    f2=>supplptr(i4b(7):)
    f4=>supplptr(i4b(9):)
    psi=>orbptr

    fock1=>scratchptr(          1: 1*mxsize8)
    fock2=>scratchptr(1*mxsize8+1: 2*mxsize8)    
    wk1 =>scratchptr( 2*mxsize8+1: 3*mxsize8)
    wk2 =>scratchptr( 3*mxsize8+1: 4*mxsize8)

    iborb=i1b(iorb)
    ibpot=i2b(iorb)

    call dcopy (mxsize,f0,ione,fock1,ione)
    call daxpy (mxsize,ee(iorb,iorb),f1,ione,fock1,ione)
    
    if (mm(iorb).ne.0) then
       ! e enters the expression with minus sign which is already incorporated in e
       w=dble(mm(iorb)*mm(iorb))
       call daxpy (mxsize,w,e,ione,fock1,ione)
    endif

    ! in the asymptotic region E-V should be negative for orbitals to
    ! converge; it may be not the case when external electric field is
    ! too strong or the practical infinity is too small
    if (checkev(nni,mxnmu,fock1).lt.zero) then
       write(*,'(" ... Warning! E-V>0 in tail region for orbital",i2,1x,a8,a1)') &
            iorn(iorb),bond(iorb),gusym(iorb)
    endif

#ifdef PRINT
! print= 95: checking (E-V) values for orbital
    if (iprint(95).ne.0) then
       write(*,'("... Checking (E-V) values for orbital",i2,1x,a8,a1)') &
            iorn(iorb),bond(iorb),gusym(iorb)
       call checkasymp (nni,mxnmu,fock1,ione,ione,incrni,incrmu)
    endif
#endif
    
    call zeroArray(mxsize,fock2)
    if (nel.eq.1) return

    call zeroArray(mxsize,wk1)

    do iorb1=firstOrb,lastOrb
       iborb1=i1b(iorb1)
       kex=iorb+(iorb1-1)*norb
       ipc=iorb1+(iorb-1)*norb
       if (iorb.le.iorb1) then
          idexp=iorb+iorb1*(iorb1-1)/2
       else
          idexp=iorb1+iorb*(iorb-1)/2
       endif
       
       ibpot1=i2b(iorb1)
       ibexp=i3b(idexp)
       coo=occ(iorb1)

       if (iorb.eq.iorb1) coo=coo-1.0_PREC
       call daxpy (mxsize,coo,excp(ibpot1:),ione,wk1,ione)

       if (iorb1.ne.iorb)  then
          call dcopy (mxsize,excp(ibexp:),ione,wk2,ione)
          call dscal (mxsize,gec(kex),wk2,ione)
          if (ilc(idexp).gt.1) then
             call daxpy (mxsize,gec(kex+norb*norb),excp(ibexp+mxsize:),ione,wk2,ione)
          endif

          if (abs(ee(iorb1,iorb))>epsilon(zero)) then
             call daxpy (mxsize,ee(iorb1,iorb),f4,ione,wk2,ione)
          endif

          call prod (mxsize,psi(iborb1:),wk2)          
          call add  (mxsize,wk2,fock2)
       else
          if ((mm(iorb).gt.0).and.(ilc(idexp).gt.0)) then
             call dcopy (mxsize,excp(ibexp:),ione,wk2,ione)
             call prod  (mxsize,psi(iborb1:),wk2)
             call daxpy(mxsize,gec(kex),wk2,ione,fock2,ione)
          endif
       endif

#ifdef DEBUG
! debug= 92: fockHF: ee(iorb1,iorb)
       if (idebug(92).ne.0) then
          write(*,1000) iorn(iorb),bond(iorb),gusym(iorb),&
               iorn(iorb1),bond(iorb1),gusym(iorb1),ee(iorb1,iorb)
1000      format(4x,i2,1x,a5,1x,a1,1x,i2,1x,a5,1x,a1,1x,4f8.3)
       endif
#endif
    enddo

#ifdef DEBUG
! debug= 93: fockHF: Coulomb potential contribution
    if (idebug(93).ne.0) then
       write(*,*) 'fock: Coulomb potential contribution'
       call pmtx (nni,nmu(1),wk1,ione,ione,incrni,incrmu)
    endif
#endif

#ifdef DEBUG    
! debug= 94: fockHF: Exchange potential contribution 
    if (idebug(94).ne.0) then
       write(*,*) 'fock: Exchange potential contribution'
       call pmtx (nni,nmu(1),fock2,ione,ione,incrni,incrmu)
    endif
#endif

    call dcopy (mxsize,wk1,ione,coulombptr(i2b(iorb):),ione)    
    call dcopy (mxsize,fock2,ione,exchangeptr(i2b(iorb):),ione)

    ! multiply Coulomb and exchange potentials by f2    
    call prod (mxsize,f2,wk1)
    call prod (mxsize,f2,fock2)
    call add  (mxsize,wk1,fock1)

  end subroutine fockHF



  ! ### fock4TED ###
  !
  subroutine fock4TED (iorb)
    use commons
    use params
    use discrete
    use scfshr
    use sharedMemory
    use utils
    use blas
    implicit none
    integer (KIND=IPREC) :: i,ibexp,iborb,iborb1,ibpot,ibpot1,idexp,&
         iorb,iorb1,ipc,kex,istart,istop
    real (PREC) :: coo,w
    real (PREC), dimension(:), pointer :: psi,excp,e,f0,f1,f2,f4,&
         fock1,fock2,wk1,wk2

    e=>supplptr(i4b(4):)
    excp=>exchptr
    f0=>supplptr(i4b(5):)
    f1=>supplptr(i4b(6):)
    f2=>supplptr(i4b(7):)
    f4=>supplptr(i4b(9):)
    psi=>orbptr

    fock1=>scratchptr(          1: 1*mxsize8)
    fock2=>scratchptr(1*mxsize8+1: 2*mxsize8)    
    wk1 =>scratchptr( 2*mxsize8+1: 3*mxsize8)
    wk2 =>scratchptr( 3*mxsize8+1: 4*mxsize8)

    if (loffDiagLM) then
       firstOrb=iorb
       lastOrb=norb
    else
       firstOrb=iorb
       lastOrb=iorb
    endif
    
    iborb=i1b(iorb)
    ibpot=i2b(iorb)

    call dcopy (mxsize,f0,ione,fock1,ione)
    call daxpy (mxsize,ee(iorb,iorb),f1,ione,fock1,ione)
    
    if (mm(iorb).ne.0) then
       ! e enters the expression with minus sign which is already incorporated in e
       w=dble(mm(iorb)*mm(iorb))
       call daxpy (mxsize,w,e,ione,fock1,ione)
    endif

    ! in the asymptotic region E-V should be negative for orbitals to
    ! converge; it may be not the case when external electric field is
    ! too strong or the practical infinity is too small
    if (checkev(nni,mxnmu,fock1).lt.zero) then
       write(*,'(" ... Warning! E-V>0 in tail region for orbital",i2,1x,a8,a1)') &
            iorn(iorb),bond(iorb),gusym(iorb)
    endif

#ifdef PRINT
! print= 95: checking (E-V) values for orbital
    if (iprint(95).ne.0) then
       write(*,'("... Checking (E-V) values for orbital",i2,1x,a8,a1)') &
            iorn(iorb),bond(iorb),gusym(iorb)
       call checkasymp (nni,mxnmu,fock1,ione,ione,incrni,incrmu)
    endif
#endif

    call zeroArray(mxsize,fock2)
    if (nel.eq.1) return
    
    call zeroArray(mxsize,wk1)
    do iorb1=firstOrb,lastOrb
       iborb1=i1b(iorb1)
       kex=iorb+(iorb1-1)*norb
       ipc=iorb1+(iorb-1)*norb
       if (iorb.le.iorb1) then
          idexp=iorb+iorb1*(iorb1-1)/2
       else
          idexp=iorb1+iorb*(iorb-1)/2
       endif
       
       ibpot1=i2b(iorb1)
       ibexp=i3b(idexp)
       coo=occ(iorb1)
       
       if (iorb.eq.iorb1) then
          coo=coo-1.0_PREC
          call daxpy (mxsize,coo,excp(ibpot1:),ione,wk1,ione)
       endif
       
       if (iorb1.ne.iorb)  then
          if (abs(ee(iorb1,iorb))>epsilon(zero)) then
             call daxpy (mxsize,ee(iorb1,iorb),f4,ione,wk2,ione)
             call prod (mxsize,psi(iborb1:),wk2)          
             call add  (mxsize,wk2,fock2)
          endif
       else
          if ((mm(iorb).gt.0).and.(ilc(idexp).gt.0)) then
             call dcopy (mxsize,excp(ibexp:),ione,wk2,ione)
             call prod  (mxsize,psi(iborb1:),wk2)
             call daxpy(mxsize,gec(kex),wk2,ione,fock2,ione)
          endif
       endif

#ifdef DEBUG
! debug= 92: fockHF: ee(iorb1,iorb)
       if (idebug(92).ne.0) then
          write(*,1000) iorn(iorb),bond(iorb),gusym(iorb),&
               iorn(iorb1),bond(iorb1),gusym(iorb1),ee(iorb1,iorb)
1000      format(4x,i2,1x,a5,1x,a1,1x,i2,1x,a5,1x,a1,1x,4f8.3)
       endif
#endif
    enddo

#ifdef DEBUG
! debug= 93: fockHF: Coulomb potential contribution
    if (idebug(93).ne.0) then
       write(*,*) 'fock: Coulomb potential contribution'
       call pmtx (nni,nmu(1),wk1,ione,ione,incrni,incrmu)
    endif
#endif

#ifdef DEBUG    
! debug= 94: fockHF: Exchange potential contribution 
    if (idebug(94).ne.0) then
       write(*,*) 'fock: Exchange potential contribution'
       call pmtx (nni,nmu(1),fock2,ione,ione,incrni,incrmu)
    endif
#endif

    call dcopy (mxsize,wk1,ione,coulombptr(i2b(iorb):),ione)    
    call dcopy (mxsize,fock2,ione,exchangeptr(i2b(iorb):),ione)

    ! multiply Coulomb and exchange potentials by f2    
    call prod (mxsize,f2,wk1)
    call prod (mxsize,f2,fock2)
    call add  (mxsize,wk1,fock1)

  end subroutine fock4TED



  ! ### fockExch ###
  !
  !     Calculates the exchange part of the Fock potential for a given orbital.
  !
  subroutine fockExch (iorb)
    use commons
    use params
    use discrete
    use scfshr
    use sharedMemory
    use utils
    use blas

    implicit none
    integer (KIND=IPREC) :: i,ibexp,iborb,iborb1,ibpot,ibpot1,idexp,iorb,iorb1,ipc,kex
    real (PREC) :: coo,w
    !real (PREC), dimension(*) :: psi,pot,excp,e,f0,f1,f2,f4,fock1,fock2,wk1,wk2

    real (PREC), dimension(:), pointer :: psi,excp,e,f0,f1,f2,f4,&
         fock1,fock2,wk2
    print *,"no fockExch"
    return
    
    e=>supplptr(i4b(4):)
    excp=>exchptr
    f0=>supplptr(i4b(5):)
    f1=>supplptr(i4b(6):)
    f2=>supplptr(i4b(7):)
    f4=>supplptr(i4b(9):)
    psi=>orbptr

    fock1=>scratchptr(          1: 1*mxsize8)
    fock2=>scratchptr(1*mxsize8+1: 2*mxsize8)    
    !wk1 =>scratchptr( 2*mxsize8+1: 3*mxsize8)
    wk2 =>scratchptr( 3*mxsize8+1: 4*mxsize8)

    iborb=i1b(iorb)
    ibpot=i2b(iorb)

    do i=1,mxsize
       fock2(i)=0.0_PREC
    enddo

    if (nel.eq.1) return

    do iorb1=1,norb
       iborb1=i1b(iorb1)
       kex=iorb+(iorb1-1)*norb
       ipc=iorb1+(iorb-1)*norb
       
       if (iorb.le.iorb1) then
          idexp=iorb+iorb1*(iorb1-1)/2
       else
          idexp=iorb1+iorb*(iorb-1)/2
       endif
       
       ibexp=i3b(idexp)

       if (iorb1.ne.iorb)  then
          call dcopy (mxsize,excp(ibexp:),ione,wk2,ione)
          call dscal (mxsize,gec(kex),wk2,ione)
          if (ilc(idexp).gt.1) then
             call daxpy (mxsize,gec(kex+norb*norb),excp(ibexp+mxsize:),ione,wk2,ione)
          endif

          if (abs(ee(iorb1,iorb))>epsilon(zero)) then
             call daxpy (mxsize,ee(iorb1,iorb),f4,ione,wk2,ione)
          endif
          
          call prod (mxsize,psi(ibpot1:),wk2)
          call add  (mxsize,wk2,fock2)
       else
          if ((mm(iorb).gt.0).and.(ilc(idexp).gt.0)) then
             call dcopy (mxsize,excp(ibexp:),ione,wk2,ione)
             call prod  (mxsize,psi(iborb1:),wk2)
             call dscal (mxsize,gec(kex),wk2,ione)
             call add (mxsize,wk2,fock2)
          endif
       endif

#ifdef DEBUG
! debug= 92: fockExch: ee(iorb1,iorb)
       if (idebug(92).ne.0) then
          write(*,1000) iorn(iorb),bond(iorb),gusym(iorb),&
               iorn(iorb1),bond(iorb1),gusym(iorb1),ee(iorb1,iorb)
1000      format(4x,i2,1x,a5,1x,a1,1x,i2,1x,a5,1x,a1,1x,4f8.3)
       endif
#endif
    enddo

    ! multiply exchange potentials by f2
    call prod (mxsize,f2,fock2)

    if (lxcHyb) then
       call dscal (mxsize,alphaf,fock2,ione)
    endif

    return
  end subroutine fockExch
  
#ifdef LIBXC
  ! ### fockLXCunpol ###
  !
  ! Calculates exchange and correlation potentials via calls to 
  ! libxcf90 and libxc libraries (unpolarized version).
  !
  subroutine fockLXCunpol(iorb)
    use params
    use detect
    use discrete
    use memory
    use scfshr
    use commons
    use blas
    use elocc
    use dftvxc
    use nabla
    use sharedMemory
    use utils

    use, intrinsic :: iso_c_binding
    use xc_f90_lib_m
    use libxc_funcs_m    
    
    implicit none
    integer (KIND=IPREC) :: i,ii,iborb,iborb1,ibpot,ibpot1,iorb,iorb1,ipc,isiorb1
    integer (KIND=IPREC) :: ibexp,idexp,kex,n,nan

    real (PREC) :: const13,coo,ocdown,ocup,tmpf,w
    !real (PREC),dimension(*) :: psi,pot,excp,e,f0,f1,f2,f4,fock1,fock2
    !rhot,wk3,wk4,wk5,wk6,wk7,wk8,wk9,wk10,wk11,wk12,wk13

    real (PREC), dimension(:), pointer :: psi,excp,e,f0,f1,f2,f4,fock1,fock2
   
    !real (PREC), target :: work(length7) 
    real (PREC), dimension(:), pointer :: wk3,wk4,wk5,wk6,wk7,wk8,wk9,&
         wk10,wk11,wk12,wk13,work

    real (PREC), dimension(:), pointer :: rhoup,rhodown,rho,sigma,vrho,vsigma

    type(xc_f90_func_info_t) :: xc_info
    type(xc_f90_func_t) :: xc_func
    type(xc_f90_func_info_t) :: info
    integer(c_int) :: vmajor, vminor, vmicro, family_id, func_id, kind_id, err
    integer(c_size_t) :: size 
    character(len=120) :: name, kind1, family, ref
    real*8 :: xc_omega,xc_alpha,xc_beta
    integer (KIND=IPREC) :: allocstat

    e=>supplptr(i4b(4):)
    excp=>exchptr
    f0=>supplptr(i4b(5):)
    f1=>supplptr(i4b(6):)
    f2=>supplptr(i4b(7):)
    f4=>supplptr(i4b(9):)
    psi=>orbptr
    !work=>scratch4lxcptr
    work=>scratchptr

    fock1   => work(           1: 1*mxsize8)
    ! fock2==rhs in calling routine
    fock2   => work( 1*mxsize8+1: 2*mxsize8)
    rhoup   => work( 2*mxsize +1: 3*mxsize)
    rhodown => work( 3*mxsize +1: 4*mxsize)    
    rho     => work( 4*mxsize +1: 6*mxsize)
    vrho    => work( 6*mxsize +1: 8*mxsize)    
    sigma   => work( 9*mxsize +1:12*mxsize)
    vsigma  => work(12*mxsize +1:15*mxsize)
    wk3     => work(15*mxsize +1:16*mxsize8)
    wk4     => work(16*mxsize8+1:17*mxsize8)    
    wk5     => work(17*mxsize8+1:18*mxsize8)    
    wk6     => work(18*mxsize8+1:19*mxsize8)
    wk7     => work(19*mxsize8+1:20*mxsize8)
    wk8     => work(21*mxsize8+1:22*mxsize8)
    wk9     => work(22*mxsize8+1:23*mxsize8)
    wk10    => work(23*mxsize8+1:24*mxsize8)    
    wk11    => work(24*mxsize8+1:25*mxsize8)
    wk12    => work(25*mxsize8+1:26*mxsize8)
    wk13    => work(26*mxsize8+1:27*mxsize8)    

    iborb=i1b(iorb)
    ibpot=i2b(iorb)

    call zeroArray(mxsize,fock1)
    call zeroArray(mxsize,fock2)

    if (nel==1) goto 100

    call zeroArray(mxsize,rho)
    !call zeroArray(mxsize,vrho)

    do iorb1=1,norb
       iborb1 =i1b (iorb1)
       call prodas (mxsize,occ(iorb1),psi(iborb1:),psi(iborb1:),rho)
    enddo

    ! choose a given functional from libxc library     
    size=mxsize

    ! XC_UNPOLARIZED version
    do n=1,lxcFuncs
       call xc_f90_func_init(xc_func, lxcFuncs2use(n), XC_UNPOLARIZED,  err)
       call xc_f90_func_set_dens_threshold(xc_func, dens_threshold)
       
       ! It turns out that for all the functionals used xc_f90_func_init returns err=0
       ! so the following piece of code is redundant
       ! if (err>0) then
       !    write(*,'("Error: xc_f90_func_init failed with err=",i5)') err
       !    stop "fockLXC"
       ! endif
       
       info=xc_f90_func_get_info(xc_func)
       
       select case (xc_f90_func_info_get_family(info))
       case(XC_FAMILY_LDA)
          call xc_f90_lda_vxc(xc_func, size, rho, vrho)
       case(XC_FAMILY_HYB_LDA)
          call xc_f90_lda_vxc(xc_func, size, rho, vrho)
       case(XC_FAMILY_GGA)
          ! calculate nabla rho nabla rho
          
          ! sigma = \naba \rho \nabla \rho: contracted gradients of the density
          
          ! vrho: first partial derivative of the energy per unit
          ! volume in terms of the density
          
          ! =vsigma: first partial derivative of the energy per unit
          ! volume in terms of sigma
          call nfng (rho,rho,wk3,wk4,wk5,wk6,wk7,wk8,wk9,sigma)
          
          !call xc_f90_gga_vxc(xc_func, size, rho, sigma, wk12, wk3 )
          !call xc_f90_gga_vxc(xc_func, size, rho, sigma, vrho, vsigma )
          call xc_f90_gga_vxc(xc_func, size, rho, sigma, vrho, vsigma )
          ! do i=1,mxsize,1
          !      write(*,'(i6,4e16.8)') i, rho(i), sigma(i), vrho(i), vsigma(i) 
          ! enddo
          ! stop "upol"
          
          !call xc_hyb_cam_coef(xc_func, xc_omega, xc_alpha, xc_beta)
          !print *,xc_omega,xc_alpha,xc_beta
          !stop "test"
          ! It turns out that for for some GGA functionals xc_f90_gga_exc
          ! return NaN(s) value(s), and therefore the following piece of
          ! code must be used to detect them and fix via interpolation
          
          if (ldetectNaN) then
             nan=detectNaN(mxsize,wk3)
             if (nan>0) then
                write(*,'("Error: NaN detected at",i5)') nan
             endif
          endif
          if (lfixNaN) then          
             call fixNaN(nni,mxnmu,wk3)
          endif
          
          ! do i=1,5
          !    write(*,'(i6,4e16.8)') i, rho(i),sigma(i),vrho(i),vsigma(i)
          ! enddo
          
          !stop "npol"
          
          ! v_xc=vrho - 2 (nabla vsigma  nabla n + vsigma nabla^2 rho)
          ! vrho==wk12
          ! vsigma==wk3 
          ! wk11==(nabla vsignma  nabla n)
          call nfng (vsigma,rho,wk4,wk5,wk6,wk7,wk8,wk9,wk10,wk11)
          
          ! wk8==nabla^2 rho
          call n2f (rho,wk5,wk6,wk7,wk8)
          
          
          ! wk8==vsigma nabla^2 rho
          call prod (mxsize,vsigma,wk8)
          
          ! do i=1,10,1
          !    write(*,'(i6,4e16.8)') i, rho(i),sigma(i),vrho(i),vsigma(i)
          ! enddo
          
          ! wk8==(nabla vsignma  nabla n) + vsigma nabla^2 n
          call add (mxsize,wk11,wk8)
          
          ! wk12==vrho - 2 (nabla vsignma  nabla n) + vsigma nabla^2 rho
          call daxpy (mxsize,-two,wk8,ione,vrho,ione)
          
          ! do i=1,mxsize,1
          !    write(*,'(i6,4e16.8)') i, rho(i), vrho(i)
          ! enddo
          !stop "npol"
          
          
       case(XC_FAMILY_MGGA)
          ! calculate nabla rho nabla rho
          stop "fockLXC"
          
          ! double *rho, double *sigma, double *lapl, double *tau,
          ! double *vrho, double *vsigma, double *vlapl, double *vtau)
          ! sigma = \naba \rho \nabla \rho: contracted gradients of the density
          
          ! vrho: first partial derivative of the energy per unit
          ! volume in terms of the density
          
          ! =vsigma: first partial derivative of the energy per unit
          ! volume in terms of sigma
          call nfng (rho,rho,wk3,wk4,wk5,wk6,wk7,wk8,wk9,sigma)
          
          !call xc_f90_gga_vxc(xc_func, size, rho, sigma, wk12, wk3 )
          !call xc_f90_gga_vxc(xc_func, size, rho, sigma, vrho, vsigma )
          
          ! call xc_f90_mgga_vxc(xc_func, size, rho,  sigma,  lapl,  tau,&
          !                                    vrho, vsigma, vlapl, vtau )
          
          ! return NaN(s) value(s), and therefore the following piece of
          ! code must be used to detect them and fix via interpolation
          
          if (ldetectNaN) then
             nan=detectNaN(mxsize,wk3)
             if (nan>0) then
                write(*,'("Error: NaN detected at",i5)') nan
             endif
          endif
          if (lfixNaN) then          
             call fixNaN(nni,mxnmu,wk3)
          endif
          
          ! do i=1,5
          !    write(*,'(i6,4e16.8)') i, rho(i),sigma(i),vrho(i),vsigma(i)
          ! enddo
          
          !stop "npol"
          
          ! v_xc=vrho - 2 (nabla vsigma  nabla n + vsigma nabla^2 rho)
          ! vrho==wk12
          ! vsigma==wk3 
          ! wk11==(nabla vsignma  nabla n)
          call nfng (vsigma,rho,wk4,wk5,wk6,wk7,wk8,wk9,wk10,wk11)
          
          ! wk8==nabla^2 rho
          call n2f (rho,wk5,wk6,wk7,wk8)
          
          ! wk8==vsigma nabla^2 rho
          call prod (mxsize,vsigma,wk8)
          
          ! do i=1,10,1
          !    write(*,'(i6,4e16.8)') i, rho(i),sigma(i),vrho(i),vsigma(i)
          ! enddo
          
          ! wk8==(nabla vsignma  nabla n) + vsigma nabla^2 n
          call add (mxsize,wk11,wk8)
          
          ! wk12==vrho - 2 (nabla vsignma  nabla n) + vsigma nabla^2 rho
          call daxpy (mxsize,-two,wk8,ione,vrho,ione)
          
          ! do i=1,mxsize,1
          !    write(*,'(i6,4e16.8)') i, rho(i), vrho(i)
          ! enddo
          !stop "npol"
          
       case(XC_FAMILY_HYB_GGA)
          call nfng (rho,rho,wk3,wk4,wk5,wk6,wk7,wk8,wk9,sigma)
          call xc_f90_gga_vxc(xc_func, size, rho, sigma, vrho, vsigma )
          ! if (ldetectNaN) then
          !    nan=detectNaN(mxsize,wk3)
          !    if (nan>0) then
          !       write(*,'("Error: NaN detected at",i5)') nan
          !    endif
          ! endif
          
          ! if (lfixNaN) then          
          !    call fixNaN(nni,mxnmu,wk3)
          ! endif
          
          ! call nfng (vsigma,rho,wk4,wk5,wk6,wk7,wk8,wk9,wk10,wk11)
          ! call n2f (rho,wk5,wk6,wk7,wk8)
          ! call prod (mxsize,vsigma,wk8)
          ! call add (mxsize,wk11,wk8)
          ! call daxpy (mxsize,-two,wk8,ione,vrho,ione)
          
          ! v_xc=vrho - 2 (nabla vsigma  nabla n + vsigma nabla^2 rho)
          ! vrho==wk12
          ! vsigma==wk3 
          ! wk11==(nabla vsignma  nabla n)
          call nfng (vsigma,rho,wk4,wk5,wk6,wk7,wk8,wk9,wk10,wk11)
          
          ! wk8==nabla^2 rho
          call n2f (rho,wk5,wk6,wk7,wk8)
          
          ! wk8==vsigma nabla^2 rho
          call prod (mxsize,vsigma,wk8)
          
          ! wk8==(nabla vsignma  nabla n) + vsigma nabla^2 n
          call add (mxsize,wk11,wk8)
          
          ! wk12==vrho - 2 (nabla vsignma  nabla n) + vsigma nabla^2 rho
          call daxpy (mxsize,-two,wk8,ione,vrho,ione)
          
       case default
          write(*,'("Error! Unsupported family of libxc functionals.")')
          stop 'fockLXC'
       end select
       call xc_f90_func_end(xc_func)
       
       call add (mxsize,vrho,fock1)
      
    enddo
    ! take care of F4 factor (r\xi/2) which transforms the potential into \tilde form
    call multf4(fock1)
    !alphaf=zero
    ! let's simulate HF case. no LXC functional
    if (idebug(999)==1) then
       call zeroArray(mxsize,fock1)    
       alphaf=one
    endif
    
    ! store the local exchange/correlation potential in extra portion
    ! of pot array (to be used by EaDFT and EabDFT)
    !call dcopy (mxsize,fock1,ione,pot(length2-2*mxsize+1),ione)
    !call dcopy (mxsize,fock1,ione,excp(length2-2*mxsize+1:),ione)
    call dcopy (mxsize,fock1,ione,coulombptr(i2b(iorb):),ione)    

    call zeroArray(mxsize,wk3)       
    call zeroArray(mxsize,fock2)    

    do iorb1=1,norb
       if (inDFT(iorb1).eq.0) cycle
       iborb1=i1b(iorb1)
       ibpot1=i2b(iorb1)
       !ipc=iorb1+norb*(iorb-1)
       
       ! add the Coulomb potential which in DFT also includes the
       ! contribution from the orbital
       call daxpy (mxsize,occ(iorb1),excp(ibpot1:),ione,wk3,ione)
       
       ! add contributions from off-diagonal Lagrange multipliers.
       if (iorb1.ne.iorb.and.abs(ee(iorb1,iorb))>epsilon(zero)) then
          do i=1,mxsize
             fock2(i)=fock2(i)+ee(iorb1,iorb)*f4(i)*psi(iborb1+i-1)
          enddo
       endif
    enddo

    ! add Coulomb potential to fock1
    call add (mxsize,wk3,fock1)
    call dcopy (mxsize,fock1,ione,coulombptr(i2b(iorb):),ione)    

    ! do not add wk3 to pot(length2-2*mxsize+1) since this
    ! contribution is taken care of by EaDFT/EabDFT itself 
    
    ! multiply fock2 by f2 and save for EaDFT

    call prod (mxsize,f2,fock1)
    
    !call dcopy (mxsize,fock2,ione,pot(length2-mxsize+1),ione)
    call dcopy (mxsize,fock2,ione,excp(length2-mxsize+1:),ione)
    call dcopy (mxsize,fock2,ione,exchangeptr(ibpot:),ione)    
    
    call prod (mxsize,f2,fock2) 
    ! no hybrid functionals OK
    
    !call zeroArray(mxsize,fock2)

    ! two-electron contributions OK
    if (lxcHyb) then
       ! In case of hybrid functionals one has to add the alphaf fraction of
       ! the exact HF exchange. This must also include q/2<i|K(i)|i> part.
       call zeroArray(mxsize,wk3)        
       ibpot=i2b(iorb)
       ! in the local exchange approximation the Coulomb potential also includes the
       ! contribution from the orbital
       !call dcopy (mxsize,pot(ibpot),ione,wk3,ione)

       ! add real exchange contributions
       ! substract a part of exchange contributions 
       call zeroArray(mxsize,wk4)

       ! add contributions to the exchange potential
       do iorb1=1,norb
          iborb1=i1b(iorb1)
          kex=iorb+norb*(iorb1-1)
          ipc=iorb1+norb*(iorb-1)
          if (iorb.le.iorb1) then
             idexp=iorb+iorb1*(iorb1-1)/2
          else
             idexp=iorb1+iorb*(iorb-1)/2
          endif

          ibpot1=i2b(iorb1)
          ibexp=i3b(idexp)

          coo=occ(iorb1)
          if (iorb1==iorb) coo=coo-one
          coo=alphaf*coo
          
          if (iorb1.ne.iorb) then
             call dcopy (mxsize,excp(ibexp:),ione,wk3,ione)
             call dscal (mxsize,gec(kex),wk3,ione)
             if (ilc(idexp).gt.1) then
                call daxpy (mxsize,gec(kex+norb*norb),excp(ibexp+mxsize:),ione,wk3,ione)
             endif
             call prod (mxsize,psi(iborb1:),wk3)
             call add  (mxsize,wk3,wk4)
          else
             if ((mm(iorb).gt.0).and.(ilc(idexp).gt.0)) then
                call dcopy (mxsize,excp(ibexp:),ione,wk3,ione)
                call prod  (mxsize,psi(iborb1:),wk3)
                call daxpy (mxsize,gec(kex),wk3,ione,wk4,ione)
             endif
          endif
       enddo

       ! add exchange potential of the form <i|K(i)|i> to the local exchange one
       ! call dscal(mxsize,-two*alphaf,wk3,ione)

       ! add an extra (-two*alphaf) portion of exchange energy 
       !call dscal (mxsize,-two*alphaf,wk4,ione)
       call dscal (mxsize,alphaf,wk4,ione)

       call add (mxsize,wk4,excp(length2-mxsize+1:))
       call add (mxsize,wk4,exchangeptr(ibpot:))    

       ! multiply coulomb/exchange potentials by f2
       call prod (mxsize,f2,wk4)
       call add (mxsize,wk4,fock2)
    endif
    
100 continue

    ! OK
    ! add one-electron contributions
    call dcopy (mxsize,f0,ione,wk3,ione)
    call daxpy (mxsize,ee(iorb,iorb),f1,ione,wk3,ione)
    
    if (mm(iorb).ne.0) then
       ! e enters the expression with minus sign which is already incorporated in e
       w=dble(mm(iorb)*mm(iorb))
       call daxpy (mxsize,w,e,ione,wk3,ione)
    endif
    
    ! add one-electron contributions to Coulomb and local exchange/correlation one
    call add (mxsize,wk3,fock1)
    
  end subroutine fockLXCunpol

  subroutine fockLXCpol(iorb)
    use params
    use detect
    use discrete
    use memory
    use scfshr
    use commons
    use blas
    use elocc
    use dftvxc
    use nabla
    use sharedMemory
    use utils

    use, intrinsic :: iso_c_binding
    use xc_f90_lib_m
    use libxc_funcs_m    
    
    implicit none
    integer (KIND=IPREC) :: i,ii,iborb,iborb1,ibpot,ibpot1,iorb,iorb1,ipc,isiorb1
    integer (KIND=IPREC) :: ibexp,idexp,kex,n,nan

    real (PREC) :: const13,ocdown,ocup,tmpf,w
    !real (PREC),dimension(*) :: psi,pot,excp,e,f0,f1,f2,f4,fock1,fock2
    !rhot,wk3,wk4,wk5,wk6,wk7,wk8,wk9,wk10,wk11,wk12,wk13

    real (PREC), dimension(:), pointer :: psi,excp,e,f0,f1,f2,f4,fock1,fock2
   
    !real (PREC), target :: work(length7) 
    real (PREC), dimension(:), pointer :: wk3,wk4,wk5,wk6,wk7,wk8,wk9,&
         wk10,wk11,wk12,wk13,work

    real (PREC), dimension(:), pointer :: rhoup,rhodown,rho,sigma,vrho,vsigma

    type(xc_f90_func_info_t) :: xc_info
    type(xc_f90_func_t) :: xc_func
    type(xc_f90_func_info_t) :: info
    integer(c_int) :: vmajor, vminor, vmicro, family_id, func_id, kind_id, err
    integer(c_size_t) :: size 
    character(len=120) :: name, kind1, family, ref
    real*8 :: xc_omega,xc_alpha,xc_beta
    integer (KIND=IPREC) :: allocstat

    e=>supplptr(i4b(4):)
    excp=>exchptr
    f0=>supplptr(i4b(5):)
    f1=>supplptr(i4b(6):)
    f2=>supplptr(i4b(7):)
    f4=>supplptr(i4b(9):)
    psi=>orbptr
    !work=>scratch4lxcptr
    work=>scratchptr

    fock1   => work(           1: 1*mxsize8)
    ! fock2==rhs in calling routine
    fock2   => work( 1*mxsize8+1: 2*mxsize8)
    rhoup   => work( 2*mxsize +1: 3*mxsize)
    rhodown => work( 3*mxsize +1: 4*mxsize)    
    rho     => work( 4*mxsize +1: 6*mxsize)
    vrho    => work( 6*mxsize +1: 8*mxsize)    
    sigma   => work( 9*mxsize +1:12*mxsize)
    vsigma  => work(12*mxsize +1:15*mxsize)
    wk3     => work(15*mxsize +1:16*mxsize8)
    wk4     => work(16*mxsize8+1:17*mxsize8)    
    wk5     => work(17*mxsize8+1:18*mxsize8)    
    wk6     => work(18*mxsize8+1:19*mxsize8)
    wk7     => work(19*mxsize8+1:20*mxsize8)
    wk8     => work(21*mxsize8+1:22*mxsize8)
    wk9     => work(22*mxsize8+1:23*mxsize8)
    wk10    => work(23*mxsize8+1:24*mxsize8)    
    wk11    => work(24*mxsize8+1:25*mxsize8)
    wk12    => work(25*mxsize8+1:26*mxsize8)
    wk13    => work(26*mxsize8+1:27*mxsize8)    

    iborb=i1b(iorb)
    ibpot=i2b(iorb)

    call zeroArray(mxsize,fock1)
    !call zeroArray(mxsize,fock2)
    
    if (nel==1) goto 100

    call zeroArray(mxsize,rho)
    !call zeroArray(mxsize,vrho)

    do iorb1=1,norb
       iborb1 =i1b (iorb1)
       call prodas (mxsize,occ(iorb1),psi(iborb1:),psi(iborb1:),rho)
    enddo

    ! choose a given functional from libxc library     
    size=mxsize

    ! XC_UNPOLARIZED branch goes first

    ! XC_POLARIZED version
    call zeroArray(mxsize,rhoup)
    call zeroArray(mxsize,rhodown)
    call zeroArray(2*mxsize,vrho)              
    call zeroArray(3*mxsize,vsigma)       
    
    do iorb1=norb,1,-1
       call exocc (iorb1,ocup,ocdown)
       call prodas (mxsize,ocup,  psi(i1b(iorb1):),psi(i1b(iorb1):),rhoup)
       call prodas (mxsize,ocdown,psi(i1b(iorb1):),psi(i1b(iorb1):),rhodown)
    enddo
    
    ! Libxc DFT functionals are implemented in C and expect 2*size rho
    ! array with spin-up and spin-down densities packed row-wise
    do i=1,mxsize
       rho(2*i-1)=rhoup(i)
       rho(2*i  )=rhodown(i)
    enddo
    
    call nfng (rhoup,rhoup,wk3,wk4,wk5,wk6,wk7,wk8,wk9,wk10)
    ! wk11=sigma[2]
    call nfng (rhoup,rhodown,wk3,wk4,wk5,wk6,wk7,wk8,wk9,wk11)
    ! wk12=sigma[3]
    call nfng (rhodown,rhodown,wk3,wk4,wk5,wk6,wk7,wk8,wk9,wk12)             
    
    do i=1,mxsize
       sigma(3*i-2)=wk10(i)
       sigma(3*i-1)=wk11(i)
       sigma(3*i  )=wk12(i)
    enddo
    
    do n=1,lxcFuncs
       call xc_f90_func_init(xc_func, lxcFuncs2use(n), XC_POLARIZED, err)
       call xc_f90_func_set_dens_threshold(xc_func, dens_threshold)
       info=xc_f90_func_get_info(xc_func)
       
       select case (xc_f90_func_info_get_family(info))
       case(XC_FAMILY_LDA)
          call xc_f90_lda_vxc(xc_func, size, rho, vrho)
          
          do i=1,mxsize
             wk13(i)=vrho(2*i-1)
          enddo
          
          call dcopy (mxsize,wk13,ione,fock1,ione)
          
       case(XC_FAMILY_GGA)
          ! calculate nabla rho nabla rho
          
          ! wk10=\naba \rho \nabla \rho: contracted gradients of the density
          
          ! wk12=vrho: first partial derivative of the energy per unit
          ! volume in terms of the density
          
          ! wk3 =vsigma: first partial derivative of the energy per unit
          ! volume in terms of sigma
          ! wk10=sigma[1]
          
          call xc_f90_gga_vxc(xc_func, size, rho, sigma, vrho, vsigma )
          
          ! It turns out that for for some GGA functionals xc_f90_gga_exc
          ! return NaN(s) value(s), and therefore the following piece of
          ! code must be used to detect them and fix via interpolation
          if (ldetectNaN) then
             nan=detectNaN(mxsize,wk3)
             if (nan>0) then
                write(*,'("Error: NaN detected in (XC_FAMILY_GGA) xc_f90_gga_vxc at",i5)') nan
             endif
          endif
          if (lfixNaN) then          
             call fixNaN(nni,mxnmu,wk3)
          endif
          
          ! w2=vrho
          ! w3(1:mxsize)=vsigma^{up,up}
          ! w3(mxsize+1,2*mxsize)=vsigma^{up,down}
          
          ! add vrho (up) to fock1
          do i=1,mxsize
             wk13(i)=vrho(2*i-1)
             !if (mod(i,1000)==0) print *,i,w4(i)
          enddo
          
          ! do i=1,10,1
          !    write(*,'(i6,4e16.8)') i, rho(2*i-1), sigma(3*i-2), vrho(2*i-1),vsigma(3*i-2) 
          ! enddo
          ! ! ! ! stop "pol"
          
          ! do i=1,5,1
          !    write(*,'(i6,4e16.8)') i, rho(2*i), sigma(3*i-1), vrho(2*i),vsigma(3*i-1) 
          ! enddo
          ! ! ! stop "pol"
          !call dcopy (mxsize,wk13,ione,fock1,ione)
          call daxpy  (mxsize,one,wk13,ione,fock1,ione)
          
          ! ! v_xc=vrho - 2 (nabla vsigma  nabla n + vsigma nabla^2 rho)
          ! ! wk11==(nabla vsignma  nabla n)
          ! ! extract vsigma (up,up)
          do i=1,mxsize
             wk12(i)=vsigma(3*i-2)
             !    !if (mod(i,1000)==0) print *,i,w4(i)
          enddo
          
          ! ! do i=1,mxsize,1
          ! !    write(*,'(i6,5e16.8)') i, rho(i), wk10(i), wk13(i),&
          ! !         wk12(i)
          ! ! enddo
          
          ! !stop "pol"
          
          call nfng (wk12,rhoup,wk4,wk5,wk6,wk7,wk8,wk9,wk10,wk11)
          
          ! ! wk8==nabla^2 rho
          call n2f (rhoup,wk5,wk6,wk7,wk8)
          
          ! wk8==vsigma nabla^2 rho
          call prod (mxsize,wk12,wk8)
          
          ! wk8==(nabla vsignma  nabla n) + vsigma nabla^2 rho
          call add (mxsize,wk11,wk8)
          
          ! ! fock1==vrho - 2 (nabla vsignma  nabla n) + vsigma nabla^2 rho
          call daxpy (mxsize,-two,wk8,ione,fock1,ione)
          
          ! extract vsigma (up,down)
          do i=1,mxsize
             wk12(i)=vsigma(3*i-1)
          enddo
          
          ! wk11==(nabla vsignma'  nabla n')
          call nfng (wk12,rhodown,wk4,wk5,wk6,wk7,wk8,wk9,wk10,wk11)
          
          ! wk8==nabla^2 rho
          call n2f (rhodown,wk5,wk6,wk7,wk8)
          
          ! wk8==vsigma nabla^2 rho
          call prod (mxsize,wk12,wk8)
          
          ! wk8==(nabla vsignma  nabla n) + vsigma nabla^2 rho
          call add (mxsize,wk11,wk8)
          call daxpy (mxsize,-one,wk8,ione,fock1,ione)
          
       case default
          write(*,'("Error! Unsupported family of libxc functionals.")')
          stop 'fockLXC'
       end select
       call xc_f90_func_end(xc_func)
    enddo

    ! take care of F4 factor (r\xi/2) which transforms the potential into \tilde form
    call multf4(fock1)

    !alphaf=zero
    ! let's simulate HF case. no LXC functional
    if (idebug(999)==1) then
       call zeroArray(mxsize,fock1)    
       alphaf=one
    endif
    
    ! store the local exchange/correlation potential in extra portion
    ! of pot array (to be used by EaDFT and EabDFT)
    !call dcopy (mxsize,fock1,ione,pot(length2-2*mxsize+1),ione)
    !call dcopy (mxsize,fock1,ione,excp(length2-2*mxsize+1:),ione)
    call dcopy (mxsize,fock1,ione,coulombptr(i2b(iorb):),ione)

    call zeroArray(mxsize,wk3)       
    call zeroArray(mxsize,fock2)    

    do iorb1=1,norb
       if (inDFT(iorb1).eq.0) cycle
       iborb1=i1b(iorb1)
       ibpot1=i2b(iorb1)
       !ipc=iorb1+norb*(iorb-1)
       
       ! add the Coulomb potential which in DFT also includes the
       ! contribution from the orbital
       !call daxpy (mxsize,occ(iorb1),pot(ibpot1),ione,wk3,ione)
       call daxpy (mxsize,occ(iorb1),excp(ibpot1:),ione,wk3,ione)
       
       ! add contributions from off-diagonal Lagrange multipliers.
       if (iorb.ne.iorb1.and.abs(ee(iorb1,iorb))>epsilon(zero)) then
          do i=1,mxsize
             fock2(i)=fock2(i)+ee(iorb1,iorb)*f4(i)*psi(iborb1+i-1)
          enddo
       endif
    enddo
    ! add Coulomb potential to fock1
    call add (mxsize,wk3,fock1)
    call dcopy (mxsize,fock1,ione,coulombptr(i2b(iorb):),ione)
    
    ! do not add wk3 to pot(length2-2*mxsize+1) since this
    ! contribution is taken care of by EaDFT/EabDFT itself 
    
    ! multiply fock2 by f2 and save for EaDFT

    call prod (mxsize,f2,fock1)
    
    !call dcopy (mxsize,fock2,ione,pot(length2-mxsize+1),ione)
    call dcopy (mxsize,fock2,ione,excp(length2-mxsize+1:),ione)
    call dcopy (mxsize,fock2,ione,exchangeptr(ibpot:),ione)
    
    call prod (mxsize,f2,fock2) 
    
    if (lxcHyb) then
       write(*,'("Warning! Hybrid polarized functionals not supported yet.")')
       stop "fockLXCpol"

       call zeroArray(mxsize,wk3)        
       ! add contributions from Coulomb and off-diagonal Lagrange multipliers.
       ibpot=i2b(iorb)
       ! in the local exchange approximation the Coulomb potential also includes the
       ! contribution from the orbital
       !call dcopy (mxsize,pot(ibpot),ione,wk3,ione)
       call dcopy (mxsize,excp(ibpot:),ione,wk3,ione)
       
       ! add exchange potential of the form <i|K(i)|i> to the local exchange one
       !call dscal(mxsize,-two*alphaf,wk3,ione)
       call dscal(mxsize,-alphaf,wk3,ione)
       !!! call add (mxsize,wk3,excp(length2-2*mxsize+1:))
       call prod (mxsize,f2,wk3)
       call add (mxsize,wk3,fock1)
       
       ! substract a part of exchange contributions 
       if (norb>1) then
          call zeroArray(mxsize,wk4)
          ! add contributions from exchange potentials times -two*alphaf
          do iorb1=1,norb
             iborb1=i1b(iorb1)
             kex=iorb+norb*(iorb1-1)
             ipc=iorb1+norb*(iorb-1)
             
             if (iorb.le.iorb1) then
                idexp=iorb+iorb1*(iorb1-1)/2
             else
                idexp=iorb1+iorb*(iorb-1)/2
             endif
             
             ibpot1=i2b(iorb1)
             ibexp=i3b(idexp)
             if (iorb1.ne.iorb) then
                call dcopy (mxsize,excp(ibexp:),ione,wk3,ione)
                call dscal (mxsize,gec(kex),wk3,ione)
                if (ilc(idexp).gt.1) then
                   call daxpy (mxsize,gec(kex+norb*norb),excp(ibexp+mxsize:),ione,wk3,ione)
                endif
                
                call prod (mxsize,psi(ibpot1:),wk3)
                call add  (mxsize,wk3,wk4)
             else
                if ((mm(iorb).gt.0).and.(ilc(idexp).gt.0)) then
                   call dcopy (mxsize,excp(ibexp:),ione,wk3,ione)
                   call prod  (mxsize,psi(iborb1:),wk3)
                   call daxpy (mxsize,gec(kex),wk3,ione,wk4,ione)
                endif
             endif
          enddo
          
          ! add an extra (-alphaf) portion of exchange energy 
          call dscal (mxsize,-alphaf,wk4,ione)
          call add (mxsize,wk4,excp(length2-mxsize+1:))
          call add (mxsize,wk4,exchangeptr(ibpot:))

          ! multiply coulomb/exchange potentials by f2
          call prod (mxsize,f2,wk4)
          call add (mxsize,wk4,fock2)
       endif
    endif
 
100 continue
    
    ! add one-electron contributions
    call dcopy (mxsize,f0,ione,wk3,ione)
    call daxpy (mxsize,ee(iorb,iorb),f1,ione,wk3,ione)
    
    if (mm(iorb).ne.0) then
       ! e enters the expression with minus sign which is already incorporated in e
       w=dble(mm(iorb)*mm(iorb))
       call daxpy (mxsize,w,e,ione,wk3,ione)
    endif
    
    ! add one-electron contributions to Coulomb and local exchange/correlation one
    call add (mxsize,wk3,fock1)

  end subroutine fockLXCpol
#endif  

  ! ### fockDFT ###
  !
  !     Calculates the exchange potentials as the local Slater approximation.
  !
  subroutine fockDFT(iorb)
    use params
    use discrete
    use memory
    use scfshr
    use commons
    use blas
    use elocc
    use dftvxc
    use nabla
    use sharedMemory
    use utils
    
    implicit none

    
    integer (KIND=IPREC) :: i,iborb,iborb1,ibpot,ibpot1,iorb,iorb1,ipc
    real (PREC) :: const13,ocdown,ocup,w
    parameter (const13=1.0_PREC/3.0_PREC)
    real (PREC), dimension(:), pointer :: psi,excp,e,f0,f1,f2,f4,&
              fock1,fock2,wk2,wk3,wk4,wk5,wk6,wk7,wk8,wk9,wk10,wk11,wk12,wk13

    e=>supplptr(i4b(4):)
    excp=>exchptr
    f0=>supplptr(i4b(5):)
    f1=>supplptr(i4b(6):)
    f2=>supplptr(i4b(7):)
    f4=>supplptr(i4b(9):)
    psi=>orbptr

    ! fock2==rhs in calling routine
    fock1=>scratchptr(          1: 1*mxsize8)
    fock2=>scratchptr(1*mxsize8+1: 2*mxsize8)    
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

    ! contributions from one-electron terms
    iborb=i1b(iorb)
    ibpot=i2b(iorb)

    call zeroArray(mxsize,fock1)
    call zeroArray(mxsize,fock2)

    ! store wk0 in fock2 for a time being
    !    call dcopy (mxsize,wk0,ione,fock2,ione)
    ! use fock2 as wk0

    ! contributions from the local exchange approximation

    ! if (idebug(707)==1.or.idebug(807)==1) then
    !    call prodas (mxsize,one, psi(iborb),psi(iborb),fock2)
    !    !call test4n2f (psi(iborb),wk3,wk4,wk5,wk6,wk7)
    !    !call test4n2f (fock2,wk3,wk4,wk5,wk6,wk7)
    !    call printn2f (fock2,wk3,wk4,wk5,wk6,wk7)
    !    stop "fockDFT"
    ! endif

    !if (nel.gt.1) then
    if (nel.ge.1) then
       call zeroArray(mxsize,wk2)
       call zeroArray(mxsize,fock2)
       
       ! store wk0 in fock2 for a time being call dcopy (mxsize,wk0,ione,fock2,ione) use
       ! fock2 as wk0

       ! contributions from the local exchange approximation
       select case (idftex) 
       case (1) ! LDA
          do iorb1=1,norb
             if (inDFT(iorb1).eq.0) cycle

             iborb1=i1b (iorb1)
             call exocc (iorb1,ocup,ocdown)
             call prodas (mxsize,ocup,psi(iborb1:),psi(iborb1:),wk2)
          enddo
          do i=1,mxsize
             wk2(i)=(wk2(i))**const13
          enddo
          call dscal(mxsize,fdftpot(alphaf),wk2,ione)
          
          ! multiply the local exchange potential by f4 to make it commensurate with the
          ! Coulomb potential
          call multf4(wk2)

          call dcopy (mxsize,wk2,ione,excp(length2-2*mxsize+1:),ione)
          call dcopy (mxsize,wk2,ione,coulombptr(ibpot:),ione)
       case (2) ! B88
          ! call fbe88(psi,f4,wk8,wk9,wk10,wk11,wk12,wk13,excp,fock2,wk2,wk3,wk4,wk5,wk6,wk7)
          call fbe88(psi,f4,wk8,wk9,wk10,wk11,wk12,wk13,excp(length2-2*mxsize+1:),&
               fock2,wk2,wk3,wk4,wk5,wk6,wk7)
       case default
          if (idftex/=0) then
             write(*,'(/"Warning: unsupported exchange potential")')
             stop 'fockDFT'
          endif
       end select

       ! store the local exchange potential in a separate array
       call dcopy (mxsize,wk2,ione,fock1,ione)

       ! add contributions from correlations potentials
       call zeroArray(mxsize,wk2)
       select case (idftcorr)
       case (1) ! LYP
          !call flypcs(psi,f4,wk8,wk9,wk10,wk11,wk12,wk13,excp,fock2,wk2,wk3,wk4,wk5,wk6,wk7)
          call flypcs(psi,f4,wk8,wk9,wk10,wk11,wk12,wk13,excp(length2-2*mxsize+1:),&
               fock2,wk2,wk3,wk4,wk5,wk6,wk7)
       case (2) ! VWN
          call fvwncs(psi,f4,wk8,wk9,wk10,wk11,wk12,wk13,excp(length2-2*mxsize+1:),&
               fock2,wk2,wk3,wk4,wk5,wk6,wk7)
       case default
          if (idftcorr/=0) then
             write(*,'(/"Warning: unsupported correlation potential")')
             stop 'fockDFT'
          endif
       end select
       ! add correlation contributions to the local exchange ones (f4 factor is already included)
       call add (mxsize,wk2,fock1)
    endif
    
    ! add local exchange and correlations contributions ones to one-electron contribution
    ! call dcopy (mxsize,wk2,ione,fock1,ione)
    ! call dcopy (mxsize,fock2,ione,wk0,ione)
    
    call zeroArray(mxsize,wk2)
    call zeroArray(mxsize,fock2)
    
    if (nel.gt.1) then
       ! add contributions from Coulomb and off-diagonal Lagrange multipliers.
       
       do iorb1=1,norb
          if (inDFT(iorb1).eq.0) cycle
          iborb1=i1b(iorb1)
          ibpot1=i2b(iorb1)
          !          ipc=iorb1+norb*(iorb-1)
          ! in the local exchange approximation the Coulomb potential also includes the
          ! contribution from the orbital
          call daxpy (mxsize,occ(iorb1),excp(ibpot1:),ione,wk2,ione)
          if (iorb.ne.iorb1.and.abs(ee(iorb1,iorb))>epsilon(zero)) then
             do i=1,mxsize
                fock2(i)=fock2(i)+ee(iorb1,iorb)*f4(i)*psi(iborb1+i-1)
             enddo
          endif
       enddo
       ! store the local exchange potential in exch array as its not used in HFS/DFT (to
       ! be used by EaDFT and EabDFT) at the end of excp array (important when scmc is on)
       !call dcopy (mxsize,fock1,ione,excp(length3-mxsize),ione)
       call add (mxsize,wk2,fock1)
       call dcopy (mxsize,fock1,ione,excp(length2-2*mxsize+1:),ione)
       call dcopy (mxsize,fock1,ione,coulombptr(ibpot:),ione)       
       ! add the coulomb potential to the local exchange one


       ! multiply coulomb/exchange potentials and off-diagonal Lagrange
       ! multipliers by f2
       call prod (mxsize,f2,fock1)

       call dcopy (mxsize,fock2,ione,excp(length2-mxsize+1:),ione)
       call dcopy (mxsize,fock2,ione,exchangeptr(ibpot:),ione)              

       call prod (mxsize,f2,fock2)
       ! nel.gt.1
    endif

    call dcopy (mxsize,f0,ione,wk2,ione)
    call daxpy (mxsize,ee(iorb,iorb),f1,ione,wk2,ione)

    if (mm(iorb).ne.0) then
       !      e enter the expression with minus sign which is already incorporated in e
       w=dble(mm(iorb)*mm(iorb))
       call daxpy (mxsize,w,e,ione,wk2,ione)
    endif
    
    ! add Coulomb contributions to one-electron one (containing local
    ! exchange and correlation contributions)
    
    call add (mxsize,wk2,fock1)

  end subroutine fockDFT
  
end module fock
