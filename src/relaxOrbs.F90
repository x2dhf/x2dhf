! SPDX-License-Identifier: GPL-2.0-or-later

! Copyright (C) 1996       Leif Laaksonen, Dage Sundholm               
! Copyright (C) 1996-2024  Jacek Kobus 

module relaxOrbs
  implicit none
contains

  ! ### orbSOR ###
  !
  !     Evaluates the Fock potential for a given orbital, sets up the
  !     right-hand side of the Poisson equation for that orbital and
  !     performs a few SOR iterations.
  !
  subroutine orbSOR (iorb)
    use params
    use discrete
    use scfshr
    use solver
    use commons
    use fock
    use inout
    use sormcsor
    use sharedMemory
    use utils
    
    implicit none
    integer (KIND=IPREC) :: i,iborb,iborb1,ii,ij,ioffs1,ioffst,iorb,isym,itr1,itr2,j,nmut
    
    !integer (KIND=IPREC),dimension(*) :: cw_sor
    !real (PREC), dimension(*) :: psi,pot,excp,b,d,e,f0,f1,f2,f4,wgt2,lhs,rhs, &
    !     wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7,wk8,wk9,fock1,fock2
    !real (PREC), dimension(*) :: work

    integer (KIND=IPREC),dimension(:), pointer :: cw_sor    
    real (PREC), dimension(:), pointer :: psi,b,d,f1,fock1,lhs,rhs,wk1,wk2,wk3

    real (PREC), dimension(:), pointer :: wk4    
    integer (KIND=IPREC) :: numax,mumax
    real (PREC) :: fcurr,fprev,fmax

    
    b=>supplptr(i4b(1):)
    d=>supplptr(i4b(3):)    
    f1=>supplptr(i4b(6):)
    psi=>orbptr
    cw_sor=>sorptr    

    ! rhs==fock2 in fockLXC routine
    fock1=>scratchptr(           1: 1*mxsize8)
    rhs  =>scratchptr( 1*mxsize8+1: 2*mxsize8)    
    wk1  =>scratchptr( 2*mxsize8+1: 3*mxsize8)
    wk2  =>scratchptr( 3*mxsize8+1: 4*mxsize8)
    wk3  =>scratchptr( 4*mxsize8+1: 5*mxsize8)
    lhs  =>scratchptr( 5*mxsize8+1: 6*mxsize8)
    wk4  =>scratchptr( 6*mxsize8+1: 7*mxsize8)    

    ! prepare right and left-hand side of the Poisson equation include
    ! the diagonal part of the differentiation operator in lhs

    !if (ifixorb(iorb)==1) return

    if (HF.or.OED) call fockHF (iorb)
    if (TED) call fockTED (iorb)

    ! use self-coded exchange-correlation functionals
    if (DFT.or.HFS.or.SCMC) call fockDFT(iorb)

#ifdef LIBXC
    ! use exchange-correlation functionals from libxc
    if (LXC) then
       if (lxcPolar) then
          call fockLXCpol(iorb)
       else
          call fockLXCunpol(iorb)
       endif
    endif
#endif    

    iborb=i1b(iorb)
    omega=ovforb
    omega1=1.0_PREC-omega
    isym=isymOrb(iorb)
    
    maxsor1=maxsororb(1)    
    maxsor2=maxsororb(2)
    
    nsor4orb=nsor4orb+maxsor1*maxsor2
    
    !call zeroArray(mxsize8,wk2)
    do itr1=1,maxsor1
       do i=1,mxsize
          lhs(i)=fock1(i)+diag
       enddo

       ! Initialize edecay array needed to evaluate orbital values
       ! in the asymptotic region after relaxation. For every
       ! 2<=nu<=nni-1 orbital values at 2<=mu<nmu-4 are
       ! relaxed. Upon relaxation the values in the unrelaxed tail
       ! region are calculated using the known exponential decay
       ! behaviour of the function
       
       do i=mxnmu-4,mxnmu
          ii=(i-1)*nni
          do j=1,nni
             ij=ii+j
             if (abs(f1(ij)).le.precis) then
                wk1(ij)=0.0_PREC
             else
                wk1(ij)=fock1(ij)/f1(ij)
             endif
          enddo
       enddo
       
       call orbAsymptGet (wk3,wk1)
       call putin (nni,mxnmu,isym,psi(iborb:),wk2)
       wk4=wk2
       call sor (isym,wk2,lhs,rhs,b,d,          &
            cw_sor(iadext(1):),cw_sor(iadnor(1):),&
            cw_sor(iadex1(1):),cw_sor(iadex2(1):),&
            cw_sor(iadex3(1):))
       call putout (nni,mxnmu,psi(iborb:),wk2)
       call orbAsymptSet (psi(iborb:),wk3)

    fmax=zero
    numax=0
    mumax=0
    do i=1,mxsize8
    !   if (mod(i,5000)==0) print *,i,wk2(i),wk3(i)
       if (abs(wk2(i)-wk4(i))>fmax) then
          fmax=abs(wk2(i)-wk4(i))
          numax=indx1nu(i)
          mumax=indx1mu(i)
       endif
    enddo
    !print *,"orbSOR:    ",numax,mumax,fmax

    enddo


  end subroutine orbSOR

  ! ### orbMCSOR ###
  !
  !     Evaluates the Fock potential for a given orbital, sets up the
  !     right-hand side of the Poisson equation for that orbital and
  !     performs a few MCSOR iterations.
  !
  subroutine orbMCSOR (iorb)
    use commons
    use discrete
    use fock
    use sormcsor
    use params
    use inout
    use scfshr
    use sharedMemory
    use solver
    implicit none
    integer (KIND=IPREC) :: i,iborb,ii,ij,iorb,isym,itr1,itr2,j,nmut
    
    integer (KIND=IPREC),dimension(:), pointer :: cw_sor    
    real (PREC), dimension(:), pointer :: psi,b,d,f1,fock1,lhs,rhs,wk1,wk2,wk3

    b=>supplptr(i4b(1):)
    d=>supplptr(i4b(3):)    
    f1=>supplptr(i4b(6):)
    psi=>orbptr
    cw_sor=>sorptr

    ! rhs==fock2 in fockLXC routine
    fock1=>scratchptr(           1: 1*mxsize8)
    rhs  =>scratchptr( 1*mxsize8+1: 2*mxsize8)    
    wk1  =>scratchptr( 2*mxsize8+1: 3*mxsize8)
    wk2  =>scratchptr( 3*mxsize8+1: 4*mxsize8)
    wk3  =>scratchptr( 4*mxsize8+1: 5*mxsize8)
    lhs  =>scratchptr( 5*mxsize8+1: 6*mxsize8)

    if (ifixorb(iorb)==1) return

    ! prepare right and left-hand side of the Poisson equation
    ! include the diagonal part of the differentiation operator in lhs

    if (HF.or.OED) call fockHF (iorb)
    if (TED) call fockTED (iorb)
    
    ! use self-coded exchange-correlation functionals
    if (DFT.or.HFS.or.SCMC) call fockDFT(iorb)

#ifdef LIBXC
    ! use exchange-correlation functionals from libxc
    if (LXC) then
       if (lxcPolar) then
          call fockLXCpol(iorb)
       else
          call fockLXCunpol(iorb)
       endif
    endif
#endif    

    iborb=i1b(iorb)
    omega=ovforb
    omega1=1.0_PREC-omega
    isym=isymOrb(iorb)

    maxsor1=maxsororb(1)    
    maxsor2=maxsororb(2)
    nsor4orb=nsor4orb+maxsor1*maxsor2
    
    do itr1=1,maxsor1
       do i=1,mxsize
          lhs(i)=fock1(i)+diag
       enddo
       
       ! initialize aorb array needed to evaluate orbital values in
       ! the asymptotic region

       do i=mxnmu-4,mxnmu
          ii=(i-1)*nni
          do j=1,nni
             ij=ii+j
             if (abs(f1(ij)).le.precis) then
                wk1(ij)=0.0_PREC
             else
                wk1(ij)=fock1(ij)/f1(ij)
             endif
          enddo
       enddo
       
       call orbAsymptGet (wk3,wk1)
       call putin (nni,mxnmu,isym,psi(iborb:),wk2)
       call mcsor(isym,wk2,lhs,rhs,b,d,         &
            cw_sor(iadext(1):),cw_sor(iadnor(1):),&
            cw_sor(iadex1(1):),cw_sor(iadex2(1):),&
            cw_sor(iadex3(1):))
       call putout (nni,mxnmu,psi(iborb:),wk2)
       call orbAsymptSet (psi(iborb:),wk3)
    enddo
    
  end subroutine orbMCSOR

#if ( defined PTHREAD || defined TPOOL )
  
  ! ### orbMCSORPT ###
  !
  !     Evaluates the Fock potential for a given orbital, sets up the
  !     right-hand side of the Poisson equation for that orbital and
  !     performs a few MCSOR iterations.
  !
  subroutine orbMCSORPT (iorb)
    use commons
    use detect
    use discrete
    use fock
    use sormcsor
    use params
    use inout
    use scfshr
    use sharedMemory
    use solver
    
    implicit none
    integer (KIND=IPREC) :: i,iborb,iborb1,ii,ij,ioffs1,ioffst,iorb,isym,itr1,itr2,j,&
         nmut,nan
    integer (KIND=IPREC) :: mxnmuc,mxsizec,mxsize8c,ngrid6ac,ngrid6bc,ngrid7c,&
         nnu1c,nnu2c,nnu3c,nnu4c,nnu5c,isstartc,isstopc,maxsor1c,maxsor2c
    real (PREC) omegac,omega1c

    integer (KIND=IPREC),dimension(:), pointer :: cw_sor    
    real (PREC), dimension(:), pointer :: psi,b,d,f1,fock1,lhs,rhs,wk1,wk2,wk3

    ! interface
    !    subroutine mcsorptdrv(cisym,cnthreads,cb,cd,cwk2,clhs,crhs,&
    !         cw_sor1,cw_sor2,cw_sor3,cw_sor4,cw_sor5) bind(c, name='mcsorptdrv_')
    !      use, intrinsic :: iso_c_binding
    !      integer(c_int), intent(in) :: cisym,cnthreads
    !      integer(c_INT), intent(in), dimension(:) :: cw_sor1,cw_sor2,cw_sor3,cw_sor4,cw_sor5
    !      real (kind=8), intent(in), dimension(:) :: cb,cd,clhs,crhs
    !      real (kind=8), intent(out), dimension(:) :: cwk2
    !    end subroutine mcsorptdrv
    ! end interface

    !external mcsorptdrv

#ifdef PTHREAD    
    interface
       subroutine mcsor_pthread (isym,nthreads,wk2,lhs,rhs) bind(c,name='mcsor_pthread_')
         use, intrinsic :: iso_c_binding, only: c_int,c_double,c_ptr
         integer(c_int), intent(in)  :: isym,nthreads
         real (c_double) :: wk2(*),lhs(*),rhs(*)         
       end subroutine mcsor_pthread
    end interface
#endif

#ifdef TPOOL    
    interface
       subroutine mcsor_tpool (isym,nthreads,wk2,lhs,rhs) bind(c,name='mcsor_tpool_')
         use, intrinsic :: iso_c_binding, only: c_int,c_double,c_ptr
         integer(c_int), intent(in)  :: isym,nthreads
         !real(c_ptr) :: wk2(1),lhs(1),rhs(1)
         real (c_double) :: wk2(*),lhs(*),rhs(*)         
       end subroutine mcsor_tpool
    end interface
#endif
    
!    external mcsor_pthread,mcsor_tpool
    
    common /c_interface_17/ mxnmuc,mxsizec,ngrid6ac,ngrid6bc,ngrid7c,&
         nnu1c,nnu2c,nnu3c,nnu4c,nnu5c,isstartc,isstopc,maxsor1c,maxsor2c

    common /c_interface_45/ omegac,omega1c

    b=>supplptr(i4b(1):)
    d=>supplptr(i4b(3):)    
    f1=>supplptr(i4b(6):)
    psi=>orbptr
    cw_sor=>sorptr

    ! rhs==fock2 in fockLXC routine    
    fock1=>scratchptr(           1: 1*mxsize8)
    rhs  =>scratchptr( 1*mxsize8+1: 2*mxsize8)    
    wk1  =>scratchptr( 2*mxsize8+1: 3*mxsize8)
    wk2  =>scratchptr( 3*mxsize8+1: 4*mxsize8)
    wk3  =>scratchptr( 4*mxsize8+1: 5*mxsize8)
    lhs  =>scratchptr( 5*mxsize8+1: 6*mxsize8)

    if (ifixorb(iorb)==1) return

    if (HF.or.OED) call fockHF (iorb)
    if (TED) call fockTED (iorb)
    
    ! use self-coded exchange-correlation functionals
    if (DFT.or.HFS.or.SCMC) call fockDFT(iorb)

    
    ! prepare right and left-hand side of the Poisson equation
    ! include the diagonal part of the differentiation operator in lhs
#ifdef LIBXC           
    if (LXC) then
       if (lxcPolar) then
          call fockLXCpol(iorb)
       else
          call fockLXCunpol(iorb)
       endif
    endif
#endif    

    isym=isymOrb(iorb)
    iborb=i1b(iorb)

    omega=ovforb
    maxsor1=maxsororb(1)    
    maxsor2=maxsororb(2)
    nsor4orb=nsor4orb+maxsor1*maxsor2
    
    ! maxsor1c=maxsororb(1)    
    ! maxsor2c=maxsororb(2)
    ! omega1c=1.0_PREC-omegac

    do itr1=1,maxsor1
       do i=1,mxsize
          lhs(i)=fock1(i)+diag
          !if (mod(i,100)==0) write(*,'(i5,e16.8)') i,lhs(i)
       enddo
       ! initialize aorb array needed to evaluate orbital values in the asymptotic region
       do i=mxnmu-4,mxnmu
          ii=(i-1)*nni
          do j=1,nni
             ij=ii+j
             if (abs(f1(ij)).le.precis) then
                wk1(ij)=0.0_PREC
             else
                wk1(ij)=lhs(ij)/f1(ij)
             endif
          enddo
       enddo

       if (ldetectNaN) then
          nan=detectNaN(mxsize,psi(iborb:))
          if (nan>0) then
             write(*,'("Error: NaN detected in orbMCSOR-a: psi at",i5)') nan
          endif
       endif

       call orbAsymptGet (wk3,wk1)
       call putin (nni,mxnmu,isym,psi(iborb:),wk2)
       
#if ( defined PTHREAD )
#ifdef TRACE
       write(*,'("TRACE:  orbMCSORPT/iorb ",i3)') iorb
#endif
       call mcsor_pthread(isym,nthreads,wk2,lhs,rhs)
#elif ( defined TPOOL )
       call mcsor_tpool(isym,nthreads,wk2,lhs,rhs)
#endif       
       call putout (nni,mxnmu,psi(iborb:),wk2)
       call orbAsymptSet (psi(iborb:),wk3)

       if (ldetectNaN) then
          nan=detectNaN(mxsize,psi(iborb:))
          if (nan>0) then
             write(*,'("Error: NaN detected in orbMCSOR-b: psi at",i5)') nan
          endif
       endif
    enddo
  end subroutine orbMCSORPT
#endif

  ! ### orbAsymptGet ###
  !
  !   Initializes array edecay which is used by orbAsymptSet to calculate
  !   boundary values of a given orbital at practical infinity.
  !
  subroutine orbAsymptGet (edecay,fa)
    use params
    use discrete
    use commons

    implicit none
    integer (KIND=IPREC) :: i,j,jj
    real (PREC) :: raiq,raiq1
    real (PREC), dimension(*) :: fa
    real (PREC), dimension(nni,4) :: edecay

    jj=0
    do j=mxnmu-3,mxnmu
       jj=jj+1
       do i=1,nni
          raiq1=sqrt(vxisq(j-1)+vetasq(i)-1.0_PREC)
          raiq =sqrt(vxisq(j)+vetasq(i)-1.0_PREC)
          edecay(i,jj)=(abs(fa(i+(j-1)*nni)/fa(i+(j-1)*nni-nni)))**0.250_PREC&
               *exp(sqrt(abs(fa(i+(j-1)*nni-nni)))*(raiq1-raiq))
       enddo
    enddo

  end subroutine orbAsymptGet

  ! ### orbAsymptSet ###
  !
  !     Recalculates asymptotic values of a given orbital at the practical
  !     infinity using exponential decay values prepared by calling orbAsymptDet.
  !
  subroutine orbAsymptSet (psi,edecay)
    use params
    use discrete
    use commons

    implicit none
    integer (KIND=IPREC) :: i,j,jj
    real (PREC), dimension(*) :: psi
    real (PREC), dimension(nni,4) :: edecay

    jj=0
    do j=mxnmu-3,mxnmu
       jj=jj+1
       do i=1,nni
          psi(i+(j-1)*nni)=psi(i+(j-1)*nni-nni)*edecay(i,jj)
       enddo
    enddo

  end subroutine orbAsymptSet
  
end module relaxOrbs
