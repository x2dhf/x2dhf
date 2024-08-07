! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2023  Jacek Kobus 

module relaxPots
  implicit none
contains
  subroutine coulExchSOR (iorb)
    use blas
    use commons
    use discrete
    use inout
    use params
    use sormcsor
    use scfshr
    use sharedMemory
    use solver
    use utils

#ifdef OPENMP    
    use omp_lib
#endif
    implicit none
    integer (KIND=IPREC) :: i,ib1,ib2,ibexp,idel,idel1,ifirst,in1,in2,ioffs1,&
         iorb,iorb1,isym,iswtch,itr1,itr2,nexchpot
    integer (KIND=IPREC),dimension(:), pointer :: cw_sor
    real (PREC), dimension(:), pointer :: psi,excp,bpot,d,e,f3,g
    
#ifdef OPENMP
    real (PREC), dimension(:), allocatable :: lhs,rhs,wk2
#else
    real (PREC), dimension(:), pointer :: lhs,rhs,wk2
    lhs =>scratchptr(          1:  1*mxsize8)
    rhs =>scratchptr(   mxsize8+1: 2*mxsize8)
    wk2 =>scratchptr( 2*mxsize8+1: 3*mxsize8)
#endif

    d=>supplptr(i4b(3):)
    e=>supplptr(i4b(4):)
    excp=>exchptr
    bpot=>supplptr(i4b(2):)
    f3=>supplptr(i4b(8):)
    g=>supplptr(i4b(10):)
    psi=>orbptr

    cw_sor=>sorptr
    
    omega=ovfexch
    omega1=1.0_PREC-omega
    maxsor1=maxsorpot(1)    
    maxsor2=maxsorpot(2)
    
    !$OMP PARALLEL NUM_THREADS(nexchpots(iorb)) DEFAULT(SHARED) PRIVATE(nexchpot,isym,i,ib1,ib2,in1,in2,ibexp,idel,itr1,lhs,rhs,wk2)
#ifdef OPENMP

    i=OMP_get_thread_num()
    nexchpot=i+1
#else
    do nexchpot=1,nexchpots(iorb)
#endif
#ifdef OPENMP
       allocate(lhs(mxsize))
       allocate(rhs(mxsize))
       allocate(wk2(mxsize8))
#endif    
       in1=ins1(iorb,nexchpot)
       in2=ins2(iorb,nexchpot)
       
       ib1=i1b(in1)
       ib2=i1b(in2)
       
       !ipc=ipcs(iorb,nexchpot)
       ibexp=ibexcp(iorb,nexchpot)
       
       idel=deltam4pot(iorb,nexchpot)
       isym=isyms(iorb,nexchpot)
       
       ! prepare the right-hand side of the Poisson equation
       
       call prod2 (mxsize,psi(ib1:),psi(ib2:),wk2)
       call prod2 (mxsize,wk2,g,rhs)
       
       ! prepare left-hand side of the poisson equation include the diagonal part of
       ! the differentiation operator in lhs

       if (idel==0) then
          do i=1,mxsize
             lhs(i)=f3(i)
          enddo
       else
          do i=1,mxsize
             lhs(i)=f3(i)+dble(idel*idel)*e(i)
          enddo
       endif
       
       do itr1=1,maxsor1
          call putin (nni,mxnmu,isym,excp(ibexp:),wk2)
          call sor (isym,wk2,lhs,rhs,bpot,d,          &
               cw_sor(iadext(1):),cw_sor(iadnor(1):), &
               cw_sor(iadex1(1):),cw_sor(iadex2(1):), &
               cw_sor(iadex3(1):) )
          call putout (nni,mxnmu,excp(ibexp:),wk2)
       enddo
       
#ifdef OPENMP
       deallocate(lhs)
       deallocate(rhs)
       deallocate(wk2)
#else
    enddo
#endif
    !$OMP END PARALLEL  
  end subroutine coulExchSOR

#if ( defined PTHREAD || defined TPOOL )
  subroutine coulExchSORPT (iorb)
    use commons
    use discrete
    use inout
    use params
    use sormcsor
    use scfshr
    use sharedMemory
    use solver
    use utils
    use, intrinsic :: iso_c_binding
    implicit none
    integer (KIND=IPREC) :: i,ib1,ib2,ibexp,idel,idel1,ifirst,in1,in2,ioffs1, &
         iorb,iorb1,isym,iswtch,itr1,itr2,nexchpot
    integer (KIND=IPREC),dimension(:), pointer :: cw_sor
    real (PREC), dimension(:), pointer :: psi,excp,bpot,d,e,f3,g
    
    integer (KIND=IPREC) :: mxnmuc,mxsizec,mxsize8c,ngrid6ac,ngrid6bc,ngrid7c,&
         nnu1c,nnu2c,nnu3c,nnu4c,nnu5c,isstartc,isstopc,maxsor1c,maxsor2c
    real (PREC) omegac,omega1c
    common /c_interface_17/ mxnmuc,mxsizec,ngrid6ac,ngrid6bc,ngrid7c,&
         nnu1c,nnu2c,nnu3c,nnu4c,nnu5c,isstartc,isstopc,maxsor1c,maxsor2c
    common /c_interface_45/ omegac,omega1c

    integer (KIND=IPREC) :: iorbc,isymc,nthreadsc
    common /c_interface_46/ iorbc,isymc,nthreadsc
    !integer (KIND=IPREC),pointer :: cw_sor4p
    !common /c_interface_47/ cw_sor4p

#ifdef PTHREAD    
    interface
       subroutine coulExch_pthread () bind(c,name='coulExch_pthread_')
         use, intrinsic :: iso_c_binding, only: c_int
       end subroutine coulExch_pthread
    end interface
#endif
    
#ifdef TPOOL
    interface
       subroutine coulExch_tpool () bind(c,name='coulExch_tpool_')
         use, intrinsic :: iso_c_binding, only: c_int
       end subroutine coulExch_tpool
    end interface
#endif
  
    d=>supplptr(i4b(3):)
    e=>supplptr(i4b(4):)
    excp=>exchptr
    bpot=>supplptr(i4b(2):)
    f3=>supplptr(i4b(8):)
    g=>supplptr(i4b(10):)
    psi=>orbptr
    cw_sor=>sorptr

    iorbc=iorb
    maxsor1c=maxsorpot(1)    
    maxsor2c=maxsorpot(2)

    omegac=ovfexch
    omega1c=1.0_PREC-omegac
    nthreadsc=nexchpots(iorb)

#ifdef TRACE    
    write(*,'("TRACE:  coulExchSORPT/iorb, nexchpots, maxpots", 4i5)') iorb,nexchpots(iorb),maxpots
#endif
    if (nexchpots(iorb)==0) return

#ifdef PTHREAD
    call coulExch_pthread ()
#endif

#ifdef TPOOL
    call coulExch_tpool ()
#endif
  end subroutine coulExchSORPT
#endif

  subroutine coulExchMCSOR (iorb)
    use params
    use discrete
    use scfshr
    use solver
    use commons
    use sharedMemory
#ifdef OPENMP
    use sormcsor
# else
    use sormcsor
#endif
    use inout
    use utils
#ifdef OPENMP
    use omp_lib
#endif

    implicit none
    integer (KIND=IPREC) :: i,ib1,ib2,ibexp,idel,in1,in2,ioffs1,&
         iorb,iorb1,isym,itr1,nexchpot

    integer (KIND=IPREC),dimension(:), pointer :: cw_sor
    real (PREC), dimension(:), pointer :: psi,excp,bpot,d,e,f3,g

#ifdef OPENMP
    real (PREC), dimension(:), allocatable :: lhs,rhs,wk2
#else
    real (PREC), dimension(:), pointer :: lhs,rhs,wk2
#endif

    d=>supplptr(i4b(3):)
    e=>supplptr(i4b(4):)
    excp=>exchptr
    bpot=>supplptr(i4b(2):)
    f3=>supplptr(i4b(8):)
    g=>supplptr(i4b(10):)
    psi=>orbptr
    cw_sor=>sorptr
    
    omega=ovfexch
    omega1=1.0_PREC-omega
    maxsor1=maxsorpot(1)    
    maxsor2=maxsorpot(2)

    ! It turns out that the parallel region below cannot be nested due to a heavy overhead
    ! (with 'mcsor 4' the execution time of Ar on 181/40 grid increases
    ! tenfold). Therefore SOR routine must be used insted of MCSOR one.
    ! call omp_set_num_threads(nexchpots(iorb))

    !$OMP PARALLEL NUM_THREADS(nexchpots(iorb)) DEFAULT(SHARED) PRIVATE(nexchpot,isym,i,ib1,ib2,in1,in2,ibexp,idel,itr1,lhs,rhs,wk2)
#ifdef OPENMP
    i=OMP_get_thread_num()
    nexchpot=i+1
#else
    do nexchpot=1,nexchpots(iorb)
#endif

#ifdef OPENMP
       allocate(lhs(mxsize))
       allocate(rhs(mxsize))
       allocate(wk2(mxsize8))
#else
       lhs =>scratchptr(          1:   mxsize8)
       rhs =>scratchptr(   mxsize8+1: 2*mxsize8)
       wk2 =>scratchptr( 2*mxsize8+1: 3*mxsize8)
#endif    

       in1=ins1(iorb,nexchpot)
       in2=ins2(iorb,nexchpot)
       
       ib1=i1b(in1)
       ib2=i1b(in2)
       
       ibexp=ibexcp(iorb,nexchpot)

       idel=deltam4pot(iorb,nexchpot)
       isym=isyms(iorb,nexchpot)

       do i=1,mxsize
          rhs(i)=psi(ib1+i-1)*psi(ib2+i-1)*g(i)
       enddo

       if (idel==0) then
          do i=1,mxsize
             lhs(i)=f3(i)
          enddo
       else
          do i=1,mxsize
             lhs(i)=f3(i)+dble(idel*idel)*e(i)
          enddo
       endif

       do itr1=1,maxsor1
          ! prepare left-hand side of the poisson equation include the diagonal part of
          ! the differentiation operator in lhs
          
          !call zeroArray(mxsize8,wk2)
          !wk2=0
          
          call putin (nni,mxnmu,isym,excp(ibexp:),wk2)
#ifdef OPENMP
          call mcsor (isym,wk2,lhs,rhs,bpot,d,       &
               cw_sor(iadext(1):),cw_sor(iadnor(1):),&
               cw_sor(iadex1(1):),cw_sor(iadex2(1):),&
               cw_sor(iadex3(1):))
#else
          call sor (isym,wk2,lhs,rhs,bpot,d,         &   
               cw_sor(iadext(1):),cw_sor(iadnor(1):),&
               cw_sor(iadex1(1):),cw_sor(iadex2(1):),&
               cw_sor(iadex3(1):))
#endif
          call putout (nni,mxnmu,excp(ibexp:),wk2)
       enddo

#ifdef OPENMP
       deallocate(lhs)
       deallocate(rhs)
       deallocate(wk2)
#else
    enddo
#endif
    !$OMP END PARALLEL 

  end subroutine coulExchMCSOR

#if ( defined PTHREAD || defined TPOOL )
  subroutine coulExchMCSORPT (iorb)
    use commons
    use discrete
    use sormcsor
    use params
    use inout
    use scfshr
    use sharedmemory
    use solver
    use utils

    implicit none
    integer (KIND=IPREC) :: i,ib1,ib2,ibexp,idel,ifirst,in1,in2, &
         iorb,iorb1,isym,itr1,nexchpot
    integer (KIND=IPREC),dimension(:), pointer :: cw_sor
    real (PREC), dimension(:), pointer :: psi,excp,bpot,d,e,f3,g
    real (PREC), dimension(:), pointer :: lhs,rhs,wk2
    
    real (PREC) omegac,omega1c
    integer (KIND=IPREC) :: mxnmuc,mxsizec,mxsize8c,ngrid6ac,ngrid6bc,ngrid7c,&
         nnu1c,nnu2c,nnu3c,nnu4c,nnu5c,isstartc,isstopc,maxsor1c,maxsor2c

    common /c_interface_17/ mxnmuc,mxsizec,ngrid6ac,ngrid6bc,ngrid7c,&
         nnu1c,nnu2c,nnu3c,nnu4c,nnu5c,isstartc,isstopc,maxsor1c,maxsor2c
    common /c_interface_45/ omegac,omega1c

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
    
    !external mcsor_pthread,mcsor_tpool

    d=>supplptr(i4b(3):)
    e=>supplptr(i4b(4):)
    excp=>exchptr
    bpot=>supplptr(i4b(2):)
    f3=>supplptr(i4b(8):)
    g=>supplptr(i4b(10):)
    psi=>orbptr

    lhs =>scratchptr(          1:   mxsize8)
    rhs =>scratchptr(   mxsize8+1: 2*mxsize8)
    wk2 =>scratchptr( 2*mxsize8+1: 3*mxsize8)
    
    cw_sor=>sorptr

    omegac=ovfexch
    omega1c=1.0_PREC-omegac

    maxsor1=maxsorpot(1)
    maxsor2c=maxsorpot(2)    

    do nexchpot=1,nexchpots(iorb)
       in1=ins1(iorb,nexchpot)
       in2=ins2(iorb,nexchpot)
       
       ib1=i1b(in1)
       ib2=i1b(in2)

       ibexp=ibexcp(iorb,nexchpot)

       idel=deltam4pot(iorb,nexchpot)
       isym=isyms(iorb,nexchpot)

       in1=ins1(iorb,nexchpot)
       in2=ins2(iorb,nexchpot)
       
       ib1=i1b(in1)
       ib2=i1b(in2)
       
       ibexp=ibexcp(iorb,nexchpot)
       
       idel=deltam4pot(iorb,nexchpot)
       isym=isyms(iorb,nexchpot)
       !write(*,'(10i10)') iorb,nexchpot,in1,in2,ibexp,idel,isym,ib1,ib2
       
       call prod2 (mxsize,psi(ib1:),psi(ib2:),wk2)
       call prod2 (mxsize,wk2,g,rhs)

       ! prepare left-hand side of the poisson equation include the diagonal part of
       ! the differentiation operator in lhs
       if (idel==0) then
          do i=1,mxsize
             lhs(i)=f3(i)
          enddo
       else
          do i=1,mxsize
             lhs(i)=f3(i)+dble(idel*idel)*e(i)
          enddo
       endif
       
       do itr1=1,maxsor1       
          call putin (nni,mxnmu,isym,excp(ibexp:),wk2)
#if ( defined PTHREAD )
          call mcsor_pthread(isym,nthreads,wk2,lhs,rhs)
#elif ( defined TPOOL )
          call mcsor_tpool(isym,nthreads,wk2,lhs,rhs)
#endif
          call putout (nni,mxnmu,excp(ibexp:),wk2)
       enddo
    enddo
  end subroutine coulExchMCSORPT
#endif

end module relaxPots
