! SPDX-License-Identifier: GPL-2.0-or-later

! Copyright (C) 1996       Leif Laaksonen, Dage Sundholm               
! Copyright (C) 1996-2024  Jacek Kobus 

! ### x2dhf ###
!
!     This is the main routine of the 2DHF program.
!

PROGRAM  x2dhf
  use params
  use commons
  use dateTime  
  use initCBAllocArrays
  use initOrbitalsPotentials
  use initVariables
  use inputInterface
  use memory
  use printInitData
  use summary
  use doSCF
  use scfshr
  use sharedMemory
  use solver
  use utils
  
#ifdef TPOOL  
!  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR, C_INT
  USE, INTRINSIC :: ISO_C_BINDING
#endif
  implicit none  

  character*80 :: currDateTime

  integer (KIND=IPREC) :: i4,i8,iorb1,iorb2,tpoolStart
  integer (KIND=IPREC8) :: n
  integer (KIND=IPREC) :: ni,mu,no,nons,ni8,mu8,ngrid,ngrid8,size
  integer (KIND=IPREC) :: count1,count2,countRate
  integer (KIND=IPREC) :: max_num_threads_coulexch
  real (PREC) :: trelaxReal
  data i4,i8/4,8/

#if defined TPOOL  
  interface
     subroutine tpoolStop4mcsor ( cnthreads ) bind(c,name='tpoolStop4mcsor_')
       use, intrinsic :: iso_c_binding, only: c_int
       integer(c_int), intent(in)  :: cnthreads
     end subroutine tpoolStop4mcsor

     subroutine tpoolStop4pots ( cnthreads4pots ) bind(c,name='tpoolStop4pots_')
       use, intrinsic :: iso_c_binding, only: c_int
       integer(c_int), intent(in)  :: cnthreads4pots
     end subroutine tpoolStop4pots
  end interface

  type(c_ptr) :: barrier4mcsor
  common/tpool_pointer2/ barrier4mcsor
  integer (c_int) :: cnthreads4pots,cnthreads4mcsor
#endif  

  ! set the number of additional arrays
  nsuppl=14
  nsctch=14
#ifdef LIBXC  
  nsctch  = 30
#endif
 
  call separator
  call printBanner
  call separator

  call getDateTime(currDateTime)
  call getRealTime(count1,countRate)
  
  ! determine and set precision constants
  call setPrecision
  
  ! set default values of various parameters
  call setDefaults

  ! no and nons denote the number of all the orbitals and number
  ! of non sigma type orbitals, respectively. ni and mu correspond to
  ! the number of grid points in ni and mu variables.
  
  ! read input data and echo input parameters
  ! open input/output files
  
  call inputData(ni,mu,no,nons)
  firstOrb=1
  lastOrb=norb

  ! set array dimensions
  ni8=ni+8
  mu8=mu+8
  ngrid=ni*mu
  ngrid8=ni8*mu8
  
  length1=no*ngrid
  ! length2   = no*ngrid
  ! extra space is needed to store fock1 and fock2 
  length2= (no+2)*ngrid

  ! Let's extend the space originally reserved for storing exchange
  ! potentials to store also the Coulomb potentials. The original excp area
  ! is shifted rightward by no*ngrid (no*mxsize) elements.

  length3=no*ngrid + ( no*(no+1)/2+nons*(nons+1)/2 )*ngrid8

  ! ngrid8 is used below to allow excp to be used as a working array in fockDFT
  
  ! HFS+DFT
  if (HFS.or.DFT) length3 = no*ngrid+3*ngrid8
  
  ! SCMC calculates exchange potentials (as in HF case) and needs extra space to
  ! store the local exchange potential at the end of excp array
  if (SCMC) then
     length3=no*ngrid+( no*(no+1)/2+nons*(nons+1)/2 + 1 )*ngrid8
  endif

  if (lcoulexch) then
     length3= length3+length2
  endif

  length4   = nsuppl*ngrid
  length5   = nsctch*ngrid8
  length6   = 2*ngrid8+4*ni*(mu+16)+2*maxgrids*(ni+mu)
  lexchrecl = ni*mu*8

  !length7   = 24*ngrid8
  length7   = 30*ngrid8
  length8   = 8*ngrid
  size=lengthfp
  
  ! dynamic allocation of memory acccording to the current case
  allocate(orbptr(length1))
  allocate(exchptr(length3))
  allocate(supplptr(length4))
  allocate(scratchptr(length5))
  allocate(sorptr(length6))
  allocate(scratch4lxcptr(length7))
  allocate(coulombptr(length1))
  allocate(exchangeptr(length1))
  allocate(legendreptr(length8))

  
  ! zeroise arrays
  call zeroArray8(length1,orbptr)
  !if (.not.lcoulexch) call zeroArray8(length2,cw_coul_add)
  !if (.not.lcoulexch) call zeroArray8(length2,coulptr)
  call zeroArray8(length3,exchptr)
  call zeroArray8(length4,supplptr)
  call zeroArray8(length5,scratchptr)
  call zeroArray8(length7,scratch4lxcptr)
  call zeroArray8(length8,legendreptr)    
  
  do n=1,length6
     sorptr(n)=0
  enddo
  
  ! off-diagonal Lagrange multipliers must be zero (see fock)
  do iorb1=1,norb
     do iorb2=1,norb
        ee(iorb1,iorb2)=zero
     enddo
  enddo

  ! initialize arrays of variable length and check input data
  ! prepare environment for Poisson's equation solving routines
  ! prepare arrays determining ordering of mesh points

  call initArrays


#if ( defined PTHREAD )
  write(*,'("   PTHREAD: ",i2," extra threads (at most) will be used for ",&
       "relaxing Coulomb/exchange potentials")') maxpots     
  if (lmcsorpt) then
     write(*,'("   PTHREAD: ",i2," extra threads will be used for ",&
          "parallelisation of MCSORPT routine")') nthreads
  endif
#endif

 
#if ( defined TPOOL )  
  ! nthreads determine how many threads are to be used when the
  ! parallelised version of the MCSOR method is used.
  if (lmcsorpt) then
     cnthreads4mcsor=nthreads
  endif

  ! maxpots is equal to the maximum number of exchange potentials (plus a
  ! Coulomb one) to be relaxed for the "worst" orbital can be relaxed
  ! during the SCF process. Therefore a pool of maxpots threads can be
  ! created
  cnthreads4pots=maxpots
#endif

  call printCase
  call initOrbPot 
  call prepSCF
  call SCF 
  call printResults 

  deallocate(orbptr)
  deallocate(exchptr)
  deallocate(supplptr)
  deallocate(scratchptr)
  deallocate(scratch4lxcptr)
  deallocate(sorptr)

  ! record the wall time needed for the program run
  call separator
  write(*,'(a15,a80)') '       start:  ',currDateTime
  call getDateTime(currDateTime)
  call getRealTime (count2,countRate)
  write(*,'(a15,a80)') '        stop:  ',currDateTime

  if (count2>count1) then
     trelaxReal=dble(count2-count1)/dble(countRate)
  else
     !trelaxReal=dble(count2+(2**31-1-count1))/dble(countRate)
     trelaxReal=dble(count2+(2147483647-count1))/dble(countRate)
  endif
  write(*,'("  start-stop:  ",f23.2/)') trelaxReal    

#if defined TPOOL
  ! the thread pool(s) can now be destroyed
  if (lcoulexch.and..not.lfixexch) then
     call tpoolStop4pots(cnthreads4pots)
  endif

  ! if (lmcsorpt) then
  !    call tpoolStop4mcsor(cnthreads4mcsor)
  ! endif
  write(*,*)
#endif
end program x2dhf

