! ### x2dhf ###
!
!     This is the main routine of the 2DHF program.

PROGRAM  x2dhf
  use params
  use memory
  use scf
  use commons8

  use doSCF_m
  use separator_m
  use printBanner_m
  use getDateTime_m
  use setPrecision_m
  use setDefaults_m
  use inputData_m
  use zeroArray_m
  use initArrays_m
  use printCase_m
  use initOrbPot_m
  use prepSCF_m
  use printResults_m


  implicit none

  character*80 :: datetime
!      character*1, dimension(80) ::  str
!      character*80 :: datetimex
!      equivalence(datetimex,str(1))

      integer :: i4,i8,n,next,error
      integer :: ni,mu,no,nons,ni8,mu8,ngrid,ngrid8,size

      data i4,i8/4,8/

      integer, dimension(:), allocatable :: cw_sor_add
      real (PREC), dimension(:), allocatable :: cw_orb_add,cw_coul_add,cw_exch_add,cw_suppl_add,cw_sctch_add

 100  continue
      next=0
!     set number of additional arrays
      nsuppl=14
      nsctch=14

      call separator
      call printBanner
      call separator

      call  getDateTime(datetime)
!      dtstart=datetime

!     determine and set precision constants
      call setPrecision

!     set default values of various parameters
      call setDefaults

!     no and nons denote the number of all the orbitals and number
!     of non sigma type orbitals, respectively. ni and mu correspond to
!     the number of grid points in ni and mu variables.

!     read input data and echo input parameters
!     open input/output files

      call inputData(ni,mu,no,nons)

!     set array dimensions

      ni8=ni+8
      mu8=mu+8
      ngrid=ni*mu
      ngrid8    = ni8*mu8

      length1   = no*ngrid
      length2   = no*ngrid

      if     (iform.eq.1.or.iform.eq.3) then
         length3=( no*(no+1)/2+nons*(nons+1)/2 )*ngrid8
      elseif (iform.eq.0.or.iform.eq.2) then
         length3 = (no+nons)*ngrid8
      endif

!     ngrid8 is used below to allow excp to be used as a working array in fockDFT

!     OED
      if (imethod.eq.2 ) length3 =( no*(no+1)/2+nons*(nons+1)/2 )*ngrid8

!     HFS+DFT
      if (imethod.eq.3.or.imethod.eq.4) length3 = 2*ngrid8

!     SCMC calculates exchange potentials (as in HF case) and needs extra space to
!     store the local exchange potential at the end of excp array
      if (imethod.eq.5) then
         length3=( no*(no+1)/2+nons*(nons+1)/2 + 1 )*ngrid8
      endif

      length4   = nsuppl*ngrid
      length5   = nsctch*ngrid8
      length6   = 2*ngrid8+4*ni*(mu+16)+2*maxgrids*(ni+mu)
      lexchrecl = ni*mu*8

!     dynamic allocation of memory acccording to the current case

      size=lengthfp

      allocate(cw_orb_add(length1),stat=error)
      ! if (stat.ne.0) then
      !    print*,"x2dhf: Error! Couldn't allocate memory for array cw_orb_add, dim=",length1
      ! else
      !    write(*,'("  allocated",f6.2,"MB for cw_orb_add function")') length1*size/dble(1024*1024)
      ! endif

      allocate(cw_coul_add(length2),stat=error)

      ! if (stat.ne.0) then
      !    print*,"x2dhf: Error! Couldn't allocate memory for array cw_coul_add, dim=",length2
      ! else
      !    write(*,'("  allocated",f6.2,"MB for cw_coul_add function")') length2*size/dble(1024*1024)
      ! endif

      allocate(cw_exch_add(length3),stat=error)

      ! if (stat.ne.0) then
      !    print*,"x2dhf: Error! Couldn't allocate memory for array cw_exch_add, dim=",length3
      ! else
      !    write(*,'("  allocated",f6.2,"MB for cw_exch_add function")') length3*size/dble(1024*1024)
      ! endif


      allocate(cw_suppl_add(length4),stat=error)
      ! if (stat.ne.0) then
      !    print*,"x2dhf: Error! Couldn't allocate memory for array cw_suppl_add, dim=",length4
      !    stop 'x2dhf'
      ! else
      !    write(*,'("  allocated",f6.2,"MB for supplementary arrays")') length4*size/dble(1024*1024)
      ! endif

      allocate(cw_sctch_add(length5),stat=error)
      ! if (stat.ne.0) then
      !    print*,"x2dhf: Error! Couldn't allocate memory for array cw_sctch_add, dim=",length5
      !    stop 'x2dhf'
      ! else
      !    write(*,'("  allocated",f6.2,"MB for scratch arrays")') length5*size/dble(1024*1024)
      ! endif


      allocate(cw_sor_add(length6),stat=error)
      ! if (stat.ne.0) then
      !    print*,"x2dhf: Error! Couldn't allocate memory for array cw_sor_add, dim=",length6
      !    stop 'x2dhf'
      ! else
      !    write(*,'("  allocated",f6.2,"MB for cw_sor arrays")') length6*size/dble(1024*1024)
      ! endif

!     zeroise arrays
      call zeroArray(length1,cw_orb_add)
      call zeroArray(length2,cw_coul_add)
      call zeroArray(length3,cw_exch_add)
      call zeroArray(length4,cw_suppl_add)
      call zeroArray(length5,cw_sctch_add)

      do n=1,length6
         cw_sor_add(n)=0
      enddo

      n=maxorb*maxorb
      call zeroArray(n,engo)

      !     initialize arrays within common blocks,
      !     initialize arrays of variable length and check input data
      !     prepare environment for Poisson's equation solving routines
      !     prepare arrays determining ordering of mesh points

      call initArrays (cw_suppl_add,cw_sor_add)

      !     print out data defining the case

       call printCase

       !     initialize orbitals, Coulomb and exchange potentials

       call initOrbPot (cw_orb_add,cw_coul_add,cw_exch_add,cw_suppl_add,cw_sctch_add)

       !     prepare SCF process: normalize/orthogonalize orbitals, calculate
       !     Lagrange multipliers prepare data for boundary value evaluation

       call prepSCF (cw_sor_add,cw_orb_add,cw_coul_add,cw_exch_add,cw_suppl_add,cw_sctch_add)

! c  perform scf
       call  doSCF (cw_sor_add,cw_orb_add,cw_coul_add,cw_exch_add,cw_suppl_add,cw_sctch_add)

       call printResults (cw_orb_add,cw_coul_add,cw_exch_add,cw_suppl_add,cw_sctch_add)

      deallocate(cw_orb_add,stat=error)
      if (error.ne.0) then
         print*,'error in deallocating array cw_sctch'
      endif

      deallocate(cw_coul_add,stat=error)
      if (error.ne.0) then
         print*,'error in deallocating array cw_sctch'
      endif

      deallocate(cw_exch_add,stat=error)
      if (error.ne.0) then
         print*,'error in deallocating array cw_sctch'
      endif

      deallocate(cw_suppl_add,stat=error)
      if (error.ne.0) then
         print*,'error in deallocating array cw_sctch'
      endif

      deallocate(cw_sctch_add,stat=error)
      if (error.ne.0) then
         print*,'error in deallocating array cw_sctch'
      endif

      deallocate(cw_sor_add,stat=error)
      if (error.ne.0) then
         print*,'error in deallocating array cw_sor'
      endif

! record the wall time needed for the program run
!       datetimex=dtstart
!       write(*,'(//,1x,a7,80a1)') 'start: ',(str(i),i=1,24)
       write(*,'(//,1x,a7,a80)') 'start: ',datetime
       call  getDateTime(datetime)
!       datetimex=datetime
!       write(*,'(1x,a7,80a1)') 'stop:  ',(str(i),i=1,24)
       write(*,'(1x,a7,a80)') 'stop:  ',datetime
       call separator

      if (next.eq.1) then
         goto 100
      else
         stop 'x2dhf'
      endif

end program x2dhf
