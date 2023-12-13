! SPDX-License-Identifier: GPL-2.0-or-later

! Copyright (C) 1996-2023  Jacek Kobus 

module initVariables
  implicit none
contains
  ! ### setDefaults ###
  !
  !     This routine sets default values of various (scalar and array)
  !     variables.
  !
  subroutine setDefaults
    use params
    use discrete
    use scfshr
    use solver
    use commons

    implicit none

    integer (KIND=IPREC) :: i,j,lip_order

    common /sap4interp/ lip_order

    ! lip_order<=3 provides safe interpolation; higher values lead to +/-Infinity at some grid points
    ! (they might disaper if the internucler separation is changed just a bit; see lookup_array (sap.F90).
    ! 
    ! lip_oder=4 produces a single infinity for H2
    lip_order=3

    ! It is assumed that 2dhf_input.dat file is present
    ldatPresent=.true.
    lomegaz=.true.

    ! default data input format depends on compilation options
    lengthintin=lengthint
    lengthfpin=lengthfp

    ! set format of floating-point numbers
    if (lengthfpin.eq.8) then
       formfp=formfp64
    else
       formfp=formfp128
    endif

    ! printout according to the precision used
    if (precis.gt.1.e-20_PREC) then
       iprint16=0
    else
       iprint16=1
    endif

    lplot=.false.        
    iplot=1
    istep=2
    plotThreshold=1.0e-10_PREC
    ltail=.false.    
    homoIncl=.false.
    lcoulexch=.true.
    openShell=.false.
    ! for a moment LXC functionals can be used for closed-shell systems only
    lxcPolar=.false.
    !lxcPolar=.true.
    ldetectNaN=.false.
    lfixNaN=.false.
    !dens_threshold=1e-48_PREC
    !dens_threshold=1e-18_PREC
    dens_threshold=1e-10_PREC
    !dens_threshold=1e-5_PREC
    lorbmcsor=.false.
    lpotmcsor=.false.
    lmcsorpt=.false.
    nthreads= 1
    lmmoments=.false.

    ! imethod=1 -- hf
    ! imethod=2 -- oed
    ! imethod=3 -- hfs
    ! imethod=4 -- dft
    ! imethod=5 -- scmc

    ! HF method is chosen as a default (card METHOD is not compulsory any more)
    imethod=1
    DFT=.false.
    HF=.false.
    HFinput=.false.    
    HFS=.false.
    LXC=.false.
    lxcHyb=.false.    
    OED=.false.
    TED=.false.    
    SCMC=.false.    

    ! Has the execuitable been compiled with the LXC library support?
    lxcSupp=.false.
#ifdef LIBXC
    lxcSupp=.true.
#endif

    lpotCoulomb=.true.
    lpotCoul2=.false.
    lpotCoul3=.false.            
    lpotFermi=.false.        
    lpotGauss=.false.
    lpotGSZ=.false.
    lpotGSZG=.false.
    lpotHarm2=.false.    
    lpotHarm3=.false.
    lpotHooke=.false.
    lpotKH=.false.    
    lpotSAP=.false.
    lspherium=.false.    
    lextracule=.false.
    lintracule=.false.
    hooke=half
    
    initAddData=.false.
    !ini=1
    linitFuncsHydrogen=.false.
    !ini=2
    linitFuncsGauss=.false.
    !ini=3
    linitFuncsGaussC=.false.
    !ini=4
    initFuncsOED=.false.
    !ini=5
    linitFuncsOld=.false.
    !ini=6
    linitFuncsNoexch=.false.
    !ini=11
    linitFuncsQRHF=.false.
    !ini=12
    linitFuncsLDA=.false.
    !ini=13
    linitFuncsHF=.false.
    !ini=22
    linitFuncsMolcas=.false.
    !ini=55
    linitFuncsNodat=.false.
    
    ! The results of test calculations for KrH+ on 331x554/200 grid suggest
    ! the following number of threads for MCSOR (mcsor-o and mcsor-ce);
    ! multi-threaded version of coulExchSOR used. 2.8 speedup in relaxation of
    ! orbitals and potentials.
#ifdef OPENMP
    openmp=.true.
    nthreads= 8
#else
    openmp=.false.
#endif

#ifdef PTHREAD
    pthread=.true.
    mcsorpt=.true.
    nthreads= 6    
    nthreads4coulexch = 44        
    lpthreadpoolq=.false.
#else
    pthread=.false.
#endif

#if defined TPOOL
    tpool=.true.
    mcsorpt=.true.
    nthreads= 6
    nthreads4coulexch = 44
    lpthreadpoolq=.true.
#else
    tpool=.false.
#endif    
    
    ! if |z1-z2|<homolevl then the case is treated as a homonuclear one
    ! (ihomon=2)

    homolevl=1.0e-6_PREC
    !ihomon=0
    lhomonucl=.false.
    lbreakCi=.false.
    iinterp=0

    cutorb =1.e-6_PREC
    cutcoul=1.e-6_PREC
    cutexch=1.e-6_PREC
    iftail  =13

    ! orbEnergyIncThld - if demax(1) (doSCF routine) is greater than this
    !           parameter the divergence in the SCF process is detected and the
    !           program stops

    orbEnergyIncThld = 1.0e4_PREC

    ! default canonicalization coefficients

    icanon=1

    ! multipol default value; the order of multipole expansion cannot
    ! exceed maxmpole-1 (and mpole maxmpole)

    mpole=4

    ! default number of grids and the default ordering of mesh points
    ngrids=1
    lversion=.false.
    altSweeps=.false.
    meshOrdering='middle'
    
    ! scf default values
    maxSCF=1000
    saveScfData=20
    ienterm=10
    inoterm=10
    verboseLevel=2

    ! check total energy uncertainty due to orbital norms not being equal 1:
    lcheckTotEnergy=.true.
    !lcheckTotEnergy=.false.

    lfastscf=.false.    
    lfastexch=.true.
    lfixorb=.false.
    lfixcoul=.false.
    lfixexch=.false.

    lfixorb4init=.false.    
    ienterm4init=0
    
    ! and DFT approach is switched off
    idftex=0
    idftcorr=0
    alphaf = two/three
    
    ! initials orbitals and potentials
    inclorb=1
    ini4=0
    inunFormatted   =.true.
    outunFormatted  =.true.
    inFormatI32   =.false.
    inFormatI64   =.false.
    inFormatR128  =.false.
    outFormatI32  =.false.
    outFormatI64  =.false.
    outFormatR128 =.false.

    ! sor default value
    maxsor1 = 1
    maxsor2 =10
    maxsor3 =10
    maxsororb(1)=1
    maxsororb(2)=10    
    maxsorpot(1)=1
    maxsorpot(2)=10    
    lsor    = 1
    iomega  = 0
    nsor4orb=0
    nsor4pot=0
    ireset=-1
    ovforb=-one
    ovfcoul=-one

    ! scaling factor for omega value for orbitals and potentials
    ! previous default value decreased due to "orbpot hf|lda" tests that suggested that
    ! the smaller value increases the convergence rate
    ! suitable for large grids
    omegasf=0.980_PREC
    ! suitable for smaller grids    
    !omegasf=0.920_PREC    

    ! values based on tests performed on FH 157x229/45 (used when omegaopt label is used)
    omegasfOrb=0.979_PREC
    omegasfPot=0.996_PREC


    
    ! skip nscf2skip iterations before starting the search for saturation;
    ! examine the last nenlast and nnolast iterations in case of orbital
    ! energy and normalization, respectively
    nscf2skip=600
    nenlast=20
    nnolast=20
    ! The normalization threshold is used to stop SCF iterations on par
    ! with the orbital energy threshold. If either of the thresholds is
    ! reached in nscfextra consecutive iterations the SCF process
    ! is stopped.
    ! In case fastSCF option is on nscfExtra parameters is used to stop
    ! relaxing particular orbitals that have reached one of the threshold
    ! for nscfExtra consequitive iterations.
    maxresets=5
    nscfExtra=3
    ! ipot default value (point nuclei)
    ipot=0
    sflagra=1.0_PREC
    dflagra=0.0_PREC
    
    loffDiagLM=.false.
    do i=1,norb
       do j=1,norb
          nlm(i,j)=izero
          nlm(j,i)=izero
          offDiagLM(i,j)=.false.
          offDiagLM(j,i)=.false.
       enddo
    enddo
    
    ! Select a way in which off-diagonal Lagrange multipliers are
    ! calculated (use lm1 and lm2 labels to choose two alternative schemes
    ! (see EabHF|DFT).
    lmtype=0

    lfefield=.false.
    ffield =0.0_PREC
    fgrad  =0.0_PREC
    zcutoff=1.e5_PREC
    iharm2xy=0

    ! default order of interpolation polynomials when changing grids
    ! for orbitals and potentials

    iord_nu_orb=5
    iord_nu_coul=5
    iord_nu_exch=5

    iord_mu_orb=5
    iord_mu_coul=5
    iord_mu_exch=5

    co1lda=one
    co2lda=one
    ! array iscforder is used to change the defaukt order of relaxation
    ! of orbitals during scf process
    do i=1,maxorb
       iscforder(i)=0
    enddo

    iout4kinpot=0
    
    incrni=10
    incrmu=10

    do i=1,maxorb
       gusym(i)  ='        '
       bond(i) ='        '
    enddo
    do i=1,4*maxorb
       spin(i) ='        '
    enddo

    do i=1,2*maxorb
       lagraon(i)=0
    enddo


    do i=1,1000
       idebug(i)=0
       iprint(i)=0
    enddo

    ! If recalcMMfactor>0 multipole moment expansion coefficients are recalculated
    ! every time deltaee, i.e. maximum error in orbital energy, changes by
    ! this factor.
    ! If recalcMMfactor<0 these coefficients are not calculated

    recalcMMfactor=1.150_PREC

    !     set maximum allowed unit number
    call setmaxunit()

  end subroutine setDefaults

  ! ### setmaxunit ###
  !
  !     Determines and sets maximum unit number allowed
  !
  subroutine setmaxunit
    use commons
    implicit none

    integer (KIND=IPREC) :: i,iunit
    do i=99,1000
       iunit=i
       open(iunit,status='scratch',form='unformatted',err=100)
       close(iunit)
    enddo
    maxunit=1000
    return
100 maxunit=i-1

  end subroutine setmaxunit

  ! ### setPrecision ###
  !
  !     Calculates floating-point precision and integer/real variable
  !     lengths
  !
  subroutine setPrecision
    use commons
    implicit none
    integer (KIND=IPREC) :: i,i4tmp1,i8tmp1,i8tmp2
    real (PREC) :: o

    precis=epsilon(one)

    ! set default lengths of integer and real constants and variables used

    i8tmp1=2**30
    i8tmp2=100000*i8tmp1
    i4tmp1=i8tmp2
    if (i4tmp1.eq.0) then
       lengthint=4
    else
       lengthint=8
    endif

    if (precis>1e-20_PREC) then
       lengthfp =8
    else
       lengthfp =16
    endif

  end subroutine setPrecision
  
end module initVariables
