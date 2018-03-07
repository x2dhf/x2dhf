! ### setDefaults ###
!
!     This routine sets default values of various (scalar and array)
!     variables.

module setDefaults_m
  implicit none
contains
  subroutine setDefaults
    use params
    use discret
    use scf
    use solver
    use commons8
    use setmaxunit_m

    implicit none

    integer :: i

    !     no=isum=norb
    !     iexd=maxorb
    !     bond=orbsym

    inpform=0
    ioutform=0
    inout32=0
    inout64=0
    inout128=0

    iout4dd=0

    !     2dhf_input.dat is assumed
    idat=0


    !     default data input format depends on compilation options
    lengthintin=lengthint
    lengthfpin=lengthfp

    !     set format of floating-point numbers
    if (lengthfpin.eq.8) then
       formfp=formfp64
    else
       formfp=formfp128
    endif

    !     printout according to the precision used
    if (precis.gt.1.e-20_PREC) then
       iprint16=0
    else
       iprint16=1
    endif

    !  initialization of some variables

    !     if |z1-z2|<homolevl then the case is treated as a homonuclear one
    !     (ihomon=2)

    homolevl=1.0e-6_PREC
    ihomon=0
    ibreak=0
    iinterp=0

    cutorb =1.e-6_PREC
    cutcoul=1.e-6_PREC
    cutexch=1.e-6_PREC
    iftail  =13

    !     diver  -  if demax(1) in proce is greater than diver
    !               divergence in the scf process is detected and
    !               the program stops

    diver = 1.e20_PREC

    !     default canonicalization coefficients

    icanon=1

    !     multipol default value; the order of multipole expansion cannot
    !     exceed maxmpole-1 (and mpole maxmpole)

    mpole=4

    !     default number of grids and the default ordering of mesh points
    ngrids=1
    iorder(ngrids)=2
    ialtsweeps=0

    !     scf default values

    maxscf  = 1000
    nobckup =  20
    ienterm =  10
    inoterm =  10
    iprtlev =   2

    !     fix default values

    iexlorb  = 0
    iexlcoul = 0
    iexlexp  = 2

    exlorb  = 0.0_PREC
    exlcoul = 0.0_PREC
    exlexp  = 2.0_PREC


    !     HF method is chosen as a default (card METHOD is not compulsory any more)
    imethod=1
    !     and DFT approach is switched off
    idft=0
    idftex=0
    idftcorr=0
    alphaf = two/three
    iscmc=0


    ! initials orbitals and potentials
    inclorb = 1
    iform=3
    ini4=0


    !     sor default value
    maxsor1 = 1
    maxsor2 =10
    maxsor3 =10
    lsor    = 1
    ipoiss  = 1
    iomega  = 0
    !
    ovforb(1)=-one
    ovfcoul(1)=-one

    !     scaling factor for omega value for orbitals and potentials
    omegasf=0.980_PREC

    !     values based on test performed on FH 157x229/45
    omegasfOrb=0.979_PREC
    omegasfPot=0.996_PREC


    !     skip nscf2skip iterations before starting the search for
    !     saturation; examine the last nenlast and nnolast iterations in case
    !     of orbital energy and normalization, respectively

    nscf2skip=600
    nenlast=20
    nnolast=20

    !     ipot default value (point nuclei)

    ipot=0
    ifermi=0

    !     plot default value (no data exported for making plots)

    iplot=0

    ilagra=0
    sflagra=1.0_PREC
    dflagra=0.0_PREC

    ifefield=0
    ffield =0.0_PREC
    fgrad  =0.0_PREC
    zcutoff=1.e5_PREC
    iharm2xy=0

    !     default order of interpolation polynomials when changing grids
    !     for orbitals and potentials

    iord_nu_orb=5
    iord_nu_coul=5
    iord_nu_exch=5

    iord_mu_orb=5
    iord_mu_coul=5
    iord_mu_exch=5

    !     array iscforder is used to change the defaukt order of relaxation
    !     of orbitals during scf process

    do i=1,maxorb
       iscforder(i)=0
    enddo

    incrni=10
    incrmu=10

    do i=1,maxorb
       gut(i)  ='        '
       bond(i) ='        '
    enddo
    do i=1,4*maxorb
       spin(i) ='        '
    enddo

    do i=1,2*maxorb
       lagraon(i)=0
    enddo


    do i=1,1000
       idbg(i)=0
       iprint(i)=0
    enddo

    !     if facmul>0 multipole moment expansion coefficients are
    !     recalculated every time denmax, i.e. maximum error in orbital
    !     energy, changes by the factor facmul if facmul<0 these
    !     coefficients are not calculated

    facmul=1.150_PREC

    !     set maximum allowed unit number
    call setmaxunit()

  end subroutine setDefaults
end module setDefaults_m
