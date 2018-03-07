! ### initCBlocks ###
!
!     Initializes arrays within common blocks and prepares grids.

module initCBlocks_m
  implicit none
contains
  subroutine initCBlocks
    use params
    use discret
    use scf
    use solver
    use commons8

    use setOmega_m
    use setOmegaOrb_m
    use setOmegaPot_m
    use setOmegaPot2_m

    implicit none

    character*8 :: sigma,pi,delta,phi,spac,ger,uger

    integer :: i,ial,ib,iba,ibo,ie,idi,idu,ig,ipi,ips,ipsu,ipu,isi,isu

    real (PREC) :: heta,hxibeg,hxiend,rinfig,xmi0

    data sigma/'sigma'/,pi/'pi'/,delta/'delta'/,phi/'phi'/
    data spac/' '/,ger/'g'/,uger/'u'/

    !     homonuclear case is treated as a heteronuclear one; if ihomon=2, i.e
    !     if HOMO label is used, the g/u symmetry is imposed explicitely

    !     ihomon=0 --  heteronuclear case
    !     ihomon=1 --  homonuclear case |z1-z2|<1.d-6
    !     ihomon=2 --  homonuclear case with forced symmetry

    !     check first if ngrids does not exceed maxgrids

    if (ngrids.gt.3) then
       write(iout6,*) 'number of grids can not exceed',maxgrids
       stop 'initCBlocks'
    endif

    !     When counting the maximum number of points in mu variable do not forget
    !     that the first point of the next grid is the last on of the previous grid.

    !     Arrays ibmu and iemu contain the first and last addresses of
    !     particular grid within the combined grid

    ibmu(1)=1
    iemu(1)=nmu(1)
    mxnmu=nmu(1)
    ioffs(1)=0
    do ig=2,ngrids
       mxnmu=mxnmu+nmu(ig)-1
       ibmu(ig)=iemu(ig-1)
       iemu(ig)=mxnmu
       ioffs(ig)=nni*(iemu(ig-1)-1)
    enddo

    !     determine the size of each grid

    do ig=1,ngrids
       ngsize(ig)=nni*nmu(ig)
    enddo

    !     total no of grid points

    mxsize=nni*mxnmu
    if ( norb.gt.maxorb ) then
       write(*,*) 'Error: number of orbitals cannot exceed',maxorb
       stop 'initCBlocks'
    endif

    nel=0
    do i=1,norb
       nn(i)=mgx(1,i)
       ll(i)=mgx(3,i)
       mm(i)=mgx(3,i)
       iocc(i)=NINT(occ(i)+.000001_PREC)
       nel=nel+iocc(i)

       if (mod(ll(i),2).eq.0) then
          isymOrb(i)= 1
       else
          isymOrb(i)=-1
       endif
    enddo

    !     hni - step in ni variable
    !     hmu - step in mu variable

    !     determine step size in ni variable

    hni=pii/dble(nni-1)

    !     determine step size in mu variable for each grid

    !     when ngrids=1 hmu can be calculated from the practical infinity
    !     variable if its value is provided in the input;

    if (ngrids.eq.1) then
       xmi0=2._PREC*rgrid(1)/r
       xmi0=log(xmi0+sqrt(xmi0*xmi0-1._PREC))
       hmu(1)=xmi0/dble(mxnmu-1)
       rgrid(1)=hmu(1)
    else
       do ig=1,ngrids
          if(ig.eq.1) then
             hmu(ig)=rgrid(ig)
          else
             hmu(ig)=hmu(1)*rgrid(ig)
          endif
       enddo
    endif

    !    initialize mu and xi arryas

    vmu(1) =0._PREC
    vxi(1) =cosh(vmu(1))
    vxisq(1)=vxi(1)*vxi(1)
    vxi1(1)=sqrt(vxi(1)*vxi(1)-1._PREC)
    vxi2(1)=0.0_PREC

    ib=2

    do ig=1,ngrids
       ie=iemu(ig)
       do i=ib,ie
          vmu(i)=vmu(i-1)+hmu(ig)
          vxi(i)=cosh(vmu(i))
          vxisq(i)=vxi(i)*vxi(i)

          vxi1(i)=sqrt(vxi(i)*vxi(i)-1._PREC)
          vxi2(i)=vxi(i)/vxi1(i)
          !           enddo over i
       enddo
       ib=ie+1
       !        enddo over ig
    enddo

    !    initialize ni and eta arrays

    do i=1,nni
       vni(i)=dble((i-1))*hni
       veta(i)=cos(vni(i))
       vetasq(i)=veta(i)*veta(i)

       !     sqrt(1-veta*veta) and  veta/sqrt(1-veta*veta)

       veta1(i)=sqrt(1._PREC-veta(i)*veta(i))
       if (veta1(i).lt.precis) then
          veta2(i)=0.0
       else
          veta2(i)=veta(i)/veta1(i)
       endif
    enddo

    !     determine rinf

    if (iprint(180).ne.0) then
       write(iout6,*)
       write(iout6,*) 'subgrid   nni  nmu     hni     hmu      rb    heta     hxib     hxie'
    endif

    rinf=r*vxi(iemu(ngrids))/2._PREC
    heta=abs(veta(nni)-veta(nni-1))
    ib=1
    do ig=1,ngrids
       rinfig=r*vxi(iemu(ig))/2._PREC
       ie=iemu(ig)
       hxibeg=vxi(ib+1)-vxi(ib)
       hxiend=vxi(ie)-vxi(ie-1)
       if (iprint(180).ne.0) write(iout6,1000) ig,nni,nmu(ig),hni,hmu(ig),rinfig,heta,hxibeg,hxiend
       ib=ie
    enddo
01000 format(4x,i2,3x,i4,2x,i4,2f9.5,f8.3,3f9.5)

    !     ige(i)=1 ungerade
    !     ige(i)=2 gerade or heteronuclear case
    !     if break is on

    if (ibreak.eq.1) then
       do i=1,norb
          gut(i)=spac
       enddo
    endif

    !     give a number within a symmetry

    isu=0
    isi=0
    ipu=0
    ipi=0
    idu=0
    idi=0
    ipsu=0
    ips=0

    do ibo=1,norb
       iba=norb+1-ibo
       if (bond(iba).eq.sigma) then
          if (gut(iba).eq.ger) then
             isi=isi+1
             iorn(iba)=isi
          endif

          if (gut(iba).eq.spac) then
             isi=isi+1
             iorn(iba)=isi
          endif

          if (gut(iba).eq.uger) then
             isu=isu+1
             iorn(iba)=isu
          endif
          goto 100
       endif

       if (bond(iba).eq.pi) then
          if (gut(iba).eq.ger) then
             ipi=ipi+1
             iorn(iba)=ipi
          endif

          if (gut(iba).eq.spac) then
             ipi=ipi+1
             iorn(iba)=ipi
          endif

          if (gut(iba).eq.uger) then
             ipu=ipu+1
             iorn(iba)=ipu
          endif
          goto 100
       endif

       if (bond(iba).eq.delta) then
          if (gut(iba).eq.ger) then
             idi=idi+1
             iorn(iba)=idi
          endif

          if (gut(iba).eq.spac) then
             idi=idi+1
             iorn(iba)=idi
          endif

          if (gut(iba).eq.uger) then
             idu=idu+1
             iorn(iba)=idu
          endif
          goto 100
       endif

       if (bond(iba).eq.phi) then
          if (gut(iba).eq.ger) then
             ips=ips+1
             iorn(iba)=ips
          endif

          if (gut(iba).eq.spac) then
             ips=ips+1
             iorn(iba)=ips
          endif

          if (gut(iba).eq.uger) then
             ipsu=ipsu+1
             iorn(iba)=ipsu
          endif
       endif
100    continue
    enddo

    !     ihsym = 1 - symmetry g
    !     ihsym =-1 - symmetry u
    do i=1,norb
       ihomo(i)=1
       !        if (gut(i).eq.'u') ihomo(i)=-1
       if (gut(i).eq.'u'.and.orbsym(i).eq.sigma) ihomo(i)=-1
       if (gut(i).eq.'u'.and.orbsym(i).eq.delta) ihomo(i)=-1
       if (gut(i).eq.'g'.and.orbsym(i).eq.pi)    ihomo(i)=-1
       if (gut(i).eq.'g'.and.orbsym(i).eq.phi)   ihomo(i)=-1
    enddo

    do ial=1,norb
       ige(ial)=2
       if (gut(ial).eq.'u') ige(ial)=1
    enddo

    !     set (more or less) optimal value of omega for orbitals ovforb(1)< 0
    !     by scaling the optimal value for potentials

    !     in case of convergence problems decrease this value explicitely or
    !     change the omega scaling factor accordingly

    if (ovforb(1).lt.0._PREC) then
       ovforb(1)=omegasf*setOmega()
    endif
    !     set optimal value of omega for potentials if ovfcoul(1)< 0
    if (ovfcoul(1).lt.0._PREC) then
       ovfcoul(1)=setOmega()
       ovfexch(1)=ovfcoul(1)
    endif

    if (iomega.eq.2) then
       ovforb(1)=omegasfOrb*setOmegaOrb()
       ovfcoul(1)=omegasfPot*setOmegaPot()
       ovfexch(1)=ovfcoul(1)
    endif

    if (iomega.eq.3) then
       ovfcoul(1)=omegasfPot*setOmegaPot2()
       ovforb(1)=omegasfOrb*setOmegaPot2()
       ovfexch(1)=ovfcoul(1)
    endif

    !     initialize array calp of coefficients of associated Legendre
    !     functions needed by mulex and multi routines

    !      data c/1.25d-01,       5.d-01,         5.59016994d-01,
    !             4.33012702d-01, 1.2247448710_PREC, 3.95284708d-01,
    !             1.369306394_PREC,  6.12372436d-01, 1.479019946_PREC,
    !             5.59016994d-01, 5.22912517d-01, 7.07106781d-01/

    calp( 1)=0.125_PREC
    calp( 2)=0.5_PREC
    calp( 3)=sqrt( 5._PREC/ 16._PREC)
    calp( 4)=sqrt( 3._PREC/ 16._PREC)
    calp( 5)=sqrt( 3._PREC/  2._PREC)
    calp( 6)=sqrt( 5._PREC/ 32._PREC)
    calp( 7)=sqrt(15._PREC/  8._PREC)
    calp( 8)=sqrt( 3._PREC/  8._PREC)
    calp( 9)=sqrt(35._PREC/ 16._PREC)
    calp(10)=sqrt( 5._PREC/ 16._PREC)
    calp(11)=sqrt(35._PREC/128._PREC)
    calp(12)=sqrt( 1._PREC/  2._PREC)

    !     initialize extrapolation coefficients for even functions
    !     bondary values for odd functions must be all zero

    exeven(1)= 210._PREC/126._PREC
    exeven(2)=-120._PREC/126._PREC
    exeven(3)=  45._PREC/126._PREC
    exeven(4)= -10._PREC/126._PREC
    exeven(5)=   1._PREC/126._PREC


  end subroutine initCBlocks
end module initCBlocks_m
