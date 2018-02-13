! ***************************************************************************
! *                                                                         *
! *   Copyright (C) 1997-2010 Jacek Kobus <jkob@fizyka.umk.pl>              *
! *   Copyright (C) 2018      Susi Lehtola                                  *
! *                                                                         *
! *   This program is free software; you can redistribute it and/or modify  *
! *   it under the terms of the GNU General Public License version 2 as     *
! *   published by the Free Software Foundation.                            *
! *                                                                         *
! ***************************************************************************
! ### prepGauss ###
!
!     This routine reads the output from the GAUSSSIAN94/98 program
!     (gaussian.out and gaussian.pun files) to determine parameters
!     of the basis functions (n, l, m, exponents) and coefficients of
!     molecular orbitals.
!     The basis can contain functions defined on centres A and B and
!     also at the bond centre.

module prepGauss
contains

  subroutine make_mapping(map,invmap,iZ,x,y,z,ncen)
    use swap
    use sort
    implicit none
    ! Gaussian center index mapping
    integer, intent(out) :: map(3), invmap(3)
    ! Nuclear charges
    integer, intent(out) :: iZ(3)
    ! x, y and z coordinates
    double precision, intent(out) :: x(3), y(3), z(3)
    ! Number of centers
    integer, intent(out) :: ncen

    ! Helper for i/o
    character(80) :: line
    ! Input center number
    integer :: icen(3)
    ! Atomic number
    integer :: iZ0(3)
    ! Atom type
    integer :: itype(3)
    ! X, Y and Z coordinates
    double precision :: x0(3), y0(3), z0(3)
    ! Is atom a dummy atom?
    logical :: idum(3)

    ! Loop index, number of centers, read status, index of dummy atom, i/o status
    integer :: i, j, stat
    ! Gaussian94 format?
    logical :: g94

    ! By default, assume G03 format
    g94 = .false.
    ! Initialize itype
    itype = 0

    ! Read four lines from input, which contain the separators
    do i=1,4
       read (7,'(A)') line
       !write (*,*) 'Skipping line ' // line
    end do

    ! Now, read the atoms
    ncen=3
    do i=1,3
       read (7,'(A)') line
       !write (*,*) 'Got line ' // line
       ! Check we haven't gone too far
       if(line .eq. ' ---------------------------------------------------------------------' .or. &
          line .eq. ' ----------------------------------------------------------') then
          ncen=i-1
          exit
       end if
       ! Parse contents of line
       if( g94 ) then
          ! g94 doesn't have atom type
          read (line,*,iostat=stat) icen(i), iZ0(i), x0(i), y0(i), z0(i)
          if(stat.ne.0) stop 'Parsing error ' // line
       else
          read (line,*,iostat=stat) icen(i), iZ0(i), itype(i), x0(i), y0(i), z0(i)
          if(stat .ne. 0) then
             ! Try switching to g94 mode
             g94=.true.
             read (line,*,iostat=stat) icen(i), iZ0(i), x0(i), y0(i), z0(i)
             if(stat.ne.0) stop 'Parsing error ' // line
          end if
       end if
    end do
    if(ncen .lt. 2) then
       stop 'Not enough atoms!'
    end if

    ! Initialize mapping
    do i=1,3
       map(i)=i
       invmap(i)=i
    end do
    ! Construct mapping by sorting in increasing z
    do i=1,ncen
       do j=1,i-1
          if(z0(map(i)) .lt. z0(map(j))) then
             call iswap(map(i),map(j))
          end if
       end do
    end do
    ! Construct inverse mapping
    do i=1,ncen
       invmap(map(i))=i
    end do

    ! Figure out dummy atom index
    idum=.false.
    if(g94) then
       do i=1,ncen
          idum(i)=iZ0(i).eq.0
       end do
    else
       do i=1,ncen
          idum(i)=itype(i).ne.0
       end do
    end if

    ! Dummy atoms have charge 0
    do i=1,ncen
       if(idum(i)) iZ0(i)=0
    end do

    ! Give output in x2dhf order
    call rsort(x0,x,map)
    call rsort(y0,y,map)
    call rsort(z0,z,map)
    call isort(iZ0,iZ,map)

  end subroutine make_mapping

  subroutine read_basis(nprim,nexpon,npbasis)
    implicit none
    integer, intent(out) :: nprim(3), nexpon, npbasis
    integer :: ibc, ibp, istop, ncen

    !     start extracting data from the output
    !     ibc counts contracted basis functions
    !     ibp counts primitive basis functions
    !     n1prim - number of primitive gaussian function at the centre A
    !     n3prim - number of primitive gaussian function at the centre B
    !     n2prim - number of primitive gaussian function at the bond centre

    npbasis=0
    ibc=0
    ibp=0
    nprim=0
    ncen=1

    ! Read basis on 1st center
    call rexponents(ibc,ibp,istop)
    nprim(1)=ibp
    if (istop.eq.0) then
       ncen=ncen+1
       ! Read basis on 2nd center
       call rexponents(ibc,ibp,istop)
       nprim(2)=ibp-sum(nprim)
       if (istop.eq.0) then
          ncen=ncen+1
          ! Read basis on 3rd center
          call rexponents(ibc,ibp,istop)
          nprim(3)=ibp-sum(nprim)
       endif
    endif

    ! Store number of functions and primitives
    nexpon=ibc
    npbasis=ibp
  end subroutine read_basis

  subroutine prepare_Gaussian
    use discret, only : z1, z2, r
    use params
    use commons8
    implicit none
    integer :: i,ib,ibc,ibp,icount,ifbo,jfbo,ilines,iorb,iprint0,iprint1,iprint2, &
         l1,m1,m1abs,nexpon,nfborb

    integer, dimension(maxorb) :: ifbord,ifdord

    real (PREC) :: d1
    real (PREC), dimension(maxbasis) :: ccoeff
    ! We can have up to 2x more orbitals in Gaussian than in the 2d
    ! calculation because non-sigma orbitals appear in twos.
    real (PREC), dimension(2*maxorb,maxbasis) :: excoeff
    real (PREC), external :: factor, factor2

    ! Orbital character
    real (PREC), dimension(0:3) :: ochar
    real (PREC) :: ocharmax
    integer :: ichar, icharmax

    character*46 :: matchstr

    ! Gaussian center index
    integer :: map(3), invmap(3)
    ! Nuclear charges
    integer :: iZ(3)
    ! x, y and z coordinates
    double precision :: x(3), y(3), z(3)
    ! Number of centers
    integer :: ncen

    ! Number of primitives on the centers
    integer :: nprim(3)
    ! Running index
    integer :: ibf
    ! Target atom
    integer :: itgt

    ! Basis function start and end
    integer :: bfstart(3), bfend(3)
    ! Found atoms and basis?
    logical :: atoms_found, basis_found

    !     nforb  - number of finite basis (fb) set orbitals
    !     ifbord - ordering of fb orbitals (0=sigma, 1=pi, 2=delta, etc)
    !              the order is determined from the gaussian.out file
    !     ifdord - ordering of finite difference (fd) orbitals
    !              the order of fd orbitals is the same as in the input data
    !              (phi orbitals first, then delta, etc)

    iprint0=0
    iprint1=0
    iprint2=0

    if (iprint(553).ne.0) iprint0=1
    if (iprint(554).ne.0) iprint1=1
    if (iprint(555).ne.0) iprint2=1

    ! Open the GaussianXY output file
    open(7,file='gaussian.out', status='old',form='formatted')

    atoms_found = .false.
    basis_found = .false.

    do ilines=1,100000
       read(7,1001,end=990,err=910) matchstr

       ! Gaussian 94
       if (matchstr.eq.'                  Z-Matrix orientation:' .or. &
            matchstr.eq.'                   Standard orientation:') then
          call make_mapping(map,invmap,iZ,x,y,z,ncen)
          atoms_found = .true.
       end if
       ! Gaussian 03/09/16
       if (matchstr.eq.'                         Standard orientation:' .or. &
            matchstr.eq.'                          Input orientation:') then
          call make_mapping(map,invmap,iZ,x,y,z,ncen)
          atoms_found = .true.
       end if
       if (matchstr.eq.' Basis set in the form of general basis input:' .or. &
            matchstr.eq.' AO basis set in the form of general basis inp') then
          call read_basis(nprim,nexpon,npbasis)
          basis_found = .true.
          exit
       end if
    enddo
00910 continue
    if (.not. atoms_found) then
       stop 'Could not find system in GAUSSIANxy output file.'
    end if
    if (.not. basis_found) then
       stop 'Could not find basis in GAUSSIANxy output file. Did you remember to specify gfinput?'
    end if

    close(7)

    ! Check atoms are on z axis
    do i=1,ncen
       if(abs(x(i)) .gt. epsilon(x(i)) .or. abs(y(i)) .gt. epsilon(y(i))) then
          stop 'Atoms are not on z axis in Gaussian output'
       end if
    end do

    ! If we have three atoms, then the middle one must be a dummy
    if(ncen==3 .and. iZ(2).ne.0) then
       stop 'Bond-center atom must be a dummy!'
    end if
    ! Check charges are correct
    if(abs(z1-iZ(1)).gt.precis) then
       write (*,'(A,F4.1,A,I3)') 'z1=',z1,', iZ(1)=',iZ(1)
       stop "Gaussian calculation is inconsistent with x2dhf for atom on the left"
    end if
    if(abs(z2-iZ(ncen)).gt.precis) then
       write (*,'(A,F4.1,A,I1,A,I3)') 'z2=',z2,', iZ(',ncen,')=',iZ(ncen)
       stop "Gaussian calculation is inconsistent with x2dhf for atom on the right"
    end if

    ! Check bond length is correct
    if( abs((z(ncen)-z(1))/bohr2ang - r) .gt. 1e-6*r ) then
       write (*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       write (*,*) '! Warning: bond distance mismatch detected in Gaussian calculation !'
       write (*,'(A,F11.6,A,F11.6,A)') ' !  bond lengths: Gaussian ',(z(ncen)-z(1))/bohr2ang,' a.u., x2dhf ',r,' a.u. !'
       write (*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    end if

    ! Form basis set table
    bfstart=0
    bfend=0
    ibf=1

    ! In the gaussian output, we have the functions on centers 1, 2,
    ! and 3 in that order. In x2dhf this corresponds to icen(1),
    ! icen(2) and icen(3).
    do i=1,ncen
       if(nprim(i).gt.0) then
          ! Target center is
          itgt=invmap(i)
          ! In the case of two centers, we must have just the lhs and
          ! rhs atom.
          if(ncen.eq.2 .and. itgt.eq.2) itgt=3
          bfstart(itgt)=ibf
          bfend(itgt)=ibf+nprim(i)-1
          ibf=ibf+nprim(i)
       end if
    end do

    write (*,'(A,I3,A,I3)') ' Center on the left   houses basis functions ',bfstart(1),'-',bfend(1)
    if((bfend(2)-bfstart(2)).gt.0) then
       write (*,'(A,I3,A,I3)') ' Center on the middle houses basis functions ',bfstart(2),'-',bfend(2)
    end if
    write (*,'(A,I3,A,I3)') ' Center on the right  houses basis functions ',bfstart(3),'-',bfend(3)

    open(7,file='gaussian.pun', status='old',form='formatted')

    !  read(7,1001,end=991,err=992) matchstr
    !    read(7,'(a8)',end=991,err=992) matchstr
    read(7,*) matchstr

    !     start extracting basis set expansion coefficients from
    !     gaussian94.pun file

    !     determine number of molecular orbitals as defined in the GAUSSIAN9x
    !     program (one nonsigma finite difference orbital corresponds
    !     to two GAUSSIAN9x ones

    nfborb=0
    do iorb=1,norb
       ifdord(iorb)=mm(iorb)
       ! write (*,*) 'fd orbital ',iorb,' m=',mm(iorb)
       if (mm(iorb).eq.0) then
          nfborb=nfborb+1
       else
          nfborb=nfborb+2
       endif
    enddo

    write(*,1142) norb,nfborb

    ! Retrieve expansion coefficients of the contracted gaussians
    if (iprint1.ne.0) then
       print *,'no. of basis function  no. of exp. coeff.'
       do ibp=1,npbasis
          ibc=ixref(ibp)
       enddo
    endif

    do ifbo=1,nfborb
       ! Read the MO coefficients
       read(7,1011,end=991,err=992)
       do ibc=1,nexpon,5
          read(7,1012) (ccoeff(ibc+i),i=0,4)
       enddo

       ! Transform expansion coefficients of the contracted gaussians
       ! into coefficients of primitive gaussians
       do ibp=1,npbasis
          ibc=ixref(ibp)
          excoeff(ifbo,ibp)=ccoeff(ibc)*coeff(ibp)
       enddo
    enddo

    ! Assign centers of basis functions
    icgau=-1
    do ib=1,npbasis
       do i=1,3
          if(ib .ge. bfstart(i) .and. ib .le. bfend(i)) icgau(ib)=i
       end do
       if(icgau(ib) .eq. -1) then
          write (*,*) 'Could not assign center for basis function ',ib
          stop
       end if
    end do

    ! Determine symmetry of GAUSSIANxy molecular orbitals
    do ifbo=1,nfborb
       ! Reset orbital character
       ochar = 0.0
       ! Compute orbital character checksum
       do ib=1,npbasis
          ochar(abs(mprim(ib))) = ochar(abs(mprim(ib))) + excoeff(ifbo,ib)**2
       end do
       ! Determine orbital symmetry from maximum
       icharmax = 0
       ocharmax = ochar(icharmax)
       do ichar=1,3
          if(ochar(ichar) > ocharmax) then
             icharmax=ichar
             ocharmax=ochar(icharmax)
          end if
       end do
       ifbord(ifbo)=icharmax
       !write (*,'(A,I2,A,4F8.4)') 'fb orbital ',ifbo,' symmetry is ',ochar
    end do

    ! Initialize memory
    do iorb=1,norb
       do ib=1,npbasis
          primcoef(iorb,ib)=0.0_PREC
       end do
    end do

    ! Associate finite basis orbitals with the corresponding fd ones
    do iorb=1,norb
       do ifbo=nfborb,1,-1
          icount=0
          ! Check if orbital characters match
          if (ifdord(iorb).eq.ifbord(ifbo)) then
             !write (*,'(A5,A,I3,A,I3)') bond(iorb),' orbital ',iorb,' corresponds to Gaussian orbital ',ifbo
             do ib=1,npbasis
                primcoef(iorb,ib)=excoeff(ifbo,ib)
             enddo
             ! If this is not a sigma orbital, we need to zero out the degenerate orbital as well.
             if (ifbord(ifbo).gt.0) then
                do jfbo=ifbo-1,1,-1
                   if(ifbord(ifbo) .eq. ifbord(jfbo)) then
                      ifbord(jfbo)=-1
                      exit
                   end if
                end do
             end if
             ! Reset orbital characters so that they'll be skipped next iteration
             ifbord(ifbo)=-1
             ifdord(iorb)=-2
          endif
       enddo
    enddo

    ! Calculate normalization factors
    do ib=1,npbasis
       d1=primexp(ib)
       l1=lprim  (ib)
       m1=mprim  (ib)
       m1abs=abs(m1)
       ! Normalization for the radial part
       fngau2(ib)=( d1**(2*l1+3) * 2.0_PREC**(4*l1+7)/pii/(factor2(2*l1+1))**2)**0.250_PREC
       ! Normalization for the angular part (spherical harmonic)
       shngau(ib)=(-1)**m1 * sqrt( ((2*l1+1)*factor(l1-m1abs)) / (4.0_PREC*pii*factor(l1+m1abs)) )
    enddo

    return
00990 print *,'PREPGAUSS: end of file encountered when reading gaussian.out'
    stop
00991 print *,'PREPGAUSS: end of file encountered when reading gaussian.pun'
    stop
00992 print *,'PREPGAUSS: error encountered when reading gaussian.pun'
    stop "prepGauss"

01001 format(a46)
01011 format(15x,e15.8)
01012 format(5e15.8)
    !01050 format(5i5,e15.8)
01142 format(15x,i4,' finite difference orbitals',/15x,i4,' finite basis set orbitals'/)
  end subroutine prepare_Gaussian
end module prepGauss
