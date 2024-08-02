! SPDX-License-Identifier: GPL-2.0-or-later

! Copyright (C) 1996-2024  Jacek Kobus 
! Copyright (C) 2018       Susi Lehtola

module prepGauss
  use params, only : IPREC, PREC
  implicit none
contains

  ! ### rexponents ###
  !
  !     Reads exponents and contraction coefficients from
  !     a GAUSSIAN94 output obtained with gfinput keyword
  !
  subroutine rexponents(ibc,ib,istop)
    use params
    use commons

    implicit none
    integer (KIND=IPREC) :: i,ib,ibc,icent,istop,j,k,maxcf,ncf

    parameter (maxcf=100)

    real (PREC), dimension(maxcf) :: expon(maxcf)
    real (PREC), dimension(maxcf,7) :: dcoef

    character(80) :: line
    character(2) :: symlab

    integer (KIND=IPREC) :: stat

    istop=0

    ! Read the centre number
    read (7,'(A)') line
    if (trim(line).eq.'') then
       ! Check if the file ended
       istop=1
       return
    end if

    read (line,*,iostat=stat) icent
    if(stat.ne.0) then
       istop=0
       return
    end if
    if (icent.eq.0) then
       istop=0
       return
    elseif (icent.eq.-9999) then
       return
    endif

    !     read symmetry label and the number of exponents
    do i=1,maxbasis
       read (7,'(A)') line

       if (trim(line).ne.' ****') then
          read (line,*) symlab, ncf
          if (ncf.gt.maxcf) then
             print *,'r_exponents: too many primitive functions per single contracted gaussian; increase maxcf'
             stop
          endif

          if (trim(symlab).eq.'S') then
             do j=1,ncf
                read (7,'(A)') line
                read (line, *) expon(j), dcoef(j,1)
             enddo
             !              define s-type gaussians
             ibc=ibc+1
             do j=1,ncf
                ib=ib+1
                if (ib.gt.maxbasis) goto 990
                coeff(ib)=dcoef(j,1)
                primexp(ib)=expon(j)
                ixref(ib)=ibc
                lprim(ib)=0
                mprim(ib)=0
             enddo
          elseif (trim(symlab).eq.'SP') then
             do j=1,ncf
                read (7,'(A)') line
                read (line,*) expon(j), dcoef(j,1), dcoef(j,2)
             enddo
             !              define s-type gaussians
             ibc=ibc+1
             do j=1,ncf
                ib=ib+1
                if (ib.gt.maxbasis) goto 990
                coeff(ib)=dcoef(j,1)
                primexp(ib)=expon(j)
                ixref(ib)=ibc
                lprim(ib)=0
                mprim(ib)=0
             enddo
             !              define p-type orbitals (x,y,z)
             do k=1,3
                ibc=ibc+1
                do j=1,ncf
                   ib=ib+1
                   if (ib.gt.maxbasis) goto 990
                   coeff(ib)=dcoef(j,2)
                   primexp(ib)=expon(j)
                   ixref(ib)=ibc
                   lprim(ib)=1
                   if (k.eq.1) mprim(ib)=+1
                   if (k.eq.2) mprim(ib)=-1
                   if (k.eq.3) mprim(ib)= 0
                enddo
             enddo
          elseif (trim(symlab).eq.'P') then
             do j=1,ncf
                read (7,'(A)') line
                read (line, *) expon(j), dcoef(j,2)
             enddo
             !              define p-type orbitals (x,y,z)
             do k=1,3
                ibc=ibc+1
                do j=1,ncf
                   ib=ib+1
                   if (ib.gt.maxbasis) goto 990
                   coeff(ib)=dcoef(j,2)
                   primexp(ib)=expon(j)
                   ixref(ib)=ibc
                   lprim(ib)=1
                   if (k.eq.1) mprim(ib)=+1
                   if (k.eq.2) mprim(ib)=-1
                   if (k.eq.3) mprim(ib)= 0
                enddo
             enddo

             !              define d-type orbitals (d0,d1,d-1,d2,d-2)
          elseif (trim(symlab).eq.'D') then
             do j=1,ncf
                read (7,'(A)') line
                read (line, *) expon(j), dcoef(j,3)
             enddo

             do k=1,5
                ibc=ibc+1
                do j=1,ncf
                   ib=ib+1
                   if (ib.gt.maxbasis) goto 990
                   coeff(ib)=dcoef(j,3)
                   primexp(ib)=expon(j)
                   ixref(ib)=ibc
                   lprim(ib)=2
                   if (k.eq.1) mprim(ib)= 0
                   if (k.eq.2) mprim(ib)=+1
                   if (k.eq.3) mprim(ib)=-1
                   if (k.eq.4) mprim(ib)=+2
                   if (k.eq.5) mprim(ib)=-2
                enddo
             enddo
             !              define f-type orbitals (f0,f1,f-1,f2,f-2,f3,f-3)`
          elseif (trim(symlab).eq.'F') then
             do j=1,ncf
                read (7,'(A)') line
                read (line, *) expon(j), dcoef(j,4)
             enddo
             do k=1,7
                ibc=ibc+1
                do j=1,ncf
                   ib=ib+1
                   if (ib.gt.maxbasis) goto 990
                   coeff(ib)=dcoef(j,4)
                   primexp(ib)=expon(j)
                   ixref(ib)=ibc
                   lprim(ib)=3
                   if (k.eq.1) mprim(ib)= 0
                   if (k.eq.2) mprim(ib)=+1
                   if (k.eq.3) mprim(ib)=-1
                   if (k.eq.4) mprim(ib)=+2
                   if (k.eq.5) mprim(ib)=-2
                   if (k.eq.6) mprim(ib)=+3
                   if (k.eq.7) mprim(ib)=-3
                enddo
             enddo
             !              define g-type orbitals (g0,g1,g-1,g2,g-2,g3,g-3,g4,g-4)`
          elseif (trim(symlab).eq.'G') then
             do j=1,ncf
                read (7,'(A)') line
                read (line, *) expon(j), dcoef(j,5)
             enddo
             do k=1,9
                ibc=ibc+1
                do j=1,ncf
                   ib=ib+1
                   if (ib.gt.maxbasis) goto 990
                   coeff(ib)=dcoef(j,5)
                   primexp(ib)=expon(j)
                   ixref(ib)=ibc
                   lprim(ib)=4
                   if (k.eq.1) mprim(ib)= 0
                   if (k.eq.2) mprim(ib)=+1
                   if (k.eq.3) mprim(ib)=-1
                   if (k.eq.4) mprim(ib)=+2
                   if (k.eq.5) mprim(ib)=-2
                   if (k.eq.6) mprim(ib)=+3
                   if (k.eq.7) mprim(ib)=-3
                   if (k.eq.8) mprim(ib)=+4
                   if (k.eq.9) mprim(ib)=-4
                enddo
             enddo
             !              define h-type orbitals (h0,h1,h-1,h2,h-2,h3,h-3,h4,h-4,h5,h-5)`
          elseif (trim(symlab).eq.'H') then
             do j=1,ncf
                read (7,'(A)') line
                read (line, *) expon(j), dcoef(j,6)
             enddo
             do k=1,11
                ibc=ibc+1
                do j=1,ncf
                   ib=ib+1
                   if (ib.gt.maxbasis) goto 990
                   coeff(ib)=dcoef(j,6)
                   primexp(ib)=expon(j)
                   ixref(ib)=ibc
                   lprim(ib)=4
                   if (k.eq.1) mprim(ib)= 0
                   if (k.eq.2) mprim(ib)=+1
                   if (k.eq.3) mprim(ib)=-1
                   if (k.eq.4) mprim(ib)=+2
                   if (k.eq.5) mprim(ib)=-2
                   if (k.eq.6) mprim(ib)=+3
                   if (k.eq.7) mprim(ib)=-3
                   if (k.eq.8) mprim(ib)=+4
                   if (k.eq.9) mprim(ib)=-4
                   if (k.eq.10) mprim(ib)=+5
                   if (k.eq.11) mprim(ib)=-5
                enddo
             enddo
             !              define i-type orbitals (i0,i1,i-1,i2,i-2,i3,i-3,i4,i-4,i5,i-5,i6,i-6)`
          elseif (trim(symlab).eq.'I') then
             do j=1,ncf
                read (7,'(A)') line
                read (line, *) expon(j), dcoef(j,7)
             enddo
             do k=1,13
                ibc=ibc+1
                do j=1,ncf
                   ib=ib+1
                   if (ib.gt.maxbasis) goto 990
                   coeff(ib)=dcoef(j,7)
                   primexp(ib)=expon(j)
                   ixref(ib)=ibc
                   lprim(ib)=4
                   if (k.eq.1) mprim(ib)= 0
                   if (k.eq.2) mprim(ib)=+1
                   if (k.eq.3) mprim(ib)=-1
                   if (k.eq.4) mprim(ib)=+2
                   if (k.eq.5) mprim(ib)=-2
                   if (k.eq.6) mprim(ib)=+3
                   if (k.eq.7) mprim(ib)=-3
                   if (k.eq.8) mprim(ib)=+4
                   if (k.eq.9) mprim(ib)=-4
                   if (k.eq.10) mprim(ib)=+5
                   if (k.eq.11) mprim(ib)=-5
                   if (k.eq.12) mprim(ib)=+6
                   if (k.eq.13) mprim(ib)=-6
                enddo
             enddo
          else
             print *,'r_exponents: basis functions higher than i are not allowed'
             stop
          endif
       else
          istop=0
          return
       endif
    enddo
    istop=1
    return

    !00904 stop 'r_exponents: error encountered when reading gauss94.out'
    !00950 stop 'r_exponents: end of gauss94.out file encountered'
00990 stop 'r_exponents: too many basis functions; increase maxbasis'
    !01000 format(a3,2x,i2)
  end subroutine rexponents
  
  ! ### prepGauss ###
  !
  !     This routine reads the output from the GAUSSSIAN94/98 program
  !     (gaussian.out and gaussian.pun files) to determine parameters
  !     of the basis functions (n, l, m, exponents) and coefficients of
  !     molecular orbitals.
  !     The basis can contain functions defined on centres A and B and
  !     also at the bond centre.
  !
  subroutine makeMapping(map,invmap,iZ,x,y,z,ncen)
    use sort
    use utils
    implicit none
    ! Gaussian center index mapping
    integer (KIND=IPREC),intent(out) :: map(3), invmap(3)
    ! Nuclear charges
    integer (KIND=IPREC),intent(out) :: iZ(3)
    ! x, y and z coordinates
    double precision, intent(out) :: x(3), y(3), z(3)
    ! Number of centers
    integer (KIND=IPREC),intent(out) :: ncen

    ! Helper for i/o
    character(80) :: line
    ! Input center number
    integer (KIND=IPREC) :: icen(3)
    ! Atomic number
    integer (KIND=IPREC) :: iZ0(3)
    ! Atom type
    integer (KIND=IPREC) :: itype(3)
    ! X, Y and Z coordinates
    double precision :: x0(3), y0(3), z0(3)
    ! Is atom a dummy atom?
    logical :: idum(3)

    ! Loop index, number of centers, read status, index of dummy atom, i/o status
    integer (KIND=IPREC) :: i, j, stat
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
          if(stat.ne.0) then
             write (*,*) 'Error parsing ', line
             stop
          end if
       else
          read (line,*,iostat=stat) icen(i), iZ0(i), itype(i), x0(i), y0(i), z0(i)
          if(stat .ne. 0) then
             ! Try switching to g94 mode
             g94=.true.
             read (line,*,iostat=stat) icen(i), iZ0(i), x0(i), y0(i), z0(i)
             if(stat.ne.0) then
                write (*,*) 'Error parsing ', line
                stop
             end if
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

  end subroutine makeMapping

  subroutine readBasis(nprim,nexpon,npbasis)
    implicit none
    integer (KIND=IPREC),intent(out) :: nprim(3), nexpon, npbasis
    integer (KIND=IPREC) :: ibc, ibp, istop, ncen

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
  end subroutine readBasis

  subroutine prepareGaussian
    use discrete, only : z1, z2, r
    use params
    use commons
    use utils

    implicit none
    integer (KIND=IPREC) :: i,ib,ibc,ibp,icount,ifbo,jfbo,ilines,iorb,&
         l1,m1,m1abs,nexpon,nfborb

    integer (KIND=IPREC),dimension(maxorb) :: ifbord,ifdord

    real (PREC) :: d1
    real (PREC), dimension(maxbasis) :: ccoeff
    ! We can have up to 2x more orbitals in Gaussian than in the 2d
    ! calculation because non-sigma orbitals appear in twos.
    !real (PREC), dimension(2*maxorb,maxbasis) :: 
    real (PREC), dimension(:,:), allocatable :: excoeff

    ! Orbital character
    real (PREC), dimension(0:3) :: ochar
    real (PREC) :: ocharmax
    integer (KIND=IPREC) :: ichar, icharmax

    character*46 :: matchstr

    ! Gaussian center index
    integer (KIND=IPREC) :: map(3), invmap(3)
    ! Nuclear charges
    integer (KIND=IPREC) :: iZ(3)
    ! x, y and z coordinates
    double precision :: x(3), y(3), z(3)
    ! Number of centers
    integer (KIND=IPREC) :: ncen

    ! Number of primitives on the centers
    integer (KIND=IPREC) :: nprim(3)
    ! Running index
    integer (KIND=IPREC) :: ibf
    ! Target atom
    integer (KIND=IPREC) :: itgt

    ! Basis function start and end
    integer (KIND=IPREC) :: bfstart(3), bfend(3)
    ! Found atoms and basis?
    logical :: atoms_found, basis_found

    !     nforb  - number of finite basis (fb) set orbitals
    !     ifbord - ordering of fb orbitals (0=sigma, 1=pi, 2=delta, etc)
    !              the order is determined from the gaussian.out file
    !     ifdord - ordering of finite difference (fd) orbitals
    !              the order of fd orbitals is the same as in the input data
    !              (phi orbitals first, then delta, etc)

    allocate (excoeff(2*maxorb,maxbasis))
    
    ! Open the GaussianXY output file
    open(7,file='gaussian.out', status='old',form='formatted')

    atoms_found = .false.
    basis_found = .false.

    do ilines=1,100000
       read(7,1001,end=990,err=910) matchstr

       ! Gaussian 94
       if (matchstr.eq.'                  Z-Matrix orientation:' .or. &
            matchstr.eq.'                   Standard orientation:') then
          call makeMapping(map,invmap,iZ,x,y,z,ncen)
          atoms_found = .true.
       end if
       ! Gaussian 03/09/16
       if (matchstr.eq.'                         Standard orientation:' .or. &
            matchstr.eq.'                          Input orientation:') then
          call makeMapping(map,invmap,iZ,x,y,z,ncen)
          atoms_found = .true.
       end if
       if (matchstr.eq.' Basis set in the form of general basis input:' .or. &
            matchstr.eq.' AO basis set in the form of general basis inp') then
          call readBasis(nprim,nexpon,npbasis)
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
       write (*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       write (*,*) "Gaussian calculation is inconsistent with x2dhf for atom on the left"
       write (*,'(A,F4.1,A,I3)') ' z1=',z1,', iZ(1)=',iZ(1)
       write (*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    end if
    if(abs(z2-iZ(ncen)).gt.precis) then
       write (*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       write (*,*) "Gaussian calculation is inconsistent with x2dhf for atom on the right"
       write (*,'(A,F4.1,A,I1,A,I3)') ' z2=',z2,', iZ(',ncen,')=',iZ(ncen)
       write (*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
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

    ! ! Retrieve expansion coefficients of the contracted gaussians
    ! if (iprint1.ne.0) then
    !    print *,'no. of basis function  no. of exp. coeff.'
    !    do ibp=1,npbasis
    !       ibc=ixref(ibp)
    !    enddo
    ! endif

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
       ochar = 0.0_PREC
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
       fngau2(ib)=( d1**(2*l1+3)*&
            2.0_PREC**(4*l1+7)/pii/(factor2(l1+1))**2)**0.250_PREC
       ! Normalization for the angular part (spherical harmonic)
       shngau(ib)=(-1.0_PREC)**m1*&
            sqrt( ((2.0_PREC*l1+1.0_PREC)*factor(l1-m1abs)) / (4.0_PREC*pii*factor(l1+m1abs)) )
    enddo

    deallocate (excoeff)
    
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
  end subroutine prepareGaussian
end module prepGauss
