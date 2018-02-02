! ***************************************************************************
! *                                                                         *
! *   Copyright (C) 1997-2010 Jacek Kobus <jkob@fizyka.umk.pl>              *
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

subroutine prepGauss
  use params
  use commons8

  implicit none
  integer :: i,ib,ibc,ibp,icent,icount,ifbo,jfbo,ilines,iorb,ipbasis,iprint0,iprint1,iprint2,istop, &
       l1,m1,m1abs,n1prim,n2prim,n3prim,nexpon,nfborb

  integer, dimension(60) :: ifbord,ifdord

  real (PREC) :: d1
  real (PREC), dimension(maxbasis) :: ccoeff
  real (PREC), dimension(0:60,0:maxbasis) :: excoeff
  real (PREC), external :: factor,factor2

  ! Orbital character
  real (PREC), dimension(0:3) :: ochar
  real (PREC) :: ocharmax
  integer :: ichar, icharmax

  character*46 :: matchstr

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

  do ilines=1,1000
     read(7,1001,end=990,err=910) matchstr
     ! Gaussian 94
     if (matchstr.eq.' Basis set in the form of general basis input:') goto 900
     ! Gaussian 03
     if (matchstr.eq.' AO basis set in the form of general basis inp') goto 900
  enddo
  write(*,*) 'prepGauss: '
  write(*,*) ' "AO basis set in the form of general basis input not found in GAUSSIANxy output file'
  stop 'prepGauss'
00900 continue

  ipbasis=1
  npbasis=0

  !     start extracting data from the output

  !     ibc counts contracted basis functions
  !     ibp counts primitive basis functions
  !     n1prim - number of primitive gaussian function at the centre A
  !     n3prim - number of primitive gaussian function at the centre B
  !     n2prim - number of primitive gaussian function at the bond centre

  ibc=0
  ibp=0
  n1prim=0
  n2prim=0
  n3prim=0
  call rexponents(ibc,ibp,istop)
  n1prim=ibp
  if (istop.eq.0) then
     call rexponents(ibc,ibp,istop)
     n3prim=ibp-n1prim
     if (istop.eq.0) then
        icent=2
        call rexponents(ibc,ibp,istop)
        n2prim=ibp-n3prim-n1prim
     endif
  endif

00910 continue

  nexpon=ibc
  npbasis=ibp
  write(*,1140) nexpon,npbasis,n1prim,n2prim,n3prim

  close(7)

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

  ! Prepare expansion coefficients of a molecular orbital in the basis
  ! set of primitive gaussian functions
  do ib=1,npbasis
     if (ib.le.n1prim) icent=1
     if (ib.gt.n1prim.and.ib.le.(n1prim+n2prim)) icent=2
     if (ib.gt.(n1prim+n2prim).and.ib.le.(n1prim+n2prim+n3prim)) icent=3
     icgau(ib)=icent
     !     write(*,1050) ib,ixref(ib),lprim(ib),mprim(ib),icgau(ib),primexp(ib)
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
     ! write (*,*) 'fb orbital ',ifbo,' symmetry is ',icharmax
  end do

  ! Initialize memory
  do iorb=1,norb
     do ib=1,npbasis
        primcoef(iorb,ib)=0.0_PREC
     end do
  end do

  ! Associate finite basis orbitals with the corresponding fd ones
  do iorb=1,norb
     do ifbo=1,nfborb
        icount=0
        ! Check if orbital characters match
        if (ifdord(iorb).eq.ifbord(ifbo)) then
           !           write (*,'(A,I3,A,I3)') 'Orbital ',iorb,' corresponds to Gaussian orbital ',ifbo
           do ib=1,npbasis
              primcoef(iorb,ib)=excoeff(ifbo,ib)
           enddo
           ! If this is not a sigma orbital, we need to zero out the degenerate orbital as well.
           if (ifbord(ifbo).gt.0) then
              do jfbo=ifbo+1,nfborb
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
01140 format(/15x,i4,' exponents ',                         &
       /15x,i4,' primitive basis functions ',      &
       /15x,i4,' primitives on centre 1'           &
       /15x,i4,' primitives on centre 2'           &
       /15x,i4,' primitives on centre 3'/)
01142 format(15x,i4,' finite difference orbitals',/15x,i4,' finite basis set orbitals'/)
end subroutine prepGauss
