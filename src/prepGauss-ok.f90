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
  integer :: i,ib,ibc,ibp,icent,icount,ifbo,ilines,iorb,ipbasis,iprint0,iprint1,iprint2,istop, &
       l1,m1,m1abs,n1prim,n2prim,n3prim,nexpon,nfborb

  integer, dimension(60) :: ifbord,ifdord

  real :: d1,oe,symthresh
  real (PREC), dimension(maxbasis) :: ccoeff
  real (PREC), dimension(0:60,0:maxbasis) :: excoeff
  real (PREC), external :: factor,factor2

  character*46 :: matchstr

  !     nforb  - number of finite basis (fb) set orbitals
  !     ifbord - ordering of fb orbitals (0=sigma, 1=pi, 2=delta, etc)
  !              the order is determined from the gaussian94.output file
  !     ifdord - ordering of finite difference (fd) orbitals
  !              the order of fd orbitals is the same as in the input data
  !              (phi orbitals first, then delta, etc)

  iprint0=0
  iprint1=0
  iprint2=0

  if (iprint(553).ne.0) iprint0=1
  if (iprint(554).ne.0) iprint1=1
  if (iprint(555).ne.0) iprint2=1

  !  opening a file with the corresponding GAUSSSIAN94/98 output

  open(7,file='gaussian.out', status='old',form='formatted')

  do ilines=1,1000
     read(7,1001,end=990,err=910) matchstr
     if (matchstr.eq.' Basis set in the form of general basis input:') goto 900
  enddo
  write(*,*) 'prepGauss: '
  write(*,*) ' "Basis set in the form of general basis input: not found in GAUSSIAN9x output file'
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
  read(7,1001,end=991,err=992)

  !     start extracting basis set expansion coefficients from
  !     gaussian94.pun file

  !     determine number of molecular orbitals as defined in the GAUSSIAN9x
  !     program (one nonsigma finite difference orbital corresponds
  !     to two GAUSSIAN9x ones

  nfborb=0
  do iorb=1,norb
     ifdord(iorb)=mm(iorb)
     if (mm(iorb).eq.0) then
        nfborb=nfborb+1
     else
        nfborb=nfborb+2
     endif
  enddo
  write(*,1142) norb,nfborb

  !     retrieve expansion coefficients of the contracted gaussians

  if (iprint1.ne.0) then
     print *,'no. of basis function  no. of exp. coeff.'
     do ibp=1,npbasis
        ibc=ixref(ibp)
     enddo
  endif

  do ifbo=1,nfborb
     read(7,1011,end=991,err=992) oe
     do ibc=1,nexpon,5
        read(7,1012,end=991,err=992) (ccoeff(ibc+i),i=0,4)
     enddo

     !        transform expansion coefficients of the contracted gaussians
     !        into coefficients of primitive gaussians

     do ibp=1,npbasis
        ibc=ixref(ibp)
        excoeff(ifbo,ibp)=ccoeff(ibc)*coeff(ibp)
        !           write(*,'(3i5,3e16.6,i5)') ifbo,ibp,ibc,ccoeff(ibc),
        !           &           coeff(ibp),excoeff(ifbo,ibp),mprim(ibp)
     enddo
  enddo

  !     prepare expansion coefficients of a molecular orbital in the basis
  !     set of primitive gaussian functions

  do ib=1,npbasis
     if (ib.le.n1prim) icent=1
     if (ib.gt.n1prim.and.ib.le.(n1prim+n2prim)) icent=2
     if (ib.gt.(n1prim+n2prim).and.ib.le.(n1prim+n2prim+n3prim)) icent=3
     icgau(ib)=icent
     !        write(*,1050) ib,ixref(ib),lprim(ib),mprim(ib),icgau(ib),
     !        &              primexp(ib)
  enddo

  !     determine symmetry of GAUSSIAN9x molecular orbitals

  symthresh=0.000010_PREC
  icount=0
  do ifbo=nfborb,1,-1
     do ib=npbasis,1,-1
        if (abs(excoeff(ifbo,ib)).gt.symthresh) then
           if(icount.ne.ifbo) then
              icount=ifbo
              if (iprint2.ne.0) print *,'ifbo,ib,mprim(ib),ibc ',ifbo,ib,mprim(ib),ixref(ib)
              ifbord(ifbo)=abs(mprim(ib))
           endif
        endif
     enddo
  enddo

  !     fb orbitals are being associated with the corresponding fd ones

  do iorb=1,norb
     do ib=1,npbasis
        primcoef(iorb,ib)=0.0_PREC
     enddo
  enddo

  do iorb=1,norb
     do ifbo=nfborb,1,-1
        icount=0
        if (ifdord(iorb).eq.ifbord(ifbo)) then
           do ib=1,npbasis
              primcoef(iorb,ib)=excoeff(ifbo,ib)
           enddo
           if (ifbord(ifbo).gt.0) ifbord(ifbo-1)=-1
           ifbord(ifbo)=-1
           ifdord(iorb)=-2
        endif
     enddo
  enddo

  !     normalization factor for a spherical harmonic Gaussian-type functions

  do ib=1,npbasis
     d1=primexp(ib)
     l1=lprim  (ib)
     m1=mprim  (ib)
     m1abs=abs(m1)
     fngau2(ib)=( d1**dble(2*l1+3) * 2.0_PREC**dble(4*l1+7)/pii/(factor2(2*l1+1))**2)**0.250_PREC

     !        normalization factor for spherical harmonics

     shngau(ib)=(-1.0_PREC)**dble((m1+m1abs)/2) /sqrt(4.0_PREC*pii)* &
          sqrt((2*l1+1)*factor(l1-m1abs)/factor(l1+m1abs))
  enddo

  return
00990 print *,'PREPGAUSSIAN: end of file encountered when reading gaussian.out'
  stop
00991 print *,'PREPGAUSSIAN: end of file encountered when reading gaussian.pun'
  stop
00992 print *,'PREPGAUSSIAN: error encountered when reading gaussian.pun'

01001 format(a46)
01011 format(15x,e15.8)
01012 format(5e15.8)
01050 format(5i5,e15.8)
01140 format(/15x,i4,' exponents ',                         &
                /15x,i4,' primitive basis functions ',      &
                /15x,i4,' primitives on centre 1'           &
                /15x,i4,' primitives on centre 2'           &
                /15x,i4,' primitives on centre 3'/)
01142 format(15x,i4,' finite difference orbitals',/15x,i4,' finite basis set orbitals'/)
end subroutine prepGauss

