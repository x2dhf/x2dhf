! ***************************************************************************
! *                                                                         *
! *   Copyright (C) 1996 Leif Laaksonen, Dage Sundholm                      *
! *   Copyright (C) 1996-2010 Jacek Kobus <jkob@fizyka.umk.pl>              *
! *                                                                         *
! *   This program is free software; you can redistribute it and/or modify  *
! *   it under the terms of the GNU General Public License version 2 as     *
! *   published by the Free Software Foundation.                            *
! *                                                                         *
! ***************************************************************************
! ### writeDisk4dd ###

!     Writes orbitals, potenials, Lagrange multipliers (diagonal
!     and off-diagonal) and multipole expansion coefficients to a disk
!     file in either formatted or unformatted form

subroutine writeDisk4dd (psi,pot,f0,f4,wk0,wk1,wk2,wk3)
  use params
  use discret
  use scf
  use solver
  use commons8

  implicit none
  integer :: i4,i1beg,ig,iorb,ngrid,isym4nu,isym4mu,isymmetry

  real (PREC), allocatable :: r8mxsize(:),r8mxsize1(:)
  real (PREC), dimension(*) :: psi,pot,f0,f4,wk0,wk1,wk2,wk3

  open(9999,file='out4dd.dat', status='replace',form='formatted')

!  call getDateTime(datetime)
  write(9999,'(a80)') header
!  write(9999,'(a80)') datetime
  write(9999,'(" nnu nmu ")')
  write(9999,formint) nni,nmu(1)
  write(9999,'(" r  rinf")')
  write(9999,formfp64) r,rinf
  write(9999,'(" z1 z2 ")')
  write(9999,formfp64) z1,z2
  write(9999,'(" norb nel  ")')
  write(9999,formint) norb,nel


  iorb=1
  write(9999,'(" ordering of orbitals ")')
  write(9999,'(1x,i3,1x,a8,1x,a1)') (iorn(iorb),bond(iorb),gut(iorb),iorb=1,norb)

  write(9999,'(" orbital energies: ")')
  do iorb=1,norb
     write(9999,formfp64) eng(iorb)
  enddo

  allocate(r8mxsize(nni*mxnmu))
  allocate(r8mxsize1(nni*mxnmu))

  do iorb=1,norb
     i1beg=i1b(iorb)
     ngrid=i1si(iorb)

     do i4=1,ngrid
        r8mxsize(i4)=pot(i1beg+i4-1)
     enddo

     if (imethod.eq.2.or.ini.eq.4) then
        do i4=1,ngrid
           r8mxsize(i4)=zero
        enddo
     endif

     ! multiply \tilde{V}_C by 2/(R\xi) to get V_C

     !     isymmetry=isymOrb(iorb)
     ! orbital 1sigma is even, its derivatives are odd functions
     isymmetry=-1

     ig=1
     call pot2pot(r8mxsize,f4)

     write(9999,'(/" orbital=",i3,1x,a8,1x,a1," Coulomb potential" )') iorn(iorb),bond(iorb),gut(iorb)
     call prtmatcw(nni,mxnmu,r8mxsize,9999)

     !     write(6,'(/" orbital=",i3,1x,a8,1x,a1," Coulomb potential" )') iorn(iorb),bond(iorb),gut(iorb)
!     write(9999,'(/" orbital=",i3,1x,a8,1x,a1," Coulomb potential" )') iorn(iorb),bond(iorb),gut(iorb)
     ! call prtmatrw(nni,mxnmu,r8mxsize,9999)

     isym4nu=1
     isym4mu=1
     call putin1 (nni,mxnmu,isym4nu,isym4mu,r8mxsize,wk3)
     call diff1nu (mxnmu,wk3,wk0,wk1,wk2)
     call putout (nni,mxnmu,r8mxsize,wk0)
     r8mxsize1=r8mxsize

     write(9999,'(/" orbital=",i3,1x,a8,1x,a1," / Coulomb potential nu derivative" )') iorn(iorb),bond(iorb),gut(iorb)
     call prtmatcw(nni,mxnmu,r8mxsize,9999)
!     write(6,'(/" orbital=",i3,1x,a8,1x,a1," / Coulomb potential nu derivative" )') iorn(iorb),bond(iorb),gut(iorb)
!     call prtmatrw(nni,mxnmu,r8mxsize,9999)

     call diff1mu (mxnmu,wk3,wk2)
     call putout (nni,mxnmu,r8mxsize,wk2)

     write(9999,'(/" orbital=",i3,1x,a8,1x,a1," / Coulomb potential mu derivative" )') iorn(iorb),bond(iorb),gut(iorb)
     call prtmatcw(nni,mxnmu,r8mxsize,9999)
!     write(6,'(/" orbital=",i3,1x,a8,1x,a1," / Coulomb potential mu derivative" )') iorn(iorb),bond(iorb),gut(iorb)
!     call prtmatrw(nni,mxnmu,r8mxsize,9999)

     isym4nu=-1
     isym4mu=1
     call putin1 (nni,mxnmu,isym4nu,isym4mu,r8mxsize1,wk3)

     call diff1mu (mxnmu,wk3,wk2)
     call putout (nni,mxnmu,r8mxsize,wk2)

     write(9999,'(/" orbital=",i3,1x,a8,1x,a1," / Coulomb potential nu,mu derivative" )') iorn(iorb),bond(iorb),gut(iorb)
     call prtmatcw(nni,mxnmu,r8mxsize,9999)
!     write(6,'(/" orbital=",i3,1x,a8,1x,a1," / Coulomb potential nu,mu derivative" )') iorn(iorb),bond(iorb),gut(iorb)
!     call prtmatrw(nni,mxnmu,r8mxsize,9999)

     do i4=1,ngrid
        r8mxsize(i4)=psi(i1beg+i4-1)
     enddo

!1000 format(i5,i4,1x,a8,a1)
     write(9999,'(/" orbital=",i3,1x,a8,1x,a1," orbital function")') iorn(iorb),bond(iorb),gut(iorb)
     call prtmatcw(nni,mxnmu,r8mxsize,9999)
!     write(6,'(/" orbital=",i3,1x,a8,1x,a1," orbital function")') iorn(iorb),bond(iorb),gut(iorb)
!     call prtmatrw(nni,mxnmu,r8mxsize,9999)

     ! printout suitable for comparison with JM data
     ! call prtmatcw1(nni,mxnmu,r8mxsize,9999)

     isym4nu=1
     isym4mu=1
     call putin1 (nni,mxnmu,isym4nu,isym4mu,psi(i1beg),wk3)
     call diff1nu (mxnmu,wk3,wk0,wk1,wk2)
     call putout (nni,mxnmu,r8mxsize,wk0)
     r8mxsize1=r8mxsize

     write(9999,'(/" orbital=",i3,1x,a8,1x,a1," / nu derivative" )') iorn(iorb),bond(iorb),gut(iorb)
     call prtmatcw(nni,mxnmu,r8mxsize,9999)
!     write(6,'(/" orbital=",i3,1x,a8,1x,a1," / nu derivative" )') iorn(iorb),bond(iorb),gut(iorb)
!     call prtmatrw(nni,mxnmu,r8mxsize,9999)

     call diff1mu (mxnmu,wk3,wk2)
     call putout (nni,mxnmu,r8mxsize,wk2)

     write(9999,'(/" orbital=",i3,1x,a8,1x,a1," / mu derivative" )') iorn(iorb),bond(iorb),gut(iorb)
     call prtmatcw(nni,mxnmu,r8mxsize,9999)
!     write(6,'(/" orbital=",i3,1x,a8,1x,a1," / mu derivative" )') iorn(iorb),bond(iorb),gut(iorb)
!     call prtmatrw(nni,mxnmu,r8mxsize,9999)

     isym4nu=-1
     isym4mu=1
     call putin1 (nni,mxnmu,isym4nu,isym4mu,r8mxsize1,wk3)
     call diff1mu (mxnmu,wk3,wk2)
     call putout (nni,mxnmu,r8mxsize,wk2)

     write(9999,'(/" orbital=",i3,1x,a8,1x,a1," / nu,mu derivative" )') iorn(iorb),bond(iorb),gut(iorb)
     call prtmatcw(nni,mxnmu,r8mxsize,9999)
!     write(6,'(/" orbital=",i3,1x,a8,1x,a1," / nu,mu derivative" )') iorn(iorb),bond(iorb),gut(iorb)
!     call prtmatrw(nni,mxnmu,r8mxsize,9999)
  enddo

  do i4=1,ngrid
     r8mxsize(i4)=f0(i4)
  enddo

  write(9999,'(/" f0=")')
  call prtmatcw(nni,mxnmu,r8mxsize,9999)
  isym4nu=1
  isym4mu=1
  call putin1 (nni,mxnmu,isym4nu,isym4mu,r8mxsize,wk3)
  call diff1nu (mxnmu,wk3,wk0,wk1,wk2)
  call putout (nni,mxnmu,r8mxsize,wk0)
  r8mxsize1=r8mxsize

  call checkd1nu(nni,mxnmu,r8mxsize)

  write(9999,'(/" f0=",i3,1x,a8,1x,a1," / nu derivative" )') iorn(iorb),bond(iorb),gut(iorb)
  call prtmatcw(nni,mxnmu,r8mxsize,9999)

  call diff1mu (mxnmu,wk3,wk2)
  call putout (nni,mxnmu,r8mxsize,wk2)

  call checkd1mu(nni,mxnmu,r8mxsize)

  write(9999,'(/" f0=",i3,1x,a8,1x,a1," / mu derivative" )') iorn(iorb),bond(iorb),gut(iorb)
  call prtmatcw(nni,mxnmu,r8mxsize,9999)

  isym4nu=-1
  isym4mu=1
  call putin1 (nni,mxnmu,isym4nu,isym4mu,r8mxsize1,wk3)
  call diff1mu (mxnmu,wk3,wk2)
  call putout (nni,mxnmu,r8mxsize,wk2)

  write(9999,'(/" f0=",i3,1x,a8,1x,a1," / nu,mu derivative" )') iorn(iorb),bond(iorb),gut(iorb)
  call prtmatcw(nni,mxnmu,r8mxsize,9999)

  close(9999)

end subroutine writeDisk4dd


! ***************************************************************************
! *                                                                         *
! *   Copyright (C) 1996 Leif Laaksonen, Dage Sundholm                      *
! *   Copyright (C) 1996-2010 Jacek Kobus <jkob@fizyka.umk.pl>              *
! *                                                                         *
! *   This program is free software; you can redistribute it and/or modify  *
! *   it under the terms of the GNU General Public License version 2 as     *
! *   published by the Free Software Foundation.                            *
! *                                                                         *
! ***************************************************************************
! ### writeDisk4dd ###

!     Transforms \tilde{V}_C into V_C (F4=R\xi/2)

subroutine pot2pot (wk0,f4)
!  use params
  use discret
!  use scf
!  use commons8

  implicit none
  integer :: inu,imu
  real (PREC),  dimension(nni,mxnmu) :: wk0,f4


  do inu=1,nni
     do imu=1,mxnmu
        wk0(inu,imu)=wk0(inu,imu)/f4(inu,imu)
     enddo
  enddo
end subroutine pot2pot


subroutine setBVpot(fun)
  use params
  use discret
  use solver
  use commons8

  implicit none
  integer :: j,ij,nnu1,nnu2,nnu3,nnu4,nnu5
  real (PREC), dimension(*) :: fun
  real (PREC) :: fun1

  nnu1=nni
  nnu2=nnu1+nnu1
  nnu3=nnu2+nnu1
  nnu4=nnu3+nnu1
  nnu5=nnu4+nnu1

  if (isym.eq.1) then
     !       if (ifill.eq.1.or.ifill.eq.2) then
!      do i=1,ngrd7
!         ij=indx7(i)
!         fun1=fun(ij)
!         fun(ij)=exeven(1)*fun(ij+nnu1)+exeven(2)*fun(ij+nnu2)+   &
!              exeven(3)*fun(ij+nnu3)+exeven(4)*fun(ij+nnu4)+exeven(5)*fun(ij+nnu5)
!         print *,'  1 ',i,ij,fun1,fun(ij),fun1-fun(ij)
! !        print *,'  1 ',fun(ij+nnu1),fun(ij+nnu2),fun(ij+nnu3),fun(ij+nnu4)
!      enddo

     do ij=2,nnu1-1
        fun1=fun(ij)
        fun(ij)=exeven(1)*fun(ij+nnu1)+exeven(2)*fun(ij+nnu2)+   &
             exeven(3)*fun(ij+nnu3)+exeven(4)*fun(ij+nnu4)+exeven(5)*fun(ij+nnu5)
!         print *,'  1 ',ij,fun1,fun(ij),fun1-fun(ij)
! !        print *,'  1 ',fun(ij+nnu1),fun(ij+nnu2),fun(ij+nnu3),fun(ij+nnu4)
      enddo


     do j=1,mxnmu
        ij=(j-1)*nnu1+1
        fun1=fun(ij)
        fun(ij)=exeven(1)*fun(ij+1)+exeven(2)*fun(ij+2)+  &
             exeven(3)*fun(ij+3)+exeven(4)*fun(ij+4)+exeven(5)*fun(ij+5)
!        print *,'  2 ',ij,fun1,fun(ij),fun1-fun(ij)
!        write(*,'(" 2",i4,6e15.6)') ij,fun(ij),fun(ij+1),fun(ij+2),fun(ij+3),fun(ij+4),fun(ij+5)
     enddo


     do j=1,mxnmu
        ij=j*nnu1
        fun1=fun(ij)
        fun(ij)=exeven(1)*fun(ij-1)+exeven(2)*fun(ij-2)+  &
             exeven(3)*fun(ij-3)+exeven(4)*fun(ij-4)+exeven(5)*fun(ij-5)
!        print *,'  3 ',ij,fun1,fun(ij),fun1-fun(ij)
     enddo
!     stop 'writeDisk4dd'

  else
     do ij=2,nnu1-1
        fun(ij)=0.0_PREC
      enddo

     do j=1,mxnmu
        ij=(j-1)*nnu1+1
        fun(ij)=0.0_PREC
     enddo

     do j=1,mxnmu
        ij=j*nnu1
        fun(ij)=0.0_PREC
     enddo
  endif

end subroutine setBVpot

subroutine setBVgen(isym4nu,isym4mu,fun)
  use params
  use discret
  use solver
  use commons8

  implicit none
  integer :: j,ij,isym4nu,isym4mu,nnu1,nnu2,nnu3,nnu4,nnu5
  real (PREC), dimension(*) :: fun
  real (PREC) :: fun1

  nnu1=nni
  nnu2=nnu1+nnu1
  nnu3=nnu2+nnu1
  nnu4=nnu3+nnu1
  nnu5=nnu4+nnu1

  if (isym4nu.eq.1) then
     do ij=2,nnu1-1
        fun1=fun(ij)
        fun(ij)=exeven(1)*fun(ij+nnu1)+exeven(2)*fun(ij+nnu2)+   &
             exeven(3)*fun(ij+nnu3)+exeven(4)*fun(ij+nnu4)+exeven(5)*fun(ij+nnu5)
        !         print *,'  1 ',ij,fun1,fun(ij),fun1-fun(ij)
        ! !        print *,'  1 ',fun(ij+nnu1),fun(ij+nnu2),fun(ij+nnu3),fun(ij+nnu4)
     enddo
  else
     do ij=2,nnu1-1
        fun(ij)=0.0_PREC
     enddo
  end if

  if (isym4mu.eq.1) then
     do j=1,mxnmu
        ij=(j-1)*nnu1+1
        fun1=fun(ij)
        fun(ij)=exeven(1)*fun(ij+1)+exeven(2)*fun(ij+2)+  &
             exeven(3)*fun(ij+3)+exeven(4)*fun(ij+4)+exeven(5)*fun(ij+5)
        !        print *,'  2 ',ij,fun1,fun(ij),fun1-fun(ij)
        !        write(*,'(" 2",i4,6e15.6)') ij,fun(ij),fun(ij+1),fun(ij+2),fun(ij+3),fun(ij+4),fun(ij+5)
     enddo

     do j=1,mxnmu
        ij=j*nnu1
        fun1=fun(ij)
        fun(ij)=exeven(1)*fun(ij-1)+exeven(2)*fun(ij-2)+  &
             exeven(3)*fun(ij-3)+exeven(4)*fun(ij-4)+exeven(5)*fun(ij-5)
        !        print *,'  3 ',ij,fun1,fun(ij),fun1-fun(ij)
     enddo
  else
     do j=1,mxnmu
        ij=(j-1)*nnu1+1
        fun(ij)=0.0_PREC
     enddo

     do j=1,mxnmu
        ij=j*nnu1
        fun(ij)=0.0_PREC
     enddo
  endif

end subroutine setBVgen


subroutine checkd1nu (m,n,a)
  use params
  use discret

  implicit none
  integer :: im,in,ioutmat,m,n
  integer :: immax,inmax
  real (PREC), dimension(m,n) :: a
  real (PREC) :: diff,vn

  diff=0.0_PREC
  immax=0
  inmax=0

  do im=1,m
     do in=1,n
        vn=-veta1(im)*(z2-z1)*r
!        write(ioutmat,'("checkd1nu: ",2i5,3e15.6)') im,in,vn,a(im,in),abs(vn-a(im,in))
        if (abs(vn-a(im,in)).gt.diff) then
           diff=abs(vn-a(im,in))
           immax=im
           inmax=in
        endif
     enddo
  enddo

  write(ioutmat,'("checkd1nu: ",2i5,e15.6)') immax,inmax,diff

end subroutine checkd1nu


subroutine checkd1mu (m,n,a)
  use params
  use discret

  implicit none
  integer :: im,in,ioutmat,m,n
  integer :: immax,inmax
  real (PREC), dimension(m,n) :: a
  real (PREC) :: diff,vn

  diff=0.0_PREC
  immax=0
  inmax=0

!   do im=1,m
!      do in=n-3,n
!          a(im,in)=a(im,n-4)
! !        a(im,in)=0.0_PREC
!      enddo
!   enddo

  do im=1,m
     do in=1,n
        vn=vxi1(in)*(z1+z2)*r
        write(ioutmat,'("checkd1mu: ",2i5,3e15.6)') im,in,vn,a(im,in),abs(vn-a(im,in))
        if (abs(vn-a(im,in)).gt.diff) then
           diff=abs(vn-a(im,in))
           immax=im
           inmax=in
        endif
     enddo
  enddo

  write(ioutmat,'("checkd1m: ",2i5,e15.6)') immax,inmax,diff

end subroutine checkd1mu

