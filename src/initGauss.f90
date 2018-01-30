! ***************************************************************************
! *                                                                         *
! *   Copyright (C) 1996-2010 Jacek Kobus <jkob@fizyka.umk.pl>              *
! *                                                                         *
! *   This program is free software; you can redistribute it and/or modify  *
! *   it under the terms of the GNU General Public License version 2 as     *
! *   published by the Free Software Foundation.                            *
! *                                                                         *
! ***************************************************************************
! ### initGauss ###

!     Molecular orbitals are initialized as a linear combination of
!     Gaussian orbitals. Parameters of the orbitals and the coeeficients
!     are provided by the GAUSSIAN program.

!     Coulomb (HF) potentials are initialized as a linear combination of
!     -ez1/r1 and -ez2/r2. For the initialization of exchange potentials
!     see routine tfpot.

subroutine initGauss (psi,pot,excp,f2,f4,wgt2,wk0)
  use params
  use discret
  use commons8
  use evalGauss

  implicit none
  integer :: i,igauss,igp,igrid,imu,in,inioff,iorb,ipb,ishift,m1, &
       ngrid,norbt
  real (PREC) :: c1,xnorm
  real (PREC), dimension(*) :: psi,pot,excp,f2,f4,wgt2,wk0
  real (PREC), dimension(:,:), allocatable :: bf
  integer :: it, jt

  !     Initialization of molecular orbitals
  if (ini.eq.4.) then
     norbt=1
  else
     norbt=norb
  endif

  ! Evaluate basis functions
  allocate(bf(nni*mxnmu,npbasis))
  call evaluate(bf)

  ! Evaluate overlap matrix over gaussian basis functions
  igauss=0
  if (idbg(562).ne.0) igauss=1
  if (igauss.ne.0) then
     write (*,*) 'Test finite basis function overlaps on grid'
     do it=1,npbasis
        !        do jt=1,it
        do jt=it,it
           xnorm=0.0
           do i=1,nni*mxnmu
              xnorm=xnorm+bf(i,it)*bf(i,jt)*wgt2(i)*f4(i)
           end do
           write (*,'(A,I3,A,I3,A,ES14.7)') 'S(',it,',',jt,') = ',xnorm
        end do
     end do
  end if

  write(*,1114)
  do iorb=1,norbt
     ishift=i1b(iorb)-1
     ngrid= i1si(iorb)
     igp=ishift
     do i=1,ngrid
        psi(ishift+i)=0.0_PREC
     enddo

     do ipb=1,npbasis
        m1 = abs(mprim(ipb))
        ! Skip functions that have different m value
        if (m1.ne.mm(iorb)) cycle

        c1 =primcoef(iorb,ipb)

        if (abs(c1).gt.0.0_PREC) then
           ! Loop over grid
           do imu=1,mxnmu
              inioff=(imu-1)*nni
              do in=1,nni
                 igrid=inioff+in
                 igp=ishift+igrid
                 psi(igp)=psi(igp)+c1*bf(igrid,ipb)
              end do
           end do
        end if
     end do
     call norm94 (iorb,psi,f4,wgt2,wk0,xnorm)
     write (*,1115) iorn(iorb),bond(iorb),gut(iorb),xnorm
  end do

  deallocate(bf)

  !  if (iprt.eq.1) write (*,1115) iorn(iorb),bond(iorb),gut(iorb),xnorm

  write(*,*)

  if (idbg(560).ne.0) stop 'inigauss'

  !     initialize Coulomb and exchange potentials

  call initPot(psi,pot,excp,f2,f4,wk0)

1114 format(/1x,'    orbital        norm      ')
1115 format(1x,i3,1x,a8,a1,e20.12)

end subroutine initGauss
