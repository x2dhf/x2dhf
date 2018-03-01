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

  ! Values of contracted basis functions
  real (PREC), dimension(:,:), allocatable :: cbf
  ! m values of contracted basis functions
  integer, dimension(:), allocatable :: cm
  ! Contracted basis function overlap
  real (PREC), dimension(:,:), allocatable :: cS

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
     write (*,*) 'Test finite basis primitive function overlaps on grid'
     do it=1,npbasis
        do jt=1,it
           ! Skip different m values
           if(mprim(it) .ne. mprim(jt)) cycle
           xnorm=0.0
           do i=1,nni*mxnmu
              xnorm=xnorm+bf(i,it)*bf(i,jt)*wgt2(i)*f4(i)
           end do
!           if(abs(xnorm) .ge. 1e-6) then
              write (*,'(A,I3,A,I3,A,ES14.7)') 'S(',it,',',jt,') = ',xnorm
!           end if
        end do
     end do

     write (*,*) 'Test finite basis contracted function overlaps on grid'
     ! Evaluate contracted basis function values on grid
     allocate(cbf(nni*mxnmu,ixref(npbasis)))
     allocate(cm(ixref(npbasis)))
     allocate(cS(ixref(npbasis),ixref(npbasis)))
     do it=1,ixref(npbasis)
        do i=1,nni*mxnmu
           cbf(i,it)=0.0
        end do
        do igp=1,npbasis
           if(ixref(igp).eq.it) then
              ! m value
              cm(it)=mprim(igp)
              do i=1,nni*mxnmu
                 cbf(i,it)=cbf(i,it) + coeff(igp)*bf(i,igp)
              end do
           end if
        end do
     end do
     do it=1,ixref(npbasis)
        do jt=1,it
           ! Skip different m blocks since those are ignored in code
           xnorm=0.0
           if(cm(it) .eq. cm(jt)) then
              ! Calculate norm
              do i=1,nni*mxnmu
                 xnorm=xnorm+cbf(i,it)*cbf(i,jt)*wgt2(i)*f4(i)
              end do
           end if
           if(abs(xnorm) <= 1e-10) xnorm=0.0
           cS(it,jt)=xnorm
           cS(jt,it)=xnorm
           write (*,'(1X,E12.6)',advance='no') cS(it,jt)
        end do
        write (*,*) ''
     end do
     deallocate(cbf)
     deallocate(cm)
     deallocate(cS)
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
!     write (*,*) 'Orbital ',iorb,' overlap'
     call norm94 (iorb,psi,f4,wgt2,wk0,xnorm)
     write (*,1115) iorn(iorb),bond(iorb),gut(iorb),xnorm
  end do

  deallocate(bf)

  !  if (iprt.eq.1) write (*,1115) iorn(iorb),bond(iorb),gut(iorb),xnorm

  write(*,*)

  if (idbg(560).ne.0) stop 'inigauss'

  !     initialize Coulomb and exchange potentials

  call initPot(psi,pot,excp,f4,wk0)

1114 format(/1x,'    orbital        norm      ')
1115 format(1x,i3,1x,a8,a1,e20.12)

end subroutine initGauss
      
