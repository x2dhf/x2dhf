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
! ### diffmu ###

!     This routine calculates 

!     (\frac{\partial^2}{\partial\mu^2} + 
!           b(\ni,\mu) \frac{\partial}{\partial \mu}) f(\ni,\mu)

!     Function f has been imersed in the array f(nni+8,nmu+8) in order
!     to calculated derivatives in all the grid points.  Originally the
!     routine was used for a single grid of constatnt step size
!     hmu. Accordingly dmu array contained the first- and second-order
!     derivative coefficients (taken from the 8th-order Sterling
!     interpolation formula) multiplied by the b array.

!     To make the routine work in the multigrid case (ngrids.ne.1) 
!     the values of dmu(k,imu) for 

!     imu=iemu(1)-3 ... iemu(1)+3
!     imu=iemu(2)-3 ... iemu(2)+3
!      .
!     imu=iemu(ngrids-1)-3 ... iemu(ngrids-1)+3

!     must be prepared with derivative coefficients which are based on
!     other interpolation formulae taking into account different grid
!     density to the left and right of the grid boundaries. See prepfix
!     for detailes.

!     This routine calculates 

!     {\partial^2 / \partial\mu^2 + 
!           b(\ni,\mu) \partial / \partial \mu} f(\ni,\mu)

!     Function f has been imersed in the array f(nni+8,nmu+8) in order
!     to calculated derivatives in all the grid points.
!     Originally the routine was used for a single grid of constatnt step
!     size hmu. Accordingly dmu array contained the first- and second-order
!     derivative coefficients (taken from the 8th-order Sterling
!     interpolation formula) multiplied by the B array. 

!     To make the routine work in the multigrid case (ngrids.ne.1!) 
!     the values of dmu(k,imu) for 

!     imu=iemu(1)-3 ... iemu(1)+3
!     imu=iemu(2)-3 ... iemu(2)+3
!     .
!     imu=iemu(ngrids-1)-3 ... iemu(ngrids-1)+3

!     must be prepared with derivative coefficients which are based on
!     other interpolation formula taking into account different grid
!     density to the left and right of the grid boundaries. See prepfix
!     for detailes.

subroutine diffmu (n,f,fd)
  use params
  use discret
  use commons8

  implicit none
  integer :: j,n,n9,nni8

  real (PREC), dimension(nni+8,*) :: f,fd

  data n9/9/


  !     The following loop runs now till j=mxnmu so that the derivatives are
  !     also defined in the tail region, i.e. for (i,mxnmu-3),...,(i,mxnmu)
  !     To this end the fill routine provides extra values for
  !     (i,mxnmu+1),...,(i,mxnmu+4) points.
  
  
  nni8=nni+8
  do j=1,n
     !         call mxv(f(1,j),nni8,dmu(1,j),n9,fd(1,j+4))
     call gemv (nni8,n9,f(1,j),dmu(1,j),fd(1,j+4))
  enddo

end subroutine diffmu



subroutine diff1mu (n,f,fd)
  use params
  use discret
  use commons8

  implicit none
  integer :: j,n,n9,nni8

  real (PREC), dimension(nni+8,*) :: f,fd

  data n9/9/

  !     The following loop runs now till j=mxnmu so that the derivatives are
  !     also defined in the tail region, i.e. for (i,mxnmu-3),...,(i,mxnmu)
  !     To this end the fill routine provides extra values for
  !     (i,mxnmu+1),...,(i,mxnmu+4) points.
  
  nni8=nni+8
  do j=1,n
     !         call mxv(f(1,j),nni8,d1mu(1,j),n9,fd(1,j+4))
     call gemv (nni8,n9,f(1,j),d1mu(1,j),fd(1,j+4))
  enddo

end subroutine diff1mu
