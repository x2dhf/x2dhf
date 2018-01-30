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
! ### diffnu ###

!    This routine calculates 

!     (\frac{\partial^2}{\partial\ni^2} + 
!           d(\ni,\mu) \frac{\partial}{\partial \ni}) f(\ni,\mu)

!     The function is imersed in the array f(nni+8,nmut+8)
!     To allow for the differentiation over ni variable in the form
!     of matrix times vector the array f has to be transposed first.
!     Results of differentiation are stored directly column-wise in 
!     another matrix which means doing back tranposition before 
!     returning to a calling routine.

!     This routine differentiates a function over ni variable.
!     the function is imersed in the array f(nni+6,nmut+6)

subroutine diffnu (n,f,fd,wtran1,wtran2)
  use params
  use discret
  use commons8

  implicit none
  integer :: j,n,n8,nni8,n9

  real (PREC), dimension(nni+8,n+8) :: f,fd
  real (PREC), dimension(n+8,nni+8) :: wtran1,wtran2

  data n9/9/

  !     To allow for the differentiation over ni variable in the form
  !     of matrix times vector the array f has to be transposed first.
  !     Results of differentiation will be stored directly
  !     column-wise in another matrix which means doing back tranposition
  !     before returning to calling routine.  
  
  nni8=nni+8
  n8=n+8
  call gmtran (f,wtran1,nni8,n8)

  !     FC3, ifort: error while passing 9 to a subroutine 
  
  do j=1,nni
     !         call mxv (wtran1(1,j),n8,dni(1,j),n9,wtran2(1,j+4))
     call gemv (n8,n9,wtran1(1,j),dni(1,j),wtran2(1,j+4))
  enddo
  
  call gmtran(wtran2,fd,n8,nni8)
  
end subroutine diffnu

subroutine diff1nu (n,f,fd,wtran1,wtran2)
  use params
  use discret
  use commons8

  implicit none
  integer :: j,n,n8,nni8,n9

  real (PREC), dimension(nni+8,n+8) :: f,fd
  real (PREC), dimension(n+8,nni+8) :: wtran1,wtran2

  data n9/9/

!     To allow for the differentiation over ni variable in the form
!     of matrix times vector the array f has to be transposed first.
!     Results of differentiation will be stored directly
!     column-wise in another matrix which means doing back tranposition
!     before returning to calling routine.  

  nni8=nni+8
  n8=n+8
  
  call gmtran (f,wtran1,nni8,n8)
  
  !     FC3, ifort: error while passing 9 to a subroutine 

  do j=1,nni
     !         call mxv (wtran1(1,j),n8,d1ni(1,j),n9,wtran2(1,j+4))
     call gemv (n8,n9,wtran1(1,j),d1ni(1,j),wtran2(1,j+4))
     
  enddo
  
  call gmtran(wtran2,fd,n8,nni8)
  
end subroutine diff1nu


