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
! ### reada8 ###
!
!     Reads the matrix A from a disk file in an unformatted form

subroutine reada8 (iunit,ndim,a,ierr)
  use params

  implicit none
  integer :: ierr,iunit,ndim
  real (PREC), dimension(ndim) :: a

  read (iunit,err=1000) a
  ierr=0
  return
1000 ierr=1
end subroutine reada8
