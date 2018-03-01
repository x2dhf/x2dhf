! ***************************************************************************
! *                                                                         *
! *   Copyright (C) 1996 Leif Laaksonen, Dage Sundholm                      *
! *   Copyright (C) 1996,2010 Jacek Kobus <jkob@fizyka.umk.pl>              *
! *                                                                         *     
! *   This program is free software; you can redistribute it and/or modify  *
! *   it under the terms of the GNU General Public License version 2 as     *
! *   published by the Free Software Foundation.                            *
! *                                                                         *
! ***************************************************************************
! ### writea32 ###	
!
!     Writes a to a disk file in an unformatted form

subroutine writea32 (iunit,ndim,a,ierr)
 use params

  implicit none

  integer :: ierr,iunit,ndim
  real (PREC), dimension(ndim) :: a
  
  write (iunit,err=1000) a
  ierr=0
  return
  
1000 ierr=1

end subroutine writea32
