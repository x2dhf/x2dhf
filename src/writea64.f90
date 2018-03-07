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
! ### writea64 ###
!
!     Writes a to a disk file in an unformatted form

module writea64_m
  implicit none
contains
  subroutine writea64 (iunit,ndim,a,ierr)
    use params

    implicit none

    integer*8 :: ierr,iunit,ndim
    real (PREC), dimension(ndim) :: a

    write (iunit,err=1000) a
    ierr=0
    return

1000 ierr=1
    return
  end subroutine writea64
end module writea64_m
