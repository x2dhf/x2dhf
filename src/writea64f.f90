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
! ### writea64f ###
!
!     Writes a to a disk file in a formatted form
!
module writea64f_m
  implicit none
contains
  subroutine writea64f (iunit,ndim,a,ierr)
    use params
    use commons8

    implicit none

    integer*8 :: iunit, ndim, ierr
    real (PREC), dimension(ndim) :: a

    write (iunit,formfp64,err=1000) a
    ierr=0
    return

1000 ierr=1
  end subroutine writea64f
end module writea64f_m
