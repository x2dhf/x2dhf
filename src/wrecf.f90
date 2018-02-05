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
! ### wrecf ###

! FIXME

subroutine wrecf(irec,ngrid,axyz)
  use params
  use commons8

  implicit none
  integer :: irec,iunit,ngrid
  real (PREC), dimension(ngrid) :: axyz

  iunit=irec+30
  if (iunit.gt.maxunit) then
     write(*,910)
     stop "wrecf"
910  format('wrec: maximum unit number exceeded (see User''s Guide)')
  endif

  open(iunit,status='old',form='formatted')
  rewind(iunit)
  write (iunit,formfp,err=900) axyz
  close(iunit)
  return
900 write(iout6,*) 'error detected when writing exchange potential',irec
  stop 'wrecf'
end subroutine wrecf
