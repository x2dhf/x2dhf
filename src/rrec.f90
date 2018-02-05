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
! ### rrec ###

! FIXME

subroutine rrec(irec,ngrid,axyz)
  use params
  use commons8

  implicit none
  integer :: irec,ngrid,iunit
  real (PREC), dimension(ngrid) :: axyz

  iunit=irec+30
  if (iunit.gt.maxunit) then
     write(*,910)
     stop "rrec"
910  format('rrec: maximum unit number exceeded (see User''s Guide)')
  endif
  open(iunit,status='old',form='unformatted')
  rewind(iunit)
  read  (iunit) axyz
  close(iunit)

end subroutine rrec
