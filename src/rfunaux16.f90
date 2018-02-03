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
! ### rfunaux16 ###

!     Reads functions from a disk file in an unformatted form

subroutine rfunaux16 (s,t)
  use params
  use commons8

  implicit none
  integer :: i,ioffset
  real (PREC16), dimension(*) :: s
  real (PREC16), dimension(*) :: t

  ioffset=0
  do i=1,60
     if (inhyd(i).eq.0) then
        t(i)=s(i-ioffset)
     else
        ioffset=ioffset+1
     endif
  enddo

end subroutine rfunaux16

