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
! ### factor ###

!     Calculates the factorial of n.

function factor(n)
  use params
  implicit none
  integer :: i,n
  real (PREC) :: factor
  
  factor=1.0_PREC
  if (n.eq.0) return
  do i=1,n
     factor=dble(i)*factor
  enddo

end function factor
