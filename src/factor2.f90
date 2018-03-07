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
! ### factor2 ###

!     Calculates n!!

module factor2_m
  implicit none
contains
  function factor2(n)
    use params
    implicit none
    integer :: i,n
    real (PREC) :: factor2

    factor2=1.0_PREC
    do i=1,n,2
       factor2=dble(i)*factor2
    enddo

  end function factor2
end module factor2_m
