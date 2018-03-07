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
! ### pmtx ###

!     Prints two-dimensional array in tabular row-wise form.

module pmtx_m
  implicit none
contains
  subroutine pmtx (n1,n2,a,n1st,n2st,incrn1,incrn2)

    use params

    implicit none
    integer :: i,j,incrn1,incrn2,n1,n2,n1st,n2st

    real (PREC), dimension(n1,n2) :: a

    !      print *,n1,n2,n1st,n2st
    write(6,1000) (j, j=n2st,n2,incrn2)
    do i=n1st,n1,incrn1
       write(6,1010) i, (a(i,j),j=n2st,n2,incrn2)
    enddo

01000 format(4i25)
01010 format(/,1x,i4,4e25.16,/(5x,4e25.16))

  end subroutine pmtx


  subroutine pmtxi (n1,n2,ia,n1st,n2st,incrn1,incrn2)

    use params

    implicit none
    integer :: i,j,incrn1,incrn2,n1,n2,n1st,n2st
    integer, dimension(n1,n2) :: ia

    !      print *,n1,n2,n1st,n2st
    write(6,1000) (j, j=n2st,n2,incrn2)
    do i=n1st,n1,incrn1
       write(6,1010) i, (ia(i,j),j=n2st,n2,incrn2)
    enddo

01000 format(4i25)
01010 format(1x,i4,10i5,/(5x,10i5))

  end subroutine pmtxi
end module pmtx_m
