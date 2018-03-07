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
! ### gmtran ###

!     Transposes an array A and stores it in ATR.

module gmtran_m
  implicit none
contains
  subroutine gmtran(a,atr,n,m)
    use params

    implicit none
    integer :: i,ij,ir,j,n,m
    real (PREC), dimension(*) :: a,atr

    ir=0
    do i=1,n
       ij=i-n
       do j=1,m
          ij=ij+n
          ir=ir+1
          atr(ir)=a(ij)
       enddo
    enddo

  end subroutine gmtran
end module gmtran_m
