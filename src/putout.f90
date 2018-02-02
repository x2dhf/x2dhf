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
! ### putout ###
!
!     Reverses the action of putin and putin[2-4]

subroutine putout (nni,nmi,fun,work)
  use params

  implicit none
  integer :: i,j,nni,nmi

  real (PREC), dimension(nni,nmi) :: fun
  real (PREC), dimension(nni+8,nmi+8) :: work

!     refill the interior of work array

  do i=1,nmi
     do j=1,nni
        fun(j,i)=work(j+4,i+4)
     enddo
  enddo

end subroutine putout
