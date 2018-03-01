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
! ### putout34 ###
!
!     Reverses the action of putin and putout[2-4]

subroutine putout34 (nnit,nmit,fun,work)
  use params
  use discret
  use solver
  use commons8

  implicit none
  integer :: i,j,nnit,nmit

  real (PREC), dimension(nnit,*) :: fun
  real (PREC), dimension(nnit+8,nmit+8) :: work

  !     refill the interior of work array
  
  do i=1,nmit
     do j=1,nnit
        fun(j,muoffs+i)=work(j+4,i+4)
     enddo
  enddo

end subroutine putout34
