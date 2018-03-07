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
! ### exocc ###
!
!     Calculates the number of alpha and beta electrons for a given orbital.

module exocc_m
  implicit none
contains
  subroutine exocc (iorb,ocup,ocdown)
    use params
    use commons8

    implicit none
    integer :: i,iend,iorb,ip,ipe
    real (PREC) :: ocdown,ocup

    character*8 ::  alpha,beta
    data alpha,beta /'+','-'/

    ocup  =0.0_PREC
    ocdown=0.0_PREC

    ipe =mgx(6,iorb)
    ip=4*(iorb-1)

    if(ipe.eq.0) then
       iend=2
    else
       iend=4
    endif

    do i=1,iend
       if(spin(ip+i).ne.alpha.and.spin(ip+i).ne.beta) goto 10
       if(spin(ip+i).eq.alpha) ocup  =ocup+1.0_PREC
       if(spin(ip+i).eq.beta ) ocdown=ocdown+1.0_PREC
00010  continue
    enddo

  end subroutine exocc
end module exocc_m
