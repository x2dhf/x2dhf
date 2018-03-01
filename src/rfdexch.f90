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
! ### rfdexch ###
!
!     Reads from a disk file exchange potentials involving orbital iorb

subroutine rfdexch (iorb1,excp)
  use params
  use commons8

  implicit none
  integer :: i3beg,iadd,iorb1,iorb2,irec,k,ngrid
  real (PREC), dimension(*) :: excp

  iadd=0
  do iorb2=1,i3nexcp(iorb1)
     k=i3breck(iorb1,iorb2)
     i3beg=i3brec(iorb1,iorb2)
     if (ilc(k).eq.1) then
        i3b(k)=i3beg
     elseif (ilc(k).eq.2) then
        if (iadd.eq.0) then
           i3b(k)=i3beg
           iadd=1
        else
           iadd=0
        endif
     endif

     irec=i3xpair(iorb1,iorb2)
     ngrid=i3si(k)

     call rrec(irec,ngrid,excp(i3beg))
  enddo
  return

!900 write(iout6,*) 'error detected when reading exchange potential',iorb1,iorb2,k,irec
!  stop 'rfdexch'
!1050 format(/1x,'... writing functions to disk ...'//)
!1070 format(//1x,'error! can not write data to disk'//)
end subroutine rfdexch

