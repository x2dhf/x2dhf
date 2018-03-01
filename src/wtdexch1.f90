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
! ### wtdexch1 ###
!
!     Writes exchange potentials as separate files

subroutine wtdexch1 (excp)
  use params
  use commons8

  implicit none
  integer :: i3beg,iorb1,iorb2,irec,k,ngrid
  real (PREC), dimension(*) :: excp

  do iorb1=1,norb
     do iorb2=iorb1,norb
        k=i3xk(iorb1,iorb2)
        if (ilc(k).ne.0) then
           i3beg=i3b(k)
           irec=i3xrec1(k)
           ngrid=i3si(k) 
           
           call wrec(irec,ngrid,excp(i3beg))
           
           if (ilc(k).eq.2) then
              i3beg=i3b(k)+ngrid
              irec=i3xrec2(k)
              call wrec(irec,ngrid,excp(i3beg))
           endif
        endif
     enddo
  enddo
  return
  
  write(iout6,*) 'error detected when writing exchange potential',iorb1,iorb2,k,irec
  stop 'wtdexch1'
end subroutine wtdexch1

