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
! ### wtdexch ###
!
!     Writes exchange potentials involving orbital iorb as separate files

module wtdexch_m
  implicit none
contains
  subroutine wtdexch (iorb1,excp)
    use params
    use commons8
    use wrec_m

    implicit none
    integer :: i3beg,iorb1,iorb2,irec,itemp,k,ngrid
    real (PREC), dimension(*) :: excp

    do iorb2=1,i3nexcp(iorb1)
       k=i3breck(iorb1,iorb2)
       itemp=iwexch(k)
       if (iwexch(k).ne.0) then
          if (iwexch(k).eq.-1) iwexch(k)= 0
          if (iwexch(k).eq.-2) iwexch(k)=-1
          i3beg=i3brec(iorb1,iorb2)
          irec=i3xpair(iorb1,iorb2)
          ngrid=i3si(k)
          call wrec(irec,ngrid,excp(i3beg))
       endif
    enddo

  end subroutine wtdexch
end module wtdexch_m
