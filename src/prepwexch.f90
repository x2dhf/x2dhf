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
! c ### prepwexch ###

! c   Prepears data for writing exchange potentials involving a given
! c   orbital (see wtdexch).

module prepwexch_m
  implicit none
contains
  subroutine prepwexch
    use params
    use commons8

    implicit none
    integer :: iorb1,iorb2,k

    do iorb1=1,norb
       do iorb2=iorb1,norb
          k=i3xk(iorb2,iorb1)
          if (k.ne.0) iwexch(k)=k
       enddo
    enddo

  end subroutine prepwexch
end module prepwexch_m
