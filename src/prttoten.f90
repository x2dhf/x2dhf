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
! ### prttoten ####

!     Prints total energy and its components

subroutine prttoten
  use params
  use commons8

  implicit none

  if (iprint16.eq.0) then
     write(iout6,6110)  etot
     write(iout6,6100)  evt

     !        if (islat.eq.0) write(iout6,6120)  virrat
     write(iout6,6120)  virrat
     if (ifermi.eq.2.and.iprint(590).ne.0) then
        write(iout6,6130)  etotFN
     endif
  else
     write(iout6,7110)  etot
     write(iout6,7100)  evt

     !        if (islat.eq.0) write(iout6,7120)  virrat
     write(iout6,7120)  virrat
     if (ifermi.eq.2.and.iprint(590).ne.0) then
        write(iout6,7130)  etotFN
     endif
  endif
  return

6100 format(1x,'total electronic energy: ',e28.16)
6110 format(1x,'total energy:            ',e28.16)
6120 format(1x,'virial ratio:            ',e28.16)
6130 format(1x,'total energy + FN:       ',e28.16)

7100 format(1x,'total electronic energy: ',e44.32)
7110 format(1x,'total energy:            ',e44.32)
7120 format(1x,'virial ratio:            ',e44.32)
7130 format(1x,'total energy + FN:       ',e44.32)

end subroutine prttoten

