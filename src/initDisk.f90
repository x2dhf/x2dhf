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
! ### initDisk ###

!     Controls retrieving orbitals and potentials from a disk file

subroutine initDisk (cw_orb,cw_coul,cw_exch,cw_sctch)
  use params
  use commons8

  implicit none
  integer :: lengthint_curr,lengthfp_curr,norb_p

  real (PREC), dimension(*) :: cw_orb,cw_coul,cw_exch,cw_sctch

  !     input files can be in i32, i64 and r128 formats
  !     save the current formats
  lengthint_curr=lengthint
  lengthfp_curr=lengthfp

  !     and determine ones used in input files
  lengthint=lengthintin
  lengthfp=lengthfpin

  call rheader(norb_p)
  if (iinterp.eq.0) then
     call rfun (norb_p,cw_orb,cw_coul,cw_exch,cw_sctch(i5b( 1)),cw_sctch(i5b( 2)))
  else
     call grid_old
     call rfun_int (norb_p,cw_orb,cw_coul,cw_exch,cw_sctch(i5b( 1)),cw_sctch(i5b( 2)),cw_sctch(i5b(3)))
  endif

  lengthint=lengthint_curr
  lengthfp=lengthfp_curr

end subroutine initDisk

