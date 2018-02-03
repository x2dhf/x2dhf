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
! ### setCi ###

!     Determines Ci symmetry of a given orbital (or set of orbitals) and
!     uses it to replace the values in the [pi/2,pi] region by these
!     from the [0,pi/2] one.

subroutine setCi (iorb,psi)
  use params
  use discret
  use commons8

  implicit none
  integer :: i,ihsym,ibeg,iorb,iorb1,iorb2,n2,nmut
  real (PREC), dimension(*) :: psi

!   ihsym = 1 - symmetry g
!   ihsym =-1 - symmetry u

  if (iorb.eq.0) then
     iorb1=1
     iorb2=norb
  else
     iorb1=iorb
     iorb2=iorb
  endif

  n2=nni/2+1

  do i=iorb1,iorb2
     ibeg = i1b(i)
     nmut = i1mu(i)
     ihsym=ihomo(i)
     call setCiOrb (n2,nmut,psi(ibeg),ihsym)
     if (iprint(55).ne.0) then
        write(*,*) 'homo: ihomo(iorb),ihsym',iorn(i),bond(i),gut(i),ihomo(i),ihsym
     endif
  enddo
end subroutine setCi




