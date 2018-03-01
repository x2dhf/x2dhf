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
! ### checkSym ###

!     Checks Ci symmetry of all orbitals

subroutine checkSym(psi)
  use params
  use commons8

  implicit none
  integer :: iorb,ibeg,ihsym,nmut

  real (PREC), dimension(*) :: psi

  character*8 :: sigma,pi,delta,phi


  data sigma/'sigma'/,pi/'pi'/,delta/'delta'/,phi/'phi'/

  write(*,1000)
  do iorb=1,norb
     ibeg = i1b(iorb)
     nmut = i1mu(iorb)
     call checkOrbSym(nmut,psi(ibeg),ihsym)
     if (ihsym.eq.1) then
        if (bond(iorb).eq.sigma.or.bond(iorb).eq.delta) then
           write(*,1005) iorn(iorb),bond(iorb),gut(iorb)
        else
           write(*,1010) iorn(iorb),bond(iorb),gut(iorb)
        endif
     endif
     if (ihsym.eq.-1) then
        if(bond(iorb).eq.sigma.or.bond(iorb).eq.delta) then
           write(*,1010) iorn(iorb),bond(iorb),gut(iorb)
        else
           write(*,1005) iorn(iorb),bond(iorb),gut(iorb)
        endif
     endif
  enddo
01000 format(/,' checking symmetry of orbitals:'/,'          required    actual ')
01005 format(1x,i3,1x,a8,1x,a1,10x,'g')
01010 format(1x,i3,1x,a8,1x,a1,10x,'u')

end subroutine checkSym



