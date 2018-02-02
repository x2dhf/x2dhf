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
! ### checkPot ###
!
!     Checks contributions of multipole moments to Coulomb and exchange
!     potentials at the practical infinity.
!
!     Determines asymptotic (boundary) values of Coulomb and exchange
!     potentials from the multipole expansion of a given order.

subroutine checkPot (pot,excp)
  use params
  use commons8

  implicit none
  integer :: iax,ibeg,idel,ido,iorb,iorb1,iorb2,ipc

  real (PREC), dimension(*) :: pot,excp

  !     write(*,*) 'asympt: for two different nonsigma orbitals',
  !     & 	  'the case (+-) is now also included'

  write(*,*)
  write(*,*) 'Checking multipole expansion for Coulomb potentials'
  do iorb=norb,1,-1
     ibeg = i1b (iorb)
     call coulAsymptCont(iorb,pot(ibeg))
  enddo

  if (iform.eq.0.or.iform.eq.2) return
  if (islat.eq.1) return

  write(*,*)
  write(*,*) 'Checking multipole expansion for exchange potentials'
  do iorb1=1,norb
     do iorb2=iorb1,norb
        if (iorb1.eq.iorb2.and.mgx(6,iorb1).eq.0 ) goto 10
        if ((iorb1.eq.iorb2).and.(ilc(iorb1*(iorb1+1)/2).lt.1)) goto 10

        !           orbitals in increasing order

        ipc=iorb1+iorb2*(iorb2-1)/2
        iax=i3b(ipc)

        idel=abs(mgx(6,iorb1)-mgx(6,iorb2))
        if (iorb1.eq.iorb2) idel=2*mgx(6,iorb1)

        ido=0
1234    ido=ido+1
        if (ido.eq.2) then
           idel=mgx(6,iorb2)+mgx(6,iorb1)
           ipc=ipc+norb*(norb+1)/2
           iax=iax+i3si(ipc)
        endif

        write(*,1000) iorn(iorb1),bond(iorb1),gut(iorb1),iorn(iorb2),bond(iorb2),gut(iorb2)
        write(*,*) '------ ilc(ipc) ',iorb1,iorb2,'...',idel,'...',ipc,ilc(ipc)
1000    format(/i4,1x,a8,a1,3x,i4,1x,a8,a1,3x)

        call exchAsymptCont (idel,ipc,excp(iax))

        if (ilc(ipc).eq.2.and.ido.eq.1) go to 1234
10      continue
     enddo
  enddo

end subroutine checkPot
