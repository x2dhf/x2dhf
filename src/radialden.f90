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
! ### radialden ###

!     Calculates total radial densities relative to centres A (z-R/2)
!     and B (z+R/2) along the internuclear axis (from A to -\infty or B
!     to +\infty).

subroutine radialden(psi,wk0)
  use params
  use discret
  use commons8

  implicit none
  integer :: i,iborb,iorb,iunit,ngorb
  real (PREC) :: coo
  real (PREC), dimension(*) :: psi,wk0

  ngorb=i1si(norb)
  call zeroArray(ngorb,wk0)

  do iorb=1,norb
     iborb=i1b(iorb)
     coo=occ(iorb)
     do i=1,ngorb
        wk0(i)=wk0(i)+coo*psi(iborb+i-1)*psi(iborb+i-1)
     enddo
  enddo

  iunit=99
  if (iprint(110).eq.1) then

     !        print radial density relative to centre A along the
     !        internuclear axis (-R_{\infty}<=z<=-R/2)

     open(iunit,file='density-A',status='replace',form='formatted')
     call prtdenA(nni,i1mu(1),wk0,iunit)
     close(iunit)
  endif

  if (iprint(111).eq.1) then

     !        print radial density relative to centre B along the
     !        internuclear axis (R/2<=z<=R_{\infty})

     open(iunit,file='density-B',status='replace',form='formatted')
     call prtdenB(ione,i1mu(1),wk0,iunit)
     close(iunit)
  endif

end subroutine radialden
