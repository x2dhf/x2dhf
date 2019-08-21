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
! ### slaterp ###

!     Calculates exchange potentials according to the local Slater
!     approximation

module slaterp_m
  implicit none
contains
  subroutine slaterp(psi,excp,f4)
    use params
    use discret
    use commons8
    use util
    use zeroArray_m

    implicit none
    integer :: i,iborb1,iorb1,ngorb1

    real (PREC) :: const13,coo,xa
    real (PREC), dimension(*) :: psi,excp,f4

    parameter (const13=1.0_PREC/3.0_PREC)

    if (nel.eq.1) return

    call zeroArray(mxsize,excp)

    !   contributions from the local slater exchange

    do iorb1=1,norb
       iborb1=i1b(iorb1)
       ngorb1=i1si(iorb1)
       coo=occ(iorb1)
       call prodas(ngorb1,coo,psi(iorb1),psi(iorb1),excp)
    enddo

    !     multiply exchange potential by f4 to make it commensurate with
    !     Coulomb potential

    xa=-three/two*alphaf*(three/pii)**const13

    do i=1,mxsize
       excp(i)=xa*f4(i)*(excp(i))**const13
    enddo

  end subroutine slaterp
end module slaterp_m
