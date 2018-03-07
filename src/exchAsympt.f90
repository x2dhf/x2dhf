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
! ### exchAsympt ###

!     Evaluates asymptotic values of a given exchange potential from the
!     multipole expansion.

!     excptail array is used to provide asymptotic values for 'immersed'
!     exchange potentials (see fill).

!     For odd values of q in P(k,q) (here q=idel) there should be no
!     (-1)**q factor in eq. (19). That is why this factor is missing in
!     this routine.

module exchAsympt_m
  implicit none
contains
  subroutine exchAsympt (idel,ipc,excp)
    use params
    use discret
    use commons8
    use vexch_m

    implicit none
    integer :: i,idel,ipc,itt,j,kk

    real (PREC), dimension(*) :: excp
    real (PREC), dimension(maxmpole) :: pe

    do j=mxnmu-3,mxnmu
       itt=(j-1)*nni
       do i=1,nni
          kk=i+itt
          call vexch(i,j,idel,ipc,pe)
          excp(kk)= r2*vxi(j)*pe(mpole)
       enddo
    enddo
  end subroutine exchAsympt
end module exchAsympt_m
