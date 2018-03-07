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
! ### exchAsymptCont ###

!     evaluates multipole moment contributions to the boundary values of
!     exchange for at imu=mxnmu, inu=(nn1-1)/2 (Pi/2)

module exchAsymptCont_m
  implicit none
contains
  subroutine exchAsymptCont (idel,ipc,excp)
    use params
    use discret
    use commons8
    use vexch_m

    implicit none
    integer :: i,idel,ipc,j,m

    real (PREC), dimension(*) ::  excp
    real (PREC), dimension(maxmpole) :: excptmp,pe

    j=mxnmu
    i=(nni-1)/2
    call vexch(i,j,idel,ipc,pe)
    do m=1,mpole
       excptmp(m)= r2*vxi(j)*pe(m)
    enddo
    write(*,1000) excptmp(1),(excptmp(m)-excptmp(m-1),m=2,mpole),excptmp(mpole)
1000 format(4e16.8/5e16.8)

  end subroutine exchAsymptCont
end module exchAsymptCont_m
