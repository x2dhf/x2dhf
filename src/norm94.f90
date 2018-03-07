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
! ### norm94 ###
!
! FIXME

module norm94_m
  implicit none
contains
  subroutine norm94 (iorb,psi,f4,wgt2,wk0,xnorm)
    use params
    use commons8
    use util
    use blas_m

    implicit none
    integer :: i,ibeg,iorb,ngrid
    real (PREC) :: xnorm
    real (PREC), dimension(*) :: psi,f4,wgt2,wk0

    ibeg = i1b (iorb)
    ngrid= i1si(iorb)
    do i=1,ngrid
       wk0(i)=psi(ibeg-1+i)*psi(ibeg-1+i)
    enddo

    call prod  (ngrid,f4,wk0)

    xnorm=dot (ngrid,wgt2,ione,wk0,ione)
    area(iorb)=xnorm

  end subroutine norm94
end module norm94_m
