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
! ### orbAsymptSet ###
!
!     Recalculates asymptotic values of a given orbital at the practical
!     infinity using exponential decay values prepared by calling orbAsymptDet.
!
subroutine orbAsymptSet (nmi,psi,edecay)
  use params
  use discret
  use commons8

  implicit none
  integer :: i,j,jj,itt,kk,nmi
  real (PREC), dimension(*) :: psi
  real (PREC), dimension(nni,4) :: edecay

  jj=0
  do j=nmi-3,nmi
     itt=(j-1)*nni
     jj=jj+1
     do i=1,nni
        kk=i+itt
        psi(kk)=psi(kk-nni)*edecay(i,jj)
     enddo
  enddo

end subroutine orbAsymptSet


