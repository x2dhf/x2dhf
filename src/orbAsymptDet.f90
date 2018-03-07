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
! ### orbAsymptDet ###

!   Initializes array edecay which is used by orbAsymptSet to calculate
!    boundary values of a given orbital at practical infinity.

module orbAsymptDet_m
  implicit none
contains
  subroutine orbAsymptDet (nmi,iorb,edecay,fa)
    use params
    use discret
    use commons8

    implicit none
    integer :: i,iorb,j,jj,itt,kk,nmi
    real (PREC) :: raiq,raiq1
    real (PREC), dimension(*) :: fa
    real (PREC), dimension(nni,4) :: edecay

    jj=0
    do j=nmi-3,nmi
       jj=jj+1
       itt=(j-1)*nni
       do i=1,nni
          kk=i+itt
          raiq1=sqrt(vxisq(j-1)+vetasq(i)-1.0_PREC)
          raiq =sqrt(vxisq(j)+vetasq(i)-1.0_PREC)
          !        edecay(i,jj)=(abs(fa(kk-nni)/fa(kk)))**0.250_PREC*exp(sqrt(abs(fa(kk-nni)))*(raiq1-raiq))
          edecay(i,jj)=(abs(fa(kk)/fa(kk-nni)))**0.250_PREC*exp(sqrt(abs(fa(kk-nni)))*(raiq1-raiq))
       enddo
    enddo

  end subroutine orbAsymptDet
end module orbAsymptDet_m
