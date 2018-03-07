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
! ### coulAsympt ###
!
!     mpole multipole moments (stored in cmulti) and associated Legendre
!     functions are used to calculate for each i=1,nni 4 boundary values
!     of Coulomb potential for a given orbital ($\tilde{V}^a_C$)),
!     i.e. its values for j=mxnmu-3, mxnmu-2, mxnmu-1,mxnmu

module coulAsympt_m
  implicit none
contains
  subroutine coulAsympt(iorb,pot)
    use params
    use discret
    use scf
    use commons8
    use vcoul_m

    implicit none
    integer :: i,iorb,itt,j,kk
    real (PREC) :: costh,rr,rr1
    real (PREC), dimension(*) :: pot

    do j=mxnmu-3,mxnmu
       itt=(j-1)*nni
       do i=1,nni
          kk=i+itt
          rr=sqrt(vxisq(j)+vetasq(i)-1.0_PREC)
          rr1=1.0_PREC/(rr*r2)
          costh=veta(i)*vxi(j)/rr
          pot(kk)=r2*vxi(j)*(rr1+vcoul(iorb,i,j,costh))
       enddo
    enddo
  end subroutine coulAsympt
end module coulAsympt_m

