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
! ### mpoleMom ###
!
!     Recalculates multipole moment expansion coefficients every time
!     demax(1), i.e. maximum error in orbital energy, is reduced by
!     facmul. Multipole moments used to calculate asymptotic values of
!     Coulomb and exchange potentials are stored in cmulti and
!     exc(di|qu|oc|he|5-8) arrays, respectively.

!     Coefficients and then asymptotic values are recalculated only for
!     orbitals which undergo relaxation, i.e. those being touched
!     (itouch=1).

subroutine mpoleMom (cw_orb,cw_suppl,cw_sctch)
  use params
  use commons8

  implicit none
  integer :: iorb1,iorb2
  real (PREC) :: time1,time2
  real (PREC), dimension(*) :: cw_orb,cw_suppl,cw_sctch

  if (nel.eq.1) return
  
  call getCpuTime(time1)
  call coulMom (cw_orb,cw_suppl(i4b(9)),cw_suppl(i4b(14)),cw_sctch(i5b( 1)),cw_sctch(i5b( 2)),cw_sctch(i5b( 3)),   &
       cw_sctch(i5b( 4)),cw_sctch(i5b( 5)),cw_sctch(i5b( 6)),cw_sctch(i5b( 7)),cw_sctch(i5b( 8)),cw_sctch(i5b( 9)),&
       cw_sctch(i5b(10)))

  if (nel.gt.1.and.imethod.eq.1) then
     do  iorb1=1,norb
        do  iorb2=iorb1,norb
           call exchMom (iorb1,iorb2,cw_orb,cw_suppl(i4b(9)),cw_suppl(i4b(14)),cw_sctch(i5b( 1)),cw_sctch(i5b( 2)),          &
                cw_sctch(i5b( 3)),cw_sctch(i5b( 4)),cw_sctch(i5b( 5)),cw_sctch(i5b( 6)),cw_sctch(i5b( 7)),cw_sctch(i5b( 8)), &
                cw_sctch(i5b( 9)),cw_sctch(i5b(10)))
        enddo
     enddo
  endif
  call getCpuTime (time2)
  tmomen =tmomen + (time2-time1)
  if (iprtlev.ne.3) then 
     write(iout6,*) '... multipole moment expansion coefficients recalculated ...'
  endif
  
end subroutine mpoleMom


