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
! ## plegendg ###
!     Evaluates and returns a value of the associated Lagendre polynomial
!     P_{lm}(\cos \theta)

module plegendg_m
  implicit none
contains
  function plegendg(l,m,costh)
    use params
    use commons8
    use factor_m
    use hypg2f1_m

    implicit none
    integer :: l,m
    real (PREC) :: plegendg
    real (PREC) :: ct,ct1,costh,fn0

    ct=abs(one-costh*costh)

    if (m.eq.0) then
       plegendg=hypg2f1(-l,l+1,ione,(one-costh)/2.0_PREC)
    elseif (m.ne.0.and.ct.lt.precis) then
       plegendg=zero
    else
       fn0=(-one**dble(m))*factor(l+m)/factor(l-m)/factor(m)/(two**dble(m))
       ct1=fn0*ct**(dble(m)/two)
       plegendg=hypg2f1(m-l,m+l+1,m+1,(one-costh)/two)*ct1
    endif

  end function plegendg
end module plegendg_m
