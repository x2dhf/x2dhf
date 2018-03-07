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
! ### plaguer ###
!    Evaluates and returns a value of the Laguerre polynomial

module plaguer_m
  implicit none
contains
  function plaguer(n,l,cz,cr)
    use params
    use hypg1f1_m

    implicit none
    integer :: n,l,lp,np
    real (PREC) :: plaguer
    real (PREC) :: cz,cr,x

    np=n-l-1
    lp=2*l+2
    x=cz*cr/dble(n)
    plaguer=cr**dble(l)*exp(-x)*hypg1f1(-np,lp,x+x)

  end function plaguer
end module plaguer_m
