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
! ### hypg2f1 ###

!     Evaluates and returns a value of the generalized hypergeometric
!     function _2F_1 = 1+abx/c+a(a+1)b(b+1)x^2/(2!(c(c+1))+...

module hypg2f1_m
  implicit none
contains
  function hypg2f1(na,nb,nc,x)
    use params

    implicit none
    integer :: i,na,nb,nc
    real (PREC) :: hypg2f1
    real (PREC) a,b,c,t,x

    hypg2f1=1.0_PREC
    if (na.eq.0) return
    a=dble(na)
    b=dble(nb)
    c=dble(nc)
    t=a*b*x/c
    hypg2f1=hypg2f1+t
    if (-na.eq.1) return
    do i=2,-na
       t=t*x*(a+dble(i-1))*(b+dble(i-1))/(c+dble(i-1))/dble(i)
       hypg2f1=hypg2f1+t
    enddo

  end function hypg2f1
end module hypg2f1_m

