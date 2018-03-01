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
! ### hypg1f1 ###

!     Evaluates and returns a value of the confluent hypergeometric
!     function of the first kind _1F_1 = 1+ax/b+a(a+1)x^2/(2!b(b+1))+...

function hypg1f1(np,lp,x)
  use params

  implicit none
  integer :: i,lp,np
  real (PREC) :: hypg1f1
  real (PREC) :: t,x

  hypg1f1=1.0_PREC
  if (-np.eq.0) return
  t=dble(np)*x/dble(lp)
  hypg1f1=hypg1f1+t
  if (-np.eq.1) return
  do i=2,-np
     t=t*x*dble(np+i-1)/dble(lp+i-1)/dble(i)
     hypg1f1=hypg1f1+t
  enddo
end function hypg1f1
