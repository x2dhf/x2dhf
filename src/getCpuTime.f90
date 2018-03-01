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
! ### getCpuTime ###
!
!     Returns cpu time in seconds

subroutine getCpuTime(t)
  use params

  implicit none
  real (PREC) t

  !      real*4 time(2),result
  !      real*4 a,b
  !      real*4 second

  !     usage:
  !     etime - Convex, Sun, Intel (Linux)
  !     second - Cray
  !     mclock - RS6000 (current process time in 1/100 sec.
  !     t=dble(second())
  !     t=dble(mclock())/100.0_PREC

  !      call etime(time,result)
  !      t=(time(1))

  !     cpu_time works for g77, gfortran, ifort
  call cpu_time (t)

end subroutine getCpuTime


