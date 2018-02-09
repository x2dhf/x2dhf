! ***************************************************************************
! *                                                                         *
! *   Copyright (C) 2018 Susi Lehtola                                       *
! *                                                                         *
! *   This program is free software; you can redistribute it and/or modify  *
! *   it under the terms of the GNU General Public License version 2 as     *
! *   published by the Free Software Foundation.                            *
! *                                                                         *
! ***************************************************************************
module swap
  use params, only : PREC
  implicit none
contains

  subroutine iswap(x,y)
    integer, intent(inout) :: x, y
    integer :: t

    t=x
    x=y
    y=t
  end subroutine iswap

  subroutine rswap(x,y)
    real(PREC), intent(inout) :: x, y
    real(PREC) :: t

    t=x
    x=y
    y=t
  end subroutine rswap
end module swap
