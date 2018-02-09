! ***************************************************************************
! *                                                                         *
! *   Copyright (C) 2018 Susi Lehtola                                       *
! *                                                                         *
! *   This program is free software; you can redistribute it and/or modify  *
! *   it under the terms of the GNU General Public License version 2 as     *
! *   published by the Free Software Foundation.                            *
! *                                                                         *
! ***************************************************************************
module sort
  use params, only : PREC
  implicit none

contains

  subroutine isort(x,y,map)
    integer, intent(in) :: x(:)
    integer, intent(out) :: y(:)
    integer, intent(in) :: map(:)
    integer :: i
    
    do i=1,size(y)
       y(i)=x(map(i))
    end do
  end subroutine isort

  subroutine rsort(x,y,map)
    real(PREC), intent(in) :: x(:)
    real(PREC), intent(out) :: y(:)
    integer, intent(in) :: map(:)
    integer :: i
    
    do i=1,size(y)
       y(i)=x(map(i))
    end do
  end subroutine rsort
end module sort

