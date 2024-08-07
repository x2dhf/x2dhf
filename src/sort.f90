! SPDX-License-Identifier: GPL-2.0-or-later

! Copyright (C) 2018       Susi Lehtola

module sort
  use params, only : IPREC, PREC
  implicit none
contains
  subroutine isort(x,y,map)
    integer (KIND=IPREC),intent(in) :: x(:)
    integer (KIND=IPREC),intent(out) :: y(:)
    integer (KIND=IPREC),intent(in) :: map(:)
    integer (KIND=IPREC) :: i

    do i=1,size(y)
       y(i)=x(map(i))
    end do
  end subroutine isort

  subroutine rsort(x,y,map)
    real(PREC), intent(in) :: x(:)
    real(PREC), intent(out) :: y(:)
    integer (KIND=IPREC),intent(in) :: map(:)
    integer (KIND=IPREC) :: i

    do i=1,size(y)
       y(i)=x(map(i))
    end do
  end subroutine rsort
  
end module sort

