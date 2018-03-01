! ***************************************************************************
! *                                                                         *
! *   Copyright (C) 2018 Susi Lehtola                                       *
! *                                                                         *
! *   This program is free software; you can redistribute it and/or modify  *
! *   it under the terms of the GNU General Public License version 2 as     *
! *   published by the Free Software Foundation.                            *
! *                                                                         *
! ***************************************************************************
module qsort
contains
  subroutine quicksort(n, data, ind)
    use params, only : PREC
    implicit none
    ! Problem size
    integer, intent(in) :: n
    ! Array values
    real(PREC), intent(in) :: data(*)
    ! Sort index
    integer, intent(out) :: ind(*)
    ! Loop index
    integer :: i

    ! Initialize index
    do i=1,n
       ind(i)=i
    end do

    ! Do the actual sort
    call quicksort_work(1, n, data, ind)
  end subroutine quicksort

  recursive subroutine quicksort_work(start, end, data, ind)
    use params, only : PREC
    use swap
    implicit none
    ! Array values
    real(PREC), intent(in) :: data(*)
    ! Sort index
    integer, intent(out) :: ind(*)
    ! Data value
    real(PREC) :: x
    ! Loop indices
    integer :: i, j
    ! Start and end of interval
    integer :: start, end

    ! Midpoint value is
    x = data(ind( (start+end) / 2 ))
    i = start
    j = end
    do
       do while (data(ind(i)) .lt. x)
          i=i+1
       end do
       do while (x .lt. data(ind(j)))
          j=j-1
       end do
       if (i .ge. j) exit
       call iswap(ind(i),ind(j))
       i=i+1
       j=j-1
    end do
    if (start .lt. i-1) call quicksort_work(start, i-1, data, ind)
    if (j+1 .lt. end)  call quicksort_work(j+1, end, data, ind)
  end subroutine quicksort_work
end module qsort
