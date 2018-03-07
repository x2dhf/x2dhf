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
! ### add ###
!
!     a vector plus a vector: dy(i)=dx(i)+dy(i), i=1..n

module util
  implicit none
contains

  subroutine  add(n,dx,dy)
    use params
    implicit none
    integer :: i,n
    real (PREC), dimension(*) :: dx,dy

    do i = 1,n
       dy(i) = dy(i) + dx(i)
    enddo

  end subroutine add

  ! ### prod ###
  !
  !     a vector times a vector: dy(i)=dx(i)*dy(i), i=1..n

  subroutine prod(n,dx,dy)
    use params
    implicit none
    integer :: i,n
    real (PREC), dimension(*) :: dx,dy

    do i = 1,n
       dy(i) = dy(i) * dx(i)
    enddo

  end subroutine prod

  ! ### prod2 ###

  !     a vector times a vector: dr(i)=dx(i)*dy(i), i=1..n

  subroutine   prod2(n,dx,dy,dr)
    use params
    implicit none
    integer :: i,n
    real (PREC), dimension(*) :: dx,dy,dr

    do i = 1,n
       dr(i) = dx(i) * dy(i)
    enddo

  end subroutine prod2

  ! ### proda ###

  !     a vector times a vector: dr(i)=dx(i)*dy(i)+r(i), i=1..n

  subroutine   proda(n,dx,dy,dr)
    use params
    implicit none
    integer :: i,n
    real (PREC), dimension(*) :: dx,dy,dr

    do i = 1,n
       dr(i) = dx(i) * dy(i)+dr(i)
    enddo

  end subroutine proda

  ! ### prodas ###

  !     a vector times a vector: dr(i)=s*dx(i)*dy(i)+dr(i), i=1..n

  subroutine   prodas(n,s,dx,dy,dr)
    use params
    implicit none
    integer :: i,n
    real (PREC) :: s
    real (PREC), dimension(*) :: dx,dy,dr

    do i = 1,n
       dr(i) = s*dx(i) * dy(i)+dr(i)
    enddo

  end subroutine prodas
end module util
