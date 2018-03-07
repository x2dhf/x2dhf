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
! ### blas ###
!
!     This is simplified replacement for BLAS routines: dcopy, daxpy, dscal,
!     ddot, mxv. Note that ix and iy are ignored and are set to 1.

module blas_m
  implicit none
contains
  function  dot(n,dx,ix,dy,iy)
    use params
    implicit none
    integer :: ix,iy,n
    real (PREC) :: dot
    real (PREC), dimension(*) :: dx,dy

#ifdef HAVE_BLAS
    DOUBLE PRECISION, EXTERNAL :: ddot

    dot=ddot(n,dx,ix,dy,iy)
#else
    integer :: i

    dot=0.0_PREC
    do i=1,n
       dot=dot+dx(i)*dy(i)
    enddo
#endif

  end function dot


  subroutine  copy(n,dx,ix,dy,iy)
    use params
    implicit none
    integer :: ix,iy,n
    real (PREC), dimension(*) :: dx,dy

#ifdef HAVE_BLAS
    call dcopy(n,dx,ix,dy,iy)
#else
    integer::i

    do i=1,n
       dy(i)=dx(i)
    enddo
#endif

  end subroutine copy

  subroutine  axpy(n,da,dx,ix,dy,iy)
    use params
    implicit none
    integer :: ix,iy,n
    real (PREC) :: da
    real (PREC), dimension(*) :: dx,dy

#ifdef HAVE_BLAS
    call daxpy(n,da,dx,ix,dy,iy)
#else
    integer :: i

    do i=1,n
       dy(i)=da*dx(i)+dy(i)
    enddo
#endif
  end subroutine axpy

  subroutine  scal(n,da,dx,ix)
    use params
    implicit none
    integer :: ix,n
    real (PREC) :: da
    real (PREC), dimension(*) :: dx

#ifdef HAVE_BLAS
    call dscal(n,da,dx,ix)
#else
    integer :: i

    do i=1,n
       dx(i)=da*dx(i)
    enddo
#endif

  end subroutine scal

  !     Multiplies the matrix DX(nr,nc) times the vector DV(nc) and stores
  !     the result in the vector DY(nr) (simplified version of DGEMV)
  !
  subroutine gemv (nr1,nc,dx,dv,dvr)
    use params
    implicit none
    integer :: nc, nr1
    real (PREC), dimension(nr1,*) :: dx
    real (PREC), dimension(*) :: dv,dvr

#ifdef HAVE_BLAS
    character :: trans = 'n'
    integer :: M, N
    double precision :: alpha = 1.0, beta = 0.0
    integer :: lda
    integer :: incx = 1, incy = 1

    ! Matrix size
    M = nr1
    N = nc
    lda = M

    call dgemv(trans, M, N, alpha, dx, lda, dv, incx, beta, dvr, incy)
#else
    integer :: inc, inr
    real (PREC) :: s

    do inr=1,nr1
       s=0.0_PREC
       do inc=1,nc
          s=s+dx(inr,inc)*dv(inc)
       enddo
       dvr(inr)=s
    enddo
#endif
  end subroutine gemv
end module blas_m
