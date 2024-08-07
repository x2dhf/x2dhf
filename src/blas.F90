! SPDX-License-Identifier: GPL-2.0-or-later

! Copyright (C) 1996       Leif Laaksonen, Dage Sundholm               
! Copyright (C) 1996-2023  Jacek Kobus 

! ### blas ###
!
!     This is simplified replacement for BLAS routines: dcopy, daxpy, dscal,
!     ddot, mxv. Note that ix and iy are ignored and are set to 1.

module blas
  implicit none
  
contains
#if !defined(BLAS)
  function ddot(n,dx,ix,dy,iy)
    use params
    implicit none
    integer (KIND=IPREC) :: i,ix,iy,n
    real (PREC) :: ddot
    real (PREC), dimension(*) :: dx,dy
    ddot=0.0_PREC
    do i=1,n
       ddot=ddot+dx(i)*dy(i)
    enddo
  end function ddot

  subroutine  dcopy(n,dx,ix,dy,iy)
    use params
    implicit none
    integer (KIND=IPREC) :: i,ix,iy,n
    real (PREC), dimension(*) :: dx,dy

    do i=1,n
        dy((i-1)*iy+1)=dx((i-1)*ix+1)
     enddo
   end subroutine dcopy
   
  subroutine  daxpy(n,da,dx,ix,dy,iy)
    use params
    implicit none
    integer (KIND=IPREC) :: i,ix,iy,n
    real (PREC) :: da
    real (PREC), dimension(*) :: dx,dy
    do i=1,n
       dy(i)=da*dx(i)+dy(i)
    enddo
  end subroutine daxpy

  subroutine  dscal(n,da,dx,ix)
    use params
    implicit none
    integer (KIND=IPREC) :: i,ix,n
    real (PREC) :: da
    real (PREC), dimension(*) :: dx
    do i=1,n
       dx(i)=da*dx(i)
    enddo
  end subroutine dscal

  ! ### dgemv ###
  !
  !     Multiplies the matrix DX(nr,nc) times the vector DV(nc) and stores
  !     the result in the vector DY(nr) (simplified version of DGEMV)
  !
  subroutine dgemv (trans,nr1,nc,alpha,dx,lda,dv,incx,beta,dvr,incy)
    use params
    implicit none
    integer (KIND=IPREC) :: nc, nr1
    real (PREC), dimension(nr1,*) :: dx
    real (PREC), dimension(*) :: dv,dvr
    character :: trans
    real (PREC) :: alpha,beta,s
    integer (KIND=IPREC) :: incx,incy,inc,inr,lda

    do inr=1,nr1
       s=0.0_PREC
       do inc=1,nc
          s=s+dx(inr,inc)*dv(inc)
       enddo
       dvr(inr)=s
    enddo
  end subroutine dgemv
#endif
  
  subroutine  dcopyi(n,dx,ix,dy,iy)
    use params
    implicit none
    integer (KIND=IPREC) :: i,ix,iy,n
    real (PREC), dimension(*) :: dx,dy
    do i=1,n
        dy((i-1)*iy+1)=-dx((i-1)*ix+1)
     enddo
  end subroutine dcopyi
end module blas
