! ***************************************************************************
! *                                                                         *
! *   Copyright (C) 2010 Jacek Kobus <jkob@fizyka.umk.pl>                   *
! *                                                                         *
! *   This program is free software; you can redistribute it and/or modify  *
! *   it under the terms of the GNU General Public License version 2 as     *
! *   published by the Free Software Foundation.                            *
! *                                                                         *
! ***************************************************************************
! ### n2f ###

!     Calculates nabla^2 f and returns it as wk3

subroutine n2f(f,wk0,wk1,wk2,wk3)
  use params
  use discret
  use solver
  use commons8

  use diffmu_m
  use diffnu_m

  implicit none
  integer :: i,ii,j
  real (PREC) :: t1
  real (PREC), dimension(*) :: f,wk0,wk1,wk2,wk3

  !   it is assumed that fuction have isym=1 symmetry
  isym=1

  call putin  (nni,mxnmu,isym,f,wk1)
  call diffnu (mxnmu,wk1,wk0,wk3,wk2)
  call putout (nni,mxnmu,wk3,wk0)
  call diffmu (mxnmu,wk1,wk2)
  call putout (nni,mxnmu,wk0,wk2)

  !   add derivatives over mu and ni
  call add (mxsize,wk0,wk3)

  !   take care of 4/[R^2(xi^2-eta^2)] factor in Laplasian
  do i=1,mxnmu
     ii=(i-1)*nni
     do j=1,nni
        t1=r*r*(vxi(i)*vxi(i)-veta(j)*veta(j))
        if (t1.lt.precis) then
           wk3(ii+j)=0.0_PREC
        else
           wk3(ii+j)=wk3(ii+j)*(four/t1)
        endif
        !         write (*,'(2i4,e15.4)') i,j, wk3(ii+j)
     enddo
  enddo

end subroutine n2f

