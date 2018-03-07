! ***************************************************************************
! *                                                                         *
! *   Copyright (C) 2010 Jacek Kobus <jkob@fizyka.umk.pl>                   *
! *                                                                         *
! *   This program is free software; you can redistribute it and/or modify  *
! *   it under the terms of the GNU General Public License version 2 as     *
! *   published by the Free Software Foundation.                            *
! *                                                                         *
! ***************************************************************************
! ### nfng ###

!     Calculates nabla f nabla g; the result is returned in wk3

subroutine nfng(f,g,fmu,fni,gmu,gni,wk0,wk1,wk2,wk3)
  use params
  use discret
  use commons8

  use diffmu_m
  use diffnu_m

  implicit none
  integer :: i,ii,isym,j
  real (PREC) :: t1,z
  real (PREC), dimension(*) ::  f,g,fmu,fni,gmu,gni,wk0,wk1,wk2,wk3

  !   it is assumed that functions have isym=1 symmetry
  isym=1

  !   calculate df/dmu and df/dni derivatives

  call putin   (nni,mxnmu,isym,f,wk3)
  call diff1nu (mxnmu,wk3,wk0,wk1,wk2)
  call putout (nni,mxnmu,fni,wk0)

  !    print *,'f: ', (f(i),i=1000,1003)
  !    print *,'fni: ',(fni(i),i=1000,1003)

  call diff1mu (mxnmu,wk3,wk2)
  call putout (nni,mxnmu,fmu,wk2)

  !    print *,'fmu: ',(fmu(i),i=1000,1003)

  !   calculate dg/dmu and dg/dni derivatives
  call putin   (nni,mxnmu,isym,g,wk3)
  call diff1nu (mxnmu,wk3,wk0,wk1,wk2)
  call putout (nni,mxnmu,gni,wk0)
  !    print *,'g: ', (f(i),i=1000,1003)
  !    print *,'gni: ',(gni(i),i=1000,1003)

  call diff1mu (mxnmu,wk3,wk2)
  call putout (nni,mxnmu,gmu,wk2)
  !    print *,'gmu: ',(gmu(i),i=1000,1003)

  !   calculate df/dmu dg/dmu term (wk0)
  call prod2(mxsize,fmu,gmu,wk0)

  do i=1,mxnmu
     ii=(i-1)*nni
     if (vxi1(i).lt.precis) then
        do j=1,nni
           wk0(ii+j)=0.0_PREC
        enddo
     else
        do j=1,nni
           z=(half*r)*vxi(i)*veta(j)
           wk0(ii+j)=wk0(ii+j)*((half*r*vxi1(i)*veta1(j)*vxi(i))**2 + (vxi(i)*z-half*r*veta(j))**2)/vxi1(i)**2
        enddo
     endif
  enddo

  !   calculate df/dni dg/dni term (wk1)
  call prod2(mxsize,fni,gni,wk1)

  do i=1,mxnmu
     ii=(i-1)*nni
     do j=1,nni
        if (veta1(j).lt.precis) then
           wk1(ii+j)=0.0_PREC
        else
           z=(half*r)*vxi(i)*veta(j)
           wk1(ii+j)=wk1(ii+j)*((half*r*vxi1(i)*veta1(j)*veta(j))**2 + (half*r*vxi(i)-veta(j)*z)**2)/veta1(j)**2
        endif
     enddo
  enddo


  !   calculate df/dmu dg/dni+ df/dni dg/dmu term (wk2)

  call prod2(mxsize,fmu,gni,wk2)
  call prod2(mxsize,fni,gmu,wk3)
  call add  (mxsize,wk2,wk3)

  do i=1,mxnmu
     ii=(i-1)*nni
     do j=1,nni
        if (veta1(j).lt.precis.or.vxi1(i).lt.precis) then
           wk3(ii+j)=0.0_PREC
        else
           z=(half*r)*vxi(i)*veta(j)
           wk3(ii+j)=wk3(ii+j)*( (half*r*vxi1(i)*veta1(j))**2*vxi(i)*veta(j)-&
                (half*r*vxi(i)-veta(j)*z)*(vxi(i)*z-half*r*veta(j)))/(vxi1(i)*veta1(j))
        endif
     enddo
  enddo

  call add  (mxsize,wk0,wk3)
  call add  (mxsize,wk1,wk3)

  do i=1,mxnmu
     ii=(i-1)*nni
     do j=1,nni
        t1=r*r*(vxi(i)*vxi(i)-veta(j)*veta(j))
        if (t1.lt.precis) then
           wk3(ii+j)=0.0_PREC
        else
           wk3(ii+j)=wk3(ii+j)*(four/t1)**2
        endif
        !          write (*,'(2i4,e15.4)') i,j, wk3(ii+j)
     enddo
  enddo

end subroutine nfng
