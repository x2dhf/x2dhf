! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 1996-2023  Jacek Kobus 
! Copyright (C) 2018       Susi Lehtola

module interpolate
  use params, only : IPREC, PREC
  integer (KIND=IPREC) :: iord,iord2,kbeg,kend
contains
  ! ### dointerp ###
  !
  !     Performs interpolations of orbitals and potentials between grids.
  !     Orbitals and potentials are treated separately.
  !
  subroutine dointerp (ic,nmuall_p,nmuall,fbefore,fafter)
    use params
    use discrete
    use commons

    implicit none
    integer (KIND=IPREC) :: i,ic,im,imu,imu_bext,in,nmuall,nmuall_p

    real (PREC), dimension(nni_p,nmuall_p)    ::  fbefore
    real (PREC), dimension(nni,nmuall)        ::  fafter
    real (PREC), dimension(:,:), allocatable  ::  fmiddle
    real (PREC) :: gtol

    real (PREC16), dimension(maxmu)           ::  vmuq

    logical :: muchange, nuchange, gridchange, rchange, rinfchange
    logical :: muinterp, nuinterp, usemiddle

    ! Relative tolerance for change in bond length or rinf
    gtol = 100.0_PREC*precis

    ! Figure out what has changed
    muchange = nmuall .ne. nmuall_p
    nuchange = nni .ne. nni_p
    gridchange = ngrids .ne. ngrids_p
    rchange = abs(r-r_p).gt.(gtol*r_p)
    rinfchange = abs(rinf-rinf_p).gt.(gtol*rinf_p)

    ! Do we need to interpolate in mu?
    muinterp=(muchange .or. gridchange .or. rchange .or. rinfchange)

    ! For nu, we only do interpolation if the size of the nu grid changes.
    nuinterp=(nuchange)

    ! Initialize memory
    do in=1,nni
       do im=1,nmuall
          fafter(in,im)=0.0_PREC
       end do
    end do

    ! Allocate helper if necessary
    if(muinterp .and. nuinterp) then
       usemiddle = .true.
       allocate(fmiddle(nni_p,nmuall))
       ! Initialize memory
       do in=1,nni_p
          do im=1,nmuall
             fmiddle(in,im)=0.0_PREC
          end do
       end do
    else
       usemiddle = .false.
    end if

    ! The previous grid contains nni_p points in \nu (\eta) direction and
    ! nmuall_pt points in \mu (\xi) direction
    if (muinterp) then
       ! Prepare data for routines evaluating coefficients of the
       ! Lagrange polynomial

       if (ic.eq.1) then
          iord=iord_mu_orb
       elseif (ic.eq.2) then
          iord=iord_mu_coul
       elseif (ic.eq.3) then
          iord=iord_mu_exch
       else
          stop "invalid argument ic"
       end if

       ! Interpolate in mu variable
       kbeg=1
       kend=iord
       iord2=iord/2
       do imu=1,nmuall_p
          vmuq(imu)=vmu_p(imu)
       end do

       if (usemiddle) then
          call dointerp_mu (nni_p,nmuall_p,nmuall,fbefore,fmiddle,vmuq)
       else
          call dointerp_mu (nni_p,nmuall_p,nmuall,fbefore,fafter,vmuq)
       end if

       ! Extrapolation (needed when R_infy is being increased seems to
       ! cause degradation of accuracy. That is why these values are
       ! taken equal (for each nu) to those corresponding to nmuall_p,
       ! i.e to the values of the last column of fbefore

       if(vmu(nmuall).gt.vmu_p(nmuall_p)) then
          imu_bext=nmuall
          do imu=1,nmuall
             if (vmu(imu).ge.vmu_p(nmuall_p)) then
                imu_bext=imu
                exit
             end if
          end do

          if(usemiddle) then
             do in=1,nni_p
                do im=imu_bext,nmuall
                   fmiddle(in,im)=fbefore(in,nmuall_p)
                end do
             end do
          else
             do in=1,nni_p
                do im=imu_bext,nmuall
                   fafter(in,im)=fbefore(in,nmuall_p)
                end do
             end do
          end if
       end if
    end if

    if(nuinterp) then
       !        interpolate in nu variable

       if (nni_p.gt.maxmu) then
          print *,"dointerp: error! array vmuq is too short"
          stop 'dointerp'
       end if

       if (ic.eq.1) then
          iord=iord_nu_orb
       elseif (ic.eq.2) then
          iord=iord_nu_coul
       elseif (ic.eq.3) then
          iord=iord_nu_exch
       else
          stop "invalid argument ic"
       end if

       ! Interpolate in nu variable
       kbeg=1
       kend=iord
       iord2=iord/2
       do i=1,nni_p
          vmuq(i)=vni_p(i)
       end do

       if (usemiddle) then
          call dointerp_nu (nmuall_p,fmiddle,fafter,vmuq)
       else
          call dointerp_nu (nmuall_p,fbefore,fafter,vmuq)
       end if
    end if

    if(usemiddle) deallocate(fmiddle)

    ! If grids are the same, just use the same values
    if(.not. muinterp .and. .not. nuinterp) then
       do in=1,nni
          do im=1,nmuall
             fafter(in,im)=fbefore(in,im)
          end do
       end do
    end if

  end subroutine dointerp

  ! ### dointerp_mu ###
  !
  !     Performs interpolations of functions in mu variable.
  !
  subroutine dointerp_mu (nnit,nmuall_p,nmuall,fbefore,fafter,vmuq)
    use params
    use discrete
    use commons

    implicit none
    integer (KIND=IPREC) :: i,idebug1,idebug2,idebug3,imu,imu_p,k, ini,&
         nmu_first,nmu_last,nnit,nmuall_p,nmuall
    real (PREC) :: rerror

    real (PREC16) xmu
    real (PREC16), dimension(kend) :: coeffq
    real (PREC16), dimension(kend,kend) :: coeffq2
    real (PREC), dimension(nnit,nmuall_p) :: fbefore
    real (PREC), dimension(nnit,nmuall) :: fafter
    real (PREC16), dimension(maxmu) :: vmuq
    rerror=2
    idebug1=0
    idebug2=0
    idebug3=0

    do imu=1,nmuall
       xmu=vmu(imu)
       do imu_p=1,nmuall_p-1
          if(vmu_p(imu_p).ge.xmu) exit
       end do

       ! Handle case at beginning or end of array
       if(imu_p .lt. iord2+1) then
          imu_p=iord2+1
       end if
       if(imu_p .gt. nmuall_p-iord2) then
          imu_p=nmuall_p-iord2
       end if

       do k=1,kend
          call lpcoeffq(imu_p,k,coeffq,vmuq)
          do i=1,kend
             coeffq2(i,k)=coeffq(i)
          enddo
       enddo

       do ini=1,nnit
          fafter(ini,imu)=0.0_PREC
          do k=1,kend
             fafter(ini,imu)=fafter(ini,imu)+fbefore(ini,imu_p-iord2-1+k)*vpoly1q(xmu,coeffq2(1,k))
          enddo
          if (idebug3.eq.1) then
             if (abs(fafter(ini,imu)-fbefore(ini,imu_p-iord2+1)).gt. &
                  abs(fbefore(ini,imu_p-iord2+1))*rerror) then
                write(*,*) 'inside'
                write(*,'(2i5,e15.3,4x,5e15.3,2i5)') ini,imu,fafter(ini,imu), &
                     (fbefore(ini,imu_p-iord2-1+k),k=1,kend),nmu_first,nmu_last
             endif
          endif
       enddo
    enddo
  end subroutine dointerp_mu

  ! ### dointerp_nu ###
  !
  !     Performs interpolations of functions in ni variable.
  !
  subroutine dointerp_nu (nmuall,fbefore,fafter,vmuq)
    use params
    use discrete
    use commons

    implicit none
    integer (KIND=IPREC) :: i,imu,ini,ini_p,k,nmuall

    real (PREC), dimension(nni_p,nmuall) :: fbefore
    real (PREC), dimension(nni,nmuall) :: fafter
    real (PREC16) xni
    real (PREC16), dimension(kend) :: coeffq
    real (PREC16), dimension(kend,kend) :: coeffq2
    real (PREC16), dimension(maxmu) ::  vmuq

    do ini=1,nni
       xni=vni(ini)
       do ini_p=1,nni_p-1
          if(vni_p(ini_p).ge.xni) exit
       enddo

       ! Handle case at beginning or end of array
       if(ini_p .lt. iord2+1) then
          ini_p=iord2+1
       end if
       if(ini_p .gt. nni_p-iord2) then
          ini_p=nni_p-iord2
       end if

       do k=1,kend
          call lpcoeffq(ini_p,k,coeffq,vmuq)
          do i=1,kend
             coeffq2(i,k)=coeffq(i)
          enddo
       enddo

       do imu=1,nmuall
          fafter(ini,imu)=0.0_PREC
          do k=1,kend
             fafter(ini,imu)=fafter(ini,imu)+fbefore(ini_p-iord2-1+k,imu)*vpoly1q(xni,coeffq2(1,k))
          enddo
       enddo
    enddo
  end subroutine dointerp_nu

  ! ### lpcoeff ###
  !
  !     This routine calculates coefficients of the (sub)Lagrange
  !     polynomial for a grid point k
  !     \prod_{i=1,i\ne k}^{9} {(r-r_{i}) \over (r_{k}-r_{i})}
  !
  !     r_{k}= r(istart+k), k=1,...,iord
  !
  subroutine lpcoeff(iord,istart,k,r,coeff)
    use params

    implicit none
    integer (KIND=IPREC) :: i,ic1,ic2,iord,istart,j,k
    real (PREC) :: c1,denom
    real (PREC), dimension(12) :: a,b
    real (PREC), dimension(9) :: coeff
    real (PREC), dimension(*) :: r

    ! calculate denominator product
    denom=1.0_PREC
    do i=1,iord
       if (i.ne.k) denom=denom*(r(istart+k)-r(istart+i))
    enddo

    do i=1,12
       a(i)=0.0_PREC
       b(i)=0.0_PREC
    enddo

    ! calculate nominator product:
    ! a(1)+a(2)*r+a(3)*r^2+...+a(9)*r^8
    ! storing coefficients in a and b
    a(1)=1.0_PREC

    ! multiply polynomial a by (r+c1), c1=-r(istart+i)
    ic2=1
    do i=1,iord
       if (i.ne.k) then
          ic2=ic2+1
          c1=-r(istart+i)
          b(1)=a(1)*c1
          do ic1=2,ic2
             b(ic1)=a(ic1)*c1+a(ic1-1)
          enddo
          b(ic2+1)=a(ic2)
          do j=1,ic2+1
             a(j)=b(j)
          enddo
       endif
    enddo

    do i=1,iord
       coeff(i)=a(i)/denom
    enddo

  end subroutine lpcoeff

  !  ### lpcoeffq ###
  !
  !      This routine calculates coefficients of the (sub)Lagrange
  !      polynomial for a grid point k
  !      \prod_{i=1,i\ne k}^{9} {(x-x_{i}) \over (x_{k}-x_{i})}
  !
  !      x_{k}= vmu(mup-5+k), k=1,...,9
  !
  subroutine lpcoeffq (mup,k,coeffq,vmuq)
    use params
    implicit none
    integer (KIND=IPREC) :: i,ib,ic1,ic2,j,k,mup
    real (PREC16) :: c1,denom
    real (PREC16), dimension(9) :: coeffq
    real (PREC16), dimension(12) :: a,b
    real (PREC16), dimension(maxmu) :: vmuq

    ! calculate denominator product
    denom=1.0_PREC16
    ib=mup-(iord/2+1)
    do i=kbeg,kend
       if (i.ne.k) denom=denom*(vmuq(ib+k)-vmuq(ib+i))
    enddo

    do i=1,12
       a(i)=0.0_PREC16
       b(i)=0.0_PREC16
    enddo

    ! calculate nominator product:
    ! a(1)+a(2)*x+a(3)*x^2+...+a(9)*x^8
    ! storing coefficients in a and b

    ! if (k.ne.kbeg) then
    ! a(1)=-vmuq(ib+kbeg)
    ! else
    ! a(1)=-vmuq(ib+kbeg+1)
    ! endif
    ! a(2)=1.0_PREC
    a(1)=1.0_PREC16

    ! multiply polynomial a by (x+c1), c1=-vmuq(ib+i)
    ic2=1
    do i=kbeg,kend
       if (i.ne.k) then
          ic2=ic2+1
          c1=-vmuq(ib+i)
          b(1)=a(1)*c1
          do ic1=2,ic2
             b(ic1)=a(ic1)*c1+a(ic1-1)
          enddo
          b(ic2+1)=a(ic2)
          do j=1,ic2+1
             a(j)=b(j)
          enddo
       endif
    enddo
    !
    do i=1,iord
       coeffq(i)=a(i)/denom
    enddo

  end subroutine lpcoeffq
  ! ### vpoly1q ###
  !
  !     This function uses the Horner scheme to calculate value of the polynomial
  !     stored in array a at a particular point
  !
  function vpoly1q (x,a)
    use params
    implicit none
    integer (KIND=IPREC) :: i
    real (PREC16) :: vpoly1q
    real (PREC16) :: x
    real (PREC16), dimension(kend) :: a

    vpoly1q=0.0_PREC
    do i=iord,2,-1
       vpoly1q=(vpoly1q+a(i))*x
    enddo
    vpoly1q=vpoly1q+a(1)
    return
  end function vpoly1q

  ! ### vpolyq ###
  !
  !     This function uses the Horner scheme to calculate value of the polynomial
  !     stored in array a at a particular point
  !
  function vpolyq (mup,a,vmuq)
    use params

    implicit none
    integer (KIND=IPREC) :: i,mup
    real (PREC16) :: vpolyq
    real (PREC16) :: x
    real (PREC16), dimension(kend) :: a
    real (PREC16), dimension(maxmu) :: vmuq

    x=vmuq(mup)
    vpolyq=0.0_PREC16
    do i=iord,2,-1
       vpolyq=(vpolyq+a(i))*x
    enddo
    vpolyq=vpolyq+a(1)

  end function vpolyq

  ! ### lpderq ###
  !
  !     This routine calculates coefficients of the first and second
  !     derivative of the polynomial stored in a
  !
  subroutine lpderq (a,a1,a2)

    use params
    use commons

    implicit none

    integer (KIND=IPREC) :: i
    real (PREC16), dimension(9) :: a,a1,a2

    do i=1,iord-1
       a1(i)=a(i+1)*dble(i)
    enddo
    a1(iord)=0.0_PREC16

    do i=1,iord-2
       a2(i)=a1(i+1)*dble(i)
    enddo
    a2(iord-1)=0.0_PREC16
    a2(iord)  =0.0_PREC16

  end subroutine lpderq


  ! ### flp ###

  !     It is assumed that function f(r) is tabulated at grid points r_i,
  !     i=1,...,n which are stored in array r and the corresponding values
  !     f_i are in array f.
  
  !     The function fr0 uses the Lagrange polynomial of a given order to
  !     calculate f(r0). The interpolation employs grid points adjacent to
  !     r0.
  
  !   iord - number of grid point used for interpolation
  !     n  - dimension of arrys r and f
  !     r  - array of abscissas
  !     f  - array of corresponding values of a given function f
  !     r0 - an abscissa for which f(r0) is calculated via 8th-order
  !          Lagrange polynomial

  function flp(iorder,n,r,f,r0)
    use params
    implicit none
    integer (KIND=IPREC) :: i,iorder,istart,k,n,nearest
    real (PREC) :: flp
    real (PREC) :: r0
    real (PREC), dimension(9) :: coeff
    real (PREC), dimension(9,9) :: coeff2
    real (PREC), dimension(n) :: r,f

    ! iorder value for 2th order (3-point) Lagrange polynomial
    ! parameter (iorder=3)

    ! iorder value for 8th order (9-point) Lagrange polynomial
    ! parameter (iorder=9)

    ! find the smallest element of arry r that is greater than
    ! r0. Whenever possible [iorder/2] points to its left and right are
    ! used to construct the interpolation polynomial (the start and the
    ! end of r array are taken care of)
    do nearest=1,n
       if (r(nearest).gt.r0) exit
    enddo

    ! calculate coefficients of iorder Lagrange polynomials (begining with
    ! istart+1) and store them in coeff2
    istart=nearest-(iorder/2+1)
    if (istart.lt.1) istart=0
    if (istart+iorder.gt.n) istart=n-iorder
    do k=1,iorder
       call lpcoeff(iorder,istart,k,r,coeff)
       do i=1,iorder
          coeff2(i,k)=coeff(i)
       enddo
    enddo

    ! evaluate the value of the Lagrange interpolation polynomial at r0
    flp=0.0_PREC
    do k=1,iorder
       flp=flp+f(istart+k)*vlpcoeff(iorder,r0,coeff2(1,k))
    enddo
    return
  end function flp

  ! ### vlpcoeff ###
  !
  !     This function uses the Horner scheme to calculate the value of a
  !     polynomial defined by coefficients stored in array coeff at a
  !     point r0
  !
  !     vlpcoeff(r0)=coeff(1)+coeff(2)*r0+...+coeff(iorder)*r0**(iorder-1)
  !
  function vlpcoeff (iorder,r0,coeff)
    use params
    
    implicit none
    integer (KIND=IPREC) :: i,iorder
    real (PREC) :: vlpcoeff
    real (PREC) :: r0
    real (PREC), dimension(*) ::coeff

    vlpcoeff=0.0_PREC
    do i=iorder,2,-1
       vlpcoeff=(vlpcoeff+coeff(i))*r0
    enddo
    vlpcoeff=vlpcoeff+coeff(1)

  end function vlpcoeff

end module interpolate
