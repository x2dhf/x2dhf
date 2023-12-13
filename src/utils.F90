! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 1996-2023  Jacek Kobus 

module utils
  use params, only : IPREC, PREC
  implicit none
contains
  ! ### add ###
  !
  !     Adds a vector to a vector: dy(i)=dx(i)+dy(i), i=1..n
  !
  subroutine  add(n,dx,dy)
    use params
    implicit none
    integer (KIND=IPREC) :: i,n
    real (PREC), dimension(*) :: dx,dy

    do i = 1,n
       dy(i) = dy(i) + dx(i)
    enddo

  end subroutine add

  ! ### prod ###
  !
  !     Multiplies a vector by a vector: dy(i)=dx(i)*dy(i), i=1..n
  !
  subroutine prod(n,dx,dy)
    use params
    implicit none
    integer (KIND=IPREC) :: i,n
    real (PREC), dimension(*) :: dx,dy

    do i = 1,n
       dy(i) = dy(i) * dx(i)
    enddo

  end subroutine prod

  ! ### prod2 ###
  !
  !     Multiplies a vector by a vector: dr(i)=dx(i)*dy(i), i=1..n
  !
  subroutine prod2(n,dx,dy,dr)
    use params
    implicit none
    integer (KIND=IPREC) :: i,n
    real (PREC), dimension(*) :: dx,dy,dr

    do i = 1,n
       dr(i) = dx(i) * dy(i)
    enddo

  end subroutine prod2

  ! ### proda ###
  !
  !     Adds a vector to a product of two vectors: dr(i)=dx(i)*dy(i)+r(i), i=1..n
  !
  subroutine   proda(n,dx,dy,dr)
    use params
    implicit none
    integer (KIND=IPREC) :: i,n
    real (PREC), dimension(*) :: dx,dy,dr

    do i = 1,n
       dr(i) = dx(i) * dy(i)+dr(i)
    enddo

  end subroutine proda

  ! ### prodas ###
  !
  !     Adds a vector to a scaled product of two vectors: dr(i)=s*dx(i)*dy(i)+r(i), i=1..n
  !
  subroutine   prodas(n,s,dx,dy,dr)
    use params
    implicit none
    integer (KIND=IPREC) :: i,n
    real (PREC) :: s
    real (PREC), dimension(*) :: dx,dy,dr

    do i = 1,n
       dr(i) = s*dx(i) * dy(i)+dr(i)
    enddo

  end subroutine prodas

  ! SPDX-License-Identifier: GPL-2.0-or-later
  !
  ! Copyright (C) 2018 Susi Lehtola                                    
  !
  subroutine quicksort(n, data, ind)
    use params, only : PREC
    implicit none
    ! Problem size
    integer (KIND=IPREC),intent(in) :: n
    ! Array values
    real(PREC), intent(in) :: data(*)
    ! Sort index
    integer (KIND=IPREC),intent(out) :: ind(*)
    ! Loop index
    integer (KIND=IPREC) :: i

    ! Initialize index
    do i=1,n
       ind(i)=i
    end do

    ! Do the actual sort
    call quicksort_work(1, n, data, ind)
  end subroutine quicksort

  recursive subroutine quicksort_work(start, end, data, ind)
    use params, only : PREC
    implicit none
    ! Array values
    real(PREC), intent(in) :: data(*)
    ! Sort index
    integer (KIND=IPREC),intent(out) :: ind(*)
    ! Data value
    real(PREC) :: x
    ! Loop indices
    integer (KIND=IPREC) :: i, j
    ! Start and end of interval
    integer (KIND=IPREC) :: start, end

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


! ### setBvals ###                                                                                                 
  subroutine setBvals(isym,inu,imu,wk0)
    use params
    use commons
    use discrete

    implicit none
    integer (KIND=IPREC) :: imu,inu,isym
    real (PREC), dimension (nni,mxnmu) :: wk0

    if (isym==1) then
       if (inu<nni/2) then
       !wk0(inu,imu)=exeven(1)*wk0(inu,imu+1)+exeven(2)*wk0(inu,imu+2)+exeven(3)*wk0(inu,imu+3)&                     
       !     +exeven(4)*wk0(inu,imu+4)+exeven(5)*wk0(inu,imu+5)                                                      
          wk0(inu,imu)=exeven(1)*wk0(inu+1,imu)+exeven(2)*wk0(inu+2,imu)+exeven(3)*wk0(inu+3,imu)&
            +exeven(4)*wk0(inu+4,imu)+exeven(5)*wk0(inu+5,imu)
       else
          wk0(inu,imu)=exeven(1)*wk0(inu-1,imu)+exeven(2)*wk0(inu-2,imu)+exeven(3)*wk0(inu-3,imu)&
               +exeven(4)*wk0(inu-4,imu)+exeven(5)*wk0(inu-5,imu)
       endif
    else
       wk0(inu,imu)=zero
    endif
  end subroutine setBvals


! ### setBvalsmu ###                                                                                                 
  subroutine setBvalsmu(wk0)
    use params
    use commons
    use discrete

    implicit none
    integer (KIND=IPREC) :: i
    real (PREC), dimension (mxnmu) :: wk0
    wk0(1)=exeven(1)*wk0(2)+exeven(2)*wk0(3)+exeven(3)*wk0(4)&
            +exeven(4)*wk0(5)+exeven(5)*wk0(6)
  end subroutine setBvalsmu

! ### setBvalsmu ###                                                                                                 
  subroutine setBvalsnu(wk0)
    use params
    use commons
    use discrete

    implicit none
    integer (KIND=IPREC) :: i
    real (PREC), dimension (nni) :: wk0
    wk0(1)=exeven(1)*wk0(2)+exeven(2)*wk0(3)+exeven(3)*wk0(4)&
            +exeven(4)*wk0(5)+exeven(5)*wk0(6)

    wk0(nni)=exeven(1)*wk0(nni-1)+exeven(2)*wk0(nni-2)+exeven(3)*wk0(nni-3)&
            +exeven(4)*wk0(nni-4)+exeven(5)*wk0(nni-5)

  end subroutine setBvalsnu

  subroutine zeroArray (n,array)
    use params
    implicit none
    integer (KIND=IPREC) :: i,n
    real (PREC), dimension(*) ::  array

    do i=1,n
       array(i)=0.0_PREC
    enddo
  end subroutine zeroArray

  subroutine zeroArray8 (n,array)
    use params
    implicit none
    integer (KIND=8) :: i,n
    real (PREC), dimension(*) ::  array
    do i=1,n
       array(i)=0.0_PREC
    enddo
  end subroutine zeroArray8

  ! ### factor ###
  !     Calculates the factorial of n.
  !
  function factor(n)
    use params
    implicit none
    integer (KIND=IPREC) :: i,n
    real (PREC) :: factor

    factor=1.0_PREC
    if (n.eq.0) return
    do i=1,n
       factor=dble(i)*factor
    enddo
  end function factor

  ! ### factor2 ###
  !
  !     Calculates n!!
  !
  function factor2(n)
    use params
    implicit none
    integer (KIND=IPREC) :: i,n
    real (PREC) :: factor2

    factor2=1.0_PREC
    do i=1,n,2
       factor2=dble(i)*factor2
    enddo

  end function factor2

  ! ### hypg1f1 ###
  !
  !     Evaluates and returns a value of the confluent hypergeometric
  !     function of the first kind _1F_1 = 1+ax/b+a(a+1)x^2/(2!b(b+1))+...
  !
  function hypg1f1(np,lp,x)
    use params

    implicit none
    integer (KIND=IPREC) :: i,lp,np
    real (PREC) :: hypg1f1
    real (PREC) :: t,x

    hypg1f1=1.0_PREC
    if (-np.eq.0) return
    t=dble(np)*x/dble(lp)
    hypg1f1=hypg1f1+t
    if (-np.eq.1) return
    do i=2,-np
       t=t*x*dble(np+i-1)/dble(lp+i-1)/dble(i)
       hypg1f1=hypg1f1+t
    enddo
  end function hypg1f1

  ! ### hypg2f1 ###
  !
  !     Evaluates and returns a value of the generalized hypergeometric
  !     function _2F_1 = 1+abx/c+a(a+1)b(b+1)x^2/(2!(c(c+1))+...
  !
  function hypg2f1(na,nb,nc,x)
    use params

    implicit none
    integer (KIND=IPREC) :: i,na,nb,nc
    real (PREC) :: hypg2f1
    real (PREC) a,b,c,t,x

    hypg2f1=1.0_PREC
    if (na.eq.0) return
    a=dble(na)
    b=dble(nb)
    c=dble(nc)
    t=a*b*x/c
    hypg2f1=hypg2f1+t
    if (-na.eq.1) return
    do i=2,-na
       t=t*x*(a+dble(i-1))*(b+dble(i-1))/(c+dble(i-1))/dble(i)
       hypg2f1=hypg2f1+t
    enddo
  end function hypg2f1

  ! ## plegendg ###
  !     Evaluates and returns a value of the associated Lagendre polynomial
  !     P_{lm}(\cos \theta)

  function plegendg(l,m,costh)
    use params
    use commons
    implicit none
    integer (KIND=IPREC) :: l,m
    real (PREC) :: plegendg
    real (PREC) :: ct,ct1,costh,fn0

    ct=abs(one-costh*costh)

    if (m.eq.0) then
       plegendg=hypg2f1(-l,l+1,ione,(one-costh)/2.0_PREC)
    elseif (m.ne.0.and.ct.lt.precis) then
       plegendg=zero
    else
       fn0=(-one**dble(m))*factor(l+m)/factor(l-m)/factor(m)/(two**dble(m))
       ct1=fn0*ct**(dble(m)/two)
       plegendg=hypg2f1(m-l,m+l+1,m+1,(one-costh)/two)*ct1
    endif

  end function plegendg


  
  ! ### pot2pot ### 
  !
  !     Transforms \tilde{V}_C into V_C (F4=R\xi/2)
  !
  subroutine pot2pot (wk0,f4)
    use discrete
    implicit none
    integer (KIND=IPREC) :: inu,imu
    real (PREC),  dimension(nni,mxnmu) :: wk0,f4
    
    do inu=1,nni
       do imu=1,mxnmu
          wk0(inu,imu)=wk0(inu,imu)/f4(inu,imu)
       enddo
    enddo
  end subroutine pot2pot
  
  
  subroutine setBVpot(fun)
    use params
    use discrete
    use solver
    use commons
    
    implicit none
    integer (KIND=IPREC) :: j,ij,isym,nnu1,nnu2,nnu3,nnu4,nnu5
    real (PREC), dimension(*) :: fun
    real (PREC) :: fun1
    
    nnu1=nni
    nnu2=nnu1+nnu1
    nnu3=nnu2+nnu1
    nnu4=nnu3+nnu1
    nnu5=nnu4+nnu1

    ! FIXME: isym=isymOrb(iorb)
    isym=1
    stop "setBVpot"

    if (isym.eq.1) then
       do ij=2,nnu1-1
          fun1=fun(ij)
          fun(ij)=exeven(1)*fun(ij+nnu1)+exeven(2)*fun(ij+nnu2)+   &
               exeven(3)*fun(ij+nnu3)+exeven(4)*fun(ij+nnu4)+exeven(5)*fun(ij+nnu5)
       enddo
       
       do j=1,mxnmu
          ij=(j-1)*nnu1+1
          fun1=fun(ij)
          fun(ij)=exeven(1)*fun(ij+1)+exeven(2)*fun(ij+2)+  &
               exeven(3)*fun(ij+3)+exeven(4)*fun(ij+4)+exeven(5)*fun(ij+5)
          !        print *,'  2 ',ij,fun1,fun(ij),fun1-fun(ij)
          !        write(*,'(" 2",i4,6e15.6)') ij,fun(ij),fun(ij+1),fun(ij+2),fun(ij+3),fun(ij+4),fun(ij+5)
       enddo
       
       
       do j=1,mxnmu
          ij=j*nnu1
          fun1=fun(ij)
          fun(ij)=exeven(1)*fun(ij-1)+exeven(2)*fun(ij-2)+  &
               exeven(3)*fun(ij-3)+exeven(4)*fun(ij-4)+exeven(5)*fun(ij-5)
          !        print *,'  3 ',ij,fun1,fun(ij),fun1-fun(ij)
       enddo
       !     stop 'writeDisk4dd'
       
    else
       do ij=2,nnu1-1
          fun(ij)=0.0_PREC
       enddo
       
       do j=1,mxnmu
          ij=(j-1)*nnu1+1
          fun(ij)=0.0_PREC
       enddo
       
       do j=1,mxnmu
          ij=j*nnu1
          fun(ij)=0.0_PREC
       enddo
    endif
    
  end subroutine setBVpot
  
  subroutine setBVgen(isym4nu,isym4mu,fun)
    use params
    use discrete
    use solver
    use commons
    
    implicit none
    integer (KIND=IPREC) :: i,j,ij,isym4nu,isym4mu,nnu1,nnu2,nnu3,nnu4,nnu5
    real (PREC), dimension(*) :: fun
    real (PREC) :: fun1
    
    nnu1=nni
    nnu2=nnu1+nnu1
    nnu3=nnu2+nnu1
    nnu4=nnu3+nnu1
    nnu5=nnu4+nnu1
    
    if (isym4nu.eq.1) then
       do ij=2,nnu1-1
          fun1=fun(ij)
          fun(ij)=exeven(1)*fun(ij+nnu1)+exeven(2)*fun(ij+nnu2)+   &
               exeven(3)*fun(ij+nnu3)+exeven(4)*fun(ij+nnu4)+exeven(5)*fun(ij+nnu5)
          !         print *,'  1 ',ij,fun1,fun(ij),fun1-fun(ij)
          ! !        print *,'  1 ',fun(ij+nnu1),fun(ij+nnu2),fun(ij+nnu3),fun(ij+nnu4)
       enddo
    else
       do ij=2,nnu1-1
          fun(ij)=0.0_PREC
       enddo
    end if
    
    if (isym4mu.eq.1) then
       do j=1,mxnmu
          ij=(j-1)*nnu1+1
          fun1=fun(ij)
          fun(ij)=exeven(1)*fun(ij+1)+exeven(2)*fun(ij+2)+  &
               exeven(3)*fun(ij+3)+exeven(4)*fun(ij+4)+exeven(5)*fun(ij+5)
          !        print *,'  2 ',ij,fun1,fun(ij),fun1-fun(ij)
          !        write(*,'(" 2",i4,6e15.6)') ij,fun(ij),fun(ij+1),fun(ij+2),fun(ij+3),fun(ij+4),fun(ij+5)
       enddo
       
       do j=1,mxnmu
          ij=j*nnu1
          fun1=fun(ij)
          fun(ij)=exeven(1)*fun(ij-1)+exeven(2)*fun(ij-2)+  &
               exeven(3)*fun(ij-3)+exeven(4)*fun(ij-4)+exeven(5)*fun(ij-5)
          !        print *,'  3 ',ij,fun1,fun(ij),fun1-fun(ij)
       enddo
    else
       do j=1,mxnmu
          ij=(j-1)*nnu1+1
          fun(ij)=0.0_PREC
       enddo
       
       do j=1,mxnmu
          ij=j*nnu1
          fun(ij)=0.0_PREC
       enddo
    endif
    
  end subroutine setBVgen
  
  
  subroutine checkd1nu (m,n,a)
    use params
    use discrete
    
    implicit none
    integer (KIND=IPREC) :: im,in,m,n
    integer (KIND=IPREC) :: immax,inmax
    real (PREC), dimension(m,n) :: a
    real (PREC) :: diff,vn
    
    diff=0.0_PREC
    immax=0
    inmax=0
    
    do im=1,m
       do in=1,n
          vn=-veta1(im)*(z2-z1)*r
          ! write(ioutmat,'("checkd1nu: ",2i5,3e15.6)') im,in,vn,a(im,in),abs(vn-a(im,in))
          if (abs(vn-a(im,in)).gt.diff) then
             diff=abs(vn-a(im,in))
             immax=im
             inmax=in
          endif
       enddo
    enddo
    
    write(*,'("checkd1nu: ",2i5,e15.6)') immax,inmax,diff
    
  end subroutine checkd1nu
  
  
  subroutine checkd1mu (m,n,a)
    use params
    use discrete
    
    implicit none
    integer (KIND=IPREC) :: im,in,ioutmat,m,n
    integer (KIND=IPREC) :: immax,inmax
    real (PREC), dimension(m,n) :: a
    real (PREC) :: diff,vn

    ioutmat=6
    diff=0.0_PREC
    immax=0
    inmax=0
    
    do im=1,m
       do in=1,n
          vn=vxi1(in)*(z1+z2)*r
          write(ioutmat,'("checkd1mu: ",2i5,3e15.6)') im,in,vn,a(im,in),abs(vn-a(im,in))
          if (abs(vn-a(im,in)).gt.diff) then
             diff=abs(vn-a(im,in))
             immax=im
             inmax=in
          endif
       enddo
    enddo
    
    write(ioutmat,'("checkd1m: ",2i5,e15.6)') immax,inmax,diff
    
  end subroutine checkd1mu

  ! ### diffmu ###
  !
  !     This routine calculates
  !
  !      (\frac{\partial^2}{\partial\mu^2} +
  !           b(\ni,\mu) \frac{\partial}{\partial \mu}) f(\ni,\mu)
  !
  !     Function f has been imersed in the array f(nni+8,nmu+8) in order
  !     to calculated derivatives in all the grid points.  Originally the
  !     routine was used for a single grid of constatnt step size
  !     hmu. Accordingly dmu array contained the first- and second-order
  !     derivative coefficients (taken from the 8th-order Sterling
  !     interpolation formula) multiplied by the b array.
  !
  !     To make the routine work in the multigrid case (ngrids.ne.1)
  !     the values of dmu(k,imu) for
  !
  !     imu=iemu(1)-3 ... iemu(1)+3
  !     imu=iemu(2)-3 ... iemu(2)+3
  !      .
  !     imu=iemu(ngrids-1)-3 ... iemu(ngrids-1)+3
  !
  !     must be prepared with derivative coefficients which are based on
  !     other interpolation formulae taking into account different grid
  !     density to the left and right of the grid boundaries. See prepfix
  !     for detailes.
  !
  !     This routine calculates
  !
  !     {\partial^2 / \partial\mu^2 +
  !           b(\ni,\mu) \partial / \partial \mu} f(\ni,\mu)
  !
  !     Function f has been imersed in the array f(nni+8,nmu+8) in order
  !     to calculated derivatives in all the grid points.
  !     Originally the routine was used for a single grid of constatnt step
  !     size hmu. Accordingly dmu array contained the first- and second-order
  !     derivative coefficients (taken from the 8th-order Sterling
  !     interpolation formula) multiplied by the B array.
  !
  !     To make the routine work in the multigrid case (ngrids.ne.1!)
  !     the values of dmu(k,imu) for
  !
  !     imu=iemu(1)-3 ... iemu(1)+3
  !     imu=iemu(2)-3 ... iemu(2)+3
  !     .
  !     imu=iemu(ngrids-1)-3 ... iemu(ngrids-1)+3
  !
  !     must be prepared with derivative coefficients which are based on
  !     other interpolation formula taking into account different grid
  !     density to the left and right of the grid boundaries. See prepfix
  !     for detailes.
  
  subroutine diffmu (n,f,fd)
    use params
    use discrete
    use commons
    use blas

    implicit none
    integer (KIND=IPREC) :: j,n,n9
    real (PREC), dimension(nni+8,*) :: f,fd
    data n9/9/
    character :: trans = 'n'
#ifdef BLAS    
    external dgemv
#endif

    ! The following loop runs now till j=mxnmu so that the derivatives are
    ! also defined in the tail region, i.e. for (i,mxnmu-3),...,(i,mxnmu)
    ! To this end the fill routine provides extra values for
    ! (i,mxnmu+1),...,(i,mxnmu+4) points.

    nni8=nni+8
    do j=1,n
       !call gemv (nni8,n9,f(1,j),dmu(1,j),fd(1,j+4))
       call dgemv(trans, nni8, n9, one, f(1,j), nni8, dmu(1,j), ione, zero, fd(1,j+4), ione)
    enddo
  end subroutine diffmu

  subroutine diff1mu (n,f,fd)
    use params
    use discrete
    use commons
    use blas

    implicit none
    integer (KIND=IPREC) :: j,n,n9

    real (PREC), dimension(nni+8,*) :: f,fd
    data n9/9/
    character :: trans = 'n'
#ifdef BLAS    
    external dgemv
#endif

    ! The following loop runs now till j=mxnmu so that the derivatives are
    ! also defined in the tail region, i.e. for (i,mxnmu-3),...,(i,mxnmu)
    ! To this end the fill routine provides extra values for
    ! (i,mxnmu+1),...,(i,mxnmu+4) points.

    nni8=nni+8
    do j=1,n
       call dgemv(trans, nni8, n9, one, f(1,j), nni8, d1mu(1,j), ione, zero, fd(1,j+4), ione)
    enddo

  end subroutine diff1mu

  ! ### diffnu ###
  !
  !    This routine calculates
  !
  !     (\frac{\partial^2}{\partial\ni^2} +
  !           d(\ni,\mu) \frac{\partial}{\partial \ni}) f(\ni,\mu)
  !
  !     The function is imersed in the array f(nni+8,nmut+8)
  !     To allow for the differentiation over ni variable in the form
  !     of matrix times vector the array f has to be transposed first.
  !     Results of differentiation are stored directly column-wise in
  !     another matrix which means doing back tranposition before
  !     returning to a calling routine.
  !
  !     This routine differentiates a function over ni variable.
  !     the function is imersed in the array f(nni+6,nmut+6)

  subroutine diffnu (n,f,fd,wtran1,wtran2)
    use params
    use discrete
    use commons

    use blas

    implicit none
    integer (KIND=IPREC) :: j,n,n8,n9

    real (PREC), dimension(nni+8,n+8) :: f,fd
    real (PREC), dimension(n+8,nni+8) :: wtran1,wtran2
    data n9/9/
    character :: trans = 'n'
#ifdef BLAS    
    external dgemv
#endif

    ! To allow for the differentiation over ni variable in the form of
    ! matrix times vector the array f has to be transposed first.  Results
    ! of differentiation will be stored directly column-wise in another
    ! matrix which means doing back tranposition before returning to
    ! calling routine.

    nni8=nni+8
    n8=n+8
    call gmtran (nni8,n8,f,wtran1)
    do j=1,nni
       call dgemv(trans, n8, n9, one, wtran1(1,j), n8, dni(1,j), ione, zero, wtran2(1,j+4), ione)
    enddo

    call gmtran(n8,nni8,wtran2,fd)

  end subroutine diffnu

  subroutine diff1nu (n,f,fd,wtran1,wtran2)
    use params
    use discrete
    use commons
    use blas

    implicit none
    integer (KIND=IPREC) :: j,n,n8,n9

    real (PREC), dimension(nni+8,n+8) :: f,fd
    real (PREC), dimension(n+8,nni+8) :: wtran1,wtran2
    data n9/9/
    character :: trans = 'n'
#ifdef BLAS    
    external dgemv
#endif
    
    ! To allow for the differentiation over ni variable in the form
    ! of matrix times vector the array f has to be transposed first.
    ! Results of differentiation will be stored directly
    ! column-wise in another matrix which means doing back tranposition
    ! before returning to calling routine.

    nni8=nni+8
    n8=n+8

    call gmtran (nni8,n8,f,wtran1)
    do j=1,nni
       call dgemv(trans, n8, n9, one, wtran1(1,j), n8,d1ni(1,j), ione, zero, wtran2(1,j+4), ione)
    enddo
    call gmtran(n8,nni8,wtran2,fd)

  end subroutine diff1nu

  ! ### gmtran ###
  !
  !     Transposes an array A and stores it in ATR.
  !
  subroutine gmtran(n,m,a,atr)
    use params

    implicit none
    integer (KIND=IPREC) :: i,ij,ir,j,n,m
    real (PREC), dimension(*) :: a,atr

    ir=0
    do i=1,n
       ij=i-n
       do j=1,m
          ij=ij+n
          ir=ir+1
          atr(ir)=a(ij)
       enddo
    enddo

  end subroutine gmtran

  subroutine iswap(x,y)
  use params, only : IPREC, PREC
  implicit none
    integer (KIND=IPREC),intent(inout) :: x, y
    integer (KIND=IPREC) :: t
    t=x
    x=y
    y=t
  end subroutine iswap

  subroutine rswap(x,y)
    use params, only : IPREC, PREC
    implicit none
    real(PREC), intent(inout) :: x, y
    real(PREC) :: t
    t=x
    x=y
    y=t
  end subroutine rswap

  ! ### multf4 ###
  !
  !     Multiplies a given array (of length mxnmu) by F4
  !
  subroutine multf4(a)
    use params
    use discrete
    !  use commons


    implicit none
    integer (KIND=IPREC) :: i,ii,j
    real (PREC), dimension(*) :: a

    do i=1,mxnmu
       ii=(i-1)*nni
       do j=1,nni
          a(ii+j)=a(ii+j)*half*r*vxi(i)
       enddo
    enddo
  end subroutine multf4

  ! ### pmtx ###
  !
  !     Prints two-dimensional array in tabular row-wise form.
  !
  subroutine pmtx (n1,n2,a,n1st,n2st,incrn1,incrn2)
    use params
    implicit none
    integer (KIND=IPREC) :: i,j,incrn1,incrn2,n1,n2,n1st,n2st

    real (PREC), dimension(n1,n2) :: a

    write(6,1000) (j, j=n2st,n2,incrn2)
    do i=n1st,n1,incrn1
       write(6,1010) i, (a(i,j),j=n2st,n2,incrn2)
    enddo

01000 format(4i25)
01010 format(/,1x,i4,4e25.16,/(5x,4e25.16))

  end subroutine pmtx

  subroutine pmtxi (n1,n2,ia,n1st,n2st,incrn1,incrn2)
    use params
    implicit none
    integer (KIND=IPREC) :: i,j,incrn1,incrn2,n1,n2,n1st,n2st
    integer (KIND=IPREC),dimension(n1,n2) :: ia

    write(6,1000) (j, j=n2st,n2,incrn2)
    do i=n1st,n1,incrn1
       write(6,1010) i, (ia(i,j),j=n2st,n2,incrn2)
    enddo

01000 format(4i25)
01010 format(1x,i4,10i5,/(5x,10i5))

  end subroutine pmtxi
end module utils
