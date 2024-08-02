! SPDX-License-Identifier: GPL-2.0-or-later

! Copyright (C) 2010-2024  Jacek Kobus 

module nabla
  implicit none
contains

  ! ### n2rho ###
  !
  !     Calculates nabla^2 rho and returns it as wk3
  !
  subroutine n2rho(psi,wk0,wk1,wk2,wk4,wk3)
    use params
    use commons
    use blas
    use discrete
    use inout
    use utils

    
    implicit none
    integer (KIND=IPREC) :: i,ii,iborb,iorb,isym,j,k
    !real (PREC), parameter :: const53=5.0_PREC/3.0_PREC
    real (PREC) :: t1
    real (PREC), dimension(*) :: psi,wk0,wk1,wk2,wk3,wk4

    ! it is assumed that the density has even symmetry
    isym=1
    call zeroArray(mxsize,wk3)    
    do iorb=1,norb
       iborb =i1b (iorb)
       call dcopy (mxsize,psi(iborb),ione,wk4,ione)
       call prod (mxsize,psi(iborb),wk4)
       call dscal (mxsize,occ(iorb),wk4,ione)

       call putin  (nni,mxnmu,isym,wk4,wk1)
       call diffnu (mxnmu,wk1,wk0,wk4,wk2)
       call putout (nni,mxnmu,wk4,wk0)
       call add (mxsize,wk4,wk3)
       
       call diffmu (mxnmu,wk1,wk2)
       call putout (nni,mxnmu,wk0,wk2)

       call add (mxsize,wk0,wk3)
    enddo
    

    ! take care of 4/[R^2(xi^2-eta^2)] factor in Laplasian
    do i=1,mxnmu
       ii=(i-1)*nni
       do j=1,nni
          t1=r*r*(vxi(i)*vxi(i)-veta(j)*veta(j))
          if (t1.lt.precis) then
             wk3(ii+j)=0.0_PREC
          else
             wk3(ii+j)=wk3(ii+j)*(four/t1)
          endif

       enddo
    enddo
    
    j=1
    !print *,"n2f before interpolation: i=1,j=1", wk3(j)
    wk3(j)=zero
    do k=1,5
       wk3(j)=wk3(j)+ wk3(j+k)*exeven(k)
    end do
    !print *,"n2f after interpolation:         ", wk3(j)

    j=nni
    !print *,"n2f before interpolation: i=1,j=nni", wk3(j)
    wk3(j)=zero
    do k=1,5
       wk3(j)=wk3(j)+ wk3(j-k)*exeven(k)
    end do
    !print *,"n2f after interpolation:         ", wk3(j)

  end subroutine n2rho

  
  ! ### n2f ###
  !
  !     Calculates nabla^2 f and returns it as wk3
  !
  subroutine n2f(f,wk0,wk1,wk2,wk3)
    use params
    use commons
    use discrete
    use utils
    use inout

    implicit none
    integer (KIND=IPREC) :: i,ii,isym,j,k
    !real (PREC), parameter :: const53=5.0_PREC/3.0_PREC
    real (PREC) :: t1
    real (PREC), dimension(*) :: f,wk0,wk1,wk2,wk3

    ! it is assumed that the function f has even symmetry

    isym=1

    call putin  (nni,mxnmu,isym,f,wk1)
    call diffnu (mxnmu,wk1,wk0,wk3,wk2)
    call putout (nni,mxnmu,wk3,wk0)
    call diffmu (mxnmu,wk1,wk2)
    call putout (nni,mxnmu,wk0,wk2)

    call add (mxsize,wk0,wk3)

    ! take care of 4/[R^2(xi^2-eta^2)] factor in Laplasian
    do i=1,mxnmu
       ii=(i-1)*nni
       do j=1,nni
          t1=r*r*(vxi(i)*vxi(i)-veta(j)*veta(j))
          if (t1.lt.precis) then
             wk3(ii+j)=0.0_PREC
          else
             wk3(ii+j)=wk3(ii+j)*(four/t1)
          endif
       enddo
    enddo
    
    j=1
    !print *,"n2f before interpolation: i=1,j=1", wk3(j)
    wk3(j)=zero
    do k=1,5
       wk3(j)=wk3(j)+ wk3(j+k)*exeven(k)
    end do
    !print *,"n2f after interpolation:         ", wk3(j)

    j=nni
    !print *,"n2f before interpolation: i=1,j=nni", wk3(j)
    wk3(j)=zero
    do k=1,5
       wk3(j)=wk3(j)+ wk3(j-k)*exeven(k)
    end do
    !print *,"n2f after interpolation:         ", wk3(j)

  end subroutine n2f


  subroutine laplace(psi,e,wk0,wk1,wk2,wk3,wk4)
    use params
    use commons
    use blas
    use discrete
    use inout
    use utils
    
    implicit none
    integer (KIND=IPREC) :: i,ii,iorb,isym,j,k
    real (PREC) :: t1,w
    real (PREC), dimension(*) :: psi,e,wk0,wk1,wk2,wk3,wk4

    ! it is assumed that the function f has even symmetry

    call zeroArray(mxsize,wk4)

    do iorb=1,norb
       isym=isymOrb(iorb)
       
       call putin  (nni,mxnmu,isym,psi(i1b(iorb)),wk1)
       call diffnu (mxnmu,wk1,wk0,wk3,wk2)
       call putout (nni,mxnmu,wk3,wk0)

       call diffmu (mxnmu,wk1,wk2)
       call putout (nni,mxnmu,wk0,wk2)

       call add (mxsize,wk0,wk3)
       
       if (mm(iorb).ne.0) then
          ! nuclear energy for non-sigma orbitals contains contribution from e term (in toten
          ! this term is correctly added to the kinetic energy contribution); e enters the
          ! expression with minus sign which is already incorporated in e
          w=dble(mm(iorb)*mm(iorb))
          call prodas (mxsize,w,e,psi(i1b(iorb)),wk3)
       endif

       call prod (mxsize,psi(i1b(iorb)),wk3)
       call dscal (mxsize,occ(iorb),wk3,ione)

       call add (mxsize,wk3,wk4)
    enddo

    ! take care of 4/[R^2(xi^2-eta^2)] factor in Laplasian and the factor 2 (from 2\phi\nabla^2 \phi term)
    do i=1,mxnmu
       ii=(i-1)*nni
       do j=1,nni
          t1=r*r*(vxi(i)*vxi(i)-veta(j)*veta(j))
          if (t1.lt.precis) then
             wk4(ii+j)=0.0_PREC
          else
             wk4(ii+j)=two*(four/t1)*wk4(ii+j)
          endif
          
       enddo
    enddo
    
    j=1
    !print *,"n2f before interpolation: i=1,j=1", wk3(j)
    wk4(j)=zero
    do k=1,5
       wk4(j)=wk4(j)+ wk4(j+k)*exeven(k)
    end do
    !print *,"n2f after interpolation:         ", wk3(j)

    j=nni
    !print *,"n2f before interpolation: i=1,j=nni", wk3(j)
    wk4(j)=zero
    do k=1,5
       wk4(j)=wk4(j)+ wk4(j-k)*exeven(k)
    end do
    !print *,"n2f after interpolation:         ", wk3(j)

  end subroutine laplace


  
  ! ### nfng ###
  !
  !     Calculates nabla f nabla g; the result is returned in wk3
  !
  subroutine nfng(f,g,fmu,fni,gmu,gni,wk0,wk1,wk2,wk3)
    use params
    use discrete
    use commons
    use utils
    use inout

    implicit none
    integer (KIND=IPREC) :: i,ii,isym,j
    real (PREC) :: t1,z
    real (PREC), dimension(*) ::  f,g,fmu,fni,gmu,gni,wk0,wk1,wk2,wk3

    !   it is assumed that functions have isym=1 symmetry (sigma-type
    !   orbitals or orbital density)
    isym=1

    !   calculate df/dni and df/dmu derivatives
    call putin   (nni,mxnmu,isym,f,wk3)
    call diff1nu (mxnmu,wk3,wk0,wk1,wk2)
    call putout (nni,mxnmu,fni,wk0)

    !    print *,'f: ', (f(i),i=1000,1003)
    !    print *,'fni: ',(fni(i),i=1000,1003)

    call diff1mu (mxnmu,wk3,wk2)
    call putout (nni,mxnmu,fmu,wk2)

    !    print *,'fmu: ',(fmu(i),i=1000,1003)

    !   calculate dg/dni and dg/dmu derivatives
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
             z=half*r*vxi(i)*veta(j)
             wk0(ii+j)=wk0(ii+j)*(  ( half*r*veta1(j)*vxi(i) )**2 &
                                  + ( ( vxi(i)*z-half*r*veta(j) )/vxi1(i) )**2 )
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
             z=half*r*vxi(i)*veta(j)
             !wk1(ii+j)=wk1(ii+j)*((half*r*vxi1(i)*veta1(j)*veta(j))**2 &
             wk1(ii+j)=wk1(ii+j)*(  ( half*r*vxi1(i)*veta(j) )**2 &
                                  + ( ( veta(j)*z-half*r*vxi(i) )/veta1(j) )**2 )
          endif
       enddo
    enddo

    ! calculate df/dmu dg/dni+ df/dni dg/dmu term (wk3)

    call prod2(mxsize,fmu,gni,wk2)
    call prod2(mxsize,fni,gmu,wk3)
    call add  (mxsize,wk2,wk3)

    do i=1,mxnmu
       ii=(i-1)*nni
       do j=1,nni
          if (veta1(j).lt.precis.or.vxi1(i).lt.precis) then
             wk3(ii+j)=0.0_PREC
          else
             z=half*r*vxi(i)*veta(j)
             wk3(ii+j)=wk3(ii+j) &
                                *( r*r/four*vxi1(i)*veta1(j)*vxi(i)*veta(j) &
                                  + ( vxi(i)*z-half*r*veta(j) )*( veta(j)*z-half*r*vxi(i) )/( vxi1(i)*veta1(j) ) )
          endif
       enddo
    enddo

    call add  (mxsize,wk0,wk3)
    call add  (mxsize,wk1,wk3)
    
    do i=1,mxnmu
       ii=(i-1)*nni
       do j=1,nni
          t1=r*r*(vxi(i)*vxi(i)-veta(j)*veta(j))/four
          if (t1.lt.precis) then
             wk3(ii+j)=0.0_PREC
          else
             wk3(ii+j)=wk3(ii+j)/(t1*t1)
          endif
          ! write (*,'(2i4,e15.4)') i,j, wk3(ii+j)
       enddo
    enddo

  end subroutine nfng

  ! ### tau ###
  !
  !     Calculates nabla f nabla g; the result is returned in wk3
  !
  subroutine tau(psi,f,g,fmu,fni,gmu,gni,wk0,wk1,wk2,wk3,wk4)
    use blas
    use params
    use discrete
    use commons
    use inout
    use utils
    
    implicit none
    integer (KIND=IPREC) :: i,ii,iorb,isym,j,m
    real (PREC) :: t1,x,y,z,thetax2,thetay2,xtilde,ytilde,theta
    real (PREC), dimension(*) ::  psi,f,g,fmu,fni,gmu,gni,wk0,wk1,wk2,wk3,wk4

    !   it is assumed that functions have isym=1 symmetry (sigma-type
    !   orbitals or orbital density)

    call zeroArray(mxsize,wk4)
    
    do iorb=1,norb
       isym=isymOrb(iorb)
       
       call dcopy (mxsize,psi(i1b(iorb)),ione,f,ione)
       call dcopy (mxsize,psi(i1b(iorb)),ione,g,ione)       
       
       !   calculate df/dni and df/dmu derivatives
       call putin   (nni,mxnmu,isym,f,wk3)
       call diff1nu (mxnmu,wk3,wk0,wk1,wk2)
       call putout (nni,mxnmu,fni,wk0)

       call diff1mu (mxnmu,wk3,wk2)
       call putout (nni,mxnmu,fmu,wk2)

       !   calculate dg/dni and dg/dmu derivatives
       call putin   (nni,mxnmu,isym,g,wk3)
       call diff1nu (mxnmu,wk3,wk0,wk1,wk2)
       call putout (nni,mxnmu,gni,wk0)

       call diff1mu (mxnmu,wk3,wk2)
       call putout (nni,mxnmu,gmu,wk2)

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
                z=half*r*vxi(i)*veta(j)
                wk0(ii+j)=wk0(ii+j)*(  ( half*r*veta1(j)*vxi(i) )**2 &
                     + ( ( vxi(i)*z-half*r*veta(j) )/vxi1(i) )**2 )
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
                z=half*r*vxi(i)*veta(j)
                !wk1(ii+j)=wk1(ii+j)*((half*r*vxi1(i)*veta1(j)*veta(j))**2 &
                wk1(ii+j)=wk1(ii+j)*(  ( half*r*vxi1(i)*veta(j) )**2 &
                     + ( ( veta(j)*z-half*r*vxi(i) )/veta1(j) )**2 )
             endif
          enddo
       enddo

    
       ! calculate df/dmu dg/dni+ df/dni dg/dmu term (wk3)
       
       call prod2(mxsize,fmu,gni,wk2)
       call prod2(mxsize,fni,gmu,wk3)
       call add  (mxsize,wk2,wk3)
       
       do i=1,mxnmu
          ii=(i-1)*nni
          do j=1,nni
             if (veta1(j).lt.precis.or.vxi1(i).lt.precis) then
                wk3(ii+j)=0.0_PREC
             else
                z=half*r*vxi(i)*veta(j)
                wk3(ii+j)=wk3(ii+j) &
                     *( r*r/four*vxi1(i)*veta1(j)*vxi(i)*veta(j) &
                     + ( vxi(i)*z-half*r*veta(j) )*( veta(j)*z-half*r*vxi(i) )/( vxi1(i)*veta1(j) ) )
             endif
          enddo
       enddo

       call add  (mxsize,wk0,wk3)
       call add  (mxsize,wk1,wk3)
       
       do i=1,mxnmu
          ii=(i-1)*nni
          do j=1,nni
             t1=r*r*(vxi(i)*vxi(i)-veta(j)*veta(j))/four
             if (t1.lt.precis) then
                wk3(ii+j)=0.0_PREC
             else
                wk3(ii+j)=wk3(ii+j)/(t1*t1)
             endif
             !write (*,'(2i4,e15.4)') i,j, wk3(ii+j)
          enddo
       enddo


       ! (\partial \theta over \partial x)^2 + (\partial \theta over \partial y)^2


       theta=pii/four
       do i=1,mxnmu
          ii=(i-1)*nni
          do j=1,nni
             if (abs(vxi1(i))>precis.and.abs(veta1(j))>precis &
                  .and. abs(one-xtilde*xtilde)>precis &
                  .and. abs(one-ytilde*ytilde)>precis ) then

                x=(r/two)*vxi1(i)*veta1(j)*cos(theta)
                y=(r/two)*vxi1(i)*veta1(j)*sin(theta)
                xtilde=cos(theta)
                ytilde=sin(theta)
                thetax2=one/(one-xtilde*xtilde)*four/r**6*&
                     ( r*r/vxi1(i)/veta1(j)-four*x*x/vxi1(i)**3/veta1(j)**3 )**2
                
                thetay2=one/(one-ytilde*ytilde)*four/r**6*&
                     ( r*r/vxi1(i)/veta1(j)-four*y*y/vxi1(i)**3/veta1(j)**3 )**2
                wk2(ii+j)=(thetax2+thetay2)*f(ii+j)*g(ii+j)*mm(iorb)*mm(iorb)
             else
                wk2(ii+j)=0.0_PREC
             endif
             !write(*,'(e16.8,2i5,e16.8)') theta,i,j,wk3(ii+j)
          enddo
       enddo
       call add (mxsize,wk2,wk3)
       
       ! theta=pii/four
       ! do i=1,mxnmu,20
       !    ii=(i-1)*nni
       !    do j=1,nni,20

       !       if (abs(vxi1(i))>precis.and.abs(veta1(j))>precis &
       !            .and. abs(one-xtilde*xtilde)>precis &
       !            .and. abs(one-ytilde*ytilde)>precis ) then
             
       !          x=(r/two)*vxi1(i)*veta1(j)*cos(theta)
       !          y=(r/two)*vxi1(i)*veta1(j)*sin(theta)
       !          xtilde=cos(theta)
       !          ytilde=sin(theta)
       !          thetax2=one/(one-xtilde*xtilde)*four/r**6*&
       !               ( r*r/vxi1(i)/veta1(j)-four*x*x/vxi1(i)**3/veta1(j)**3 )**2
                
       !          thetay2=one/(one-ytilde*ytilde)*four/r**6*&
       !               ( r*r/vxi1(i)/veta1(j)-four*y*y/vxi1(i)**3/veta1(j)**3 )**2
       !          wk3(ii+j)=thetax2+thetay2
       !       else
       !          wk3(ii+j)=0.0_PREC
       !       endif
       !       write(*,'(e16.8,2i5,e16.8)') theta,i,j,wk3(ii+j)
       !    enddo
       ! enddo

       
       ! if (mm(iorb1).ne.0) then 
       !    do i=1,mxnmu
       !       ii=(i-1)*nni
       !       do j=1,nni
       !          t1=four/(r*r))*mm(iorb1)*mm(iorb1)
       !          t2=(vxi1(i)*vxi1(i)*veta(j)*veta(j))**2 &
       !               + (two/r)*(two/r) (vxi1(i)*vxi1(i)*veta(j)*veta(j)*cos(theta)*sin(theta)-one)
       !          if (t1.lt.precis) then
       !             wk3(ii+j)=0.0_PREC
       !          else
       !             wk3(ii+j)=wk3(ii+j)/(t1*t1)
       !          endif
       !          !write (*,'(2i4,e15.4)') i,j, wk3(ii+j)
       !       enddo
       !    enddo
       ! endif

       call dscal (mxsize,occ(iorb),wk3,ione)    
       call add (mxsize,wk3,wk4)
    enddo

  end subroutine tau


  ! ### testn2f ###
  !
  !     Tests nabla^2 f
  !
  subroutine testn2f(f,wk0,wk1,wk2,wk3)
    use params
    use commons
    use discrete
    use utils

    implicit none
    integer (KIND=IPREC) :: i,ii,itest,j
    real (PREC) :: x,z
    real (PREC), dimension(*) :: f,wk0,wk1,wk2,wk3

    itest=1

    !   take care of 4/[R^2(xi^2-eta^2)] factor in Laplasian
    if (itest.eq.1) then
       do i=1,mxnmu
          ii=(i-1)*nni
          do j=1,nni
             !             t1=r*r*(vxi(i)*vxi(i)-veta(j)*veta(j))
             x=(r/2.0_PREC)*vxi1(i)*veta1(j)
             z=(r/2.0_PREC)*vxi(i)*veta(j)
             f(ii+j)=1.0_PREC
             !            write (*,'(2i4,e15.4)') i,j, wk3(ii+j)
          enddo
       enddo

       call n2f(f,wk0,wk1,wk2,wk3)
       call pmtx(nni,mxnmu,wk3,ione,ione,incrni,incrmu)

    endif

  end subroutine testn2f

  
  ! ### testnfng ###
  !
  !     Tests nabla f nabla g
  !
  subroutine testnfng(rhot,grhot,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7)
    use params
    use commons
    use discrete
    use utils

    implicit none
    integer (KIND=IPREC) :: i,ii,itest,j
    real (PREC), dimension(*) ::  wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7,rhot,grhot

    itest=1

    if (itest.eq.1) then
       do i=1,mxnmu
          ii=(i-1)*nni
          do j=1,nni
             rhot(ii+j) =0.d0
             grhot(ii+j)=0.d0
          enddo
       enddo

       call nfng (rhot,grhot,wk0,wk1,wk2,wk3,wk4,wk5,wk6,wk7)
       call pmtx(nni,mxnmu,wk7,ione,ione,incrni,incrmu)

    endif
  end subroutine testnfng
end module nabla
