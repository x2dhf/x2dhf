! SPDX-License-Identifier: GPL-2.0-or-later

! Copyright (C) 2024  Jacek Kobus 

module detect
  implicit none
contains

 function detectNaN(n,a)
    use params
    implicit none

    integer (KIND=IPREC) :: detectNaN,i,n
    real (PREC), dimension(*) :: a

    detectNaN=0
    !nnans=0
    do i=1,n
       if ( .not. isnan(a(i)) ) cycle
       !nnans=nnans+1
       detectNaN=i
       return
    enddo
  end function detectNaN

  subroutine fixNaN(nu,mu,array)
    use params
    use discrete
    implicit none

    integer (KIND=IPREC) :: i,imu,inu,mu,nu,nan
    real (PREC) :: a,b
    real (PREC), dimension(nu,mu) :: array

    do inu=1,nu
       do imu=1,mu
          if ( .not. isnan(array(inu,imu)) ) cycle
          if ( (inu+2).le.nu ) then
             a=(array(inu+1,imu)-array(inu+2,imu))/(vni(inu+1)-vni(inu+2))
             b=array(inu+1,imu)-a*vni(inu+1)
             array(inu,imu)=a*vni(i)+array(inu+1,imu)-a*vni(i+1)
#ifdef PRINT
! print=470: fixNaN 1: fixing NaN for array(inu,imu)
             if (iprint(470).ne.0) write(*,'("fixNaN 1: fixing at ",2i5,e12.4)') inu,imu,array(inu,imu)
#endif
          else
             a=(array(inu-2,imu)-array(inu-1,imu))/(vni(inu-2)-vni(inu-1))
             b=array(inu-1,imu)-a*vni(inu-1)
             array(inu,imu)=a*vni(i)+array(inu-1,imu)-a*vni(i-1)
#ifdef PRINT
! print=470: fixNaN 2: fixing NaN for array(inu,imu)
             if (iprint(470).ne.0) write(*,'("fixNaN 2: fixing at ",2i5,e12.4)') inu,imu,array(inu,imu)
#endif
          endif
       enddo
    enddo
  end subroutine fixNaN

  subroutine mapDerivs(nu,mu,a)
    use params
    use discrete
    implicit none
    integer (KIND=IPREC) :: i,j,inu,imu,nu,mu,n
    integer (KIND=IPREC) :: inu4mu,inu4nu,imu4mu,imu4nu
    real (PREC) :: deriv,derivMax
    real (PREC), dimension(nu,mu) :: a

    ! find maximum value of the first derivative with respect to mu variable

    inu4mu=0
    imu4mu=0
    derivMax=0.0_PREC
    do inu=1,nu
       do imu=1,mu-1
          deriv=abs(a(inu,imu+1)-a(inu,imu))/hmu(1)
          if (deriv>derivMax) then
             derivMax=deriv
             inu4mu=inu
             imu4mu=imu
          endif
       enddo
    enddo
    write(*,'("max df(nu,mu)/d mu ",2i5,1e12.4)') inu4mu,imu4mu,derivMax

    imu4nu=0
    inu4nu=0          
    derivMax=0.0_PREC
    do imu=1,mu
       do inu=1,nu-1
          deriv=abs(a(inu+1,imu)-a(inu,imu))/hni
          if (deriv>derivMax) then
             derivMax=deriv
             inu4nu=inu
             imu4nu=imu
          endif
       enddo
    enddo
    write(*,'("max df(nu,mu)/d nu ",2i5,1e12.4)') inu4nu,imu4nu,derivMax
    
  end subroutine mapDerivs
  
end module detect
