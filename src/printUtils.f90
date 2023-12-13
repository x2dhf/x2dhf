! SPDX-License-Identifier: GPL-2.0-or-later

! Copyright (C) 1996       Leif Laaksonen, Dage Sundholm               
! Copyright (C) 1996-2023  Jacek Kobus 

module printUtils
  implicit none
contains

    subroutine out4plot(iout,a)
    use params
    use discrete
    use commons

    implicit none
    integer (KIND=IPREC) :: iout,m,n,inu,imu
    real (PREC) :: x,z
    real (PREC), dimension(nni,mxnmu) :: a

    if (iplot==1) then
       do inu=1,nni,istep    
          do imu=1,mxnmu,istep
             z=(r/2.0_PREC)*vxi(imu)*veta(inu)
             x=(r/2.0_PREC)*vxi1(imu)*veta1(inu)
             if (abs(a(inu,imu))>plotThreshold) &             
                  write(iout,'(1Pe12.4,1Pe12.4,1Pe14.4)') z,x,a(inu,imu)
          enddo
       enddo
    endif

    if (iplot==2) then
       do inu=1,nni,istep    
          do imu=1,mxnmu,istep
             if (abs(a(inu,imu))>plotThreshold) &
                  write(iout,'(1Pe12.4,1Pe12.4,1Pe14.4)') vni(inu),vmu(imu),a(inu,imu)
          enddo
       enddo
    endif

    if (iplot==3) then
       do inu=1,nni,istep    
          do imu=1,mxnmu,istep
             if (abs(a(inu,imu))>plotThreshold) &
                  write(iout,'(1Pe12.4,1Pe12.4,1Pe14.4)') veta(inu),vxi(imu),a(inu,imu)
          enddo
       enddo
    endif
    
    
  end subroutine out4plot
  
! ### prtmatcw ###
  !
  !     Prints an array in a formatted way. Orbitals and potentials are
  !     stored in one-dimensional arrays. When they are printed as a
  !     two-dimensional ones (\nu=0,\mu_1) element coresponds to A centre
  !     and (\nu=\pi,\mu_1) --  B.
  
  subroutine prtmatcw (m,n,a,ioutmat)
    use params
    
    implicit none
    integer (KIND=IPREC) :: im,in,ioutmat,m,n
    real (PREC), dimension(m,n) :: a
    
    ! do in=1,n
    !    write(ioutmat,'("   mu =",i4)') in
    !    write(ioutmat,1000) (a(im,in),im=1,m)
    ! enddo
    
    do in=1,n
       write(ioutmat,1000) (a(im,in),im=1,m)
    enddo
    
    ! when preparing data for matlab use Ew.d or Fw.d format
01000 format(5E25.16)
    
  end subroutine prtmatcw
  
  ! ### prtmatrw ###
  !
  !     Prints an array in a formatted way. Orbitals and potentials are
  !     stored in one-dimensional arrays. When they are printed as a
  !     two-dimensional ones (\nu=0,\mu_1) element coresponds to A centre
  !     and (\nu=\pi,\mu_1) --  B.
  !
  subroutine prtmatrw (m,n,a,ioutmat)
    use params
    
    implicit none
    integer (KIND=IPREC) :: im,in,ioutmat,m,n
    real (PREC), dimension(m,n) :: a
    
    do im=1,m
       write(ioutmat,1000) (a(im,in),in=1,n)
    enddo
    
    ! when preparing data for matlab use Ew.d or Fw.d format
01000 format(5E25.16)
    ! 01000 format(5F15.6)
    
  end subroutine prtmatrw
  
end module printUtils

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! ### prtmat ###
   !
   !     Prints an array in a formatted way. Orbitals and potentials are
   !     stored in one-dimensional arrays. When they are printed as a
   !     two-dimensional ones (\nu=0,\mu_1) element coresponds to A centre
   !     and (\nu=\pi,\mu_1) --  B.
   !
   subroutine prtmat (m,n,a,ioutmat)
     use params

     implicit none
     integer (KIND=IPREC) :: im,in,ioutmat,m,n
     real (PREC), dimension(m,n) :: a

     do im=1,m
        write(ioutmat,1000) (a(im,in),in=1,n)
     enddo
     !     when preparing data for matlab use Ew.d or Fw.d format
 01000 format(500E16.8)
     !     01000 format(5F15.6)

   end subroutine prtmat

   ! ### prtmatcw ###
   !
   subroutine prtmatcw (m,n,a,ioutmat)
     use params

     implicit none
     integer (KIND=IPREC) :: im,in,ioutmat,m,n
     real (PREC), dimension(m,n) :: a

     ! do in=1,n
     !    write(ioutmat,'("   mu =",i4)') in
     !    write(ioutmat,1000) (a(im,in),im=1,m)
     ! enddo

     do in=1,n
        write(ioutmat,1000) (a(im,in),im=1,m)
     enddo

     !   do in=1,n
     ! !     write(ioutmat,'("      mu =",i4)') in
     !      do im=1,m
     !      write(ioutmat,1000) (a(im,in),im=1,m)
     !         write(ioutmat,'(d18.10)') a(im,in)
     !      enddo
     !   enddo


     !     when preparing data for matlab use Ew.d or Fw.d format
 01000 format(5E25.16)
     !     01000 format(5F15.6)

   end subroutine prtmatcw

   subroutine prtmatcw1 (m,n,a,ioutmat)
     use params

     implicit none
     integer (KIND=IPREC) :: im,in,ioutmat,m,n
     real (PREC), dimension(m,n) :: a

     do in=1,n
        write(*,'("   mu= ",i4)') in-1
        do im=1,m
           write(*,1000) a(im,in)
        enddo
     enddo
     !     when preparing data for matlab use Ew.d or Fw.d format
 01000 format(5E25.16)
     !     01000 format(5F15.6)

   end subroutine prtmatcw1

 
 !   ! ### prtmatrw ###
 !   !
 !   !     Prints an array in a formatted way. Orbitals and potentials are
 !   !     stored in one-dimensional arrays. When they are printed as a
 !   !     two-dimensional ones (\nu=0,\mu_1) element coresponds to A centre
 !   !     and (\nu=\pi,\mu_1) --  B.

 !   subroutine prtmatrw (m,n,a,ioutmat)
 !     use params

 !     implicit none
 !     integer (KIND=IPREC) :: im,in,ioutmat,m,n
 !     real (PREC), dimension(m,n) :: a

 !     do im=1,m
 !        write(ioutmat,1000) (a(im,in),in=1,n)
 !     enddo

 !     ! do in=1,n
 !     !    write(ioutmat,1000) (a(im,in),im=1,m)
 !     ! enddo


 !     !     when preparing data for matlab use Ew.d or Fw.d format
 ! 01000 format(5E25.16)
 !     !     01000 format(5F15.6)

 !   end subroutine prtmatrw


! ### prtmatcw ###
  !
  !     Prints an array in a formatted way. Orbitals and potentials are
  !     stored in one-dimensional arrays. When they are printed as a
  !     two-dimensional ones (\nu=0,\mu_1) element coresponds to A centre
  !     and (\nu=\pi,\mu_1) --  B.
  
