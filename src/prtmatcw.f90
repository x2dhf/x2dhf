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
! ### prtmatcw ###

!     Prints an array in a formatted way. Orbitals and potentials are
!     stored in one-dimensional arrays. When they are printed as a
!     two-dimensional ones (\nu=0,\mu_1) element coresponds to A centre
!     and (\nu=\pi,\mu_1) --  B. 

subroutine prtmatcw (m,n,a,ioutmat)
  use params

  implicit none
  integer :: im,in,ioutmat,m,n
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
! ### prtmatrw ###

!     Prints an array in a formatted way. Orbitals and potentials are
!     stored in one-dimensional arrays. When they are printed as a
!     two-dimensional ones (\nu=0,\mu_1) element coresponds to A centre
!     and (\nu=\pi,\mu_1) --  B. 

subroutine prtmatrw (m,n,a,ioutmat)
  use params

  implicit none
  integer :: im,in,ioutmat,m,n
  real (PREC), dimension(m,n) :: a
  
  do im=1,m
     write(ioutmat,1000) (a(im,in),in=1,n)
  enddo

  ! do in=1,n
  !    write(ioutmat,1000) (a(im,in),im=1,m)
  ! enddo


  !     when preparing data for matlab use Ew.d or Fw.d format 
01000 format(5E25.16)
  !     01000 format(5F15.6)
  
end subroutine prtmatrw


subroutine prtmatcw1 (m,n,a,ioutmat)
  use params

  implicit none
  integer :: im,in,ioutmat,m,n
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

