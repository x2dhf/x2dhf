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
! ### prtmat ###

!     Prints an array in a formatted way. Orbitals and potentials are
!     stored in one-dimensional arrays. When they are printed as a
!     two-dimensional ones (\nu=0,\mu_1) element coresponds to A centre
!     and (\nu=\pi,\mu_1) --  B. 

subroutine prtmat (m,n,a,ioutmat)
  use params

  implicit none
  integer :: im,in,ioutmat,m,n
  real (PREC), dimension(m,n) :: a
  
  do im=1,m
     write(ioutmat,1000) (a(im,in),in=1,n)
  enddo
  !     when preparing data for matlab use Ew.d or Fw.d format 
01000 format(500E16.8)
  !     01000 format(5F15.6)
  
end subroutine prtmat


