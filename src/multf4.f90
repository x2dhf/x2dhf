! ***************************************************************************
! *                                                                         *
! *   Copyright (C) 2010 Jacek Kobus <jkob@fizyka.umk.pl>                   *
! *                                                                         *
! *   This program is free software; you can redistribute it and/or modify  *
! *   it under the terms of the GNU General Public License version 2 as     *
! *   published by the Free Software Foundation.                            *
! *                                                                         *
! ***************************************************************************
! ### multf4 ###

!     Multiplies a given array (of length mxnmu) by F4

module multf4_m
  implicit none
contains
  subroutine multf4(a)
    use params
    use discret
    !  use commons8


    implicit none
    integer :: i,ii,j
    real (PREC), dimension(*) :: a

    do i=1,mxnmu
       ii=(i-1)*nni
       do j=1,nni
          a(ii+j)=a(ii+j)*half*r*vxi(i)
       enddo
    enddo
  end subroutine multf4
end module multf4_m
