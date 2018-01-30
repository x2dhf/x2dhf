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
! ### putin ###

!     Immerses FUN array into WORK and adds boundary values from symmetry

!     Immerses FUN array into WORK and provides missing values along (1,i), (nni,i) and (j,1) 
!     boundaries according to orbital symmetry isym.
!     Values along (j,nmi) boundary are set to zero as orbitals decay exponentially in this region.

subroutine putin (nni,nmi,isym,fun,work)
  use params

  implicit none
  integer :: i,j,isym,jj,nni,nmi

  real (PREC), dimension(nni,nmi) :: fun
  real (PREC), dimension(nni+8,nmi+8) :: work

  !     fill the interior of work array

  do i=1,nmi
     do j=1,nni
        work(j+4,i+4)=fun(j,i)
     enddo
  enddo

  !     isym = 1 - even symmetry, isym =-1 - odd symmetry
  !     values over i=nmi bondary are determined from the asymptotic expansion
  
  !     mu=nmu+5....nmu+8
  !     the following code is necessary since the derivatives must be 
  !     calculated up to i=nmi
  
  if (isym.eq.1) then
     
     !  mu=1...4	
     do i=2,5
        do j=1,nni
           work(j+4,6-i)= fun(j,i)
        enddo
     enddo
     
     !        ni=1...4
     do i=1,nmi
        do j=2,5
           work(6-j,i+4)= fun(j,i)
        enddo
     enddo

     !        ni=ni+4...ni+8
     do i=1,nmi
        jj=0
        do j=nni-4,nni-1
           jj=jj+1
           work(nni+9-jj,i+4)= fun(j,i)
        enddo
     enddo
  else
     do i=2,5
        do j=1,nni
           work(j+4,6-i)=-fun(j,i)
        enddo
     enddo
     
     do i=1,nmi
        do j=2,5
           work(6-j,i+4)=-fun(j,i)
        enddo
     enddo
     
     do i=1,nmi
        jj=0
        do j=nni-4,nni-1
           jj=jj+1
           work(nni+9-jj,i+4)=-fun(j,i)
        enddo
     enddo
  endif

end subroutine putin


subroutine putin1 (nni,nmi,isym4nu,isym4mu,fun,work)
  use params

  implicit none
  integer :: i,j,isym4nu,isym4mu,jj,nni,nmi

  real (PREC), dimension(nni,nmi) :: fun
  real (PREC), dimension(nni+8,nmi+8) :: work

  !     fill the interior of work array

  do i=1,nmi
     do j=1,nni
        work(j+4,i+4)=fun(j,i)
     enddo
  enddo


  do i=nmi+1,nmi+4
     do j=1,nni
        work(j+4,i+4)=fun(j,nmi)
     enddo
  enddo

  
  !     isym = 1 - even symmetry, isym =-1 - odd symmetry
  !     values over i=nmi bondary are determined from the asymptotic expansion
  
  !     mu=nmu+5....nmu+8
  !     the following code is necessary since the derivatives must be 
  !     calculated up to i=nmi
  
  if (isym4mu.eq.1) then
     !  mu=1...4	
     do i=2,5
        do j=1,nni
           work(j+4,6-i)= fun(j,i)
        enddo
     enddo
  else
     do i=2,5
        do j=1,nni
           work(j+4,6-i)=-fun(j,i)
        enddo
     enddo
  endif

  if (isym4nu.eq.1) then
     !        ni=1...4
     do i=1,nmi
        do j=2,5
           work(6-j,i+4)= fun(j,i)
        enddo
     enddo

     !        ni=ni+4...ni+8
     do i=1,nmi
        jj=0
        do j=nni-4,nni-1
           jj=jj+1
           work(nni+9-jj,i+4)= fun(j,i)
        enddo
     enddo
  else
     do i=1,nmi
        do j=2,5
           work(6-j,i+4)=-fun(j,i)
        enddo
     enddo
     
     do i=1,nmi
        jj=0
        do j=nni-4,nni-1
           jj=jj+1
           work(nni+9-jj,i+4)=-fun(j,i)
        enddo
     enddo
  endif

end subroutine putin1
