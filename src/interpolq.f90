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
! ### interpolq ###
!
!     This routine works for ngrids<=3
!
subroutine interpolq
  use params
  use discret
  use commons8
  use commons16
  use lpcoeffq_m
  use vpoly1q_m

  implicit none

  integer :: i,imu,k,mup
  real (PREC) :: hmu1,hmu2,hmu3
  real (PREC16) xmu
  real (PREC16), dimension(9) :: coeffq

!
!      dimension coeffq(9)

!   |--------------------------------------|
!   |            ifill=1                   |    ngrids=1
!   |				          |
!   |--------------------------------------|
!
!   |--------------|-----------------------|
!   |     ifill=2  | ifill=4               |    ngrids=2
!   |              |                       |
!   |--------------|-----------------------|
!
!   |--------------|-----------|-----------|
!   |	ifill=2	  | ifill=3   | ifill=4	  |    ngrids=3
!   |              |           |           |
!   |--------------|-----------|-----------|

!      To relax a function defined on nni*nmu grid using the 8th-order central
!      difference formulae the function is immersed into (nni+8)*(nmu+8) grid
!      but the valued outside the boundary necessary to perform relaxation
!      within the boundary region have to be previded. The symmetry properties
!      are used to this end along ni=1, ni=nni and mu=1 lines. However the
!      boundaries between subgrids have to be treated differently.

!      At each subgrid boundary a set of the 8th-order Lagrange interpolation
!      formulae are formed in order to be able to calculate missing outer
!      region data.

!      For each value of mup (p) a Lagrange polynomial of the 8th-order
!      employing 4 grid points to the left and right of imu (imu is included)
!      is build (x_{0} corresponds to the grid boundary mu value)

!      f^{(p)}(x) = \sum_{k=p-4}^{p+4} f(x_{k})
!                   \prod_{i=p-4,i\ne k}^{p+4} {(x-x_{i}) \over (x_{k}-x_{i})}

!      Polynomials for each k have to be calculated and stored.

!      cint?(k,i) store the interpolation coefficients; i=1,...,4 corresponds
!      to appropriate (left or right) inter(sub)grid mu values in increasing
!      order, k numbers the coefficients.
!      The interpolation formula is constructed using 9 grid points somehow
!      distributed around the point at which interpolated value is sought
!      (see below).

!      Determine the grid points which will be used for the construction of the
!      polynomial

! !!!!!!!! iord is changed

  iord=9
  iord2=iord/2
  kbeg=1
  kend=9

  if (idbg(80).ne.0) then
     iord=3
     iord2=iord/2
     kbeg=4
     kend=6
  endif

  if (ngrids.eq.2) then

     !        ifill=2
     !        boundary between 1-2

     do i=1,4
        do k=1,9
           cint2(k,i)=0.0_PREC
           cint4(k,i)=0.0_PREC
        enddo
     enddo

     hmu1=hmu(1)
     hmu2=hmu(2)

     if (hmu2.lt.hmu1) then
        write(iout6,*) '2nd subgrid should be sparser than the 1st'
        stop 'interpolq: ifill=2'
     endif

     !          for each i=1,...,4 determine the smallest imup for which
     !          |vmu(imup)-vmu(iemu(1))| > i*hmu(1)
     !          and construct polynomial around imup+1-4 using 4 grid points to the
     !          left and right

     do i=1,4
        do imu=iemu(1)+1,iemu(2)
           if ((vmu(imu)-vmu(iemu(1))).ge.dble(i)*hmu1) then
              mup=imu+1-iord2
              goto 90
           endif
        enddo
00090   continue
        iadint2(i)=mup-4
        do k=kbeg,kend
           write (*,*) 'FIXME'
           call abort
           !call lpcoeffq(mup,k,coeffq)
           xmu=(vmu(iemu(1))+dble(i)*hmu1)
           cint2(k,i)=vpoly1q(xmu,coeffq)
        enddo
     enddo

     !         ifill=4

     !          to the left of the boundary (only 3 values are actually needed but
     !          4 values would be generated)
     !          for each i=1,...,4 determine the largest imup for which
     !          |vmu(imup)-vmu(iemu(1))| > i*hmu(2)
     !          and construct polynomial around such imup using 4 grid points to the
     !          left and right

     !          boundary between 1-2

     do i=1,4
        do imu=iemu(1)-1,1,-1
           if ((vmu(iemu(1))-vmu(imu)).ge.dble(i)*hmu2) then
              mup=imu
              goto 100
           endif
        enddo
00100   continue
        iadint4(5-i)=mup-4
        do k=kbeg,kend
           write (*,*) 'FIXME'
           call abort
           !call lpcoeffq(mup,k,coeffq)
           xmu=(vmu(iemu(1))-dble(i)*hmu2)
           !                reverse order of addressing columns, cf. fill4
           cint4(k,5-i)=vpoly1q(xmu,coeffq)
        enddo
     enddo

     if (idbg(81).ne.0) then
        write(iout6,*) 'iadint2',iadint2
        write(iout6,*) 'cint2'
        do i=1,4
           write(iout6,*) i,(cint2(k,i),k=1,9)
        enddo

        write(iout6,*) 'iadint4',iadint4
        write(iout6,*) 'cint4'
        do i=1,4
           write(iout6,*) i,(cint4(k,i),k=1,9)
        enddo
     endif
     return
  endif

  !      33333333333333333333333333333333333333333333333
  !       ifill=3 (ngrids=3)

  if (ngrids.eq.3) then

     !          boundary between 1-2

     hmu1=hmu(1)
     hmu2=hmu(2)

     if (hmu2.lt.hmu1) then
        write(iout6,*) '2nd subgrid should be sparser than the 1st'
        stop 'interpolq: ifill=2'
     endif

     do i=1,4
        do k=1,9
           cint2 (k,i)=0.0_PREC
           cint3l(k,i)=0.0_PREC
           cint3r(k,i)=0.0_PREC
           cint4 (k,i)=0.0_PREC
        enddo
     enddo

     !            for each i=1,...,4 determine the smallest imup for which
     !            |vmu(imup)-vmu(iemu(1))| > i*hmu(1)
     !            and construct polynomial around imup+1-4 using 4 grid points
     !            to the left and right

     do i=1,4
        do imu=iemu(1)+1,iemu(2)
           if ((vmu(imu)-vmu(iemu(1))).ge.dble(i)*hmu1) then
              mup=imu+1-iord2
              goto 91
           endif
        enddo
00091   continue
        iadint2(i)=mup-4
        do k=kbeg,kend 
           write (*,*) 'FIXME'
           call abort
!           call lpcoeffq(mup,k,coeffq)
           xmu=(vmu(iemu(1))+dble(i)*hmu1)
           cint2(k,i)=vpoly1q(xmu,coeffq)
        enddo
     enddo

     !            ifill=4

     !            to the left of the boundary (only 3 values are actually needed but
     !            4 values would be generated)
     !            for each i=1,...,4 determine the largest imup for which
     !            |vmu(imup)-vmu(iemu(1))| > i*hmu(2)
     !            and construct polynomial around such imup using 4 grid points to the
     !            left and right

     !            boundary between 2-3

     hmu3=hmu(3)
     if (hmu3.lt.hmu2) then
        write(iout6,*) '3rd subgrid should be sparser than the 2nd'
        stop 'interpolq: ifill=3'
     endif

     do i=1,4
        do imu=iemu(2)-1,1,-1
           if ((vmu(iemu(2))-vmu(imu)).ge.dble(i)*hmu3) then
              mup=imu
              goto 101
           endif
        enddo
00101   continue
        iadint4(5-i)=mup-4
        do k=kbeg,kend
           write (*,*) 'FIXME'
           call abort
!           call lpcoeffq(mup,k,coeffq)
           xmu=(vmu(iemu(2))-dble(i)*hmu3)
           !                  reverse order of addressing columns, cf. fill4
           cint4(k,5-i)=vpoly1q(xmu,coeffq)
        enddo
     enddo

     !            to the left of the boundary (only 3 values are needed but 4 ...)

     !            for each i=1,...,4 determine the largest imup for which
     !            |vmu(imup)-vmu(iemu(1))| > i*hmu(2)
     !            and construct polynomial around such imup using 4 grid points to the
     !            left and right

     !            boundary between 1-2

     do i=1,4
        do imu=iemu(1)-1,1,-1
           if ((vmu(iemu(1))-vmu(imu)).ge.dble(i)*hmu2) then
              mup=imu
              goto 110
           endif
        enddo
00110   continue
        iadint3l(5-i)=mup-4
        do k=kbeg,kend
           write (*,*) 'FIXME'
           call abort
!           call lpcoeffq(mup,k,coeffq)
           xmu=(vmu(iemu(1))-dble(i)*hmu2)
           !                  reverse order of addressing columns, cf. fill3
           cint3l(k,5-i)=vpoly1q(xmu,coeffq)
        enddo
     enddo

     !            boundary between 2-3

     !            for each i=1,...,4 determine the smallest imup for which
     !            |vmu(imup)-vmu(iemu(2))| > i*hmu(2)
     !            and construct polynomial around imup+1-4 using 4 grid points to the
     !            left and right

     do i=1,4
        do imu=iemu(2)+1,iemu(3)
           if ((vmu(imu)-vmu(iemu(2))).ge.dble(i)*hmu2) then
              mup=imu+1-iord2
              goto 120
           endif
        enddo
00120   continue
        iadint3r(i)=mup-4
        do k=kbeg,kend
           write (*,*) 'FIXME'
           call abort
!           call lpcoeffq(mup,k,coeffq)
           xmu=(vmu(iemu(2))+dble(i)*hmu2)
           cint3r(k,i)=vpoly1q(xmu,coeffq)
        enddo
     enddo

     if (idbg(81).ne.0) then
        write(iout6,*) 'iadint2',iadint2
        write(iout6,*) 'cint2'
        do i=1,4
           write(iout6,*) i,(cint2(k,i),k=1,9)
        enddo

        write(iout6,*) 'iadint3l',iadint3l
        write(iout6,*) 'cint3l'
        do i=1,4
           write(iout6,*) i,(cint3l(k,i),k=1,9)
        enddo

        write(iout6,*) 'iadint3r',iadint3r
        write(iout6,*) 'cint3r'
        do i=1,4
           write(iout6,*) i,(cint3r(k,i),k=1,9)
        enddo

        write(iout6,*) 'iadint4',iadint4
        write(iout6,*) 'cint4'
        do i=1,4
           write(iout6,*) i,(cint4(k,i),k=1,9)
        enddo
     endif
     return
  endif

  if (ngrids.gt.3) then
     write(iout6,*) 'No more than 3 subgrids are allowd'
     stop 'interpolq'
  endif

end subroutine interpolq
