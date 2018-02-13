! ***************************************************************************
! *                                                                         *
! *   Copyright (C) 1996-2010 Jacek Kobus <jkob@fizyka.umk.pl>              *
! *   Copyright (C) 2018      Susi Lehtola                                  *
! *                                                                         *
! *   This program is free software; you can redistribute it and/or modify  *
! *   it under the terms of the GNU General Public License version 2 as     *
! *   published by the Free Software Foundation.                            *
! *                                                                         *
! ***************************************************************************
! ### rexponents ###
!
!     Reads exponents and contraction coefficients from
!     a GAUSSIAN94 output obtained with gfinput keyword
!
!     Warning!
!     Basis functions of g, h and higher symmertries are not allowed
!     Uncontracted basis set must be used
!
subroutine rexponents(ibc,ib,istop)
  use params
  use commons8

  implicit none
  integer :: i,ib,ibc,icent,istop,j,k,maxcf,ncf

  parameter (maxcf=100)

  real (PREC), dimension(maxcf) :: expon(maxcf)
  real (PREC), dimension(maxcf,7) :: dcoef

  character(80) :: line
  character(2) :: symlab

  integer :: stat

  istop=0

  ! Read the centre number
  read (7,'(A)') line
  if (trim(line).eq.'') then
     ! Check if the file ended
     istop=1
     return
  end if

  read (line,*,iostat=stat) icent
  if(stat.ne.0) then
     istop=0
     return
  end if
  if (icent.eq.0) then
     istop=0
     return
  elseif (icent.eq.-9999) then
     return
  endif

  !     read symmetry label and the number of exponents
  do i=1,maxbasis
     read (7,'(A)') line

     if (trim(line).ne.' ****') then
        read (line,*) symlab, ncf
        if (ncf.gt.maxcf) then
           print *,'r_exponents: too many primitive functions per single contracted gaussian; increase maxcf'
           stop
        endif

        if (trim(symlab).eq.'S') then
           do j=1,ncf
              read (7,'(A)') line
              read (line, *) expon(j), dcoef(j,1)
           enddo
           !              define s-type gaussians
           ibc=ibc+1
           do j=1,ncf
              ib=ib+1
              if (ib.gt.maxbasis) goto 990
              coeff(ib)=dcoef(j,1)
              primexp(ib)=expon(j)
              ixref(ib)=ibc
              lprim(ib)=0
              mprim(ib)=0
           enddo
        elseif (trim(symlab).eq.'SP') then
           do j=1,ncf
              read (7,'(A)') line
              read (line,*) expon(j), dcoef(j,1), dcoef(j,2)
           enddo
           !              define s-type gaussians
           ibc=ibc+1
           do j=1,ncf
              ib=ib+1
              if (ib.gt.maxbasis) goto 990
              coeff(ib)=dcoef(j,1)
              primexp(ib)=expon(j)
              ixref(ib)=ibc
              lprim(ib)=0
              mprim(ib)=0
           enddo
           !              define p-type orbitals (x,y,z)
           do k=1,3
              ibc=ibc+1
              do j=1,ncf
                 ib=ib+1
                 if (ib.gt.maxbasis) goto 990
                 coeff(ib)=dcoef(j,2)
                 primexp(ib)=expon(j)
                 ixref(ib)=ibc
                 lprim(ib)=1
                 if (k.eq.1) mprim(ib)=+1
                 if (k.eq.2) mprim(ib)=-1
                 if (k.eq.3) mprim(ib)= 0
              enddo
           enddo
        elseif (trim(symlab).eq.'P') then
           do j=1,ncf
              read (7,'(A)') line
              read (line, *) expon(j), dcoef(j,2)
           enddo
           !              define p-type orbitals (x,y,z)
           do k=1,3
              ibc=ibc+1
              do j=1,ncf
                 ib=ib+1
                 if (ib.gt.maxbasis) goto 990
                 coeff(ib)=dcoef(j,2)
                 primexp(ib)=expon(j)
                 ixref(ib)=ibc
                 lprim(ib)=1
                 if (k.eq.1) mprim(ib)=+1
                 if (k.eq.2) mprim(ib)=-1
                 if (k.eq.3) mprim(ib)= 0
              enddo
           enddo

           !              define d-type orbitals (d0,d1,d-1,d2,d-2)
        elseif (trim(symlab).eq.'D') then
           do j=1,ncf
              read (7,'(A)') line
              read (line, *) expon(j), dcoef(j,3)
           enddo

           do k=1,5
              ibc=ibc+1
              do j=1,ncf
                 ib=ib+1
                 if (ib.gt.maxbasis) goto 990
                 coeff(ib)=dcoef(j,3)
                 primexp(ib)=expon(j)
                 ixref(ib)=ibc
                 lprim(ib)=2
                 if (k.eq.1) mprim(ib)= 0
                 if (k.eq.2) mprim(ib)=+1
                 if (k.eq.3) mprim(ib)=-1
                 if (k.eq.4) mprim(ib)=+2
                 if (k.eq.5) mprim(ib)=-2
              enddo
           enddo
           !              define f-type orbitals (f0,f1,f-1,f2,f-2,f3,f-3)`
        elseif (trim(symlab).eq.'F') then
           do j=1,ncf
              read (7,'(A)') line
              read (line, *) expon(j), dcoef(j,4)
           enddo
           do k=1,7
              ibc=ibc+1
              do j=1,ncf
                 ib=ib+1
                 if (ib.gt.maxbasis) goto 990
                 coeff(ib)=dcoef(j,4)
                 primexp(ib)=expon(j)
                 ixref(ib)=ibc
                 lprim(ib)=3
                 if (k.eq.1) mprim(ib)= 0
                 if (k.eq.2) mprim(ib)=+1
                 if (k.eq.3) mprim(ib)=-1
                 if (k.eq.4) mprim(ib)=+2
                 if (k.eq.5) mprim(ib)=-2
                 if (k.eq.6) mprim(ib)=+3
                 if (k.eq.7) mprim(ib)=-3
              enddo
           enddo
           !              define g-type orbitals (g0,g1,g-1,g2,g-2,g3,g-3,g4,g-4)`
        elseif (symlab.eq.'g       ') then
           do j=1,ncf
              read (7,'(A)') line
              read (line, *) expon(j), dcoef(j,5)
           enddo
           do k=1,9
              ibc=ibc+1
              do j=1,ncf
                 ib=ib+1
                 if (ib.gt.maxbasis) goto 990
                 coeff(ib)=dcoef(j,5)
                 primexp(ib)=expon(j)
                 ixref(ib)=ibc
                 lprim(ib)=4
                 if (k.eq.1) mprim(ib)= 0
                 if (k.eq.2) mprim(ib)=+1
                 if (k.eq.3) mprim(ib)=-1
                 if (k.eq.4) mprim(ib)=+2
                 if (k.eq.5) mprim(ib)=-2
                 if (k.eq.6) mprim(ib)=+3
                 if (k.eq.7) mprim(ib)=-3
                 if (k.eq.8) mprim(ib)=+4
                 if (k.eq.9) mprim(ib)=-4
              enddo
           enddo
           !              define h-type orbitals (h0,h1,h-1,h2,h-2,h3,h-3,h4,h-4,h5,h-5)`
        elseif (symlab.eq.'h       ') then
           do j=1,ncf
              read (7,'(A)') line
              read (line, *) expon(j), dcoef(j,6)
           enddo
           do k=1,11
              ibc=ibc+1
              do j=1,ncf
                 ib=ib+1
                 if (ib.gt.maxbasis) goto 990
                 coeff(ib)=dcoef(j,6)
                 primexp(ib)=expon(j)
                 ixref(ib)=ibc
                 lprim(ib)=4
                 if (k.eq.1) mprim(ib)= 0
                 if (k.eq.2) mprim(ib)=+1
                 if (k.eq.3) mprim(ib)=-1
                 if (k.eq.4) mprim(ib)=+2
                 if (k.eq.5) mprim(ib)=-2
                 if (k.eq.6) mprim(ib)=+3
                 if (k.eq.7) mprim(ib)=-3
                 if (k.eq.8) mprim(ib)=+4
                 if (k.eq.9) mprim(ib)=-4
                 if (k.eq.10) mprim(ib)=+5
                 if (k.eq.11) mprim(ib)=-5
              enddo
           enddo
           !              define i-type orbitals (i0,i1,i-1,i2,i-2,i3,i-3,i4,i-4,i5,i-5,i6,i-6)`
        elseif (trim(symlab).eq.'I') then
           do j=1,ncf
              read (7,'(A)') line
              read (line, *) expon(j), dcoef(j,7)
           enddo
           do k=1,13
              ibc=ibc+1
              do j=1,ncf
                 ib=ib+1
                 if (ib.gt.maxbasis) goto 990
                 coeff(ib)=dcoef(j,7)
                 primexp(ib)=expon(j)
                 ixref(ib)=ibc
                 lprim(ib)=4
                 if (k.eq.1) mprim(ib)= 0
                 if (k.eq.2) mprim(ib)=+1
                 if (k.eq.3) mprim(ib)=-1
                 if (k.eq.4) mprim(ib)=+2
                 if (k.eq.5) mprim(ib)=-2
                 if (k.eq.6) mprim(ib)=+3
                 if (k.eq.7) mprim(ib)=-3
                 if (k.eq.8) mprim(ib)=+4
                 if (k.eq.9) mprim(ib)=-4
                 if (k.eq.10) mprim(ib)=+5
                 if (k.eq.11) mprim(ib)=-5
                 if (k.eq.12) mprim(ib)=+6
                 if (k.eq.13) mprim(ib)=-6
              enddo
           enddo
        else
           print *,'r_exponents: basis functions higher than i are not allowed'
           stop
        endif
     else
        istop=0
        return
     endif
  enddo
  istop=1
  return

!00904 stop 'r_exponents: error encountered when reading gauss94.out'
!00950 stop 'r_exponents: end of gauss94.out file encountered'
00990 stop 'r_exponents: too many basis functions; increase maxbasis'
!01000 format(a3,2x,i2)
end subroutine rexponents
