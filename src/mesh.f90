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
! ### mesh ###

!     Initializes index arrays determining the ordering of mesh points during sor process
!    
!     type of ordering (itype==iorder):
!     itype=1 -- natural column-wise 
!     itype=2 -- middle
!     itype=3 -- natural row-wise 
!     itype=4 -- reverse row-wise

subroutine mesh(imcase,ig,indx1,indx2,indx6a,indx6b,indx7)
  use params
  use discret
  use commons8

  implicit none
  integer :: i,i2,i0,iend1,iend2,ii,imcase,ig,itype,j,j2,jj,nmu8,nmut,nmux,nni8,nnix

  integer, dimension(ngrd1) :: indx1
  integer, dimension(ngrd2) :: indx2
  integer, dimension(ngrd6a) :: indx6a
  integer, dimension(ngrd6b) :: indx6b
  integer, dimension(ngrd7) :: indx7


  !   initialize variables needed in sor and mcsor

  nni1=nni+8
  nni2=nni1+nni1
  nni3=nni2+nni1
  nni4=nni3+nni1
  nni5=nni4+nni1
  
  nni8=nni+8
  nmut=nmu(ig)
  nmu8=nmut+8
  
  if (imcase.eq.1) then
     iend1=nmut
     iend2=nmut-4
  elseif (imcase.eq.2) then
     iend1=nmut+4
     iend2=nmut
  endif
  
  itype=iorder(ig)
  
  if (itype.eq.1) then
     
!        Natural column-wise ordering, bondary points are excluded,
!        the upper limit of mu variable is nmu.
       
!        ?  values for points (j,nmu+1), (j,nmu+2) and (j,nmu+3) are calculated
!        ?  from the asymptotic exspansion (the indexes refer to the larger
!        ?  grid, see routine fill and refill).
!        
!        region i
!        
!        in case of more than one grid the upper limit in not nmut-4
!        but nmut
         
     ngrd1=0
     do i=6,iend1
        ii=(i-1)*nni8
        do j=6,nni+3
           ngrd1=ngrd1+1
           !              ires2(j,i)=ngrd1
           indx1(ngrd1)=ii+j
        enddo
     enddo
     
     ngrd2=0
     do i=2,iend2
        ii=(i-1)*nni
        do j=2,nni-1
           ngrd2=ngrd2+1
!              ires1(j,i)=ngrd2
!              iresni(ngrd2)=j
!              iresmu(ngrd2)=i
!              indx2(ngrd2)=ii+j
        enddo
     enddo
     
     goto 1000
  endif
  
  if (itype.eq.2) then
     
     !        The following sequence of mesh points is a modification of
     !        Laaksonen's "middle" type of sweeps
     
     !        mu changes from mumax/2 to 1 and from mumax/2 to mumax; for
     !        every mu value nu changes from pi/2 to 0 and then from pi/2 to pi
     
     i0=0
     nmux=(nmu8-1)/2
     nnix=(nni8-1)/2
     
     do i=6,iend1
        if (i.le.nmux) then
           ii=nmux-i+6
        else
           ii=i
        endif
        i2=ii-4
        
        do j=6,nni+3
           if(j.le.nnix) then
              jj=nnix-j+6
           else
              jj=j
           endif
           j2=jj-4
           
           i0=i0+1
           indx2(i0)=(i2-1)*nni  +j2
           indx1(i0)=(ii-1)*nni8 +jj
        enddo
     enddo
     
     ngrd1=i0
     ngrd2=i0
     
     !        bondary points like in the itype=1 case
     goto 1000
  endif
  
  if (itype.eq.3) then

     !        natural row-wise ordering
     
     !        region i
     
     ngrd1=0
     do j=6,nni+3
        do i=6,iend1
           ii=(i-1)*nni8
           ngrd1=ngrd1+1
           indx1(ngrd1)=ii+j
        enddo
     enddo
     
     ngrd2=0
     do j=2,nni-1
        do i=2,iend1-4
           ii=(i-1)*nni
           ngrd2=ngrd2+1
           indx2(ngrd2)=ii+j
        enddo
     enddo
     
     !        bondary points like in the itype=1 case
     
     goto 1000
     
  endif

  if (itype.eq.4) then
     
     !        reverse ordering: mi=nmut-4 sweep in ni variable, mi=nmut-5
     !        followed by another sweep in ni etc.. it should allow the
     !        bondary conditions to propagate quicker into the solution.
     
     !        region i
     
     ngrd1=0
     do j=nni+3,6,-1
        do i=iend1,6,-1
           ii=(i-1)*nni8
           ngrd1=ngrd1+1
           indx1(ngrd1)=ii+j
        enddo
     enddo
     
     ngrd2=0
     do j=nni-1,2,-1
        do i=iend1-4,2,-1
           ii=(i-1)*nni
           ngrd2=ngrd2+1
           indx2(ngrd2)=ii+j
        enddo
     enddo
     
     !        bondary points like in the itype=1 case
     
     goto 1000
  endif
  
  if (itype.eq.5) then
     
     !        modification of the 2nd ordering (for testing purposes)
     
     i0=0
     nmux=(nmu8-1)/2
     nnix=(nni8-1)/2
     
     do j=6,nni+3
        if(j.le.nnix) then
           jj=nnix-j+6
        else
           jj=j
        endif
        j2=jj-4
        
        do i=6,iend1
           if (i.le.nmux) then
              ii=nmux-i+6
           else
              ii=i
           endif
           i2=ii-4
           i0=i0+1
           indx2(i0)=(i2-1)*nni  +j2
           indx1(i0)=(ii-1)*nni8 +jj
        enddo
     enddo
     
     ngrd1=i0
     ngrd2=i0
     
     !        bondary points like in the itype=1 case
     
     goto 1000
     
  endif
      
  if (itype.eq.6) then
     
     !        another modification of the 2nd ordering (for testing purposes)
     
     i0=0
     nmux=(nmu8-1)/2
     nnix=(nni8-1)/2
     
     do j=6,nni+3
        if(j.le.nnix) then
           jj=nnix-j+6
        else
           jj=j
        endif
        j2=jj-4
        
        do i=iend1,6,-1
           ii=i
           i2=ii-4
           i0=i0+1
           indx2(i0)=(i2-1)*nni  +j2
           indx1(i0)=(ii-1)*nni8 +jj
        enddo
     enddo
     
     
     ngrd1=i0
     ngrd2=i0
     
     !        bondary points like in the itype=1 case
     
     goto 1000
  endif
  
  
  if (itype.eq.7) then
     
     !        nu-mu sweeps
     
     !        nu changes from pi/2 to pi and then from pi/2 to 0; for every
     !        nu value mu changes from mumax/2 to 1 and from mumax/2 to mumax
     
     i0=0
     nmux=(nmu8-1)/2
     nnix=(nni8-1)/2
     
     do j=6,nni+3
        if(j.le.nnix) then
           jj=nnix-j+6
        else
           jj=j
        endif
        j2=jj-4
        
        do i=6,iend1
           if (i.le.nmux) then
              ii=nmux-i+6
           else
              ii=i
           endif
           i2=ii-4
           
           i0=i0+1
           indx2(i0)=(i2-1)*nni  +j2
           indx1(i0)=(ii-1)*nni8 +jj
        enddo
     enddo
     
     ngrd1=i0
     ngrd2=i0
     
     !        bondary points like in the itype=1 case
     
     goto 1000
     
  endif
  
  
  if (itype.eq.8) then
     
     !        mu changes from mumax to 1; for every mu value nu changes from
     !        pi/2 to pi and then from pi/2 to 0
     
     i0=0
     nmux=(nmu8-1)/2
     nnix=(nni8-1)/2
     
     do i=iend1,6,-1
        ii=i
        i2=ii-4
        do j=6,nni+3
           if(j.le.nnix) then
              jj=nnix-j+6
           else
              jj=j
           endif
           j2=jj-4
           
           i0=i0+1
           indx2(i0)=(i2-1)*nni  +j2
           indx1(i0)=(ii-1)*nni8 +jj
        enddo
     enddo
     
     ngrd1=i0
     ngrd2=i0
     
     !        bondary points like in the itype=1 case
     goto 1000
  endif
      
  if (itype.eq.9) then
     
     !        The following sequence of mesh points is a variant of
     !        Laaksonen's "middle" type of sweeps
     
     !        nu changes from pi/2 to 0 and then from pi/2 to pi.  For every
     !        nu value mu changes from mumax/2 to 1 and from mumax/2 to mumax
     
     i0=0
     nmux=(nmu8-1)/2
     nnix=(nni8-1)/2

     do j=6,nni+3
        if(j.le.nnix) then
           jj=nnix-j+6
        else
           jj=j
        endif
        j2=jj-4
        
        do i=6,iend1
           if (i.le.nmux) then
              ii=nmux-i+6
           else
              ii=i
           endif
           i2=ii-4
           
           i0=i0+1
           indx2(i0)=(i2-1)*nni  +j2
           indx1(i0)=(ii-1)*nni8 +jj
        enddo
     enddo
     
     ngrd1=i0
     ngrd2=i0
     
     !        bondary points like in the itype=1 case
     goto 1000
  endif
  
  
  write(*,*) 'Error: undefined type of ordering'
  stop 'mesh'
  return
  
01000 continue
  
  !     region vi
  !     bondary values from extrapolation
  !     values for points (j,nmut-3),...,(j,nmut) are calculated from
  !     asymptotic expansion (orbitals) or multipole expansion (potentials)
  
  i0=0
  do i=5,iend1
     i0=i0+1
     indx6a(i0)=(i-1)*nni8+5
  enddo
  ngrd6a=i0
  
  i0=0
  do i=5,iend1
     i0=i0+1
     indx6b(i0)=(i-1)*nni8+nni+4
  enddo
  ngrd6b=i0
  
  i0=0
  do j=6,nni+3
     i0=i0+1
     indx7(i0)=4*nni8+j
  enddo
  ngrd7=i0
  
  if (imcase.eq.1) then
     ingr1(1,ig)=ngrd1
     ingr1(2,ig)=ngrd6a
     ingr1(3,ig)=ngrd6b
     ingr1(4,ig)=ngrd7
  elseif (imcase.eq.2) then
     ingr2(1,ig)=ngrd1
     ingr2(2,ig)=ngrd6a
     ingr2(3,ig)=ngrd6b
     ingr2(4,ig)=ngrd7
  endif
  
  
  !     Orbitals and potentials are stored in one-dimensional arrays. When
  !     they are printed as a two-dimensional ones (\nu=0,\mu_1) element
  !     coresponds to the centre A and (\nu=\pi,\mu_1) -- B.
  
  if (iprint(40).ne.0) then
     write(*,'(a8/,(10i5))') 'indx1 ',(indx1 (i),i=1,ngrd1)
     write(*,'(a8/,(10i5))') 'indx2 ',(indx2 (i),i=1,ngrd2)
     write(*,'(a8/,(10i5))') 'indx6a',(indx6a(i),i=1,ngrd6a)
     write(*,'(a8/,(10i5))') 'indx6b',(indx6b(i),i=1,ngrd6b)
     write(*,'(a8/,(10i5))') 'indx7 ',(indx7 (i),i=1,ngrd7 )
     
     call pmtxi(nni-2,mxnmu-4,indx2,ione,ione,incrni,incrmu)
  endif
  
  if (iprint(41).ne.0) then
     write(*,'(40a/,10i5)') 'ngrd1,ngrd2,ngrd6a,ngrd6b,ngrd7',  ngrd1,ngrd2,ngrd6a,ngrd6b,ngrd7
  endif
  
end subroutine mesh


