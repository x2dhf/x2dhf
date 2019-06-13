! ***************************************************************************
! *                                                                         *
! *   Copyright (C) 2019 Jacek Kobus <jkob@fizyka.umk.pl>                   *
! *                                                                         *
! *   This program is free software; you can redistribute it and/or modify  *
! *   it under the terms of the GNU General Public License version 2 as     *
! *   published by the Free Software Foundation.                            *
! *                                                                         *
! ***************************************************************************
! ### initLDA ###

! This routine initializes molecular orbitals as linear combinations of LDA
! orbitals used by S. Lehtola to generate the SAP potentials, see sap routine
! for details.

module initLDA_m
  implicit none
contains
  subroutine initLDA (psi,pot,excp,f2,f4,wgt2,wk0)
    use params
    use discret
    use commons8
    use factor_m
    use flp_m
    use memory
    use norm94_m
    use initPot_m
    use plegendg_m
    use qsort
    use zeroArray_m
    
    implicit none

    integer ::  igp,ihf1,ihf2,ilabel,imu,in,inioff,iord,ishift,mxmax1,mxmax2, &
         nwf1,nwf2,ouf2dhf1,ouf2dhf2
    integer :: i,j,iorb,l1,m1,n1,l2,m2,n2,ns,np,nd,nf
    real (PREC) :: costh1,costh2,psi1,psi2,psi1prev,psi2prev,shn1,shn2,r1t,r2t,rr,xnorm, &
         z
    integer :: zlda1, zlda2,icen

    parameter (iord=3,ouf2dhf1=8,ouf2dhf2=9)

    integer, dimension(:), allocatable :: nhf1,lhf1,nhf2,lhf2

    real (PREC), dimension(*) :: psi,pot,excp,f2,f4,wgt2,wk0
    real (PREC), dimension(:), allocatable :: rhf1,rhf2,phf1t,phf2t
    real (PREC), dimension(:), allocatable :: ehf1,qc1,ehf2,qc2
    real (PREC), dimension(:,:), allocatable :: phf1,phf2
    integer, dimension(maxorb,2) :: lcaomap
    real (PREC), dimension(maxorb) :: eh,f4swap
    real (PREC) :: ehc,ehmax
    integer, dimension(maxorb) :: index,i4swap,tokens1,tokens2
    integer :: i1,i2,i1start,i2start

    ! Where to find the orbitals and the file name to read
    character*512 :: x2dhfdir, orbdir, filename
    logical :: file_exists

    ! Initialization of molecular orbitals
    print *,'... initializing molecular orbitals from LDA functions ...'

    ! Get the location of x2dhf
    call get_environment_variable("X2DHF_DIRECTORY", x2dhfdir)
    ! Orbital directory is
    if (len(x2dhfdir)>0) then
       orbdir = trim(x2dhfdir) // "/lda_orbitals/"
    else
       orbdir = "./"
    end if

    do iorb=1,norb
       do icen=1,2
          lcaomap(iorb,icen)=0
       enddo
    enddo
    nwf1=0
    nwf2=0
    
    ilabel=0

    ! read LDA functions for centre A
    if (z1.ne.zero) then
       zlda1=nint(z1)
       filename = trim(orbdir) // trim(element(zlda1)) // '_orbs.dat'
       ! Check if file exists
       inquire(file=filename, exist=file_exists)
       if(.not. file_exists) then
          print *,'Cannot open file ',trim(filename),'. Did you set X2DHF_DIRECTORY?'
          stop
       end if

       open(ouf2dhf1,file=filename,form='formatted',status='old')

       read(ouf2dhf1,*) mxmax1,nwf1

       ! Allocate memory
       allocate(lhf1(nwf1))
       allocate(nhf1(nwf1))
       allocate(rhf1(mxmax1))
       allocate(phf1t(mxmax1))
       allocate(ehf1(nwf1))
       allocate(qc1(nwf1))
       allocate(phf1(nwf1,mxmax1))

       read(ouf2dhf1,*) (lhf1(i),i=1,nwf1)
       read(ouf2dhf1,*) (qc1(i),i=1,nwf1)
       read(ouf2dhf1,*) (ehf1(i),i=1,nwf1)

       ! Principle quantum numbers of LDA orbitals are not included in data files. Since
       ! orbitals are ordered according to their symmetry (s-type orbitals go first, then
       ! p-type ones, and so on) their principal quantum numbers can be easily assigned.

       ns=0
       np=1
       nd=2
       nf=3
       do i=1,nwf1
          if (lhf1(i)==0) then
             ns=ns+1
             nhf1(i)=ns
          endif
          if (lhf1(i)==1) then
             np=np+1
             nhf1(i)=np
          endif
          if (lhf1(i)==2) then
             nd=nd+1
             nhf1(i)=nd
          endif
          if (lhf1(i)==3) then
             nf=nf+1
             nhf1(i)=nf
          endif
       enddo                       

       read(ouf2dhf1,*) rhf1(1),(phf1(i,1),i=1,nwf1)
       mxmax1=mxmax1-1
       do j=1,mxmax1
          read(ouf2dhf1,*) rhf1(j),(phf1(i,j),i=1,nwf1)
       enddo

       if (iprint(220).ne.0) then
          write(*,*)
          write(*,'(" Orbitals on centre Z1 (Z=",i4,"):")') zlda1
          write(*,'(13x,"n",4x,"l",9x,"e")')
          do i=1,nwf1
             write(*,'(4x,3i5,e16.6)') i,nhf1(i),lhf1(i),ehf1(i)
          enddo
       endif
       
       if (iprint(222).ne.0) then
          do j=1,mxmax1
             write(*,'(20e24.16)') rhf1(j),(phf1(i,j),i=1,nwf1)
          enddo
       endif
    endif
    close(ouf2dhf1)

    ! read LDA functions for centre B

    if (z2.ne.zero) then
       zlda2=nint(z2)
       filename = trim(orbdir) // trim(element(zlda2)) // '_orbs.dat'
       ! Check if file exists
       inquire(file=filename, exist=file_exists)
       if(.not. file_exists) then
          print *,'Cannot open file ',trim(filename),'. Did you set X2DHF_DIRECTORY?'
          stop
       end if

       open(ouf2dhf2,file=filename,form='formatted',status='old')
       read(ouf2dhf2,*) mxmax2,nwf2

       ! Allocate memory
       allocate(lhf2(nwf2))
       allocate(nhf2(nwf2))
       allocate(rhf2(mxmax2))
       allocate(phf2t(mxmax2))
       allocate(ehf2(nwf2))
       allocate(qc2(nwf2))
       allocate(phf2(nwf2,mxmax2))

       read(ouf2dhf2,*) (lhf2(i),i=1,nwf2)
       read(ouf2dhf2,*) (qc2(i),i=1,nwf2)
       read(ouf2dhf2,*) (ehf2(i),i=1,nwf2)

       ns=0
       np=1
       nd=2
       nf=3
       do i=1,nwf2                               
          if (lhf2(i)==0) then
             ns=ns+1
             nhf2(i)=ns
          endif
          if (lhf2(i)==1) then
             np=np+1
             nhf2(i)=np
          endif
          if (lhf2(i)==2) then
             nd=nd+1
             nhf2(i)=nd
          endif
          if (lhf2(i)==3) then
             nf=nf+1
             nhf2(i)=nf
          endif
       enddo

       read(ouf2dhf2,*) rhf2(1),(phf2(i,1),i=1,nwf2)
       mxmax2=mxmax2-1
       do j=1,mxmax2
          read(ouf2dhf2,*) rhf2(j),(phf2(i,j),i=1,nwf2)
       enddo
       
       if (iprint(220).ne.0) then
          write(*,*)
          write(*,'(" Orbitals on centre Z2 (Z=",i4,"):")') zlda2
          write(*,'(13x,"n",4x,"l",9x,"e")')
          do i=1,nwf2
             write(*,'(4x,3i5,e16.6)') i,nhf2(i),lhf2(i),ehf2(i)
          enddo
       endif
       
       if (iprint(222).ne.0) then
          do j=1,mxmax2-1
             write(*,'(20e24.16)') rhf2(j),(phf2(i,j),i=1,nwf2)
          enddo
       endif
       close(ouf2dhf2)
    endif
    
    do iorb=1,norb
       co1(iorb)=zero
       co2(iorb)=zero
       eh(iorb)=-1e-30_PREC
    enddo
 
    if (.not.lcaoIncl) then
       do i1=1,nwf1
          if (lhf1(i1)==0) tokens1(i1)=1
          if (lhf1(i1)==1) tokens1(i1)=1
          if (lhf1(i1)==2) tokens1(i1)=1
          if (lhf1(i1)==3) tokens1(i1)=1     
       enddo

       do i2=1,nwf2
          if (lhf2(i2)==0) tokens2(i2)=1
          if (lhf2(i2)==1) tokens2(i2)=1
          if (lhf2(i2)==2) tokens2(i2)=1
          if (lhf2(i2)==3) tokens2(i2)=1     
       enddo

       iorb=norb
       do while (iorb>=1)
          if (trim(orbsym(iorb))=='sigma') then
             ! looking for orbital at centre A first
             ! sigma orbitals

             i1start=0
             i2start=0
             if (nwf1>0) then
                do i1=1,nwf1
                   if (tokens1(i1)>0.and.ehf1(i1)<eh(iorb)) then
                      !                print *,i1,ehf1(i1)
                      eh(iorb)=ehf1(i1)
                      i1start=i1
                   endif
                enddo
                if (i1start/=0) then
                   eh(iorb)=ehf1(i1start)                   
                   mgx(1,iorb)=nhf1(i1start)
                   mgx(2,iorb)=lhf1(i1start)
                   mgx(3,iorb)=0
                   lcaomap(iorb,1)=i1start
                   co1(iorb)=one
                   co2(iorb)=zero
                   tokens1(i1start)=tokens1(i1start)-1
                endif
             endif
             
             if (nwf2>0) then
                i2start=0
                do i2=1,nwf2
                   if (tokens2(i2)>0.and.ehf2(i2)<eh(iorb)) then
                      eh(iorb)=ehf2(i2)
                      i2start=i2
                   endif
                enddo
                if (i2start/=0) then
                   eh(iorb)=ehf2(i2start)                   
                   mgx(4,iorb)=nhf2(i2start)
                   mgx(5,iorb)=lhf2(i2start)
                   mgx(6,iorb)=0
                   lcaomap(iorb,2)=i2start
                   co1(iorb)=zero
                   co2(iorb)=one
                   if (i1start>1) tokens1(i1start)=tokens1(i1start)+1
                   tokens2(i2start)=tokens2(i2start)-1
                endif
             endif
             
             if (co1(iorb)==one) then
                write(*,'("Asigma iorb,co1,co2,eh,n,l",i4,2f6.1,f12.4,2i5)') &
                     iorb,co1(iorb),co2(iorb),eh(iorb),nhf1(i1start),lhf1(i1start)
             elseif (co2(iorb)==one) then
                write(*,'("Bsigma iorb,co1,co2,eh,n,l",i4,2f6.1,f12.4,2i5)') &
                     iorb,co1(iorb),co2(iorb),eh(iorb),nhf2(i2start),lhf2(i2start)
             else
                write(*,'("?      iorb,co1,co2,eh,n,l",i4,2f6.1,f12.4,2i5)') &
                     !iorb,co1(iorb),co2(iorb),eh(iorb),nhf2(i2start),lhf2(i2start)
                     iorb,co1(iorb),co2(iorb),eh(iorb)
                write (*,*) "      continue ..."                
             endif
          endif
          iorb=iorb-1
       enddo
       
       write(*,*)
    
       do i1=1,nwf1
          if (lhf1(i1)==0) tokens1(i1)=0
          if (lhf1(i1)==1) tokens1(i1)=1
          if (lhf1(i1)==2) tokens1(i1)=1
          if (lhf1(i1)==3) tokens1(i1)=1     
       enddo
       
       do i2=1,nwf2
          if (lhf2(i2)==0) tokens2(i2)=0
          if (lhf2(i2)==1) tokens2(i2)=1
          if (lhf2(i2)==2) tokens2(i2)=1
          if (lhf2(i2)==3) tokens2(i2)=1     
       enddo

       iorb=norb
       do while (iorb>=1)
          if (trim(orbsym(iorb))=='pi') then
             ! looking for orbital at centre A first
             ! sigma orbitals
             ! pi orbitals
             i1start=0
             i2start=0
             if (nwf1>0) then
                do i1=1,nwf1
                   if (tokens1(i1)>0.and.ehf1(i1)<eh(iorb)) then
                      eh(iorb)=ehf1(i1)
                      i1start=i1
                   endif
                enddo
                if (i1start/=0) then
                   eh(iorb)=ehf1(i1start)                   
                   mgx(1,iorb)=nhf1(i1start)
                   mgx(2,iorb)=lhf1(i1start)
                   mgx(3,iorb)=1
                   lcaomap(iorb,1)=i1start
                   co1(iorb)=one
                   co2(iorb)=zero
                   tokens1(i1start)=tokens1(i1start)-1
                endif
             endif

             if (nwf2>0) then
                do i2=1,nwf2
                   if (tokens2(i2)>0.and.ehf2(i2)<eh(iorb)) then
                      eh(iorb)=ehf2(i2)
                      i2start=i2
                   endif
                enddo
                if (i2start/=0) then
                   eh(iorb)=ehf2(i2start)                   
                   mgx(4,iorb)=nhf2(i2start)
                   mgx(5,iorb)=lhf2(i2start)
                   mgx(6,iorb)=1
                   lcaomap(iorb,2)=i2start
                   co1(iorb)=zero
                   co2(iorb)=one
                   if (i1start>0) tokens1(i1start)=tokens1(i1start)+1
                   tokens2(i2start)=tokens2(i2start)-1
                endif
             endif
             
             if (co1(iorb)==one) then
                write(*,'("Api    iorb,co1,co2,eh,n,l",i4,2f6.1,f12.4,2i5)') &
                     iorb,co1(iorb),co2(iorb),eh(iorb),nhf1(i1start),lhf1(i1start)
             elseif (co2(iorb)==one) then
                write(*,'("Bpi    iorb,co1,co2,eh,n,l",i4,2f6.1,f12.4,2i5)') &
                     iorb,co1(iorb),co2(iorb),eh(iorb),nhf2(i2start),lhf2(i2start)
             else
                write(*,'("?pi     iorb,co1,co2,eh,n,l",i4,2f6.1,f12.4,2i5)') &
                     iorb,co1(iorb),co2(iorb),eh(iorb),nhf2(i2start),lhf2(i2start)
                write (*,*) "      continue ..."
             endif
          endif
          iorb=iorb-1
       enddo
       
       write(*,*)
    
       do i1=1,nwf1
          if (lhf1(i1)==0) tokens1(i1)=0
          if (lhf1(i1)==1) tokens1(i1)=0
          if (lhf1(i1)==2) tokens1(i1)=1
          if (lhf1(i1)==3) tokens1(i1)=1     
       enddo
       
       do i2=1,nwf2
          if (lhf2(i2)==0) tokens2(i2)=0
          if (lhf2(i2)==1) tokens2(i2)=0
          if (lhf2(i2)==2) tokens2(i2)=1
          if (lhf2(i2)==3) tokens2(i2)=1     
       enddo

       iorb=norb
       do while (iorb>=1)
          if (trim(orbsym(iorb))=='delta') then
             ! looking for orbital at centre A first
             i1start=0
             i2start=0
             if (nwf1>0) then
                do i1=1,nwf1
                   if (tokens1(i1)>0.and.ehf1(i1)<eh(iorb)) then
                      eh(iorb)=ehf1(i1)
                      i1start=i1
                   endif
                enddo
                if (i1start/=0) then
                   eh(iorb)=ehf1(i1start)                   
                   mgx(1,iorb)=nhf1(i1start)
                   mgx(2,iorb)=lhf1(i1start)
                   mgx(3,iorb)=2
                   lcaomap(iorb,1)=i1start
                   co1(iorb)=one
                   co2(iorb)=zero
                   tokens1(i1start)=tokens1(i1start)-1
                endif
             endif

             if (nwf2>0) then
                do i2=1,nwf2
                   if (tokens2(i2)>0.and.ehf2(i2)<eh(iorb)) then
                      eh(iorb)=ehf2(i2)
                      i2start=i2
                   endif
                enddo
                if (i2start/=0) then
                   eh(iorb)=ehf2(i2start)                   
                   mgx(4,iorb)=nhf2(i2start)
                   mgx(5,iorb)=lhf2(i2start)
                   mgx(6,iorb)=2
                   lcaomap(iorb,2)=i2start
                   co1(iorb)=zero
                   co2(iorb)=one
                   if (i1start>0) tokens1(i1start)=tokens1(i1start)+1
                   tokens2(i2start)=tokens2(i2start)-1
                endif
             endif
             if (co1(iorb)==one) then
                write(*,'("Adelta iorb,co1,co2,eh,n,l",i4,2f6.1,f12.4,2i5)') &
                     iorb,co1(iorb),co2(iorb),eh(iorb),nhf1(i1start),lhf1(i1start)
             elseif (co2(iorb)==one) then
                write(*,'("Bdelta iorb,co1,co2,eh,n,l",i4,2f6.1,f12.4,2i5)') &
                     iorb,co1(iorb),co2(iorb),eh(iorb),nhf2(i2start),lhf2(i2start)
             else
                write(*,'("?delta  iorb,co1,co2,eh,n,l",i4,2f6.1,f12.4,2i5)') &
                     iorb,co1(iorb),co2(iorb),eh(iorb),nhf2(i2start),lhf2(i2start)
                write (*,*) "      continue ..."
             endif
          endif
          iorb=iorb-1
       enddo

       write(*,*)
    
       do i1=1,nwf1
          if (lhf1(i1)==0) tokens1(i1)=0
          if (lhf1(i1)==1) tokens1(i1)=0
          if (lhf1(i1)==2) tokens1(i1)=0
          if (lhf1(i1)==3) tokens1(i1)=1     
       enddo
       
       do i2=1,nwf2
          if (lhf2(i2)==0) tokens2(i2)=0
          if (lhf2(i2)==1) tokens2(i2)=0
          if (lhf2(i2)==2) tokens2(i2)=0
          if (lhf2(i2)==3) tokens2(i2)=1     
       enddo

       iorb=norb
       do while (iorb>=1)
          if (trim(orbsym(iorb))=='phi') then
             ! looking for orbital at centre A first
             i1start=0
             i2start=0
             if (nwf1>0) then
                do i1=1,nwf1
                   if (tokens1(i1)>0.and.ehf1(i1)<eh(iorb)) then
                      eh(iorb)=ehf1(i1)
                      i1start=i1
                   endif
                enddo
                if (i1start/=0) then
                   eh(iorb)=ehf1(i1start)                   
                   mgx(1,iorb)=nhf1(i1start)
                   mgx(2,iorb)=lhf1(i1start)
                   mgx(3,iorb)=3
                   lcaomap(iorb,1)=i1start
                   co1(iorb)=one
                   co2(iorb)=zero
                   tokens1(i1start)=tokens1(i1start)-1
                endif
             endif

             if (nwf2>0) then
                do i2=1,nwf2
                   if (tokens2(i2)>0.and.ehf2(i2)<eh(iorb)) then
                      eh(iorb)=ehf2(i2)
                      i2start=i2
                   endif
                enddo
                if (i2start/=0) then
                   eh(iorb)=ehf2(i2start)                   
                   mgx(4,iorb)=nhf2(i2start)
                   mgx(5,iorb)=lhf2(i2start)
                   mgx(6,iorb)=3
                   lcaomap(iorb,2)=i2start
                   co1(iorb)=zero
                   co2(iorb)=one
                   if (i1start>0) tokens1(i1start)=tokens1(i1start)+1
                   tokens2(i2start)=tokens2(i2start)-1
                endif
             endif
             
             if (co1(iorb)==one) then
                write(*,'("Aphi   iorb,co1,co2,eh,n,l",i4,2f6.1,f12.4,2i5)') &
                     iorb,co1(iorb),co2(iorb),eh(iorb),nhf1(i1start),lhf1(i1start)
             elseif (co2(iorb)==one) then
                write(*,'("Bphi   iorb,co1,co2,eh,n,l",i4,2f6.1,f12.4,2i5)') &
                     iorb,co1(iorb),co2(iorb),eh(iorb),nhf2(i2start),lhf2(i2start)
             else
                write(*,'("?phi    iorb,co1,co2,eh,n,l",i4,2f6.1,f12.4,2i5)') &
                     iorb,co1(iorb),co2(iorb),eh(iorb),nhf2(i2start),lhf2(i2start)
                write (*,*) "      continue ..."
             endif
          endif
          iorb=iorb-1
       enddo
 endif

 ! loop over orbitals
    do iorb=1,norb
       ishift=i1b(iorb)-1
       n1=mgx(1,iorb)
       l1=mgx(2,iorb)
       m1=mgx(3,iorb)

       ! normalization factor for spherical harmonics
       shn1=(-1.0_PREC)**dble(m1)/sqrt(4.0_PREC*pii)*sqrt((2*l1+1)*factor(l1-m1)/factor(l1+m1))
       if (m1.eq.0) shn1=1.0_PREC/sqrt(4.0_PREC*pii)*sqrt((2*l1+1)*factor(l1-m1)/factor(l1+m1))

       if ((co1(iorb).ne.zero) .and. (z1 .ne. zero)) then
          ihf1=0
          do i=1,nwf1
             if (nhf1(i).eq.n1.and.lhf1(i).eq.l1) ihf1=i
          enddo
          if (ihf1.eq.0) then
             print *,"initLDA: no proper atomic orbital for centre Z1 found"
             stop 'initLDA'
          endif

          do j=1,mxmax1
             !phf1t(j)=phf1(ihf1,j)
             phf1t(j)=phf1(lcaomap(iorb,1),j)
          enddo
       endif

       n2=mgx(4,iorb)
       l2=mgx(5,iorb)
       m2=mgx(6,iorb)

       ! normalization factor for spherical harmonics
       shn2=(-1.0_PREC)**dble(m2)/sqrt(4.0_PREC*pii)*sqrt((2*l2+1)*factor(l2-m2)/factor(l2+m2))

       if ((co2(iorb).ne.zero) .and. (z2 .ne. zero)) then
          ihf2=0
          do i=1,nwf2
             if (nhf2(i).eq.n2.and.lhf2(i).eq.l2) ihf2=i
             !if (lhf2(i).eq.l2) ihf2=i
          enddo
          if (ihf2.eq.0) then
             print *,"initLDA: no proper atomic orbital for centre Z2 found"
             stop 'initLDA'
          endif

          do j=1,mxmax2
             !phf2t(j)=phf2(ihf2,j)
             phf2t(j)=phf2(lcaomap(iorb,2),j)
          enddo
       endif

       if (iprint(221).ne.0) then
          print *,'inihf: n1,l1,m1,shn1,ihf1',n1,l1,m1,shn1,ihf1
          print *,'inihf: n2,l2,m2,shn2,ihf2',n2,l2,m2,shn2,ihf2
       endif

       ! loop over grid points
       psi1prev=zero
       psi2prev=zero
       do imu=1,mxnmu
          inioff=(imu-1)*nni
          do in=1,nni
             igp=ishift+inioff+in

             ! for each grid point, e.i. for (vmu(imu),vni(ini))
             ! determine its distance |_r1| and |_r2| from the nuclei A
             ! and B and cosine of the polar angles costh1 and costh2
             ! between z axis and the vectors _r1 and _r2

             psi1=zero
             psi2=zero

             rr=(r/2.0_PREC)*sqrt(vxisq(imu)+vetasq(in)-1.0_PREC)
             z=(r/2.0_PREC)*vxi(imu)*veta(in)
             r1t=(r/2.0_PREC)*(vxi(imu)+veta(in))
             r2t=(r/2.0_PREC)*(vxi(imu)-veta(in))

             ! calculate radial part of the hydrogenic orbital centered
             ! on both the nuclei

             if (r1t.lt.precis) then
                costh1=zero
             else
                costh1=(z+r/2.0_PREC)/r1t
             endif
             !
             if (r2t.lt.precis) then
                costh2=zero
             else
                costh2=-(z-r/2.0_PREC)/r2t
             endif

             ! calculate radial part of the LDA orbital
             if(z1 .ne. zero) then
                if     (co1(iorb).ne.zero) then
                   if (r1t.ge.rhf1(mxmax1)) then
                      psi1=zero
                   elseif (r1t.ge.precis) then
                      psi1=flp(iord,mxmax1,rhf1,phf1t,r1t)*shn1*plegendg(l1,m1,costh1)
                      psi1prev=psi1
                   else
                      psi1=psi1prev
                   endif
                endif
             end if

             if(z2 .ne. zero) then
                if (co2(iorb).ne.zero) then
                   if (r2t.ge.rhf2(mxmax2)) then
                      psi2=zero
                   elseif (r2t.ge.precis) then
                      psi2=flp(iord,mxmax2,rhf2,phf2t,r2t)*shn2*plegendg(l2,m2,costh2)
                      psi2prev=psi2
                   else
                      psi2=psi2prev
                   endif
                endif
             end if
             psi(igp)=co1(iorb)*psi1+co2(iorb)*psi2
          enddo
       enddo
       if (iprint(220).ne.0) then
          if (ilabel.eq.0) then
             ilabel=1
             write (*,*)
             write (*,*) 'Normalization of LCAOs'
          endif
          call norm94 (iorb,psi,f4,wgt2,wk0,xnorm)
          !write (*,1115) iorn(iorb),bond(iorb),gut(iorb),xnorm
          if (lcaomap(iorb,1)/=0) ehc=ehf1(lcaomap(iorb,1))
          if (lcaomap(iorb,2)/=0) ehc=ehf2(lcaomap(iorb,2))

          write (*,1116) iorn(iorb),bond(iorb),gut(iorb),xnorm,&
               lcaomap(iorb,1),lcaomap(iorb,2),ehc,eh(iorb)
1115      format(i4,1x,a8,a1,3x,e22.16,2e16.2)
          !1116      format(i4,1x,a8,a1,3x,e22.16,2i6,e14.4)
1116      format(i4,1x,a8,a1,3x,e22.16,2i6,2f14.4)
       endif
    enddo
    write(*,*)
    
    ! initialize Coulomb and exchange potentials unless OED method is chosen 
    if (imethod/=2) then
       print *,'... setting exchange potential to zero...'
       call zeroArray(length3,excp)
    else
       call initPot(psi,pot,excp,f2,f4,wk0)
    endif
    
  end subroutine initLDA
end module initLDA_m
