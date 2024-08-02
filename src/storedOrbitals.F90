! SPDX-License-Identifier: GPL-2.0-or-later

! Copyright (C) 2019-2024  Jacek Kobus 

module storedOrbitals
  implicit none
contains
  subroutine initStoredOrbs4homo
    use params
    use discrete
    use commons
    use interpolate
    use memory
    use normOrtho
    use sapPotential
    use sharedMemory
    use coulombExchangePotentials
    use utils
    
    implicit none

    integer (KIND=IPREC) ::  igp,ihf1,ihf2,ilabel,imu,in,inioff,iorderl,ishift,mxmax1,mxmax2, &
         nwf1,nwf2,ouf2dhf1,ouf2dhf2
    integer (KIND=IPREC) :: i,j,iorb,l1,m1,n1,l2,m2,n2,ns,np,nd,nf
    real (PREC) :: costh1,costh2,psi1,psi2,psi1prev,psi2prev,shn1,shn2,r1t,r2t,rr,xnorm,z
    integer (KIND=IPREC) :: zlda1, zlda2,icen

    parameter (iorderl=3,ouf2dhf1=8,ouf2dhf2=9)

    integer (KIND=IPREC),dimension(:), allocatable :: nhf1,lhf1,nhf2,lhf2

    !real (PREC), dimension(*) :: psi,pot,excp,f4,wgt2,wk0
    real (PREC), dimension(:), allocatable :: rhf1,rhf2,phf1t,phf2t
    real (PREC), dimension(:), allocatable :: ehf1,qc1,ehf2,qc2,effz1,effz2
    real (PREC), dimension(:,:), allocatable :: phf1,phf2
    integer (KIND=IPREC),dimension(maxorb,2) :: lcaomap
    real (PREC), dimension(maxorb) :: eh
    real (PREC) :: ehc
    integer (KIND=IPREC),dimension(maxorb) :: tokens1,tokens2

    integer (KIND=IPREC) :: i1,i2,icount,i1start,i2start

    ! Where to find the orbitals and the file name to read
    character*512 :: x2dhfdir, orbdir, filename
    logical :: file_exists

    real (PREC), dimension(:), pointer :: f4,psi,wgt2,wk0

    !(psi,pot,excp,f4,wgt2,wk0)
    psi=>orbptr
    f4=>supplptr(i4b(9):)
    wgt2=>supplptr(i4b(14):)
    wk0=>scratchptr

    
    ! Initialization of molecular orbitals
    if (hfIncl) print *,'... initializing molecular orbitals from HF functions ...'
    if (ldaIncl) print *,'... initializing molecular orbitals from LDA functions ...'

    ! Get the location of x2dhf
    call get_environment_variable("X2DHF_DIRECTORY", x2dhfdir)
    ! Orbital directory is
    if (len(x2dhfdir)>0) then
       if (hfIncl)  orbdir = trim(x2dhfdir) // "/hf_orbitals/"
       if (ldaIncl) orbdir = trim(x2dhfdir) // "/lda_orbitals/"
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

    ! read HF/LDA functions for centre A
    if ( abs(z1)>epsilon(zero)) then
       zlda1=nint(z1)
       filename = trim(orbdir) // trim(element(zlda1)) // '_orbs.dat'
       ! Check if file exists
       inquire(file=filename, exist=file_exists)
       if(.not. file_exists) then
          print *,'Cannot open file ',trim(filename),'. Did you set X2DHF_DIRECTORY?'
          stop "initStoreOrbitals4homo"
       end if

       open(ouf2dhf1,file=filename,form='formatted',status='old')

       read(ouf2dhf1,*) mxmax1,nwf1

       ! Allocate memory
       allocate(lhf1(nwf1))
       allocate(nhf1(nwf1))
       allocate(rhf1(mxmax1))
       allocate(phf1t(mxmax1))
       allocate(ehf1(nwf1))
       allocate(effz1(nwf1))
       allocate(qc1(nwf1))
       allocate(phf1(nwf1,mxmax1))

       read(ouf2dhf1,*) (lhf1(i),i=1,nwf1)
       read(ouf2dhf1,*) (qc1(i),i=1,nwf1)
       read(ouf2dhf1,*) (ehf1(i),i=1,nwf1)

       if (hfIncl) then
          read(ouf2dhf1,*) (effz1(i),i=1,nwf1)
       endif

       
       ! Principle quantum numbers of HF orbitals are not included in data files. Since
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

       if (hfIncl) then       
          do i=1,nwf1
             do j=1,mxmax1
                phf1(i,j)=phf1(i,j)/rhf1(j)
             enddo
          enddo
       endif
       
#ifdef PRINT
! print=220: initStoredOrb4homo: Orbitals on centre Z1: i,nhf1(i),lhf1(i),ehf1(i) 
       if (iprint(220).ne.0) then
          write(*,*)
          write(*,'(" Orbitals on centre Z1 (Z=",i4,"):")') zlda1
          write(*,'(13x,"n",4x,"l",9x,"e")')
          do i=1,nwf1
             write(*,'(4x,3i5,e16.6)') i,nhf1(i),lhf1(i),ehf1(i)
          enddo
       endif
#endif
       
#ifdef PRINT
! print=222: initStoredOrb4homo: Orbitals on centre Z1: rhf1(j),(phf1(i,j),i=1,nwf1)
       if (iprint(222).ne.0) then
          do j=1,mxmax1
             write(*,'(20e24.16)') rhf1(j),(phf1(i,j),i=1,nwf1)
          enddo
       endif
#endif
    endif
    close(ouf2dhf1)

    ! read HF/LDA functions for centre B

    if ( abs(z2)>epsilon(zero)) then
       zlda2=nint(z2)
       filename = trim(orbdir) // trim(element(zlda2)) // '_orbs.dat'
       ! Check if file exists
       inquire(file=filename, exist=file_exists)
       if(.not. file_exists) then
          print *,'Cannot open file ',trim(filename),'. Did you set X2DHF_DIRECTORY?'
          stop "initStoredOrbs4homo"
       end if

       open(ouf2dhf2,file=filename,form='formatted',status='old')
       read(ouf2dhf2,*) mxmax2,nwf2

       ! Allocate memory
       allocate(lhf2(nwf2))
       allocate(nhf2(nwf2))
       allocate(rhf2(mxmax2))
       allocate(phf2t(mxmax2))
       allocate(ehf2(nwf2))
       allocate(effz2(nwf2))
       allocate(qc2(nwf2))
       allocate(phf2(nwf2,mxmax2))

       read(ouf2dhf2,*) (lhf2(i),i=1,nwf2)
       read(ouf2dhf2,*) (qc2(i),i=1,nwf2)
       read(ouf2dhf2,*) (ehf2(i),i=1,nwf2)

       if (hfIncl) then
          read(ouf2dhf2,*) (effz2(i),i=1,nwf2)
       endif
       
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

       if (hfIncl) then
          do i=1,nwf2
             do j=1,mxmax2
                phf2(i,j)=phf2(i,j)/rhf2(j)
             enddo
          enddo
       endif

#ifdef PRINT
! print=220: initStoredOrb4homo: Orbitals on centre Z2: i,nhf2(i),lhf2(i),ehf2(i) 
       if (iprint(220).ne.0) then
          write(*,*)
          write(*,'(" Orbitals on centre Z2 (Z=",i4,"):")') zlda2
          write(*,'(13x,"n",4x,"l",9x,"e")')
          do i=1,nwf2
             write(*,'(4x,3i5,e16.6)') i,nhf2(i),lhf2(i),ehf2(i)
          enddo
       endif
#endif

#ifdef PRINT
! print=222: initStoredOrb4homo: Orbitals on centre Z2: rhf2(j),(phf2(i,j),i=1,nwf2)
       if (iprint(222).ne.0) then
          do j=1,mxmax2-1
             write(*,'(20e24.16)') rhf2(j),(phf2(i,j),i=1,nwf2)
          enddo
       endif
#endif
       close(ouf2dhf2)
    endif
    
    do iorb=1,norb
       co1(iorb)=zero
       co2(iorb)=zero
       eh(iorb)=-1e-30_PREC
    enddo

    do i1=1,nwf1
       tokens1(i1)=0
    enddo
    
    do i2=1,nwf2
       tokens2(i2)=0
    enddo
    
    ! lcaoIncl==.false.

    ! atomic orbitals of s, p, d and f symmetry may contribute to sigma-type orbitals     
    do i1=1,nwf1
       if (lhf1(i1)==0) tokens1(i1)=1
       if (lhf1(i1)==1) tokens1(i1)=1
       if (lhf1(i1)==2) tokens1(i1)=1
       if (lhf1(i1)==3) tokens1(i1)=1     
    enddo
    
    icount=0

    ! the orbitals with the lowest energies are treated first
    iorb=norb
    do while (iorb>=1)
       ! the same hf orbital is used for the two consecutive 2dHF orbitals

       if (trim(orbsym(iorb))=='sigma') then
          ! looking for orbital at centre A first; sigma orbitals
          
          i1start=0
          icount=icount+1

          do i1=1,nwf1
             if (tokens1(i1)>0.and.ehf1(i1)<eh(iorb)) then
                eh(iorb)=ehf1(i1)
                i1start=i1
                ! the same hf orbital is used for the two consecutive HF orbitals
                if (mod(icount,2)==0) tokens1(i1)=tokens1(i1)-1
             endif
             !print *,"i1,tokens1(i1),ehf1(i1),icount",i1,tokens1(i1),ehf1(i1),icount
          enddo
          !print *,"iorb,eh(iorb) ",iorb,eh(iorb)
          
          eh(iorb)=ehf1(i1start)                   

          mgx(1,iorb)=nhf1(i1start)
          mgx(2,iorb)=lhf1(i1start)
          mgx(3,iorb)=0
          lcaomap(iorb,1)=i1start
          co1(iorb)=half
          eza1(iorb)=effz1(i1start)
          
          mgx(4,iorb)=nhf2(i1start)
          mgx(5,iorb)=lhf2(i1start)
          mgx(6,iorb)=0
          lcaomap(iorb,2)=i1start
          co2(iorb)=half
          eza2(iorb)=effz2(i1start)
          if (gusym(iorb)=="u") co2(iorb)=-one
          
          if (abs(co1(iorb))>epsilon(zero)) then
             write(*,'(1x,"Asigma iorb,co1,co2,eh,n,l",i4,2f6.1,f12.4,2i5)') &
                  iorb,co1(iorb),co2(iorb),eh(iorb),nhf1(i1start),lhf1(i1start)
          endif
          if (abs(co2(iorb))>epsilon(zero)) then
             write(*,'(1x,"Bsigma iorb,co1,co2,eh,n,l",i4,2f6.1,f12.4,2i5)') &
                  iorb,co1(iorb),co2(iorb),eh(iorb),nhf2(i1start),lhf2(i1start)
          endif
          if (abs(co1(iorb))<epsilon(zero) .and. abs(co2(iorb))<epsilon(zero)) then
             write(*,'(1x,"?      iorb,co1,co2,eh,n,l",i4,2f6.1,f12.4,2i5)') &
                  iorb,co1(iorb),co2(iorb),eh(iorb)
             write (*,*) "      continue ..."                
          endif
       endif
       if (mod(icount,2)==0) icount=0
       iorb=iorb-1
    enddo

    do i1=1,nwf1
       if (lhf1(i1)==0) tokens1(i1)=0
       if (lhf1(i1)==1) tokens1(i1)=1
       if (lhf1(i1)==2) tokens1(i1)=1
       if (lhf1(i1)==3) tokens1(i1)=1     
    enddo
    
    icount=0
    
    iorb=norb
    do while (iorb>=1)
       if (trim(orbsym(iorb))=='pi') then
          i1start=0
          icount=icount+1
          do i1=1,nwf1
             if (tokens1(i1)>0.and.ehf1(i1)<eh(iorb)) then
                eh(iorb)=ehf1(i1)
                i1start=i1
                ! the same hf orbital is used for the two consecutive HF orbitals
                if (mod(icount,2)==0) tokens1(i1)=tokens1(i1)-1
             endif
          enddo
          
          eh(iorb)=ehf1(i1start)                   

          mgx(1,iorb)=nhf1(i1start)
          mgx(2,iorb)=lhf1(i1start)
          mgx(3,iorb)=0
          lcaomap(iorb,1)=i1start
          co1(iorb)=half
          eza1(iorb)=effz1(i1start)
          mgx(4,iorb)=nhf2(i1start)
          mgx(5,iorb)=lhf2(i1start)
          mgx(6,iorb)=0
          lcaomap(iorb,2)=i1start
          co2(iorb)=half
          eza2(iorb)=effz2(i1start)
          if (gusym(iorb)=="u") co2(iorb)=-one
          
          if (abs(co1(iorb))>epsilon(zero)) then
             write(*,'(1x,"Api    iorb,co1,co2,eh,n,l",i4,2f6.1,f12.4,2i5)') &
                  iorb,co1(iorb),co2(iorb),eh(iorb),nhf1(i1start),lhf1(i1start)
          endif
          if (abs(co2(iorb))>epsilon(zero)) then
             write(*,'(1x,"Bpi    iorb,co1,co2,eh,n,l",i4,2f6.1,f12.4,2i5)') &
                  iorb,co1(iorb),co2(iorb),eh(iorb),nhf2(i1start),lhf2(i1start)
          endif
          if (abs(co1(iorb))<epsilon(zero).and.abs(co2(iorb))<epsilon(zero)) then
             write(*,'(1x,"?      iorb,co1,co2,eh,n,l",i4,2f6.1,f12.4,2i5)') &
                  iorb,co1(iorb),co2(iorb),eh(iorb)
             write (*,*) "      continue ..."                
          endif
       endif
       if (mod(icount,2)==0) icount=0
       iorb=iorb-1
    enddo
    
    do i1=1,nwf1
       if (lhf1(i1)==0) tokens1(i1)=0
       if (lhf1(i1)==1) tokens1(i1)=0
       if (lhf1(i1)==2) tokens1(i1)=1
       if (lhf1(i1)==3) tokens1(i1)=1     
    enddo
    
    icount=0
    
    iorb=norb
    do while (iorb>=1)
       if (trim(orbsym(iorb))=='delta') then
          i1start=0
          icount=icount+1

          do i1=1,nwf1
             if (tokens1(i1)>0.and.ehf1(i1)<=eh(iorb)) then
                eh(iorb)=ehf1(i1)
                i1start=i1
                ! the same hf orbital is used for the two consecutive HF orbitals
                if (mod(icount,2)==0) tokens1(i1)=tokens1(i1)-1
             endif
          enddo
          
          eh(iorb)=ehf1(i1start)                   

          mgx(1,iorb)=nhf1(i1start)
          mgx(2,iorb)=lhf1(i1start)
          mgx(3,iorb)=0
          lcaomap(iorb,1)=i1start
          co1(iorb)=half
          eza1(iorb)=effz1(i1start)
          
          mgx(4,iorb)=nhf2(i1start)
          mgx(5,iorb)=lhf2(i1start)
          mgx(6,iorb)=0
          lcaomap(iorb,2)=i1start
          co2(iorb)=half
          eza2(iorb)=effz2(i1start)
          if (gusym(iorb)=="u") co2(iorb)=-one
          
          if (abs(co1(iorb))>epsilon(zero)) then
             write(*,'(1x,"Adelta iorb,co1,co2,eh,n,l",i4,2f6.1,f12.4,2i5)') &
                  iorb,co1(iorb),co2(iorb),eh(iorb),nhf1(i1start),lhf1(i1start)
          endif
          if (abs(co2(iorb))>epsilon(zero)) then          
             write(*,'(1x,"Bdelta iorb,co1,co2,eh,n,l",i4,2f6.1,f12.4,2i5)') &
                  iorb,co1(iorb),co2(iorb),eh(iorb),nhf2(i1start),lhf2(i1start)
          endif
          if (abs(co1(iorb))<epsilon(zero).and.abs(co2(iorb))<epsilon(zero)) then
             write(*,'(1x,"?      iorb,co1,co2,eh,n,l",i4,2f6.1,f12.4,2i5)') &
                  iorb,co1(iorb),co2(iorb),eh(iorb)
             write (*,*) "      continue ..."                
          endif
       endif
       if (mod(icount,2)==0) icount=0
       iorb=iorb-1
    enddo

    do i1=1,nwf1
       if (lhf1(i1)==0) tokens1(i1)=0
       if (lhf1(i1)==1) tokens1(i1)=0
       if (lhf1(i1)==2) tokens1(i1)=0
       if (lhf1(i1)==3) tokens1(i1)=1     
    enddo
    
    icount=0

    iorb=norb
    do while (iorb>=1)
       if (trim(orbsym(iorb))=='phi') then
          i1start=0
          icount=icount+1

          do i1=1,nwf1
             if (tokens1(i1)>0.and.ehf1(i1)<=eh(iorb)) then
                eh(iorb)=ehf1(i1)
                i1start=i1
                ! the same hf orbital is used for the two consecutive HF orbitals
                if (mod(icount,2)==0) tokens1(i1)=tokens1(i1)-1
             endif
          enddo
          
          eh(iorb)=ehf1(i1start)                   

          mgx(1,iorb)=nhf1(i1start)
          mgx(2,iorb)=lhf1(i1start)
          mgx(3,iorb)=0
          lcaomap(iorb,1)=i1start
          co1(iorb)=half
          eza1(iorb)=effz1(i1start)
          
          mgx(4,iorb)=nhf2(i1start)
          mgx(5,iorb)=lhf2(i1start)
          mgx(6,iorb)=0
          lcaomap(iorb,2)=i1start
          co2(iorb)=half
          eza2(iorb)=effz2(i1start)
          if (gusym(iorb)=="u") co2(iorb)=-one
          
          if (abs(co1(iorb))>epsilon(zero)) then
             write(*,'(1x,"Aphi   iorb,co1,co2,eh,n,l",i4,2f6.1,f12.4,2i5)') &
                  iorb,co1(iorb),co2(iorb),eh(iorb),nhf1(i1start),lhf1(i1start)
          endif
          if (abs(co2(iorb))>epsilon(zero)) then
             write(*,'(1x,"Bphi   iorb,co1,co2,eh,n,l",i4,2f6.1,f12.4,2i5)') &
                  iorb,co1(iorb),co2(iorb),eh(iorb),nhf2(i1start),lhf2(i1start)
          endif
          if (abs(co1(iorb))<epsilon(zero).and.abs(co2(iorb))<epsilon(zero)) then
             write(*,'(1x,"?      iorb,co1,co2,eh,n,l",i4,2f6.1,f12.4,2i5)') &
                  iorb,co1(iorb),co2(iorb),eh(iorb)
             write (*,*) "      continue ..."                
          endif
       endif
       if (mod(icount,2)==0) icount=0
       iorb=iorb-1
    enddo

 ! loop over orbitals
    do iorb=1,norb
       ishift=i1b(iorb)-1
       n1=mgx(1,iorb)
       l1=mgx(2,iorb)
       m1=mgx(3,iorb)

       ! normalization factor for spherical harmonics
       shn1=(-1.0_PREC)**dble(m1)/sqrt(4.0_PREC*pii)*&
            sqrt((2.0_PREC*l1+1.0_PREC)*factor(l1-m1)/factor(l1+m1))
       if (m1.eq.0) shn1=1.0_PREC/sqrt(4.0_PREC*pii)*&
            sqrt((2.0_PREC*l1+1.0_PREC)*factor(l1-m1)/factor(l1+m1))

       if (abs(co1(iorb))>epsilon(zero) .and. abs(z1)>epsilon(zero) ) then
          ihf1=0
          do i=1,nwf1
             !print *, "i,nhf1(i),n1,lhf1(i),l1 ",i,nhf1(i),n1,lhf1(i),l1
             if (nhf1(i).eq.n1.and.lhf1(i).eq.l1) ihf1=i
          enddo
          if (ihf1.eq.0) then
             print *,"initStoredOrbs4homo: no proper atomic orbital for centre Z1 found"
             stop 'initStoredOrbs4homo'
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
       shn2=(-1.0_PREC)**dble(m2)/sqrt(4.0_PREC*pii)*&
            sqrt((2.0_PREC*l2+1.0_PREC)*factor(l2-m2)/factor(l2+m2))

       if (abs(co2(iorb))>epsilon(zero) .and. abs(z2)>epsilon(zero) ) then       
          ihf2=0
          do i=1,nwf2
             if (nhf2(i).eq.n2.and.lhf2(i).eq.l2) ihf2=i
             !if (lhf2(i).eq.l2) ihf2=i
          enddo
          if (ihf2.eq.0) then
             print *,"initStoredOrbs4homo: no proper atomic orbital for centre Z2 found"
             stop 'initStoredOrbs4homo'
          endif

          do j=1,mxmax2
             !phf2t(j)=phf2(ihf2,j)
             phf2t(j)=phf2(lcaomap(iorb,2),j)
          enddo
       endif

#ifdef PRINT
! print=221: initStoredOrb4homo: inihf: n1,l1,m1,shn1,ihf1',n1,l1,m1,shn1,ihf1
       if (iprint(221).ne.0) then
          print *,'inihf: n1,l1,m1,shn1,ihf1',n1,l1,m1,shn1,ihf1
          print *,'inihf: n2,l2,m2,shn2,ihf2',n2,l2,m2,shn2,ihf2
       endif
#endif       

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

             ! calculate radial part of the HF orbital
             if(z1 .ne. zero) then
                if     (abs(co1(iorb))>epsilon(zero)) then
                   if (r1t.ge.rhf1(mxmax1)) then
                      psi1=zero
                   elseif (r1t.ge.precis) then
                      psi1=flp(iorderl,mxmax1,rhf1,phf1t,r1t)*shn1*plegendg(l1,m1,costh1)
                      psi1prev=psi1
                   else
                      psi1=psi1prev
                   endif
                endif
             end if

             if(abs(z2)>epsilon(zero)) then
                if (abs(co2(iorb))>epsilon(zero)) then
                   if (r2t.ge.rhf2(mxmax2)) then
                      psi2=zero
                   elseif (r2t.ge.precis) then
                      psi2=flp(iorderl,mxmax2,rhf2,phf2t,r2t)*shn2*plegendg(l2,m2,costh2)
                      psi2prev=psi2
                   else
                      psi2=psi2prev
                   endif
                endif
             end if
             psi(igp)=co1lda*psi1+co2lda*psi2
          enddo
       enddo

       if (iprint(220).ne.0) then
          if (ilabel.eq.0) then
             ilabel=1
             write (*,*)
             write (*,*) 'Normalization of LCAOs'
          endif
          call norm94 (iorb,xnorm)
          !write (*,1115) iorn(iorb),bond(iorb),gusym(iorb),xnorm
          if (lcaomap(iorb,1)/=0) ehc=ehf1(lcaomap(iorb,1))
          if (lcaomap(iorb,2)/=0) ehc=ehf2(lcaomap(iorb,2))

          write (*,1116) iorn(iorb),bond(iorb),gusym(iorb),xnorm,&
               lcaomap(iorb,1),lcaomap(iorb,2),ehc,eh(iorb)
1115      format(i4,1x,a8,a1,3x,e22.16,2e16.2)
          !1116      format(i4,1x,a8,a1,3x,e22.16,2i6,e14.4)
1116      format(i4,1x,a8,a1,3x,e22.16,2i6,2f14.4)
       endif
    enddo
    write(*,*)

    ! initialize Coulomb and exchange potentials
    if (ldaIncl) then
       if (ldaSAPIncl) then
          call initCoulombSAP
          call initExchange
          return
       endif
    endif
    call initCoulomb
    call initExchange

  end subroutine initStoredOrbs4homo

  subroutine initStoredOrbs 
    use params
    use sharedMemory
    use blas
    use discrete
    use commons
    use interpolate
    use memory
    use normOrtho
    use integrals
    use sapPotential
    use coulombExchangePotentials
    use utils
    
    implicit none

    integer (KIND=IPREC) ::  igp,ihf1,ihf2,ilabel,imu,in,inioff,iorderl,ishift,mxmax1,mxmax2, &
         nwf1,nwf2,ouf2dhf1,ouf2dhf2
    integer (KIND=IPREC) :: i,j,iorb,l1,m1,n1,l2,m2,n2,ns,np,nd,nf,nwf1n
    real (PREC) :: co12,costh1,costh2,psi1,psi2,psi1prev,psi2prev,shn1,shn2,r1t,r2t,rr,xnorm,z
    integer (KIND=IPREC) :: zlda1, zlda2,icen

    real (PREC) :: sf4effz
    
    parameter (iorderl=3,ouf2dhf1=8,ouf2dhf2=9)
    parameter (sf4effz=1.0e0_PREC)

    integer (KIND=IPREC),dimension(:), allocatable :: nhf1,lhf1,nhf2,lhf2

    real (PREC), dimension(:), allocatable :: rhf1,rhf2,phf1t,phf2t
    real (PREC), dimension(:), allocatable :: ehf1,qc1,ehf2,qc2,effz1,effz2
    real (PREC), dimension(:,:), allocatable :: phf1,phf2

    real (PREC), dimension(:), allocatable :: psi4sort
    
    integer (KIND=IPREC),dimension(maxorb,2) :: lcaomap
    real (PREC), dimension(maxorb) :: eh
    real (PREC) :: ehc
    real (PREC) :: tt    
    integer (KIND=IPREC),dimension(maxorb) :: tokens1,tokens2

    real (PREC), dimension(maxorb) :: vt4ord
    integer (KIND=IPREC),dimension(maxorb) :: ind
    
    integer (KIND=IPREC) :: i1,i2,i1start,i2start

    ! Where to find the orbitals and the file name to read
    character*512 :: x2dhfdir, orbdir, filename
    logical :: file_exists

    real (PREC), dimension(:), pointer :: f4,psi,wgt2,wk0

    psi=>orbptr
    f4=>supplptr(i4b(9):)
    wgt2=>supplptr(i4b(14):)
    wk0=>scratchptr

    ! Initialization of molecular orbitals
    if (hfIncl) print *,'... initializing molecular orbitals from HF functions ...'
    if (ldaIncl) print *,'... initializing molecular orbitals from LDA functions ...'

    ! Get the location of x2dhf
    call get_environment_variable("X2DHF_DIRECTORY", x2dhfdir)
    ! Orbital directory is
    ! Orbital directory is
    if (len(x2dhfdir)>0) then
       if (hfIncl)  orbdir = trim(x2dhfdir) // "/hf_orbitals/"
       if (ldaIncl) orbdir = trim(x2dhfdir) // "/lda_orbitals/"
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

    zlda1=nint(z1)       
    zlda2=nint(z2)

    ! read HF/LDA functions for centre A
    if ( abs(z1)>epsilon(zero)) then
       filename = trim(orbdir) // trim(element(zlda1)) // '_orbs.dat'
       
       ! Check if file exists
       inquire(file=filename, exist=file_exists)
       if(.not. file_exists) then
          print *,'Cannot open file ',trim(filename),'. Did you set X2DHF_DIRECTORY?'
          stop "initStoreOrbitals"
       end if

       open(ouf2dhf1,file=filename,form='formatted',status='old')

       read(ouf2dhf1,*) mxmax1,nwf1

       ! Allocate memory
       allocate(lhf1(nwf1))
       allocate(nhf1(nwf1))
       allocate(rhf1(mxmax1))
       allocate(phf1t(mxmax1))
       allocate(ehf1(nwf1))
       allocate(effz1(nwf1))
       allocate(qc1(nwf1))
       allocate(phf1(nwf1,mxmax1))

       read(ouf2dhf1,*) (lhf1(i),i=1,nwf1)
       read(ouf2dhf1,*) (qc1(i),i=1,nwf1)
       read(ouf2dhf1,*) (ehf1(i),i=1,nwf1)

       if (hfIncl) then
          read(ouf2dhf1,*) (effz1(i),i=1,nwf1)
          do i=1,nwf1
             effz1(i)=sf4effz*effz1(i)
          enddo
       else
          do i=1,nwf1
             effz1(i)=z1
          enddo
       endif
       
       ! Principle quantum numbers of HF orbitals are not included in data files. Since
       ! orbitals are ordered according to their symmetry (s-type orbitals go first, then
       ! p-type ones, and so on) their principal quantum numbers can easily be assigned.

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
       if (hfIncl) then
          do i=1,nwf1
             do j=1,mxmax1
                phf1(i,j)=phf1(i,j)/rhf1(j)
             enddo
          enddo
       endif
       close(ouf2dhf1)
#ifdef PRINT
! print=220: initStoredOrb: Orbitals on centre Z1: i,nhf1(i),lhf1(i), ehf1(i)       
       if (iprint(220).ne.0) then
          write(*,*)
          write(*,'(" Orbitals on centre Z1 (Z=",i4,"):")') zlda1
          write(*,'(13x,"n",4x,"l",9x,"e")')
          do i=1,nwf1
             write(*,'(4x,3i5,e16.6)') i,nhf1(i),lhf1(i),ehf1(i)
          enddo
       endif
#endif

#ifdef PRINT      
! print=222: initStoredOrb: Orbitals on centre Z1: rhf1(j),(phf1(i,j),i=1,nwf1)       
       if (iprint(222).ne.0) then
          do j=1,mxmax1
             write(*,'(20e24.16)') rhf1(j),(phf1(i,j),i=1,nwf1)
          enddo
       endif
#endif
    endif

    ! read HF/LDA functions for centre B

    if (abs(z2)>epsilon(zero)) then
       filename = trim(orbdir) // trim(element(zlda2)) // '_orbs.dat'
       ! Check if file exists
       inquire(file=filename, exist=file_exists)
       if(.not. file_exists) then
          print *,'Cannot open file ',trim(filename),'. Did you set X2DHF_DIRECTORY?'
          stop "initStoreOrbitals"
       end if

       open(ouf2dhf2,file=filename,form='formatted',status='old')
       read(ouf2dhf2,*) mxmax2,nwf2

       ! Allocate memory
       allocate(lhf2(nwf2))
       allocate(nhf2(nwf2))
       allocate(rhf2(mxmax2))
       allocate(phf2t(mxmax2))
       allocate(ehf2(nwf2))
       allocate(effz2(nwf2))
       allocate(qc2(nwf2))
       allocate(phf2(nwf2,mxmax2))

       read(ouf2dhf2,*) (lhf2(i),i=1,nwf2)
       read(ouf2dhf2,*) (qc2(i),i=1,nwf2)
       read(ouf2dhf2,*) (ehf2(i),i=1,nwf2)

       if (hfIncl) then
          read(ouf2dhf2,*) (effz2(i),i=1,nwf2)
          do i=1,nwf2
             effz2(i)=sf4effz*effz2(i)
          enddo
       else
          do i=1,nwf2
             effz2(i)=z2
          enddo
       endif
       
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

       if (hfIncl) then
          do i=1,nwf2
             do j=1,mxmax2
                phf2(i,j)=phf2(i,j)/rhf2(j)
             enddo
          enddo
       endif
       
#ifdef PRINT
! print=220: initStoredOrb: Orbitals on centre Z2: i,nhf2(i),lhf2(i), ehf2(i)        
       if (iprint(220).ne.0) then
          write(*,*)
          write(*,'(" Orbitals on centre Z2 (Z=",i4,"):")') zlda2
          write(*,'(13x,"n",4x,"l",9x,"e")')
          do i=1,nwf2
             write(*,'(4x,3i5,e16.6)') i,nhf2(i),lhf2(i),ehf2(i)
          enddo
       endif
#endif
#ifdef PRINT      
! print=222: initStoredOrb: Orbitals on centre Z2: rhf2(j),(phf2(i,j),i=1,nwf2)        
       if (iprint(222).ne.0) then
          do j=1,mxmax2-1
             write(*,'(20e24.16)') rhf2(j),(phf2(i,j),i=1,nwf2)
          enddo
       endif
#endif       
       close(ouf2dhf2)
    endif
    
    do iorb=1,norb
       co1(iorb)=zero
       co2(iorb)=zero
       eh(iorb)=-1e-30_PREC
    enddo

    if (.not.lcaoIncl) then
       ! atomic orbitals of any symmetry may contribute to sigma-type orbitals 
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
                   eza1(iorb)=effz1(i1start)
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
                   eza2(iorb)=effz2(i2start)
                   if (i1start>1) tokens1(i1start)=tokens1(i1start)+1
                   tokens2(i2start)=tokens2(i2start)-1
                endif
             endif

             if (abs(co1(iorb))>epsilon(zero)) then
                write(*,'(1x,"Asigma iorb,co1,co2,eh,n,l",i4,2f6.1,f12.4,2i5)') &
                     iorb,co1(iorb),co2(iorb),eh(iorb),nhf1(i1start),lhf1(i1start)
             endif
             if (abs(co2(iorb))>epsilon(zero)) then
                write(*,'(1x,"Bsigma iorb,co1,co2,eh,n,l",i4,2f6.1,f12.4,2i5)') &
                     iorb,co1(iorb),co2(iorb),eh(iorb),nhf2(i2start),lhf2(i2start)
             endif
             if (abs(co1(iorb))<epsilon(zero) .and.abs(co2(iorb))<epsilon(zero)) then
                write(*,'(1x,"?      iorb,co1,co2,eh,n,l",i4,2f6.1,f12.4,2i5)') &
                     iorb,co1(iorb),co2(iorb),eh(iorb)
                write (*,*) "      continue ..."                
             endif
          endif
          iorb=iorb-1
       enddo

       ! atomic orbitals of any p, d, f symmetry may contribute to pi-type orbitals        
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
                   if (tokens1(i1)>0.and.ehf1(i1)<=eh(iorb)) then
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
                   eza1(iorb)=effz1(i1start)
                   tokens1(i1start)=tokens1(i1start)-1
                endif
             endif

             if (nwf2>0) then
                do i2=1,nwf2
                   if (tokens2(i2)>0.and.ehf2(i2)<=eh(iorb)) then
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
                   co2(iorb)=one
                   co1(iorb)=zero
                   if (i1start>0) tokens1(i1start)=tokens1(i1start)+1
                   tokens2(i2start)=tokens2(i2start)-1
                endif
             endif
             
             if (abs(co1(iorb))>epsilon(zero)) then
                write(*,'(1x,"Api    iorb,co1,co2,eh,n,l",i4,2f6.1,f12.4,2i5)') &
                     iorb,co1(iorb),co2(iorb),eh(iorb),nhf1(i1start),lhf1(i1start)
             endif
             if (abs(co2(iorb))>epsilon(zero)) then
                write(*,'(1x,"Bpi    iorb,co1,co2,eh,n,l",i4,2f6.1,f12.4,2i5)') &
                     iorb,co1(iorb),co2(iorb),eh(iorb),nhf2(i2start),lhf2(i2start)
             endif
             if (abs(co1(iorb))<epsilon(zero).and.abs(co2(iorb))<epsilon(zero)) then
                write(*,'(1x,"?pi     iorb,co1,co2,eh,n,l",i4,2f6.1,f12.4,2i5)') &
                     iorb,co1(iorb),co2(iorb),eh(iorb),nhf2(i2start),lhf2(i2start)
                write (*,*) "      continue ..."
             endif
          endif
          iorb=iorb-1
       enddo

       ! atomic orbitals of any d and f symmetry may contribute to delta-type orbitals               
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
                   if (tokens1(i1)>0.and.ehf1(i1)<=eh(iorb)) then
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
                   tokens1(i1start)=tokens1(i1start)-1
                endif
             endif

             if (nwf2>0) then
                do i2=1,nwf2
                   if (tokens2(i2)>0.and.ehf2(i2)<=eh(iorb)) then
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
                   co2(iorb)=one
                   co1(iorb)=zero
                   eza2(iorb)=effz2(i2start)
                   if (i1start>0) tokens1(i1start)=tokens1(i1start)+1
                   tokens2(i2start)=tokens2(i2start)-1
                endif
             endif
             if (abs(co1(iorb))>epsilon(zero)) then
                write(*,'(1x,"Adelta iorb,co1,co2,eh,n,l",i4,2f6.1,f12.4,2i5)') &
                     iorb,co1(iorb),co2(iorb),eh(iorb),nhf1(i1start),lhf1(i1start)
             endif
             if (abs(co2(iorb))>epsilon(zero)) then
                write(*,'(1x,"Bdelta iorb,co1,co2,eh,n,l",i4,2f6.1,f12.4,2i5)') &
                     iorb,co1(iorb),co2(iorb),eh(iorb),nhf2(i2start),lhf2(i2start)
             endif
             if (abs(co1(iorb))<epsilon(zero).and.abs(co2(iorb))<epsilon(zero)) then
                write(*,'(1x,"?delta  iorb,co1,co2,eh,n,l",i4,2f6.1,f12.4,2i5)') &
                     iorb,co1(iorb),co2(iorb),eh(iorb),nhf2(i2start),lhf2(i2start)
                write (*,*) "      continue ..."
             endif
          endif
          iorb=iorb-1
       enddo

       ! atomic orbitals of any f symmetry may contribute to phi-type orbitals                      
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
                   if (tokens1(i1)>0.and.ehf1(i1)<=eh(iorb)) then
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
                   eza1(iorb)=effz1(i1start)
                   tokens1(i1start)=tokens1(i1start)-1
                endif
             endif

             if (nwf2>0) then
                do i2=1,nwf2
                   if (tokens2(i2)>0.and.ehf2(i2)<=eh(iorb)) then
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
                   co2(iorb)=one
                   co1(iorb)=zero
                   eza2(iorb)=effz2(i1start)
                   if (i1start>0) tokens1(i1start)=tokens1(i1start)+1
                   tokens2(i2start)=tokens2(i2start)-1
                endif
             endif
             
             if (abs(co1(iorb))>epsilon(zero)) then
                write(*,'(1x,"Aphi   iorb,co1,co2,eh,n,l",i4,2f6.1,f12.4,2i5)') &
                     iorb,co1(iorb),co2(iorb),eh(iorb),nhf1(i1start),lhf1(i1start)
             endif
             if (abs(co2(iorb))>epsilon(zero)) then
                write(*,'(1x,"Bphi   iorb,co1,co2,eh,n,l",i4,2f6.1,f12.4,2i5)') &
                     iorb,co1(iorb),co2(iorb),eh(iorb),nhf2(i2start),lhf2(i2start)
             endif
             if (abs(co1(iorb))<epsilon(zero).and.abs(co2(iorb))<epsilon(zero)) then
                write(*,'(1x,"?phi    iorb,co1,co2,eh,n,l",i4,2f6.1,f12.4,2i5)') &
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
       shn1=(-1.0_PREC)**dble(m1)/sqrt(4.0_PREC*pii)*&
            sqrt((2.0_PREC*l1+1.0_PREC)*factor(l1-m1)/factor(l1+m1))
       if (m1.eq.0) shn1=1.0_PREC/sqrt(4.0_PREC*pii)*&
            sqrt((2.0_PREC*l1+1.0_PREC)*factor(l1-m1)/factor(l1+m1))

       if ( abs(co1(iorb))>epsilon(zero) .and. abs(z1)>epsilon(zero) ) then
          ihf1=0
          do i=1,nwf1
             if (nhf1(i).eq.n1.and.lhf1(i).eq.l1) ihf1=i
          enddo
          if (ihf1.eq.0) then
             print *,"initStoredOrbs: no proper atomic orbital for centre Z1 found"
             stop 'initStoredOrbs'
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
       shn2=(-1.0_PREC)**dble(m2)/sqrt(4.0_PREC*pii)*&
            sqrt((2.0_PREC*l2+1.0_PREC)*factor(l2-m2)/factor(l2+m2))

       if ( abs(co2(iorb))>epsilon(zero) .and. abs(z2)>epsilon(zero) ) then
          ihf2=0
          do i=1,nwf2
             if (nhf2(i).eq.n2.and.lhf2(i).eq.l2) ihf2=i
             !if (lhf2(i).eq.l2) ihf2=i
          enddo
          if (ihf2.eq.0) then
             print *,"initStoredOrbs: no proper atomic orbital for centre Z2 found"
             stop 'initStoredOrbs'
          endif

          do j=1,mxmax2
             !phf2t(j)=phf2(ihf2,j)
             phf2t(j)=phf2(lcaomap(iorb,2),j)
          enddo
       endif

#ifdef PRINT       
! print=221: initStoredOrb: n1,l1,m1,shn1,ihf1,n2,l2,m2,shn2,ihf2
       if (iprint(221).ne.0) then
          print *,'inihf: n1,l1,m1,shn1,ihf1',n1,l1,m1,shn1,ihf1
          print *,'inihf: n2,l2,m2,shn2,ihf2',n2,l2,m2,shn2,ihf2
       endif
#endif
       
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

             ! calculate radial part of the HF orbital
             if(abs(z1)>epsilon(zero)) then
                if     (abs(co1(iorb))>epsilon(zero)) then
                   if (r1t.ge.rhf1(mxmax1)) then
                      psi1=zero
                   elseif (r1t.ge.precis) then
                      psi1=flp(iorderl,mxmax1,rhf1,phf1t,r1t)*shn1*plegendg(l1,m1,costh1)
                      psi1prev=psi1
                   else
                      psi1=psi1prev
                   endif
                endif
             end if

             if(abs(z2)>epsilon(zero)) then
                if (abs(co2(iorb))>epsilon(zero)) then
                   if (r2t.ge.rhf2(mxmax2)) then
                      psi2=zero
                   elseif (r2t.ge.precis) then
                      psi2=flp(iorderl,mxmax2,rhf2,phf2t,r2t)*shn2*plegendg(l2,m2,costh2)
                      psi2prev=psi2
                   else
                      psi2=psi2prev
                   endif
                endif
             end if
             psi(igp)=co1lda*psi1+co2lda*psi2
          enddo
       enddo

#ifdef PRINT
! print=223: initStoredOrb: normalization of LCAOs
       if (iprint(223).ne.0) then
          if (ilabel.eq.0) then
             ilabel=1
             write (*,*)
             write (*,*) 'Normalization of LCAOs'
          endif
          call norm94 (iorb,xnorm)
          !write (*,1115) iorn(iorb),bond(iorb),gusym(iorb),xnorm
          if (lcaomap(iorb,1)/=0) ehc=ehf1(lcaomap(iorb,1))
          if (lcaomap(iorb,2)/=0) ehc=ehf2(lcaomap(iorb,2))

          write (*,1116) iorn(iorb),bond(iorb),gusym(iorb),xnorm,&
               lcaomap(iorb,1),lcaomap(iorb,2),ehc,eh(iorb)
1115      format(i4,1x,a8,a1,3x,e22.16,2e16.2)
          !1116      format(i4,1x,a8,a1,3x,e22.16,2i6,e14.4)
1116      format(i4,1x,a8,a1,3x,e22.16,2i6,2f14.4)
       endif
#endif
    enddo
    write(*,*)


    ! normalize mixing coefficients
    do iorb=1,norb
       co12=abs(co1(iorb))+abs(co2(iorb))
       if (abs(co12)>epsilon(zero)) then
          co1(iorb)=co1(iorb)/co12
          co2(iorb)=co2(iorb)/co12
       endif
    enddo

    ! initialize Coulomb and exchange potentials
    if (ldaIncl) then
       if (ldaSAPIncl) then
          call initCoulombSAP
          call initExchange
          return
       endif
    endif
    call initCoulomb
    call initExchange
  end subroutine initStoredOrbs

end module storedOrbitals
