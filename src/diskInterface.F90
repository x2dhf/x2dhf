! SPDX-License-Identifier: GPL-2.0-or-later

! Copyright (C) 1996       Leif Laaksonen, Dage Sundholm               
! Copyright (C) 1996-2024  Jacek Kobus 

module diskInterface
  implicit none
contains
  ! ### initDisk ###
  !
  !     Controls retrieving orbitals and potentials from a disk file
  !

  ! ### readHeader ###
  !
   subroutine readHeader(norb_p)
    use params
    use discrete
    use scfshr
    use commons

    implicit none
    character*15 version_p
    character*80 :: datetime_p,header_p
    character*1, dimension(80) :: title
    equivalence (header_p,title(1))
    integer (KIND=IPREC) :: i,ig,imax1,imax2,mismatch,nel_p,nexch_p,norb_p,err
    real (PREC) :: rd,z1d,z2d
    integer i4tmp1,i4tmp2,i4tmp3
    integer (KIND=IPREC),dimension(1) :: i4tmp
    integer (KIND=IPREC),dimension(10) :: nmu_comp

    integer*8 i8tmp1,i8tmp2,i8tmp3
    integer*8, dimension(1) :: i8tmp

    real (PREC) r8tmp1,r8tmp2
    real (PREC), dimension(1) :: r8tmp
    real (PREC), dimension(10) :: rgrid_comp

    real (PREC16) :: r16tmp,r16tmp1,r16tmp12,r16tmp2

    mismatch=0
    initAddData=.false.
    imax1=34
    imax2=24

    ! 2dhf_input.dat file is NOT available
    if (.not.ldatPresent) then
       
       ! check default length of integer variables used in a file
       ! intlen=checkintlen()
    
       read(iinp11,err=1000) header_p
       read(iinp11,err=1000) datetime_p
       
       write(iout6,1050)
       write(iout6,'(  6x,"   timestamp: ",a80)') datetime_p
       write(iout6,'(  6x,"      header: ",a80)') header_p
       
       if (lengthint.eq.4) then
          read(iinp11,err=1000) i4tmp1,i4tmp2,i4tmp
          ngrids_p=i4tmp1
          nni_p   =i4tmp2
          nmu_p=i4tmp(1)
       else
          read(iinp11,err=1000) i8tmp1,i8tmp2,i8tmp
          ngrids_p=int(i8tmp1)
          nni_p   =int(i8tmp2)
          nmu_p=int(i8tmp(1))
       endif
       
       if (lengthfp.eq.8) then
          read(iinp11,err=1000) r8tmp1,r8tmp
          r_p=r8tmp1
          rgrid_p=r8tmp(1)
       else
          read(iinp11,err=1000) r16tmp1,r16tmp
          r_p=r16tmp1
          rgrid_p=r16tmp
       endif
       
       if (lengthfp.eq.8) then
          read(iinp11,err=1000) r8tmp1,r8tmp2
          z1_p=r8tmp1
          z2_p=r8tmp2
       else
          read(iinp11,err=1000) r16tmp1,r16tmp2
          z1_p=r16tmp1
          z2_p=r16tmp2
       endif
       
       if (lengthint.eq.4) then
          read(iinp11,err=1000) i4tmp1,i4tmp2,i4tmp3
          norb_p =i4tmp1
          nel_p  =i4tmp2
          nexch_p=i4tmp3
          nexch  =nexch_p
       else
          read(iinp11,err=1000) i8tmp1,i8tmp2,i8tmp3
          norb_p =int(i8tmp1)
          nel_p  =int(i8tmp2)
          nexch_p=int(i8tmp3)
          nexch  =nexch_p
       endif
    endif

    ! 2dhf_input.dat file is available in at least 3.0 version
    rewind(iinp14)
    read(iinp14,'(a80)') header_p
    read(iinp14,'(a25$)') datetime_p
    read(iinp14,'(a55)') version_p
    if (trim(version_p)/="") lversion=.true.
    read(iinp14,'(2i5,e25.16)') nni_p,nmu_p,rinf_p
    read(iinp14,'(3e25.16)') z1_p,z2_p,r_p
    read(iinp14,'(3i5)') norb_p,nel_p,nexch_p

    rgrid_p=rinf_p
    rd =abs(r_p-r)
    z1d=abs(z1_p-z1)
    z2d=abs(z2_p-z2)
    
    write(iout6,1050)
    write(iout6,'(  6x,"   timestamp: ",a80)') datetime_p
    if (lversion) write(iout6,'(  6x,"     version: ",a15)') version_p
    write(iout6,'(  6x,"       title: ",80a1)') (title(i),i=1,80)

    if (iinterp.eq.0) then
       if (.not.lversion.and.ngrids.ne.ngrids_p ) then
          write (iout6,'("Mismatch in ngrids: ")')
          write (iout6,'("input="i2," file="i2)') ngrids,ngrids_p
          stop 'rheader'
       endif

       if (nni.ne.nni_p) then
          write (iout6,*) 'mismatch in nni: '
          write (iout6,*) 'input=',nni,' file=',nni_p
          stop 'rheader'
       endif

       if (nmu(1).ne.nmu_p) then
          write (iout6,*) 'mismatch in nmu for grid ',1,':'
          write (iout6,*) 'input=',nmu(1),' file=',nmu_p
          stop 'rheader'
       endif
    endif


    if (rd.gt.1.d-06) then
       write(iout6,'(1x,"Warning: mismatch in bond length:"," input=",f5.2,", 2dhf.dat=",f5.2)') r,r_p
    endif
    if (z1d.gt.1.d-06) then
       write(iout6,'(1x,"Warning: mismatch in Z1:"," input=",f5.2,", 2dhf.dat=",f5.2)') z1,z1_p       
    endif
    if(z2d.gt.1.d-06) then
       write(iout6,'(1x,"Warning: mismatch in Z2:"," input=",f5.2,", 2dhf.dat=",f5.2)') z2,z2_p       
    endif

    if (iinterp.eq.0) then
       if (.not.lversion.and.abs(rgrid(1)-rgrid_p).gt.1.d-6 ) then
          write (iout6,*) 'mismatch in a subgrid ',1,':'
          write (iout6,*) 'input=',rgrid(1),' file=',rgrid_p
          write (iout6,*) 'continuing with fingers crossed ...'
       endif
    endif

    if (nel.ne.nel_p) then
       write(iout6,'(1x,"Warning: mismatch in number of electrons:"," input=",i2,", 2dhf.dat=",i2)') nel,nel_p
    endif

    if (norb.ne.norb_p ) then
       write(iout6,'(1x,"Warning: mismatch in number of orbitals:"," input=",i2,", 2dhf.dat=",i2)') norb,norb_p
       if (norb/=norb_p) then       
          write (iout6,*) '... virtual orbital(s) present ...'
          ini4=norb-norb_p
       else
          stop 'readHeader'
       endif
    endif

    return

1000 continue
    write(iout6,1070)
    stop

1050 format(/1x,'... retrieving data from disk ...')
1070 format(//1x,'rheader: cannot read data from disk'//)

  end subroutine readHeader

  ! ### grid_old ###
  !
  !     Initializes grid info to enable interpolation.
  !
  subroutine grid_old
    use params
    use discrete
    use scfshr
    use commons
    implicit none

    integer (KIND=IPREC) :: i
    real (PREC) :: rinfig_p, xmi0

    mxsize_p=nni_p*nmu_p
    hni_p=pii/dble(nni_p-1)
    if (rgrid_p.ge.0.50_PREC) then
       xmi0=2.0_PREC*rgrid_p/r_p
       xmi0=log(xmi0+sqrt(xmi0*xmi0-1.0_PREC))
       hmu_p=xmi0/dble(nmu_p-1)
    else
       hmu_p=rgrid_p
    endif

    ! initialize mu and xi arryas
    vmu_p(1)=0.0_PREC
    vxi_p(1)=cosh(vmu_p(1))
    do i=2,nmu_p
       vmu_p(i)=vmu_p(i-1)+hmu_p
       vxi_p(i)=cosh(vmu_p(i))
    enddo
    ! determine rinf
    write(iout6,'(14x,"grid:",2i5,f7.2)') nni_p,nmu_p,rinf_p

    ! initialize ni and eta arrays

    do i=1,nni_p
       vni_p(i)=dble((i-1))*hni_p
       veta_p(i)=cos(vni_p(i))
       vetasq_p(i)=veta_p(i)*veta_p(i)

       ! sqrt(1-veta*veta) and  veta/sqrt(1-veta*veta)
       veta1_p(i)=sqrt(1.0_PREC-veta_p(i)*veta_p(i))
       if (veta1_p(i).lt.precis) then
          veta2_p(i)=0.0_PREC
       else
          veta2_p(i)=veta_p(i)/veta1_p(i)
       endif
    enddo

  end subroutine grid_old
  
  ! ### rfun ###
  !
  !     Reads functions from a disk file in an unformatted form
  !
  subroutine readFuncs (norbt,cw_orb,cw_coul,cw_exch,wk8,wk16)
    use params
    use commons
    use discrete
    use memory
    use scfshr
    use utils

    implicit none
    integer (KIND=IPREC) :: i,ierr,ii,ioffset,iorb1,iorb2,ipc,j,k,norbt
    integer (KIND=IPREC),dimension(10*maxorb+5*maxorb*(maxorb+1)/2) :: i4tmp
    integer*8, dimension(10*maxorb+5*maxorb*(maxorb+1)/2) :: i8tmp

    real (PREC), dimension(*) :: wk8,cw_orb,cw_coul,cw_exch
    real (PREC), dimension(maxorb) :: r8tmp1
    real (PREC), dimension(maxorb,maxorb) :: r8tmp12
    real (PREC), dimension(maxorb+(maxmpole-1)*maxorb) :: r8tmp2
    real (PREC), dimension(maxorb*(maxorb+1)) :: r8tmp3

    real (PREC16), dimension(*) :: wk16
    real (PREC16), dimension(maxorb) :: r16tmp1
    real (PREC16), dimension(maxorb,maxorb) :: r16tmp12
    real (PREC16), dimension(maxorb+(maxmpole-1)*maxorb) :: r16tmp2
    real (PREC16), dimension(maxorb*(maxorb+1)) :: r16tmp3

    if (inUnformatted) then    
       if (lengthint.eq.4) then
          read (iinp11,err=1000) i4tmp
       else
          read (iinp11,err=1000) i8tmp
       endif
    else
       read (iinp11,formint,err=1000) i4tmp
    endif

    ! if there are no virtual orbitals defined retrieve all orbitals and Coulomb functions
    ! from the corresponding disk data files ignoring inhyd flags. The flags are restored
    ! when returning from this routine.

    ! ioffset is no longer used
    ! ioffset=norb-norbt
    ioffset=0

    if (ini4.eq.0) then
       do i=1,norb
          inhyd(i)=0
       enddo
    endif

    do i=1,norb
       if (inhyd(i).eq.1) cycle

       if (inUnformatted) then
          if (lengthfp.eq.8) then
             call reada8(iinp11,mxsize,wk8,ierr)
             if (ierr.ne.0) then
                write(iout6,*) 'error detected when reading orbital',i
                !              stop 'rfun 8'
                write(iout6,*) 'continue with crossed fingers ...'
             endif
             do j=1,mxsize
                cw_orb(i1b(i)+j-1)=wk8(j)
             enddo
          endif

          if (lengthfp.eq.16) then
             call reada16(iinp11,mxsize,wk16,ierr)
             if (ierr.ne.0) then
                write(iout6,*) 'error detected when reading orbital',i
                stop 'rfun 16'
             endif
             do j=1,mxsize
                cw_orb(i1b(i)+j-1)=wk16(j)
             enddo
          endif
       else
          read (iinp11,formfp,err=1004) (wk8(ii),ii=1,mxsize)
          do j=1,mxsize
             cw_orb(i1b(i)+j-1)=wk8(j)
          enddo
       endif
    enddo

    !     retrieve the extra data from the orbital input file

    if (lengthfp.eq.8) then
       if (inUnformatted) then
          read(iinp11,end=1008,err=1010) r8tmp1
          call rfunaux(r8tmp1,orbNorm)
          if (lversion) then
             read(iinp11,end=1008,err=1010) r8tmp12
             ee=r8tmp12
          else
             read(iinp11,end=1008,err=1010) r8tmp1
             !call rfunaux(r8tmp1,eng)
             do i=1,maxorb
                ee(i,i)=r8tmp1(i)
             enddo
             
             read(iinp11,end=1008,err=1010) r8tmp3
             !call rfunaux(r8tmp3,engo)
             do i=1,maxorb
                do j=1,maxorb
                   ipc=j+(i-1)*norb
                   ee(i,j)=r8tmp3(ipc)
                enddo
             enddo
          endif
          read(iinp11,end=1008,err=1010) r8tmp2
          do i=1,maxorb+(maxmpole-1)*maxorb
             cmulti(i)=r8tmp2(i)
          enddo
          
          read(iinp11,end=1008,err=1010) r8tmp3
          do i=1,maxorb*(maxorb+1)
             excdi(i)=r8tmp3(i)
          enddo

          read(iinp11,end=1008,err=1010) r8tmp3
          do i=1,maxorb*(maxorb+1)
             excqu(i)=r8tmp3(i)
          enddo
          read(iinp11,end=1008,err=1010) r8tmp3
          do i=1,maxorb*(maxorb+1)
             excoc(i)=r8tmp3(i)
          enddo

          read(iinp11,end=1008,err=1010) r8tmp3
          do i=1,maxorb*(maxorb+1)
             exche(i)=r8tmp3(i)
          enddo

          read(iinp11,end=1008,err=1010) r8tmp3
          do i=1,maxorb*(maxorb+1)
             exc5(i)=r8tmp3(i)
          enddo

          read(iinp11,end=1008,err=1010) r8tmp3
          do i=1,maxorb*(maxorb+1)
             exc6(i)=r8tmp3(i)
          enddo

          read(iinp11,end=1008,err=1010) r8tmp3
          do i=1,maxorb*(maxorb+1)
             exc7(i)=r8tmp3(i)
          enddo
          read(iinp11,end=1008,err=1010) r8tmp3
          do i=1,maxorb*(maxorb+1)
             exc8(i)=r8tmp3(i)
          enddo
       else
          read(iinp11,formfp64,end=1008,err=1010) r8tmp1
          call rfunaux(r8tmp1,orbNorm)
          
          read(iinp11,formfp64,end=1008,err=1010) r8tmp1
          call rfunaux(r8tmp1,ee)
          
          !read(iinp11,formfp64,end=1008,err=1010) r8tmp3
          !call rfunaux(r8tmp3,engo)
          
          read(iinp11,formfp64,end=1008,err=1010) r8tmp2
          do i=1,maxorb+(maxmpole-1)*maxorb
             cmulti(i)=r8tmp2(i)
          enddo

          read(iinp11,formfp64,end=1008,err=1010) r8tmp3
          do i=1,maxorb*(maxorb+1)
             excdi(i)=r8tmp3(i)
          enddo

          read(iinp11,formfp64,end=1008,err=1010) r8tmp3
          do i=1,maxorb*(maxorb+1)
             excqu(i)=r8tmp3(i)
          enddo
          read(iinp11,formfp64,end=1008,err=1010) r8tmp3
          do i=1,maxorb*(maxorb+1)
             excoc(i)=r8tmp3(i)
          enddo
          read(iinp11,formfp64,end=1008,err=1010) r8tmp3
          do i=1,maxorb*(maxorb+1)
             exche(i)=r8tmp3(i)
          enddo
          read(iinp11,formfp64,end=1008,err=1010) r8tmp3
          do i=1,maxorb*(maxorb+1)
             exc5(i)=r8tmp3(i)
          enddo
          read(iinp11,formfp64,end=1008,err=1010) r8tmp3
          do i=1,maxorb*(maxorb+1)
             exc6(i)=r8tmp3(i)
          enddo
          read(iinp11,formfp64,end=1008,err=1010) r8tmp3
          do i=1,maxorb*(maxorb+1)
             exc7(i)=r8tmp3(i)
          enddo

          read(iinp11,formfp64,end=1008,err=1010) r8tmp3
          do i=1,maxorb*(maxorb+1)
             exc8(i)=r8tmp3(i)
          enddo
       endif
    endif


    if (lengthfp.eq.16) then
       if (inUnformatted) then
          read(iinp11,end=1008,err=1010) r16tmp1
          call rfunaux16(r16tmp1,orbNorm)

          read(iinp11,end=1008,err=1010) r16tmp12
          ee=r16tmp12

          read(iinp11,end=1008,err=1010) r16tmp2
          do i=1,maxorb+(maxmpole-1)*maxorb
             cmulti(i)=r16tmp2(i)
          enddo
          read(iinp11,end=1008,err=1010) r16tmp3
          do i=1,maxorb*(maxorb+1)/2
             excdi(i)=r16tmp3(i)
          enddo
          read(iinp11,end=1008,err=1010) r16tmp3
          do i=1,maxorb*(maxorb+1)/2
             excqu(i)=r16tmp3(i)
          enddo
          read(iinp11,end=1008,err=1010) r16tmp3
          do i=1,maxorb*(maxorb+1)/2
             excoc(i)=r16tmp3(i)
          enddo
          read(iinp11,end=1008,err=1010) r16tmp3
          do i=1,maxorb*(maxorb+1)/2
             exche(i)=r16tmp3(i)
          enddo
          read(iinp11,end=1008,err=1010) r16tmp3
          do i=1,maxorb*(maxorb+1)/2
             exc5(i)=r16tmp3(i)
          enddo
          read(iinp11,end=1008,err=1010) r16tmp3
          do i=1,maxorb*(maxorb+1)/2
             exc6(i)=r16tmp3(i)
          enddo
          read(iinp11,end=1008,err=1010) r16tmp3
          do i=1,maxorb*(maxorb+1)/2
             exc7(i)=r16tmp3(i)
          enddo
          read(iinp11,end=1008,err=1010) r16tmp3
          do i=1,maxorb*(maxorb+1)/2
             exc8(i)=r16tmp3(i)
          enddo
       else
          read(iinp11,formfp128,end=1008,err=1010) r16tmp1
          call rfunaux16(r16tmp1,orbNorm)

          read(iinp11,formfp128,end=1008,err=1010) r16tmp12
          ee=r16tmp12

          read(iinp11,formfp128,end=1008,err=1010) r16tmp2
          do i=1,maxorb+(maxmpole-1)*maxorb
             cmulti(i)=r16tmp2(i)
          enddo

          read(iinp11,formfp128,end=1008,err=1010) r16tmp3
          do i=1,maxorb*(maxorb+1)/2
             excdi(i)=r16tmp3(i)
          enddo

          read(iinp11,formfp128,end=1008,err=1010) r16tmp3
          do i=1,maxorb*(maxorb+1)/2
             excqu(i)=r16tmp3(i)
          enddo

          read(iinp11,formfp128,end=1008,err=1010) r16tmp3
          do i=1,maxorb*(maxorb+1)/2
             excoc(i)=r16tmp3(i)
          enddo

          read(iinp11,formfp128,end=1008,err=1010) r16tmp3
          do i=1,maxorb*(maxorb+1)/2
             exche(i)=r16tmp3(i)
          enddo

          read(iinp11,formfp128,end=1008,err=1010) r16tmp3
          do i=1,maxorb*(maxorb+1)/2
             exc5(i)=r16tmp3(i)
          enddo

          read(iinp11,formfp128,end=1008,err=1010) r16tmp3
          do i=1,maxorb*(maxorb+1)/2
             exc6(i)=r16tmp3(i)
          enddo

          read(iinp11,formfp128,end=1008,err=1010) r16tmp3
          do i=1,maxorb*(maxorb+1)/2
             exc7(i)=r16tmp3(i)
          enddo

          read(iinp11,formfp128,end=1008,err=1010) r16tmp3
          do i=1,maxorb*(maxorb+1)/2
             exc8(i)=r16tmp3(i)
          enddo
       endif
    endif
    rewind iinp11

500 continue
    ! read in Coulomb potentials
    do i=1,norb
       if (inhyd(i).eq.1) cycle
       if (inUnformatted) then
          if (lengthfp.eq.8) then
             call reada8(iinp12,mxsize,wk8,ierr)
             if (ierr.ne.0) then
                write(iout6,*) 'error detected when reading Coulomb potential',i
                stop 'rfun 8'
             endif
             if (lcoulexch) then
                do j=1,mxsize
                   cw_exch(i2b(i)+j-1)=wk8(j)
                enddo
             else
                do j=1,mxsize
                   cw_coul(i2b(i)+j-1)=wk8(j)
                enddo
             endif
          endif

          if (lengthfp.eq.16) then
             call reada16(iinp12,mxsize,wk16,ierr)
             if (ierr.ne.0) then
                write(iout6,*) 'error detected when reading Coulomb potential',i
                stop 'rfun 16'
             endif
             if (lcoulexch) then
                do j=1,mxsize
                   cw_exch(i2b(i)+j-1)=wk16(j)
                enddo
             else
                do j=1,mxsize
                   cw_coul(i2b(i)+j-1)=wk16(j)
                enddo
             endif
          endif
       else
          read (iinp12,formfp,err=1005) (wk8(ii),ii=1,mxsize)
          !            call readaf(iinp12,mxsize,wk8,ierr)
          if (lcoulexch) then
             do j=1,mxsize
                cw_exch(i2b(i)+j-1)=wk8(j)
             enddo
          else
             do j=1,mxsize
                cw_coul(i2b(i)+j-1)=wk8(j)
             enddo
          endif
       endif
200    continue
    enddo

    if (HFinput) then
       if (.not.linitFuncsNoexch) then
          do iorb1=1,norbt
             do iorb2=iorb1,norbt
                k=k2(iorb1,iorb2)
                if (iorb1.eq.iorb2.and.ll(iorb1).eq.0) cycle
                if (inUnformatted) then
                   if (lengthfp.eq.8) then
                      call reada8(iinp13,mxsize,wk8,ierr)
                      if (ierr.ne.0) then
                         write(iout6,*) 'error detected when reading exchange potential',iorb1,iorb2,k
                         stop 'rfun 8'
                      endif
                      do j=1,mxsize
                         cw_exch(i3b(k)+j-1)=wk8(j)
                      enddo
                   endif
                   if (lengthfp.eq.16) then
                      call reada16(iinp13,mxsize,wk16,ierr)
                      if (ierr.ne.0) then
                         write(iout6,*) 'error detected when reading exchange potential',iorb1,iorb2,k
                         stop 'rfun 16'
                      endif
                      do j=1,mxsize
                         cw_exch(i3b(k)+j-1)=wk16(j)
                      enddo
                   endif
                else
                   read (iinp13,formfp,err=1010) (wk8(j),j=1,mxsize)
                   do j=1,mxsize
                      cw_exch(i3b(k)+j-1)=wk8(j)
                   enddo
                endif
                
                if (iorb1.eq.iorb2) cycle
                if (ll(iorb1).eq.0.or.ll(iorb2).eq.0) cycle

                if (inUnformatted) then                
                   if (lengthfp.eq.8) then
                      call reada8(iinp13,mxsize,wk8,ierr)
                      if (ierr.ne.0) then
                         write(iout6,*) 'error detected when reading exchange potential',iorb1,iorb2,k
                         stop 'rfun 8'
                      endif
                      do j=1,mxsize
                         cw_exch(i3b(k)+mxsize+j-1)=wk8(j)
                      enddo
                   endif
                   if (lengthfp.eq.16) then
                      call reada16(iinp13,mxsize,wk16,ierr)
                      if (ierr.ne.0) then
                         write(iout6,*) 'error detected when reading exchange potential',iorb1,iorb2,k
                         stop 'rfun 16'
                      endif
                      do j=1,mxsize
                         cw_exch(i3b(k)+mxsize+j-1)=wk16(j)
                      enddo
                   endif
                else
                   read(iinp13,formfp,err=1010) (wk8(j),j=1,mxsize)
                   do j=1,mxsize
                      cw_exch(i3b(k)+mxsize+j-1)=wk8(j)
                   enddo
                endif
             enddo
          enddo
          rewind iinp13
       endif
       write(*,*) '... orbitals and potentials retrieved ... '
    endif

    !if (LXC.or.DFT.or.HFS.or.SCMC) then
    if (DFT.or.HFS.or.SCMC) then
       if (inUnformatted) then
          if (lengthfp.eq.8) then
             call reada8(iinp13,i3si(1),wk8,ierr)
             if (ierr.ne.0) then
                write(iout6,*) 'error detected when reading local exchange potential'
                stop 'rfun 8'
             endif
             do j=1,i3si(1)
                cw_exch(length3-mxsize+j-1)=wk8(j)
             enddo
          endif
          if (lengthfp.eq.16) then
             call reada16(iinp13,i3si(1),wk16,ierr)
             if (ierr.ne.0) then
                write(iout6,*) 'error detected when reading local exchange potential'
                stop 'rfun 16'
             endif
             do j=1,i3si(1)
                cw_exch(length3-mxsize+j-1)=wk16(j)
             enddo
          endif
       else
          read (iinp13,formfp,err=1010) (wk8(j),j=1,i3si(1))
          do j=1,i3si(1)
             cw_exch(length3-mxsize+j-1)=wk8(j)
          enddo
       endif
    endif


    if (lengthfpin.eq.8) then
       formfp=formfp64
    else
       formfp=formfp128
    endif

#ifdef PRINT    
! print=115: readFuncs: print out input orbitals
    if(iprint(115).ne.0) then
       write(*,*)
       write(*,*) 'rfun: input orbitals'
       do i=1,norb
          write(*,*) iorn(i),' ',bond(i)
          call pmtx(nni,i1mu(i),cw_orb(i1b(i)),ione,ione,incrni,incrmu)
       enddo
    endif
#endif

#ifdef PRINT    
! print=116: readFuncs: print out input Coulomb potentials
    if(iprint(116).ne.0) then
       write(*,*)
       write(*,*) 'rfun: input Coulomb potentials'
       !        print Coulomb potential
       do i=1,norb
          write(*,*) iorn(i),' ',bond(i)
          call pmtx(nni,i1mu(i),cw_coul(i2b(i)),ione,ione,incrni,incrmu)
       enddo
    endif
#endif

#ifdef PRINT    
! print=117: readFuncs: print out input exchange potentials
    if(iprint(117).ne.0.and.nexch.ge.1) then
       write(*,*)
       write(*,*) 'rfun: input exchange potentials'
       !        print exchange potentials
       do i=1,nexch
          write(*,*) 'exchange potential ',i
          call pmtx(nni,i1mu(i),cw_exch(i3b(i)),ione,ione,incrni,incrmu)
       enddo
    endif
#endif
    do i=1,norb
       inhyd(i)=inhydlcao(i)
    enddo

    return

1000 continue
    write(*,*) '... error detected when reading disk file ... '
    stop 'rfun'

1004 continue
    write(*,*) '... error detected when reading orbital function', i
    stop 'rfun'

1005 continue
    write(*,*) '... error detected when reading potential function', i
    stop 'rfun'

1008 continue
    write(*,*) '... end of file detected when reading orbital input file ...'
    write(*,*) '... data without extension retrieved ... '
    write(*,*) '... continue reading potentials ... '    
    initAddData=.true.
    goto 500
    
1010 write(*,*) '... error detected when reading exchange function ...'
    initAddData=.true.

1100 continue

  end subroutine readFuncs


  ! ### readInterpFuncs ###
  !
  !     Reads functions from disk in unformatted form and interpolates
  !     to a new grid.
  !
  subroutine readInterpFuncs (norb_p,cw_orb,cw_coul,cw_exch,wk8,wk16,cw_sctch)
    use params
    use discrete
    use scfshr
    use commons
    use interpolate

    implicit none

    character*8 arrayName
    integer (KIND=IPREC) :: i,ica,idel,ierr,ioffset,iorb1,iorb2,ipex,isym,j,k,norb_p

    integer (KIND=IPREC),dimension(maxorb) :: i1b_p,i2b_p,i1e_p,i2e_p,i1si_p,i2si_p,i1ng_p,i2ng_p,i1mu_p,i2mu_p
    integer (KIND=IPREC),dimension(maxorb*(maxorb+1)/2) :: i3b_p,i3e_p,i3si_p,i3ng_p,i3mu_p

    real (PREC), dimension(*) :: wk8,cw_orb,cw_coul,cw_exch,cw_sctch

    real (PREC16), dimension(*) :: wk16
    real (PREC), dimension(:), allocatable :: fbefore

    read (iinp11,err=1000) i1b_p,i2b_p,i3b_p,i1e_p,i2e_p,i3e_p, &
         i1si_p,i2si_p,i3si_p,i1ng_p,i2ng_p,i3ng_p,i1mu_p,i2mu_p,i3mu_p

    ioffset=norb-norb_p
    ica=1
    write(iout6,1100)
01100 format(' ... interpolating orbitals:',/,'  ',$)
01110 format(i4,$)

    allocate(fbefore(nni_p*nmu_p))

    do i=1,norb_p
       write(iout6,1110) i
       ipex=mgx(6,i)
       ipex=ipex-2*(ipex/2)
       if (ipex.eq.0) then
          isym= 1
       else
          isym=-1
       endif
       if (lengthfp.eq.8) then
          call reada8(iinp11,i1si_p(i+ioffset),wk8,ierr)
          if (ierr.ne.0) then
             write(iout6,*) 'error detected when reading orbital',i
             stop 'rfun 8'
          endif
          do j=1,i1si_p(i+ioffset)
             fbefore(j)=wk8(j)
          enddo
       endif

       if (lengthfp.eq.16) then
          call reada16(iinp11,i1si_p(i+ioffset),wk16,ierr)
          if (ierr.ne.0) then
             write(iout6,*) 'error detected when reading orbital',i
             stop 'rfun 16'
          endif
          do j=1,i1si_p(i+ioffset)
             fbefore(j)=wk16(j)
          enddo
       endif

       call dointerp (ica,nmu_p,nmu(1),fbefore,cw_orb(i1b(i+ioffset)))

       if (ierr.ne.0) then
          write(iout6,*) 'error detected when reading orbital',i
          stop 'inifun'
       endif
    enddo

    !     retrieve extra data from the orbital input file

    arrayName="orbNorm"    
    read(iinp11,end=1010,err=1012) orbNorm
    arrayName="ee"
    read(iinp11,end=1010,err=1012) ee
    arrayName="cmulti"
    read(iinp11,end=1010,err=1012) cmulti
    read(iinp11,end=1010,err=1010) excdi
    read(iinp11,end=1010,err=1010) excqu
    read(iinp11,end=1010,err=1010) excoc
    read(iinp11,end=1010,err=1010) exche
    read(iinp11,end=1010,err=1010) exc5
    read(iinp11,end=1010,err=1010) exc6
    read(iinp11,end=1010,err=1010) exc7
    read(iinp11,end=1010,err=1010) exc8

    write(iout6,*)
    rewind iinp11

    ica=2
    isym=1
    write(iout6,1102)
01102 format(' ... interpolating Coulomb potentials:',/,'  ',$)

    do i=1,norb_p
       write(iout6,1110) i
       if (lengthfp.eq.8) then
          call reada8(iinp12,i2si_p(i+ioffset),wk8,ierr)
          if (ierr.ne.0) then
             write(iout6,*) 'error detected when reading coulomb potential',i
             stop 'rfun 8'
          endif
          do j=1,i2si_p(i+ioffset)
             fbefore(j)=wk8(j)
          enddo
       endif
       if (lengthfp.eq.16) then
          call reada16(iinp12,i2si_p(i+ioffset),wk16,ierr)
          if (ierr.ne.0) then
             write(iout6,*) 'error detected when reading coulomb potential',i
             stop 'rfun 16'
          endif
          do j=1,i2si_p(i+ioffset)
             fbefore(j)=wk16(j)
          enddo
       endif

       call dointerp (ica,nmu_p,nmu(1),fbefore,cw_coul(i2b(i+ioffset)) )

       if (ierr.ne.0) then
          write(iout6,*) 'error detected when reading coulomb potential',i
          stop 'rfun_int'
       endif
    enddo
    write(iout6,*)
    rewind iinp12

    if (OED) then
       write(*,*) 'rfun_int: disk file without extension retrieved '
       write(*,*) '      '
       rewind iinp13
       deallocate(fbefore)
       return
    endif

    if (HF) then
       ica=3

       write(iout6,1104)
01104  format(' ... interpolating exchange potentials for orbitals:',/,'  ',$)
       do iorb1=1,norb_p
          write(iout6,1110) iorb1
          do iorb2=iorb1,norb_p
             k=k2(iorb1,iorb2)
             idel=abs(mgx(6,iorb1)-mgx(6,iorb2))
             if (iorb1.eq.iorb2) idel=2*mgx(6,iorb1)
             ipex=idel-2*(idel/2)
             if (ipex.eq.0) then
                isym= 1
             else
                isym=-1
             endif

             if (iorb1.eq.iorb2.and.ll(iorb1).eq.0) goto 50
             if (lengthfp.eq.8) then
                call reada8(iinp13,i3si_p(k),wk8,ierr)
                if (ierr.ne.0) then
                   write(iout6,*) 'error detected when reading exchange potential',iorb1,iorb2,k
                   stop 'rfun 8'
                endif
                do j=1,i3si_p(k)
                   fbefore(j)=wk8(j)
                enddo
             endif
             if (lengthfp.eq.16) then
                call reada16(iinp13,i3si_p(k),wk16,ierr)
                if (ierr.ne.0) then
                   write(iout6,*) 'error detected when reading exchange potential',iorb1,iorb2,k
                   stop 'rfun 16'
                endif
                do j=1,i3si_p(k)
                   fbefore(j)=wk16(j)
                enddo
             endif

             call dointerp (ica,nmu_p,nmu(1),fbefore,cw_exch(i3b(k)))
             if (ierr.ne.0) then
                write(iout6,*) 'error detected when reading exchange potential',iorb1,iorb2,k
                stop 'rfun_int'
             endif
             if (iorb1.eq.iorb2) goto 50
             if (ll(iorb1).eq.0.or.ll(iorb2).eq.0) goto 50
             if (lengthfp.eq.8) then
                call reada8(iinp13,i3si_p(k),wk8,ierr)
                if (ierr.ne.0) then
                   write(iout6,*) 'error detected when reading exchange potential',iorb1,iorb2,k
                   stop 'rfun 8'
                endif
                do j=1,i3si_p(k)
                   fbefore(j)=wk8(j)
                enddo
             endif
             if (lengthfp.eq.16) then
                call reada16(iinp13,i3si_p(k),wk16,ierr)
                if (ierr.ne.0) then
                   write(iout6,*) 'error detected when reading exchange potential',iorb1,iorb2,k
                   stop 'rfun 16'
                endif
                do j=1,i3si_p(k)
                   fbefore(j)=wk16(j)
                enddo
             endif

             call dointerp (ica,nmu_p,nmu(1),fbefore,cw_exch(i3b(k)+i3si(k)))

             if (ierr.ne.0) then
                write(iout6,*) 'error detected when reading exchange potential',iorb1,iorb2,k
                stop 'rfun_int'
             endif
00050        continue
          enddo
       enddo
       write(iout6,*)
       write(iout6,*)
       rewind iinp13
    endif

    deallocate(fbefore)
    return
01000 continue
    write(*,*) '... error detected when retrieving the disk file ... '
    stop 'rfun_int'
01010 write(*,*) '... error detected when retrieving extension of the orbital input file ...'
    write(*,*) '... disk file without extension retrieved ... '

01012 write(*,*) '... error detected when retrieving ',trim(arrayName),' of the orbital input file ...'
    write(*,*) '... disk file without extension retrieved ... '

    initAddData=.true.
  end subroutine readInterpFuncs
  
  ! ### reada8 ###
  !
  !     Reads the matrix A from a disk file in an unformatted form
  !
  subroutine reada8 (iunit,ndim,a,ierr)
    use params
    implicit none
    integer (KIND=IPREC) :: ierr,iunit,ndim
    real (PREC), dimension(ndim) :: a

    read (iunit,IOSTAT=ierr) a
    ierr=0
    return
  end subroutine reada8
  
  ! ### reada8 ###
  !
  !     Reads the matrix A from a disk file in an unformatted form
  !
  subroutine reada16 (iunit,ndim,a,ierr)
    use params
    implicit none
    integer (KIND=IPREC) :: ierr,iunit,ndim
    real (PREC16), dimension(ndim) :: a

    read (iunit,err=1000) a
    ierr=0
    return
1000 ierr=1
  end subroutine reada16

  ! ### rfunaux ###
  !
  !     Reads functions from a disk file in an unformatted form
  !
  subroutine rfunaux (s,t)
    use params
    use commons

    implicit none
    integer (KIND=IPREC) :: i,ioffset
    real (PREC), dimension(*) :: s
    real (PREC16), dimension(*) :: t

    ioffset=0
    do i=1,maxorb
       if (inhyd(i).eq.0) then
          t(i)=s(i-ioffset)
       else
          ioffset=ioffset+1
       endif
       !         print *,norb,i,ioffset,inhyd(i),t(i),s(i-ioffset)
    enddo

  end subroutine rfunaux

  ! ### rfunaux16 ###
  !
  !     Reads functions from a disk file in an unformatted form
  !
  subroutine rfunaux16 (s,t)
    use params
    use commons

    implicit none
    integer (KIND=IPREC) :: i,ioffset
    real (PREC16), dimension(*) :: s
    real (PREC16), dimension(*) :: t

    ioffset=0
    do i=1,maxorb
       if (inhyd(i).eq.0) then
          t(i)=s(i-ioffset)
       else
          ioffset=ioffset+1
       endif
    enddo

  end subroutine rfunaux16

  subroutine readFromDisk 
    use params
    use sharedMemory
    use commons
    use discrete

    implicit none
    integer (KIND=IPREC) :: lengthint_curr,lengthfp_curr
    integer (KIND=IPREC) :: norb_p

    real (PREC), dimension(:), pointer :: cw_orb,cw_exch,cw_sctch,&
         scratch1,scratch2,scratch3
    
    !     input files can be in i32, i64 and r128 formats
    !     save the current formats
    lengthint_curr=lengthint
    lengthfp_curr=lengthfp

    !     and determine ones used in input files
    lengthint=lengthintin
    lengthfp=lengthfpin

    call readHeader(norb_p)
    if (iinterp.eq.0) then
       scratch1=>scratchptr(1:mxsize)
       scratch2=>scratchptr(mxsize+1:2*mxsize)
       call readFuncs (norb_p,orbptr,exchptr,exchptr,scratch1,scratch2)
    else
       scratch1=>scratchptr(1:mxsize)
       scratch2=>scratchptr(mxsize+1:2*mxsize)
       scratch3=>scratchptr(2*mxsize+1:3*mxsize)
       call grid_old
       call readInterpFuncs (norb_p,orbptr,exchptr,exchptr,scratch1,scratch2,scratch3)
    endif

    lengthint=lengthint_curr
    lengthfp=lengthfp_curr

  end subroutine readFromDisk

! ### writeToDisk ###
!     Writes orbitals, potenials, Lagrange multipliers (diagonal
!     and off-diagonal) and multipole expansion coefficients to a disk
!     file in either formatted or unformatted form
  subroutine writeToDisk
    use params
    use discrete
    use scfshr
    use commons
    use dateTime
    use dumpDataToDisk
    use sharedMemory
    use data4II, only : i,title
    
    implicit none
    character*80 :: currDateTime

    integer (KIND=IPREC),allocatable :: i4tmp1(:)
    integer (KIND=IPREC8), allocatable :: i8tmp1(:)

    real (PREC), allocatable :: r8tmp1(:),r8tmp12(:,:),r8tmp2(:),r8tmp3(:),r8tmp4(:)
    real (PREC), allocatable :: r16tmp1(:),r16tmp12(:,:),r16tmp2(:),r16tmp3(:),r16tmp4(:)
    real (PREC), dimension(:), pointer :: cw_orb,cw_coul,cw_exch

    cw_orb=>orbptr
    cw_coul=>exchptr
    cw_exch=>exchptr

    ! fixed precision
    rewind(iout24)
    call getDateTime(currDateTime)
    !write(iout24,'(a80)') header
    write(iout24,'(75a1)') (title(i),i=1,75)
    !write(iout24,'("")') 
    !write(iout24,'(a80)') currDateTime
    write(iout24,'(a25$)') currDateTime
    write(iout24,'(a15)') version
    !write(iout24,formint) ngrids,nni,nmu
    write(iout24,'(2i5,e25.16,"   / nnu nmu Rinf")') nni,nmu,rinf
    !write(iout24,formfp64) r,rgrid
    write(iout24,'(3e25.16,"   / Z1 Z2 R")') z1,z2,r
    write(iout24,'(3i5,"    / norb nel nexch")') norb,nel,nexch

    call flush(iout24)

    if (outFormatI32.or.outFormatI64.or.outFormatR128) then
       allocate(i4tmp1(10*maxorb+5*maxorb*(maxorb+1)/2))
       allocate(i8tmp1(10*maxorb+5*maxorb*(maxorb+1)/2))
       allocate(r8tmp1(maxorb))
       allocate(r8tmp12(maxorb,maxorb))
       allocate(r8tmp2(maxorb+(maxmpole-1)*maxorb))
       allocate(r8tmp3(maxorb*(maxorb+1)))
       allocate(r8tmp4(mxsize))

       if (outFormatI32)  call wtdisk32 &
            (cw_orb,cw_coul,cw_exch,i4tmp1,r8tmp1,r8tmp12,r8tmp2,r8tmp3,r8tmp4)
       if (outFormatI64)  call wtdisk64 &
            (cw_orb,cw_coul,cw_exch,i8tmp1,r8tmp1,r8tmp12,r8tmp2,r8tmp3,r8tmp4)
       if (outFormatr128) call wtdisk128 &
            (cw_orb,cw_coul,cw_exch,i8tmp1,r8tmp1,r8tmp12,r8tmp2,r8tmp3,r8tmp4)
       
       deallocate(i4tmp1)
       deallocate(i8tmp1)
       deallocate(r8tmp1)
       deallocate(r8tmp12)
       deallocate(r8tmp2)
       deallocate(r8tmp3)
       deallocate(r8tmp4)
       return
    endif
    
    call wtdisknat(cw_orb(1:),cw_exch(1:))

  end subroutine writeToDisk
end module diskInterface
