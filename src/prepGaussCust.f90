! ***************************************************************************
! *                                                                         *
! *   Copyright (C) 1997-2010 Jacek Kobus <jkob@fizyka.umk.pl>              *
! *                                                                         *
! *   This program is free software; you can redistribute it and/or modify  *
! *   it under the terms of the GNU General Public License version 2 as     *
! *   published by the Free Software Foundation.                            *
! *                                                                         *
! ***************************************************************************
! ### prepGaussCust ###
!
!     Reads the (customized) output from the GAUSSSIAN94/98 program
!     (gausssianl.out file) c to determine parameters of the basis
!     functions (n, l, m, exponents) c and coefficients of molecular
!     orbitals.

!     The basis can contain functions defined on centres A and B and
!     also at the bond centre.

!     Last modified 7/09/02
!     This routine extracts data from the customized (long) output
!     (see ISKIP below)

subroutine prepGaussCust
  use params
  use commons8


  implicit none
  integer :: i,ib,ibc,ibpt,ibr,ibret,ibrett,ibt,icent,icentre,icount,idavid,ifbo,ilabel,ilines,iorb, &
       ipb,iprint0,iprint1,iprint2,iskip,j,k,l1,m1,m1abs,n1prim,n2prim,n3prim,nexpon,nfborb

  integer, dimension(maxorb) :: ifbord,ifdord

  character*3, dimension(maxbasis) :: conbt
  real (PREC) :: d1,excoeffm,excoeffp,expon,symthresh
  real (PREC), dimension(maxbasis) :: ibfres
  real (PREC), dimension(3,maxbasis) :: sexpon
  real (PREC), dimension(0:maxorb,0:maxbasis) :: excoeff
  real (PREC), external :: factor,factor2

  character*2 :: st
  character*9 :: matchstr

  parameter (symthresh=0.0010_PREC)

  !     nfborb - number of finite basis (fb) set orbitals
  !     ifbord - ordering of fb orbitals (0=sigma, 1=pi, 2=delta, etc)
  !              the order is determined from the gaussian output file
  !     ifdord - ordering of finite difference (fd) orbitals
  !              the order of fd orbitals is the same as in the input data
  !              (phi orbitals first, then delta, etc)
  !
  ! BF
  !      data nfborb /7/
  !      data ifbord /0,0,0,0,1,1,0,53*100/
  !      data ifdord /1,0,0,0,0,0,54*100/

  ! InF
  !          data nfborb /29/
  !          data ifbord /0,0,1,1,0,
  !     &                 0,0,1,1,0,
  !     &                 2,2,1,1,0,
  !     &                 0,1,1,0,0,
  !     &                 2,2,1,1,0,
  !     &                 0,1,1,0,100,
  !     &                 30*100/
  !     &
  !        data ifdord /2*2,6*1,13*0,39*100/

  !  opening a file with the corresponding (customized) GAUSSIAN94/98 output

  open(7,file='gaussianc.out', status='old',form='formatted')

  !     ibc counts contracted basis functions
  !     ibp counts primitive basis functions
  !     n1prim - number of primitive gaussian function at the centre A
  !     n3prim - number of primitive gaussian function at the centre B
  !     n2prim - number of primitive gaussian function at the bond centre

  iprint0=0
  iprint1=0
  iprint2=0

  if (iprint(553).ne.0) iprint0=1
  if (iprint(554).ne.0) iprint1=1
  if (iprint(555).ne.0) iprint2=1

  !     format of data should be detected automatically
  !     if Gaussian9x output contains eigenvectors in 2-column format iskip=2
  !     if Gaussian9x output contains eigenvectors in 5-column format iskip=5

  iskip=2
  !      if (idbg(560).ne.0) iskip=5

  !     determine number of molecular orbitals as defined in the GAUSSIAN9x
  !     program (one nonsigma finite difference orbital corresponds
  !     to two GAUSSIAN9x ones

  do i=0,maxorb
     do j=1,maxbasis
        excoeff(i,j)=0.0_PREC
     enddo
  enddo

  nfborb=0
  do iorb=1,norb
     ifdord(iorb)=mm(iorb)
     if (mm(iorb).eq.0) then
        nfborb=nfborb+1
     else
        nfborb=nfborb+2
     endif
  enddo

  write(*,*)
  icentre=0
  ilabel=1
  ib=0
  ibpt=0
  icent=0
  do ilines=1,100000
     !        read(7,1001,end=991,err=991) matchstr,icentt
     read(7,1001,end=991,err=991) matchstr
     !        write(*,1001) matchstr
01000 format(a9)
01001 format(a9,i5)
     if (matchstr.eq.' Centers:') then
        icent=icent+1
        icentre=icentre+1
        ibc=0
        ipb=0

        !  start extracting data from the output

00100   ib=ib+1
        ibc=ibc+1
        !           read(7,*,end=910,err=910)  st,expon
        read(7,1010,end=910,err=910)  st
        if (st.eq.'**') goto 910

01010   format(1x,a2,6x,f18.9)
01011   format(14x,f18.9)
        read(7,1011,end=910,err=910)  expon
        sexpon(icent,ibc) =expon

        if (ibpt.ge.maxbasis) then
           stop ' prepGaussCust: MAXBASIS too small'
        endif

        if (st.eq.'s '.or.st.eq.'S ') then
           ipb=ipb+1
           ibpt=ibpt+1
           primexp(ibpt)=expon
           lprim(ibpt)=0
           mprim(ibpt)=0
           icgau(ibpt)=icent
        elseif (st.eq.'p '.or.st.eq.'P ') then
           ipb=ipb+3
           do k=1,3
              ibpt=ibpt+1
              primexp(ibpt)=expon
              lprim(ibpt)=1
              icgau(ibpt)=icent
              if (k.eq.1) mprim(ibpt)=+1
              if (k.eq.2) mprim(ibpt)=-1
              if (k.eq.3) mprim(ibpt)= 0
           enddo
        elseif (st.eq.'d '.or.st.eq.'D ') then
           ipb=ipb+5
           do k=1,5
              ibpt=ibpt+1
              primexp(ibpt)=expon
              lprim(ibpt)=2
              icgau(ibpt)=icent
              if (k.eq.1) mprim(ibpt)= 0
              if (k.eq.2) mprim(ibpt)=+1
              if (k.eq.3) mprim(ibpt)=-1
              if (k.eq.4) mprim(ibpt)=+2
              if (k.eq.5) mprim(ibpt)=-2
           enddo
        elseif (st.eq.'f '.or.st.eq.'F ') then
           ipb=ipb+7
           do k=1,7
              ibpt=ibpt+1
              primexp(ibpt)=expon
              lprim(ibpt)=3
              icgau(ibpt)=icent
              if (k.eq.1) mprim(ibpt)= 0
              if (k.eq.2) mprim(ibpt)=+1
              if (k.eq.3) mprim(ibpt)=-1
              if (k.eq.4) mprim(ibpt)=+2
              if (k.eq.5) mprim(ibpt)=-2
              if (k.eq.6) mprim(ibpt)=+3
              if (k.eq.7) mprim(ibpt)=-3
           enddo
        elseif (st.eq.'g '.or.st.eq.'G ') then
           ipb=ipb+9
           do k=1,9
              ibpt=ibpt+1
              primexp(ibpt)=expon
              lprim(ibpt)=4
              icgau(ibpt)=icent
              if (k.eq.1) mprim(ibpt)= 0
              if (k.eq.2) mprim(ibpt)=+1
              if (k.eq.3) mprim(ibpt)=-1
              if (k.eq.4) mprim(ibpt)=+2
              if (k.eq.5) mprim(ibpt)=-2
              if (k.eq.6) mprim(ibpt)=+3
              if (k.eq.7) mprim(ibpt)=-3
              if (k.eq.8) mprim(ibpt)=+4
              if (k.eq.9) mprim(ibpt)=-4
           enddo
        elseif (st.eq.'h '.or.st.eq.'H ') then
           ipb=ipb+11
           do k=1,11
              ibpt=ibpt+1
              primexp(ibpt)=expon
              lprim(ibpt)=5
              icgau(ibpt)=icent
              if (k.eq. 1) mprim(ibpt)= 0
              if (k.eq. 2) mprim(ibpt)=+1
              if (k.eq. 3) mprim(ibpt)=-1
              if (k.eq. 4) mprim(ibpt)=+2
              if (k.eq. 5) mprim(ibpt)=-2
              if (k.eq. 6) mprim(ibpt)=+3
              if (k.eq. 7) mprim(ibpt)=-3
              if (k.eq. 8) mprim(ibpt)=+4
              if (k.eq. 9) mprim(ibpt)=-4
              if (k.eq.10) mprim(ibpt)=+5
              if (k.eq.11) mprim(ibpt)=-5
           enddo
        elseif (st.eq.'sp'.or.st.eq.'SP') then
           ipb=ipb+4
           ibpt=ibpt+1
           primexp(ibpt)=expon
           lprim(ibpt)=0
           mprim(ibpt)=0
           icgau(ibpt)=icent
           do k=1,3
              ibpt=ibpt+1
              primexp(ibpt)=expon
              lprim(ibpt)=1
              icgau(ibpt)=icent
              if (k.eq.1) mprim(ibpt)=+1
              if (k.eq.2) mprim(ibpt)=-1
              if (k.eq.3) mprim(ibpt)= 0
           enddo
        endif
        if (iprint0.ne.0) then
           print *,ibpt,icgau(ibpt),lprim(ibpt),mprim(ibpt),primexp(ibpt)
        endif
        goto 100

910     continue
        ib=ib-1
        if (icentre.eq.1) then
           ibret=ib
           ibrett=ibret
           n1prim=ipb
        elseif (icentre.eq.2) then
           ibret=ib-ibrett
           ibrett=ibrett+ibret
           n2prim=ipb
        elseif (icentre.eq.3) then
           ibret=ib-ibrett
           ibrett=ibrett+ibret
           n3prim=ipb
        else
           stop "Invalid center in prepGaussCust"
        endif
        write(*,1141) icent,ibret
     endif
991  continue
     if (icentre.eq.3) goto 992
  enddo
992 continue
  rewind (7)
  nexpon=ib

  npbasis=n1prim+n2prim+n3prim
  write(*,1140) nexpon,npbasis,n1prim,n2prim,n3prim

  idavid=1
  !     start extracting basis set expansion coefficientsdata from the output
  !     since David Moncrieff supplied data in various formats they are all
  !     tried one by one

900 continue
  rewind(7)
  do ilines=1,1000000
     read(7,1000,end=990,err=990) matchstr
     if (iprint1.ne.0) write(*,*) '1',matchstr,'2'
     if (matchstr.eq.'     EIGE') then
        do iorb=1,nfborb,iskip
           do ib=1,npbasis
              if (idavid.eq.1) then
                 !                    iskip.eq.2
                 read(7,1015,end=940,err=800) ibt,ibr,conbt(ib),excoeff(iorb,ib),excoeff(iorb+1,ib)
1015             format(i4,6x,i3,a3,11x,f19.16,5x,f19.16)
              elseif (idavid.eq.2) then
                 !                    iskip.eq.2
                 read(7,1115,end=940,err=801) ibt,ibr,conbt(ib),excoeff(iorb,ib),excoeff(iorb+1,ib)
                 !                    1115                format(a4,6x,i3,a3,5x,2f25.16)
1115             format(i5,6x,i3,a3,5x,2f20.16)
              elseif (idavid.eq.10) then
                 !                    iskip.eq.5
                 read(7,1016,end=940,err=810) ibt,ibr,conbt(ib),excoeff(iorb  ,ib), &
                      excoeff(iorb+1,ib),excoeff(iorb+2,ib),excoeff(iorb+3,ib),excoeff(iorb+4,ib)
1016             format(a4,6x,i3,a3,5x,5f10.5)
              elseif (idavid.eq.11) then
                 !                    iskip.eq.5
                 read(7,1016,end=940,err=811) ibt,ibr,conbt(ib),excoeff(iorb  ,ib), &
                      excoeff(iorb+1,ib),excoeff(iorb+2,ib),excoeff(iorb+3,ib),excoeff(iorb+4,ib)
!1116             format(a4,6x,i3,a3,5x,5f10.5)
              endif

              if (abs(excoeff(iorb,ib)).gt.1.0_PREC.and.iorb.le.nfborb) then
                 write(*,1130) iorb,ib,excoeff(iorb,ib)
              endif
              if (abs(excoeff(iorb+1,ib)).gt.1.0_PREC.and.iorb+1.le.nfborb) then
                 write(*,1130) iorb+1,ib,excoeff(iorb+1,ib)
              endif

              ibfres(ib)=ibr
              if (ib.le.n1prim) icent=1
              if (ib.gt.n1prim.and.ib.le.(n1prim+n2prim)) icent=2
              if (ib.gt.(n1prim+n2prim).and.ib.le.(n1prim+n2prim+n3prim)) icent=3

              if (iprint1.ne.0) then
                 if (iskip.eq.2) then
                    write(*,1013,err=990) iorb,ib,ibr,icent,conbt(ib), &
                         excoeff(iorb,ib),excoeff(iorb+1,ib),sexpon(icent,ibr)
                 elseif(iskip.eq.5) then
                    write(*,1014,err=990) iorb,ib,ibr,icent,conbt(ib), &
                         excoeff(iorb,  ib),excoeff(iorb+1,ib),excoeff(iorb+2,ib), &
                         excoeff(iorb+3,ib),excoeff(iorb+4,ib),sexpon(icent,ibr)
                 endif
              endif
           enddo

           read(7,1000,end=950,err=950) matchstr
           read(7,1000,end=950,err=950) matchstr
           read(7,1000,end=950,err=950) matchstr
        enddo
     endif
  enddo

800 continue
  idavid=2
  iskip=2
  goto 900

801 continue
  idavid=10
  iskip=5
  goto 900

810 continue
  idavid=11
  iskip=5
  goto 900

811 continue
  WRITE(*,*) 'prepGaussCust: check format of gaussianc.out file'
  stop 'prepGaussCust'

940 continue
  write(*,*) 'prepGaussCust: uncomplete gaussianc file', ib-1,ibt
  stop 'prepGaussCust'
950 continue
  if (iprint1.ne.0) write(*,*) 'warning: 950'
990 continue
  if (iprint1.ne.0) write(*,*) 'warning: 990'

  write(*,1142) norb,nfborb

  !     determine symmetry of GAUSSIAN94 molecular orbitals

  icount=0
  do ifbo=nfborb,1,-1
     do ib=npbasis,1,-1
        if (abs(excoeff(ifbo,ib)).gt.symthresh) then
           if(icount.ne.ifbo) then
              icount=ifbo
              if (iprint2.ne.0) print *,'ifbo,ib,mprim(ib) ',ifbo,ib,mprim(ib)
              ifbord(ifbo)=abs(mprim(ib))
           endif
        endif
     enddo
  enddo

  !     fb orbitals are being associated with the corresponding fd ones

  do iorb=1,norb
     do ifbo=nfborb,1,-1
        icount=0
        if (ifdord(iorb).eq.ifbord(ifbo)) then
           if (iprint2.ne.0) then
              print *,'iorb,ifbo',iorb,ifdord(iorb),ifbo,ifbord(ifbo)
           endif
           ib=0
500        continue
           ib=ib+1
           if (ib.le.npbasis) then
              if (mprim(ib).eq.0) then
                 primcoef(iorb,ib)=excoeff(ifbo,ib)
              else
                 excoeffp=excoeff(ifbo,ib)
                 excoeffm=excoeff(ifbo-1,ib)
                 primcoef(iorb,ib)=sqrt(excoeffp*excoeffp+excoeffm*excoeffm)
                 if (excoeffp.lt.0.0_PREC.or.excoeffm.lt.0.0_PREC) primcoef(iorb,ib)=-primcoef(iorb,ib)
              endif
           else
              goto 510
           endif
           goto 500
00510      continue
           if (ifbord(ifbo).gt.0) ifbord(ifbo-1)=-1
           ifbord(ifbo)=-1
           ifdord(iorb)=-2
        endif
     enddo
  enddo

  !     normalization factor for a spherical harmonic Gaussian-type functions

  do ib=1,npbasis
     d1=primexp(ib)
     l1=lprim  (ib)
     m1=mprim  (ib)
     m1abs=abs(m1)
     fngau2(ib)=( d1**dble(2*l1+3) * 2**dble(4*l1+7)/pii/(factor2(2*l1+1))**2)**0.250_PREC

     !        normalization factor for spherical harmonics

     shngau(ib)=(-1.0_PREC)**dble((m1+m1abs)/2) /sqrt(4.0_PREC*pii)* &
          sqrt((2*l1+1)*factor(l1-m1abs)/factor(l1+m1abs))
  enddo

  !     01011 format(14x,f18.9)
!01012 format(13x,a3,5x,5f10.5)
01013 format(4i5,3x,a3,2f19.9,e15.6)
01014 format(4i5,3x,a3,5f12.5,e15.6)
!01042 format(i4,i3,1x,5x,a3,5x,5f10.5)
!01052 format(2i4,5x,a3,5x,5f10.5)
!01020 format(i3,2i4,5x,f20.9,f15.5)
!01032 format(2i4,i3,a2,2x,f18.9)
!01035 format(4i4,5x,2e25.10)
01130 format(1x,'... expansion coefficient greater than 1.0 detected ...',/,1x,'    orbital:',i3, &
           ', basis function:',i4,', coefficient:',f6.3,/)
01140 format(/15x,i4,' exponents ',                   &
          /15x,i4,' primitive basis functions ',      &
          /15x,i4,' primitives on centre 1'           &
          /15x,i4,' primitives on centre 2'           &
          /15x,i4,' primitives on centre 3'/)
01141 format(6x,'centre',i2,':',i4,' exponents ')
01142 format(15x,i4,' finite difference orbitals'/15x,i4,' finite basis set orbitals'/)

end subroutine prepGaussCust



