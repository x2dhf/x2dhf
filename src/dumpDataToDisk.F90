! SPDX-License-Identifier: GPL-2.0-or-later

! Copyright (C) 1996       Leif Laaksonen, Dage Sundholm               
! Copyright (C) 1996-2023  Jacek Kobus 

module dumpDataToDisk
  implicit none
contains


  subroutine writea (iunit,ndim,a,ierr)
    use params
    implicit none
    integer (KIND=IPREC) :: ierr,iunit,ndim
    real (PREC), dimension(ndim) :: a
    write (iunit,err=1000) a
    ierr=0
    return

1000 ierr=1
  end subroutine writea

  ! ### writea32 ###
  !
  !     Writes a to a disk file in an unformatted form
  !
  subroutine writea32 (iunit,ndim,a,ierr)
    use params

    implicit none

    integer (KIND=IPREC) :: ierr,iunit,ndim
    real (PREC), dimension(ndim) :: a

    write (iunit,err=1000) a
    ierr=0
    return

1000 ierr=1

  end subroutine writea32

  ! ### writea32f ###
  !
  !     Writes a to a disk file in a formatted form
  
  subroutine writea32f (iunit,ndim,a,ierr)
    use params
    use commons

    implicit none

    integer (KIND=IPREC) :: ierr,iunit,ndim
    real (PREC), dimension(ndim) :: a
    write (iunit,formfp64,err=1000) a
    ierr=0
    return

1000 ierr=1
  end subroutine writea32f

  ! ### writea64 ###
  !
  !     Writes a to a disk file in an unformatted form
  !
  subroutine writea64 (iunit,ndim,a,ierr)
    use params

    implicit none

    integer*8 :: ierr,iunit,ndim
    real (PREC), dimension(ndim) :: a

    write (iunit,err=1000) a
    ierr=0
    return

1000 ierr=1
    return
  end subroutine writea64

  ! ### writea64f ###
  !
  !     Writes a to a disk file in a formatted form
  !
  subroutine writea64f (iunit,ndim,a,ierr)
    use params
    use commons

    implicit none

    integer*8 :: iunit, ndim, ierr
    real (PREC), dimension(ndim) :: a

    write (iunit,formfp64,err=1000) a
    ierr=0
    return

1000 ierr=1
  end subroutine writea64f

  ! ### writea128 ###
  !
  !     Writes a to a disk file in an unformatted form
  !
  subroutine writea128 (iunit,ndim,a,ierr)
    use params
    use commons

    implicit none
    integer*8 :: iunit, ndim, ierr
    real (PREC16), dimension(ndim) :: a

    write (iunit,err=1000) a
    ierr=0
    return

1000 ierr=1

  end subroutine writea128

  ! ### writea128f ###
  !
  !     Writes a to a disk file in a formatted form
  !
  subroutine writea128f (iunit,ndim,a,ierr)
    use params
    use commons
    implicit none

    integer*8 :: iunit, ndim, ierr
    real (PREC16), dimension(ndim) :: a

    write (iunit,formfp128,err=1000) a
    ierr=0
    return

1000 ierr=1
  end subroutine writea128f

  ! ### wtdisknat ###
  !
  !     Writes orbitals, potenials, Lagrange multipliers (diagonal and
  !     off-diagonal) and multipole expansion coefficients to a disk file
  !     in an unformatted form
  !
  subroutine wtdisknat (cw_orb,cw_exch)
    use params
    use discrete
    use memory
    use scfshr
    use commons
    use utils

    implicit none
    integer (KIND=IPREC) :: i,ierr,iorb1,iorb2,k

    real (PREC), dimension(*) :: cw_orb,cw_exch

    if (verboseLevel>1) then
       write(iout6,'(1x,"... saving data to disk ...")')
    endif
    
    !   write a header into the orbital output file (deprecated)
    if (outUnformatted) then
       write (iout21,err=1000) i1b,i2b,i3b,i1e,i2e,i3e,i1si,i2si,i3si,i1ng,i2ng,i3ng,i1mu,i2mu,i3mu
    else
       write (iout21,formint,err=1005) i1b,i2b,i3b,i1e,i2e,i3e,i1si,i2si,i3si,i1ng,i2ng,i3ng,i1mu,i2mu,i3mu
    endif

    !   add orbitals
    do i=1,norb
       call writea(iout21,mxsize,cw_orb(i1b(i)),ierr)
       if (ierr.ne.0) then
          write(iout6,*) 'error detected when writing orbital',i
          stop 'wtdisknat'
       endif
    enddo

    !   append the following arrays to the output orbital file
    write(iout21,err=1020) orbNorm
    write(iout21,err=1020) ee

    !write(iout21,err=1020) engo
    write(iout21,err=1020) cmulti
    write(iout21,err=1020) excdi
    write(iout21,err=1020) excqu
    write(iout21,err=1020) excoc
    write(iout21,err=1020) exche
    write(iout21,err=1020) exc5
    write(iout21,err=1020) exc6
    write(iout21,err=1020) exc7
    write(iout21,err=1020) exc8

    ! write out Coulomb potentials
    do i=1,norb
       call writea(iout22,mxsize,cw_exch(i2b(i)),ierr)
       if (ierr.ne.0) then
          write(iout6,*) 'error detected when writing coulomb potential',i
          stop 'wtdisknat'
       endif
    enddo

    ! write out exchange potentials
    if (HFinput) then
       do iorb1=1,norb
          do iorb2=iorb1,norb
             k=k2(iorb1,iorb2)
             if (iorb1.eq.iorb2.and.ll(iorb1).eq.0) cycle
             call writea(iout23,mxsize,cw_exch(i3b(k)),ierr)
             if (ierr.ne.0) then
                write(iout6,*) 'error detected when writing exchange potential',iorb1,iorb2,k
                stop 'wtdisknat'
             endif
             if (iorb1.eq.iorb2) cycle
             if (ll(iorb1).eq.0.or.ll(iorb2).eq.0) cycle
             call writea(iout23,mxsize,cw_exch(i3b(k)+mxsize),ierr)
             if (ierr.ne.0) then
                write(iout6,*) 'error detected when writing exchange potential',iorb1,iorb2,k
                stop 'wtdisknat'
             endif
          enddo
       enddo
       rewind(iout23)
    endif
    
    if (DFT.or.HFS.or.SCMC) then
       call writea(iout23,mxsize,cw_exch(length2-2*mxsize+1),ierr)
       call writea(iout23,mxsize,cw_exch(length2-mxsize+1),ierr)

       if (ierr.ne.0) then
          write(iout6,*) 'error detected when writing local exchange potential'
          stop 'wtdisknat'
       endif
    endif
    
    if (OED) then
       call writea(iout23,mxsize,cw_exch(length2-2*mxsize+1),ierr)
       if (ierr.ne.0) then
          write(iout6,*) 'error detected when writing local exchange potential'
          stop 'wtdisknat 6a'
       endif

       call writea(iout23,mxsize,cw_exch(length2-mxsize+1),ierr)
       if (ierr.ne.0) then
          write(iout6,*) 'error detected when writing local exchange potential'
          stop 'wtdisknat 6b'
       endif
    endif
    
    ! if (HF.and.lxcHyb) then
    !    call writea(iout23,mxsize,cw_exch(length2-2*mxsize+1),ierr)
    !    if (ierr.ne.0) then
    !       write(iout6,*) 'error detected when writing local exchange potential'
    !       stop 'wtdisknat 7a'
    !    endif

    !    call writea(iout23,mxsize,cw_exch(length2-mxsize+1),ierr)
    !    if (ierr.ne.0) then
    !       write(iout6,*) 'error detected when writing local exchange potential'
    !       stop 'wtdisknat 7b'
    !    endif
    ! endif

    rewind(iout21)
    rewind(iout22)
    rewind(iout23)

#ifdef PRINT
! print=121: wtdisknat: print out output orbitals
    if (iprint(121).ne.0) then
       write(*,*)
       write(*,*) 'wtdisknat: output orbitals'
       !      print orbitals
       do i=1,norb
          write(*,*) iorn(i),' ',bond(i)
          call pmtx(nni,i1mu(i),cw_orb(i1b(i)),ione,ione,incrni,incrmu)
       enddo
    endif
#endif

#ifdef PRINT
! print=122: wtdisknat: print out Coulomb potentials
    if(iprint(122).ne.0) then
       write(*,*)
       write(*,*) 'wtdisknat: output Coulomb potentials'
       !      print Coulomb potential
       do i=1,norb
          write(*,*) iorn(i),' ',bond(i)
          call pmtx(nni,i1mu(i),cw_exch(i2b(i)),ione,ione,incrni,incrmu)
       enddo
    endif
#endif
    
    !   print exchange potentials
#ifdef PRINT
! print=123: wtdisknat: print out exchange potentials
    if(iprint(123).ne.0.and.nexch.ge.1) then
       if (HFinput) then
          write(*,*)
          write(*,*) 'wtdisknat: output exchange potentials'

          do i=1,nexch
             write(*,*) 'exchange potential ',i
             call pmtx(nni,i1mu(i),cw_exch(i3b(i)),ione,ione,incrni,incrmu)
          enddo
       elseif (DFT.or.HFS.or.SCMC) then
          write(*,*)
          write(*,*) 'wtdisknat: output local exchange potential'
          call pmtx(nni,i1mu(1),cw_exch(length3-mxsize),ione,ione,incrni,incrmu)
       endif
    endif
#endif
    return

1000 continue
    write(*,*) 'wtdisknat: error encountered when writing integer arrays to disk (unformatted form)'
    stop 'wtdisknat'


1005 continue
    write(*,*) 'wtdisknat: error encountered when writing integer arrays to disk (formatted form)'
    stop 'wtdisknat'

1020 continue
    write(*,*) 'wtdisknat: error encountered when writing an extra real array to disk'
    write(*,*) '      '
    stop 'wtdisknat'

  end subroutine wtdisknat

  ! ### wtdisk32 ###
  !
  !     Writes orbitals, potenials, Lagrange multipliers and multipole
  !     expansion coefficients to a disk file in an unformatted form:
  !     integers are 4 and reals are 8-byte long.
  !
  subroutine wtdisk32 (cw_orb,cw_coul,cw_exch,i4tmp1,r8tmp1,r8tmp12,r8tmp2,r8tmp3,r8tmp4)
    use params
    use discrete
    use memory
    use scfshr
    use commons

    implicit none
    integer (KIND=IPREC) :: i,j,ierr,iorb1,iorb2,ioffset,k
    integer (KIND=IPREC),dimension(10*maxorb+5*maxorb*(maxorb+1)/2) :: i4tmp1
    real (PREC), dimension(maxorb) :: r8tmp1
    real (PREC), dimension(maxorb,maxorb) :: r8tmp12
    real (PREC), dimension(maxorb+(maxmpole-1)*maxorb) :: r8tmp2
    real (PREC), dimension(maxorb*(maxorb+1)) :: r8tmp3
    real (PREC), dimension(*) :: r8tmp4
    real (PREC), dimension(*) :: cw_orb,cw_coul,cw_exch

    if (outUnformatted) then
       write(iout6,'(/1x,"... writing functions to disk (i32) ...")')
    else
       write(iout6,'(/1x,"... writing functions to disk (i32f) ...")')
    endif
    
    ioffset=0
    do i=1,maxorb
       i4tmp1(ioffset+i)=i1b(i)
    enddo
    ioffset=ioffset+maxorb

    do i=1,maxorb
       i4tmp1(ioffset+i)=i2b(i)
    enddo
    ioffset=ioffset+maxorb

    do i=1,maxorb*(maxorb+1)/2
       i4tmp1(ioffset+i)=i3b(i)
    enddo
    ioffset=ioffset+maxorb*(maxorb+1)/2

    do i=1,maxorb
       i4tmp1(ioffset+i)=i1e(i)
    enddo
    ioffset=ioffset+maxorb

    do i=1,maxorb
       i4tmp1(ioffset+i)=i2e(i)
    enddo
    ioffset=ioffset+maxorb

    do i=1,maxorb*(maxorb+1)/2
       i4tmp1(ioffset+i)=i3e(i)
    enddo
    ioffset=ioffset+maxorb*(maxorb+1)/2

    do i=1,maxorb
       i4tmp1(ioffset+i)=i1si(i)
    enddo
    ioffset=ioffset+maxorb

    do i=1,maxorb
       i4tmp1(ioffset+i)=i2si(i)
    enddo
    ioffset=ioffset+maxorb

    do i=1,maxorb*(maxorb+1)/2
       i4tmp1(ioffset+i)=i3si(i)
    enddo
    ioffset=ioffset+maxorb*(maxorb+1)/2

    do i=1,maxorb
       i4tmp1(ioffset+i)=i1ng(i)
    enddo
    ioffset=ioffset+maxorb


    do i=1,maxorb
       i4tmp1(ioffset+i)=i2ng(i)
    enddo
    ioffset=ioffset+maxorb

    do i=1,maxorb*(maxorb+1)/2
       i4tmp1(ioffset+i)=i3ng(i)
    enddo
    ioffset=ioffset+maxorb*(maxorb+1)/2

    do i=1,maxorb
       i4tmp1(ioffset+i)=i1mu(i)
    enddo
    ioffset=ioffset+maxorb

    do i=1,maxorb
       i4tmp1(ioffset+i)=i2mu(i)
    enddo
    ioffset=ioffset+maxorb

    do i=1,maxorb*(maxorb+1)/2
       i4tmp1(ioffset+i)=i3mu(i)
    enddo
    ioffset=ioffset+maxorb*(maxorb+1)/2

    if (outUnformatted) then
       write (iout21,err=1000) (i4tmp1(i),i=1,ioffset)
    else
       write (iout21,formint,err=1000) (i4tmp1(i),i=1,ioffset)
    endif
    
    ! add orbitals
    do i=1,norb
       do j=1,mxsize
          r8tmp4(j)=cw_orb(i1b(i)+j-1)
       enddo
       if (outUnformatted) then
          call writea32(iout21,mxsize,r8tmp4,ierr)
       else
          call writea32f(iout21,mxsize,r8tmp4,ierr)
       endif
       if (ierr.ne.0) then
          write(iout6,*) 'error detected when writing orbital',i
          stop 'wtdisk32'
       endif
    enddo

    ! append the following arrays to the output orbital file
    do i=1,maxorb
       r8tmp1(i)=orbNorm(i)
    enddo
    if (outUnformatted) then
       write(iout21,err=1020) (r8tmp1(i),i=1,maxorb)
    else
       write(iout21,formfp64,err=1020) (r8tmp1(i),i=1,maxorb)
    endif
    
    do i=1,maxorb
       do j=1,maxorb
          r8tmp12(i,j)=ee(i,j)
       enddo
    enddo
    if (outUnformatted) then
       write(iout21,err=1020) ((r8tmp12(i,j),i=1,maxorb),j=1,maxorb)
    else
       write(iout21,formfp64,err=1020) ((r8tmp12(i,j),i=1,maxorb),j=1,maxorb)
    endif

    do i=1,maxorb+(maxmpole-1)*maxorb
       r8tmp2(i)=cmulti(i)
    enddo
    if (outUnformatted) then
       write(iout21,err=1020) (r8tmp2(i),i=1,maxorb+(maxmpole-1)*maxorb)
    else
       write(iout21,formfp64,err=1020) (r8tmp2(i),i=1,maxorb+(maxmpole-1)*maxorb)
    endif
   
    do i=1,maxorb*(maxorb+1)
       r8tmp3(i)=excdi(i)
    enddo
    if (outUnformatted) then
       write(iout21,err=1020) (r8tmp3(i),i=1,maxorb*(maxorb+1))
    else
       write(iout21,formfp64,err=1020) (r8tmp3(i),i=1,maxorb*(maxorb+1))
    endif

    do i=1,maxorb*(maxorb+1)
       r8tmp3(i)=excqu(i)
    enddo
    if (outUnformatted) then
       write(iout21,err=1020) (r8tmp3(i),i=1,maxorb*(maxorb+1))
    else
       write(iout21,formfp64,err=1020) (r8tmp3(i),i=1,maxorb*(maxorb+1))
    endif

    do i=1,maxorb*(maxorb+1)
       r8tmp3(i)=excoc(i)
    enddo
    if (outUnformatted) then
       write(iout21,err=1020) (r8tmp3(i),i=1,maxorb*(maxorb+1))
    else
       write(iout21,formfp64,err=1020) (r8tmp3(i),i=1,maxorb*(maxorb+1))
    endif

    do i=1,maxorb*(maxorb+1)
       r8tmp3(i)=exche(i)
    enddo
    if (outUnformatted) then
       write(iout21,err=1020) (r8tmp3(i),i=1,maxorb*(maxorb+1))
    else
       write(iout21,formfp64,err=1020) (r8tmp3(i),i=1,maxorb*(maxorb+1))
    endif

    do i=1,maxorb*(maxorb+1)
       r8tmp3(i)=exc5(i)
    enddo
    if (outUnformatted) then
       write(iout21,err=1020) (r8tmp3(i),i=1,maxorb*(maxorb+1))
    else
       write(iout21,formfp64,err=1020) (r8tmp3(i),i=1,maxorb*(maxorb+1))
    endif

    do i=1,maxorb*(maxorb+1)
       r8tmp3(i)=exc6(i)
    enddo
    if (outUnformatted) then
       write(iout21,err=1020) (r8tmp3(i),i=1,maxorb*(maxorb+1))
    else
       write(iout21,formfp64,err=1020) (r8tmp3(i),i=1,maxorb*(maxorb+1))
    endif

    do i=1,maxorb*(maxorb+1)
       r8tmp3(i)=exc7(i)
    enddo
    if (outUnformatted) then
       write(iout21,err=1020) (r8tmp3(i),i=1,maxorb*(maxorb+1))
    else
       write(iout21,formfp64,err=1020) (r8tmp3(i),i=1,maxorb*(maxorb+1))
    endif

    do i=1,maxorb*(maxorb+1)
       r8tmp3(i)=exc8(i)
    enddo
    if (outUnformatted) then
       write(iout21,err=1020) (r8tmp3(i),i=1,maxorb*(maxorb+1))
    else
       write(iout21,formfp64,err=1020) (r8tmp3(i),i=1,maxorb*(maxorb+1))
    endif

    do i=1,norb
       do j=1,mxsize
          r8tmp4(j)=cw_coul(i2b(i)+j-1)
       enddo
       if (outUnformatted) then
          call writea32(iout22,mxsize,r8tmp4,ierr)
       else
          call writea32f(iout22,mxsize,r8tmp4,ierr)
       endif
       if (ierr.ne.0) then
          write(iout6,*) 'error detected when writing coulomb potential',i
          stop 'wtdisk32'
       endif
    enddo

    !    write out exchange potentials
    if (HF) then
       do iorb1=1,norb
          do iorb2=iorb1,norb
             k=k2(iorb1,iorb2)
             if (iorb1.eq.iorb2.and.ll(iorb1).eq.0) cycle
             
             do j=1,mxsize
                r8tmp4(j)=cw_exch(i3b(k)+j-1)
             enddo
             if (outUnformatted) then
                call writea32(iout23,mxsize,r8tmp4,ierr)
             else
                call writea32f(iout23,mxsize,r8tmp4,ierr)
             endif
             if (ierr.ne.0) then
                write(iout6,*) 'error detected when writing exchange potential',iorb1,iorb2,k
                stop 'wtdisk32'
             endif
             
             if (iorb1.eq.iorb2) cycle
             if (ll(iorb1).eq.0.or.ll(iorb2).eq.0) cycle
             
             do j=1,mxsize
                r8tmp4(j)=cw_exch(i3b(k)+mxsize+j-1)
             enddo
             if (outUnformatted) then
                call writea32(iout23,mxsize,r8tmp4,ierr)
             else
                call writea32f(iout23,mxsize,r8tmp4,ierr)
             endif

             if (ierr.ne.0) then
                write(iout6,*) 'error detected when writing exchange potential',iorb1,iorb2,k
                stop 'wtdisk32'
             endif
          enddo
       enddo
       rewind(iout23)
    endif
    
    if (DFT.or.HFS.or.SCMC) then
       do j=1,mxsize
          r8tmp4(j)=cw_exch(length3-mxsize+j-1)
       enddo
       if (outUnformatted) then
          call writea32(iout23,mxsize,r8tmp4,ierr)
       else
          call writea32f(iout23,mxsize,r8tmp4,ierr)
       endif
       if (ierr.ne.0) then
          write(iout6,*) 'error detected when writing local exchange potential'
          stop 'wtdisk32'
       endif
    endif
    
    rewind(iout21)
    rewind(iout22)
    rewind(iout23)

    return

1000 continue
    write(iout6,1070)
    stop 'wtdisk32'

1020 continue
    write(*,*) 'wtdisk32: error encountered when writing an extension to orbital file'
    write(*,*) '      '
    stop 'wtdisk32'
1070 format(//1x,'error! can not write data to disk'//)
    !1080 format(//1x,'error! can not write data to disk'//)
  end subroutine wtdisk32

  ! ### wtdisk64 ###
  !
  !     Writes orbitals, potenials, Lagrange multipliers and multipole
  !     expansion coefficients to a disk file in an unformatted form:
  !     (integers and reals are 8-byte long.
  !
  subroutine wtdisk64 (cw_orb,cw_coul,cw_exch,i8tmp1,r8tmp1,r8tmp12,r8tmp2,r8tmp3,r8tmp4)
    use params
    use discrete
    use memory
    use scfshr
    use commons
    implicit none

    integer (KIND=IPREC) :: i,ierr,iorb1,iorb2,ioffset,k
    integer*8 :: ierr8,j64,mxsize64
    integer*8,dimension(10*maxorb+5*maxorb*(maxorb+1)/2) :: i8tmp1
    real (PREC), dimension(maxorb) :: r8tmp1
    real (PREC), dimension(maxorb,maxorb) :: r8tmp12    
    real (PREC), dimension(maxorb+(maxmpole-1)*maxorb) :: r8tmp2
    real (PREC), dimension(maxorb*(maxorb+1)) :: r8tmp3
    real (PREC), dimension(*) :: r8tmp4
    real (PREC), dimension(*) :: cw_orb,cw_coul,cw_exch

    mxsize64=mxsize

    if (outUnformatted) then
       write(iout6,'(/1x,"... writing functions to disk (i64) ...")')
    else
       write(iout6,'(/1x,"... writing functions to disk (i64f) ...")')
    endif
    
    ! write a header into the orbital output file
    ioffset=0
    do i=1,maxorb
       i8tmp1(ioffset+i)=i1b(i)
    enddo
    ioffset=ioffset+maxorb

    do i=1,maxorb
       i8tmp1(ioffset+i)=i2b(i)
    enddo
    ioffset=ioffset+maxorb

    do i=1,maxorb*(maxorb+1)/2
       i8tmp1(ioffset+i)=i3b(i)
    enddo
    ioffset=ioffset+maxorb*(maxorb+1)/2

    do i=1,maxorb
       i8tmp1(ioffset+i)=i1e(i)
    enddo
    ioffset=ioffset+maxorb

    do i=1,maxorb
       i8tmp1(ioffset+i)=i2e(i)
    enddo
    ioffset=ioffset+maxorb

    do i=1,maxorb*(maxorb+1)/2
       i8tmp1(ioffset+i)=i3e(i)
    enddo
    ioffset=ioffset+maxorb*(maxorb+1)/2

    do i=1,maxorb
       i8tmp1(ioffset+i)=i1si(i)
    enddo
    ioffset=ioffset+maxorb

    do i=1,maxorb
       i8tmp1(ioffset+i)=i2si(i)
    enddo
    ioffset=ioffset+maxorb

    do i=1,maxorb*(maxorb+1)/2
       i8tmp1(ioffset+i)=i3si(i)
    enddo
    ioffset=ioffset+maxorb*(maxorb+1)/2

    do i=1,maxorb
       i8tmp1(ioffset+i)=i1ng(i)
    enddo
    ioffset=ioffset+maxorb

    do i=1,maxorb
       i8tmp1(ioffset+i)=i2ng(i)
    enddo
    ioffset=ioffset+maxorb

    do i=1,maxorb*(maxorb+1)/2
       i8tmp1(ioffset+i)=i3ng(i)
    enddo
    ioffset=ioffset+maxorb*(maxorb+1)/2

    do i=1,maxorb
       i8tmp1(ioffset+i)=i1mu(i)
    enddo
    ioffset=ioffset+maxorb

    do i=1,maxorb
       i8tmp1(ioffset+i)=i2mu(i)
    enddo
    ioffset=ioffset+maxorb

    do i=1,maxorb*(maxorb+1)/2
       i8tmp1(ioffset+i)=i3mu(i)
    enddo
    ioffset=ioffset+maxorb*(maxorb+1)/2

    if (outUnformatted) then
       write (i8out21,err=1000) (i8tmp1(i),i=1,ioffset)
    else
       write (i8out21,formint,err=1000) (i8tmp1(i),i=1,ioffset)
    endif

    !   add orbitals
    do i=1,norb
       do j64=1,mxsize64
          r8tmp4(j64)=cw_orb(i1b(i)+j64-1)
       enddo
       if (outUnformatted) then
          call writea64(i8out21,mxsize64,r8tmp4,ierr8)
       else
          call writea64f(i8out21,mxsize64,r8tmp4,ierr8)
       endif
       if (ierr8.ne.0) then
          write(iout6,*) 'error detected when writing orbital',i
          stop 'wtdisk64'
       endif
    enddo

    !   append the following arrays to the output orbital file
    do i=1,maxorb
       r8tmp1(i)=orbNorm(i)
    enddo
    if (outUnformatted) then
       write(i8out21,err=1020) (r8tmp1(i),i=1,maxorb)
    else
       write(i8out21,formfp64,err=1020) (r8tmp1(i),i=1,maxorb)
    endif
    
    do i=1,maxorb
       r8tmp1(i)=ee(i,i)
    enddo
    if (outUnformatted) then
       write(i8out21,err=1020) (r8tmp1(i),i=1,maxorb)
    else
       write(i8out21,formfp64,err=1020) (r8tmp1(i),i=1,maxorb)
    endif

    do i=1,maxorb
       r8tmp3(i)=ee(i,i)
    enddo
    if (outUnformatted) then
       write(i8out21,err=1020) (r8tmp3(i),i=1,maxorb*(maxorb+1))
    else
       write(i8out21,formfp64,err=1020) (r8tmp3(i),i=1,maxorb*maxorb)
    endif

    do i=1,maxorb+(maxmpole-1)*maxorb
       r8tmp2(i)=cmulti(i)
    enddo
    if (outUnformatted) then
       write(i8out21,err=1020) (r8tmp2(i),i=1,maxorb+(maxmpole-1)*maxorb)
    else
       write(i8out21,formfp64,err=1020) (r8tmp2(i),i=1,maxorb+(maxmpole-1)*maxorb)
    endif
    
    do i=1,maxorb*(maxorb+1)
       r8tmp3(i)=excdi(i)
    enddo
    if (outUnformatted) then
       write(i8out21,err=1020) (r8tmp3(i),i=1,maxorb*(maxorb+1))
    else
       write(i8out21,formfp64,err=1020) (r8tmp3(i),i=1,maxorb*(maxorb+1))
    endif
    
    do i=1,maxorb*(maxorb+1)
       r8tmp3(i)=excqu(i)
    enddo
    if (outUnformatted) then
       write(i8out21,err=1020) (r8tmp3(i),i=1,maxorb*(maxorb+1))
    else
       write(i8out21,formfp64,err=1020) (r8tmp3(i),i=1,maxorb*(maxorb+1))
    endif

    do i=1,maxorb*(maxorb+1)
       r8tmp3(i)=excoc(i)
    enddo
    if (outUnformatted) then
       write(i8out21,err=1020) (r8tmp3(i),i=1,maxorb*(maxorb+1))
    else
       write(i8out21,formfp64,err=1020) (r8tmp3(i),i=1,maxorb*(maxorb+1))
    endif

    do i=1,maxorb*(maxorb+1)
       r8tmp3(i)=exche(i)
    enddo
    if (outUnformatted) then
       write(i8out21,err=1020) (r8tmp3(i),i=1,maxorb*(maxorb+1))
    else
       write(i8out21,formfp64,err=1020) (r8tmp3(i),i=1,maxorb*(maxorb+1))
    endif

    do i=1,maxorb*(maxorb+1)
       r8tmp3(i)=exc5(i)
    enddo
    if (outUnformatted) then
       write(i8out21,err=1020) (r8tmp3(i),i=1,maxorb*(maxorb+1))
    else
       write(i8out21,formfp64,err=1020) (r8tmp3(i),i=1,maxorb*(maxorb+1))
    endif

    do i=1,maxorb*(maxorb+1)
       r8tmp3(i)=exc6(i)
    enddo
    if (outUnformatted) then
       write(i8out21,err=1020) (r8tmp3(i),i=1,maxorb*(maxorb+1))
    else
       write(i8out21,formfp64,err=1020) (r8tmp3(i),i=1,maxorb*(maxorb+1))
    endif

    do i=1,maxorb*(maxorb+1)
       r8tmp3(i)=exc7(i)
    enddo
    if (outUnformatted) then
       write(i8out21,err=1020) (r8tmp3(i),i=1,maxorb*(maxorb+1))
    else
       write(i8out21,formfp64,err=1020) (r8tmp3(i),i=1,maxorb*(maxorb+1))
    endif

    do i=1,maxorb*(maxorb+1)
       r8tmp3(i)=exc8(i)
    enddo
    if (outUnformatted) then
       write(i8out21,err=1020) (r8tmp3(i),i=1,maxorb*(maxorb+1))
    else
       write(i8out21,formfp64,err=1020) (r8tmp3(i),i=1,maxorb*(maxorb+1))
    endif

    ! write out Coulomb potentials
    do i=1,norb
       do j64=1,mxsize64
          r8tmp4(j64)=cw_coul(i2b(i)+j64-1)
       enddo
       if (outUnformatted) then
          call writea64(i8out22,mxsize64,r8tmp4,ierr8)
       else
          call writea64f(i8out22,mxsize64,r8tmp4,ierr8)
       endif
       if (ierr8.ne.0) then
          write(iout6,*) 'error detected when writing coulomb potential',i
          stop 'wtdisk64'
       endif
    enddo

    !   write out exchange potentials
    if (HF) then
       do iorb1=1,norb
          do iorb2=iorb1,norb
             k=k2(iorb1,iorb2)
             if (iorb1.eq.iorb2.and.ll(iorb1).eq.0) cycle
             
             do j64=1,mxsize64
                r8tmp4(j64)=cw_exch(i3b(k)+j64-1)
             enddo
             if (outUnformatted) then             
                call writea64(i8out23,mxsize64,r8tmp4,ierr8)
             else
                call writea64f(i8out23,mxsize64,r8tmp4,ierr8)
             endif
             if (ierr8.ne.0) then
                write(iout6,*) 'error detected when writing exchange potential',iorb1,iorb2,k
                stop 'wtdisk64'
             endif
             if (iorb1.eq.iorb2) cycle
             if (ll(iorb1).eq.0.or.ll(iorb2).eq.0) cycle
             
             do j64=1,mxsize64
                r8tmp4(j64)=cw_exch(i3b(k)+mxsize64+j64-1)
             enddo
             if (outUnformatted) then             
                call writea64(i8out23,mxsize64,r8tmp4,ierr8)
             else
                call writea64f(i8out23,mxsize64,r8tmp4,ierr8)
             endif

             if (ierr8.ne.0) then
                write(iout6,*) 'error detected when writing exchange potential',iorb1,iorb2,k
                stop 'wtdisk64'
             endif
          enddo
       enddo
       rewind(i8out23)
    endif

    if (DFT.or.HFS.or.OED) then
       do j64=1,mxsize64
          r8tmp4(j64)=cw_exch(length3-mxsize64+j64-1)
       enddo
       if (outUnformatted) then             
          call writea64(i8out23,mxsize64,r8tmp4,ierr8)
       else
          call writea64f(i8out23,mxsize64,r8tmp4,ierr8)
       endif
       if (ierr8.ne.0) then
          write(iout6,*) 'error detected when writing local exchange potential'
          stop 'wtdisk64'
       endif
    endif
    
    rewind(i8out21)
    rewind(i8out22)
    rewind(i8out23)

    return

1000 continue
    write(iout6,1070)
    stop 'wtdisk64'

1020 continue
    write(*,*) 'wtdisk64: error encountered when writing an extension to orbital file'
    write(*,*) '      '
    stop 'wtdisk64'
1070 format(//1x,'error! can not write data to disk'//)
    ! 1080 format(//1x,'error! can not write data to disk'//)
  end subroutine wtdisk64

  ! ### wtdisk128 ###
  !
  !     Writes orbitals, potenials, Lagrange multipliers and multipole
  !     expansion coefficients to a disk file in an unformatted form
  !     (integers are 4 and reals are 8-byte long)
  !
  subroutine wtdisk128 (cw_orb,cw_coul,cw_exch,i8tmp1,r16tmp1,r16tmp12,r16tmp2,r16tmp3,r16tmp4)
    use params
    use discrete
    use memory
    use scfshr
    use commons

    implicit none
    integer (KIND=IPREC) :: i,j,iorb1,iorb2,ioffset,k
    integer*8 :: ierr8,j64,mxsize64
    integer*8,dimension(10*maxorb+5*maxorb*(maxorb+1)/2) :: i8tmp1    

    real (PREC16), dimension(maxorb) :: r16tmp1
    real (PREC16), dimension(maxorb,maxorb) :: r16tmp12    
    real (PREC16), dimension(maxorb+(maxmpole-1)*maxorb) :: r16tmp2
    real (PREC16), dimension(maxorb*(maxorb+1)) :: r16tmp3
    real (PREC16), dimension(*) :: r16tmp4
    real (PREC), dimension(*) :: cw_orb,cw_coul,cw_exch

    mxsize64=mxsize

    if (outUnformatted) then
       write(iout6,'(/1x,"... writing functions to disk (r128) ...")')
    else
       write(iout6,'(/1x,"... writing functions to disk (r128f) ...")')
   endif
   
    ioffset=0
    do i=1,maxorb
       i8tmp1(ioffset+i)=i1b(i)
    enddo
    ioffset=ioffset+maxorb

    do i=1,maxorb
       i8tmp1(ioffset+i)=i2b(i)
    enddo
    ioffset=ioffset+maxorb

    do i=1,maxorb*(maxorb+1)/2
       i8tmp1(ioffset+i)=i3b(i)
    enddo
    ioffset=ioffset+maxorb*(maxorb+1)/2

    do i=1,maxorb
       i8tmp1(ioffset+i)=i1e(i)
    enddo
    ioffset=ioffset+maxorb

    do i=1,maxorb
       i8tmp1(ioffset+i)=i2e(i)
    enddo
    ioffset=ioffset+maxorb

    do i=1,maxorb*(maxorb+1)/2
       i8tmp1(ioffset+i)=i3e(i)
    enddo
    ioffset=ioffset+maxorb*(maxorb+1)/2

    do i=1,maxorb
       i8tmp1(ioffset+i)=i1si(i)
    enddo
    ioffset=ioffset+maxorb

    do i=1,maxorb
       i8tmp1(ioffset+i)=i2si(i)
    enddo
    ioffset=ioffset+maxorb

    do i=1,maxorb*(maxorb+1)/2
       i8tmp1(ioffset+i)=i3si(i)
    enddo
    ioffset=ioffset+maxorb*(maxorb+1)/2

    do i=1,maxorb
       i8tmp1(ioffset+i)=i1ng(i)
    enddo
    ioffset=ioffset+maxorb

    do i=1,maxorb
       i8tmp1(ioffset+i)=i2ng(i)
    enddo
    ioffset=ioffset+maxorb

    do i=1,maxorb*(maxorb+1)/2
       i8tmp1(ioffset+i)=i3ng(i)
    enddo
    ioffset=ioffset+maxorb*(maxorb+1)/2

    do i=1,maxorb
       i8tmp1(ioffset+i)=i1mu(i)
    enddo
    ioffset=ioffset+maxorb

    do i=1,maxorb
       i8tmp1(ioffset+i)=i2mu(i)
    enddo
    ioffset=ioffset+maxorb

    do i=1,maxorb*(maxorb+1)/2
       i8tmp1(ioffset+i)=i3mu(i)
    enddo
    ioffset=ioffset+maxorb*(maxorb+1)/2

    if (outUnformatted) then
       write (i8out21,err=1000) (i8tmp1(i),i=1,ioffset)
    else
       write (i8out21,formint,err=1000) (i8tmp1(i),i=1,ioffset)
    endif
    
    !   add orbitals
    do i=1,norb
       do j64=1,mxsize64
          r16tmp4(j64)=cw_orb(i1b(i)+j64-1)
       enddo
       if (outUnformatted) then
          call writea128(i8out21,mxsize64,r16tmp4,ierr8)
       else
          call writea128f(i8out21,mxsize64,r16tmp4,ierr8)
       endif
       if (ierr8.ne.0) then
          write(iout6,*) 'error detected when writing orbital',i
          stop 'wtdisk128'
       endif
    enddo

    !   append the following arrays to the output orbital file
    do i=1,maxorb
       r16tmp1(i)=orbNorm(i)
    enddo
    if (outUnformatted) then
       write(i8out21,err=1020) (r16tmp1(i),i=1,maxorb)
    else
       write(i8out21,formfp128,err=1020) (r16tmp1(i),i=1,maxorb)
    endif
    
    do i=1,maxorb
       do j=1,maxorb
          r16tmp12(i,j)=ee(i,j)
       enddo
    enddo
    if (outUnformatted) then
       write(i8out21,err=1020) ((r16tmp12(i,j),i=1,maxorb),j=1,maxorb)
    else
       write(i8out21,formfp128,err=1020) ((r16tmp12(i,j),i=1,maxorb),j=1,maxorb)
    endif

    do i=1,maxorb+(maxmpole-1)*maxorb
       r16tmp2(i)=cmulti(i)
    enddo
    if (outUnformatted) then
       write(i8out21,err=1020) (r16tmp2(i),i=1,maxorb+(maxmpole-1)*maxorb)
    else
       write(i8out21,formfp128,err=1020) (r16tmp2(i),i=1,maxorb+(maxmpole-1)*maxorb)
    endif

    do i=1,maxorb*(maxorb+1)
       r16tmp3(i)=excdi(i)
    enddo
    if (outUnformatted) then
       write(i8out21,err=1020) (r16tmp3(i),i=1,maxorb*(maxorb+1))
    else
       write(i8out21,formfp128,err=1020) (r16tmp3(i),i=1,maxorb*(maxorb+1))
    endif
    
    do i=1,maxorb*(maxorb+1)
       r16tmp3(i)=excqu(i)
    enddo
    if (outUnformatted) then
       write(i8out21,err=1020) (r16tmp3(i),i=1,maxorb*(maxorb+1))
    else
       write(i8out21,formfp128,err=1020) (r16tmp3(i),i=1,maxorb*(maxorb+1))
    endif

    do i=1,maxorb*(maxorb+1)
       r16tmp3(i)=excoc(i)
    enddo
    if (outUnformatted) then
       write(i8out21,err=1020) (r16tmp3(i),i=1,maxorb*(maxorb+1))
    else
       write(i8out21,formfp128,err=1020) (r16tmp3(i),i=1,maxorb*(maxorb+1))
    endif

    do i=1,maxorb*(maxorb+1)
       r16tmp3(i)=exche(i)
    enddo
    if (outUnformatted) then
       write(i8out21,err=1020) (r16tmp3(i),i=1,maxorb*(maxorb+1))
    else
       write(i8out21,formfp128,err=1020) (r16tmp3(i),i=1,maxorb*(maxorb+1))
    endif

    do i=1,maxorb*(maxorb+1)
       r16tmp3(i)=exc5(i)
    enddo
    if (outUnformatted) then
       write(i8out21,err=1020) (r16tmp3(i),i=1,maxorb*(maxorb+1))
    else
       write(i8out21,formfp128,err=1020) (r16tmp3(i),i=1,maxorb*(maxorb+1))
    endif

    do i=1,maxorb*(maxorb+1)
       r16tmp3(i)=exc6(i)
    enddo
    if (outUnformatted) then
       write(i8out21,err=1020) (r16tmp3(i),i=1,maxorb*(maxorb+1))
    else
       write(i8out21,formfp128,err=1020) (r16tmp3(i),i=1,maxorb*(maxorb+1))
    endif

    do i=1,maxorb*(maxorb+1)
       r16tmp3(i)=exc7(i)
    enddo
    if (outUnformatted) then
       write(i8out21,err=1020) (r16tmp3(i),i=1,maxorb*(maxorb+1))
    else
       write(i8out21,formfp128,err=1020) (r16tmp3(i),i=1,maxorb*(maxorb+1))
    endif

    do i=1,maxorb*(maxorb+1)
       r16tmp3(i)=exc8(i)
    enddo
    if (outUnformatted) then
       write(i8out21,err=1020) (r16tmp3(i),i=1,maxorb*(maxorb+1))
    else
       write(i8out21,formfp128,err=1020) (r16tmp3(i),i=1,maxorb*(maxorb+1))
    endif

    !   write out Coulomb potentials
    do i=1,norb
       do j64=1,mxsize64
          r16tmp4(j64)=cw_coul(i2b(i)+j64-1)
       enddo
       if (outUnformatted) then
          call writea128(i8out22,mxsize64,r16tmp4,ierr8)
       else
          call writea128f(i8out22,mxsize64,r16tmp4,ierr8)
       endif
       if (ierr8.ne.0) then
          write(iout6,*) 'error detected when writing coulomb potential',i
          stop 'wtdisk128'
       endif
    enddo

    !   write out exchange potentials
    if (HF) then
       do iorb1=1,norb
          do iorb2=iorb1,norb
             k=k2(iorb1,iorb2)
             if (iorb1.eq.iorb2.and.ll(iorb1).eq.0) cycle
             
             do j64=1,mxsize64
                r16tmp4(j64)=cw_exch(i3b(k)+j64-1)
             enddo
             if (outUnformatted) then
                call writea128(i8out23,mxsize64,r16tmp4,ierr8)
             else
                call writea128f(i8out23,mxsize64,r16tmp4,ierr8)
             endif
             if (ierr8.ne.0) then
                write(iout6,*) 'error detected when writing exchange potential',iorb1,iorb2,k
                stop 'wtdisk128'
             endif
             
             if (iorb1.eq.iorb2) cycle
             if (ll(iorb1).eq.0.or.ll(iorb2).eq.0) cycle
             
             do j64=1,mxsize64
                r16tmp4(j64)=cw_exch(i3b(k)+i3si(k)+j64-1)
             enddo
             if (outUnformatted) then
                call writea128(i8out23,mxsize64,r16tmp4,ierr8)
             else
                call writea128f(i8out23,mxsize64,r16tmp4,ierr8)
             endif
             if (ierr8.ne.0) then
                write(iout6,*) 'error detected when writing exchange potential',iorb1,iorb2,k
                stop 'wtdisk128'
             endif
          enddo
       enddo
       rewind(i8out23)
    endif

    if (DFT.or.HFS.or.SCMC) then
       do j64=1,mxsize64
          r16tmp4(j64)=cw_exch(length3-mxsize64+j64-1)
       enddo
       if (outUnformatted) then
          call writea128(i8out23,mxsize64,r16tmp4,ierr8)
       else
          call writea128f(i8out23,mxsize64,r16tmp4,ierr8)
       endif
       if (ierr8.ne.0) then
          write(iout6,*) 'error detected when writing local exchange potential'
          stop 'wtdisk128'
       endif
    endif
    
    rewind(i8out21)
    rewind(i8out22)
    rewind(i8out23)

    return

1000 continue
    write(iout6,1070)
    stop 'wtdisk128'

1020 continue
    write(*,*) 'wtdiski128: error encountered when writing an extension to orbital file'
    write(*,*) '      '
    stop 'wtdiski128'
1070 format(//1x,'error! can not write data to disk'//)
    !1080 format(//1x,'error! can not write data to disk'//)
  end subroutine wtdisk128

end module dumpDataToDisk
