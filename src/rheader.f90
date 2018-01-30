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
! c ### rheader ###

subroutine rheader(norb_p)
  use params
  use discret
  use scf
  use commons8

  implicit none
  character*80 :: datetime,datetime_p,header_p
  integer :: i,ig,imax1,imax2,mismatch,nel_p,nexch_p,norb_p
  real (PREC) :: rd,z1d,z2d
  integer i4tmp1,i4tmp2,i4tmp3
  integer, dimension(10) :: i4tmp

  integer*8 i8tmp1,i8tmp2,i8tmp3
  integer*8, dimension(10) :: i8tmp

  real (PREC) r8tmp1,r8tmp2
  real (PREC), dimension(10) :: r8tmp

  real (PREC16) :: r16tmp1,r16tmp2
  real (PREC16), dimension(10) :: r16tmp
  
  
  mismatch=0
  idump=0
  imax1=34
  imax2=24
  
  !     check if 2dhf_input.dat is available

  if (idat.eq.0) then
     rewind(iinp14)
     read(iinp14,'(a80)') header_p
     read(iinp14,'(a80)') datetime_p
     read(iinp14,formint) ngrids_p,nni_p,nmu_p   

     read(iinp14,formfp64) r_p,rgrid_p
     read(iinp14,formfp64) z1_p,z2_p	
     read(iinp14,formint) norb_p,nel_p,nexch_p
     
     rd =abs(r_p-r)
     z1d=abs(z1_p-z1)
     z2d=abs(z2_p-z2)
     
     write(iout6,1050)
!     write(iout6,1060) (dtarr1(i),i=1,imax1)
!     write(iout6,1061) (dtarr2(i),i=1,imax2)

     write(iout6,1080) header_p
     write(iout6,1081) datetime_p


  else 
     !     check default length of integer variables used in a file
     !      intlen=checkintlen()
     
     read(iinp11,err=1000) header_p
     read(iinp11,err=1000) datetime_p
     
     write(iout6,1050)
!     write(iout6,1060) (dtarr1(i),i=1,imax1)
!     write(iout6,1061) (dtarr2(i),i=1,imax2)

     write(iout6,1080) header_p
     write(iout6,1081) datetime_p
     

     if (lengthint.eq.4) then
        read(iinp11,err=1000) i4tmp1,i4tmp2,i4tmp
        ngrids_p=i4tmp1
        nni_p   =i4tmp2 
        do i=1,10
           nmu_p(i)=i4tmp(i)
        enddo
        
     else
        read(iinp11,err=1000) i8tmp1,i8tmp2,i8tmp
        ngrids_p=i8tmp1
        nni_p   =i8tmp2 
        do i=1,10
           nmu_p(i)=i8tmp(i)
        enddo
     endif
     
     if (lengthfp.eq.8) then
        read(iinp11,err=1000) r8tmp1,r8tmp
        r_p=r8tmp1
        do i=1,10
           rgrid_p(i)=r8tmp(i)
        enddo
     else
        read(iinp11,err=1000) r16tmp1,r16tmp
        r_p=r16tmp1
        do i=1,10
           rgrid_p(i)=r16tmp(i)
        enddo
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
     
     rd =abs(r_p-r)
     z1d=abs(z1_p-z1)
     z2d=abs(z2_p-z2)
     
     if (lengthint.eq.4) then
        read(iinp11,err=1000) i4tmp1,i4tmp2,i4tmp3
        norb_p =i4tmp1
        nel_p  =i4tmp2
        nexch_p=i4tmp3
        nexch  =nexch_p
     else
        read(iinp11,err=1000) i8tmp1,i8tmp2,i8tmp3
        norb_p =i8tmp1
        nel_p  =i8tmp2
        nexch_p=i8tmp3
        nexch  =nexch_p
     endif
  endif
  
  !     test for various parameters
  
  if (iinterp.eq.0) then
     if (ngrids.ne.ngrids_p ) then
        write (iout6,*) 'mismatch in ngrids: '
        write (iout6,*) 'input=',ngrids,' file=',ngrids_p
        stop 'rheader'
     endif
     
     if (nni.ne.nni_p) then
        write (iout6,*) 'mismatch in nni: '
        write (iout6,*) 'input=',nni,' file=',nni_p
        stop 'rheader'
     endif
     
     do ig=1,ngrids   
        if (nmu(ig).ne.nmu_p(ig)) then
           write (iout6,*) 'mismatch in nmu for grid ',ig,':'
           write (iout6,*) 'input=',nmu(ig),' file=',nmu_p(ig)
           stop 'rheader'
        endif
     enddo
  endif
  
  
  if (rd.gt.1.d-06) then
     write(iout6,*) 'mismatch in bond length:'
     write(iout6,*) 'input=',r, 'disk=',r_p
     
  elseif (z1d.gt.1.d-06) then
     write(iout6,*) 'Warning! mismatch in Z1:'
     write(iout6,*) 'input=',z1, 'disk=',z1_p
     
  elseif (z2d.gt.1.d-06) then
     write(iout6,*) 'Warning! mismatch in Z2:'
     write(iout6,*) 'input=',z2, 'disk=',z2_p
  endif
  
  if (iinterp.eq.0) then
     do ig=1,ngrids   
        if (abs(rgrid(ig)-rgrid_p(ig)).gt.1.d-6 ) then
           write (iout6,*) 'mismatch in a subgrid ',ig,':'
           write (iout6,*) 'input=',rgrid(ig),' file=',rgrid_p(ig)
           write (iout6,*) 'continuing with fingers crossed ...'
        endif
     enddo
  endif
  
  
  if (nel.ne.nel_p) then
     write (iout6,*) '    Warning: mismatch in number of electrons: input=',nel,' file=',nel_p
  endif
  
  if (norb.ne.norb_p ) then
     write (iout6,*) '    Warning: mismatch in number of orbitals:  input=',norb,' file=',norb_p
     !         if (norb.gt.norb_p.and.imethod.eq.2) then 
     
     if ( norb.gt.norb_p ) then 
        write (iout6,*) '... assuming a virtual orbital is being generated ...'
        ini4=norb-norb_p
     else
        stop 'rheader'
     endif
  endif
  
  return
  
1000 continue
  write(iout6,1070)
  stop

1050 format(/1x,'... retrieving data from disk ...')	
1060 format( 1x,'    file with heading: ',80a1)
1061 format( 1x,'    created on: ',80a1)
1070 format(//1x,'rheader: cannot read data from disk'//)     
1080 format( 1x,'    file with heading: ',a80)
1081 format( 1x,'    created on: ',a80)

end subroutine rheader


