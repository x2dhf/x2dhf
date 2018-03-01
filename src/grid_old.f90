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
! ### grid_old ###

!     Initializes grid info to enable interpolation.

subroutine grid_old
  use params
  use discret
  use scf
  use commons8

  ibmu_p(1)=1
  iemu_p(1)=nmu_p(1)
  mxnmu_p=nmu_p(1)
  ioffs_p(1)=0
  do ig=2,ngrids_p
     mxnmu_p=mxnmu_p+nmu_p(ig)-1
     ibmu_p(ig)=iemu_p(ig-1)
     iemu_p(ig)=mxnmu_p
     ioffs_p(ig)=nni_p*(iemu_p(ig-1)-1)
  enddo
  
  !     determine the size of each grid
  
  do ig=1,ngrids_p
     ngsize_p(ig)=nni_p*nmu_p(ig)
  enddo
  
  !     total no of grid points 
  
  mxsize_p=nni_p*mxnmu_p
  
  !     hni - step in ni variable
  !     hmu - step in mu variable
  
  !     determine step size in ni variable 
  
  hni_p=pii/dble(nni_p-1)
  
  !     determine step size in mu variable for each grid
  
  if (ngrids_p.eq.1.and.rgrid_p(1).ge.0.50_PREC) then
     xmi0=2.0_PREC*rgrid_p(1)/r_p
     xmi0=log(xmi0+sqrt(xmi0*xmi0-1.0_PREC))
     hmu_p(1)=xmi0/dble(mxnmu_p-1)
  else
     do ig=1,ngrids_p
        if(ig.eq.1) then 
           hmu_p(ig)=rgrid_p(ig)
        else
           hmu_p(ig)=hmu_p(1)*rgrid_p(ig)
        endif
     enddo
  endif
  
  !     initialize mu and xi arryas
  
  vmu_p(1) =0.0_PREC
  vxi_p(1) =cosh(vmu_p(1))
  vxisq_p(1)=vxi_p(1)*vxi_p(1)
  vxi1_p(1)=sqrt(vxi_p(1)*vxi_p(1)-1.0_PREC)
  vxi2_p(1)=0.00_PREC
  
  ib=2
  
  do ig=1,ngrids_p
     ie=iemu_p(ig)
     do i=ib,ie
        vmu_p(i)=vmu_p(i-1)+hmu_p(ig)
        vxi_p(i)=cosh(vmu_p(i))
        vxisq_p(i)=vxi_p(i)*vxi_p(i)
        
        !           sqrt(vxi*vxi-1)  and  ssf*vxi/sqrt(vxi*vxi-1)
        
        vxi1_p(i)=sqrt(vxi_p(i)*vxi_p(i)-1.0_PREC)
        vxi2_p(i)=vxi_p(i)/vxi1_p(i)
        !           enddo over i
     enddo
     ib=ie+1
     !        enddo over ig
  enddo
  
  !     determine rinf
  
  write(iout6,*) 
  write(iout6,*) '    previous grid:'
  write(iout6,*) '    nni    nmu   hni  hmu      rb '
  rinf_p=r_p*vxi_p(iemu_p(ngrids_p))/2.0_PREC
  do ig=1,ngrids_p
     rinfig_p=r_p*vxi_p(iemu_p(ig))/2.0_PREC
     write(iout6,1000) nni_p,nmu_p(ig),hni_p,hmu_p(ig),rinfig_p
  enddo
01000 format(4x,i4,2x,i4,2f8.5,f8.3/)
  
  !    initialize ni and eta arrays
  
  do i=1,nni_p
     vni_p(i)=dble((i-1))*hni_p
     veta_p(i)=cos(vni_p(i))
     vetasq_p(i)=veta_p(i)*veta_p(i)
     
     !        sqrt(1-veta*veta) and  veta/sqrt(1-veta*veta)
     !        
     veta1_p(i)=sqrt(1.0_PREC-veta_p(i)*veta_p(i))
     if (veta1_p(i).lt.precis) then
        veta2_p(i)=0.0
     else
        veta2_p(i)=veta_p(i)/veta1_p(i)
     endif
  enddo
  
end subroutine grid_old
