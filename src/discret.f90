! SPDX-License-Identifier: GPL-2.0-or-later

! Copyright (C) 2024  Jacek Kobus 

module discrete
  use params
  integer (KIND=IPREC) :: ngrids,nni,mxnmu,mxsize,mxsize8,ngrid1,ngrid2,ngrid6a,ngrid6b,ngrid7
  integer (KIND=IPREC) :: nni1,nni2,nni3,nni4,nni5,nni8,mxnmu8,firstOrb,lastOrb
  real (PREC), dimension(4) :: dmu1,dmu2
  real (PREC), dimension(4) :: dni1,dni2
  real (PREC), dimension(maxmu,12) :: xi
  real (PREC), dimension(maxnu,6) :: eta
  real (PREC), dimension(9,maxnu) :: d1ni,dni
  real (PREC), dimension(9,maxmu) :: d1mu,dmu
  real (PREC), dimension(9,4) :: cint2,cint3l,cint3r,cint4
  real (PREC), dimension(10) :: hmu
  real (PREC) :: hmu_p  
  real (PREC), dimension(5) :: exeven,exodd
  real (PREC), dimension(maxmu) :: vmu,vxi,vxisq,vxi1,vxi2,vmu_p,vxi_p,wmu
  real (PREC), dimension(maxnu) :: vni,veta,vetasq,veta1,veta2,vni_p,veta_p,vetasq_p,veta1_p,veta2_p,wni
  real (PREC) :: r,r2,z1,z2,rinf,cutorb,cutcoul,cutexch
  real (PREC) :: hni,hni_p
end module discrete
