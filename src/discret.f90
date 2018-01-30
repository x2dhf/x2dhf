module discret
  use params
      integer :: ngrids,nni,mxnmu,mxsize,ngrd1,ngrd2,ngrd6a,ngrd6b,ngrd7,nni1,nni2,nni3,nni4,nni5,ini,ini4
      real (PREC), dimension(4,10) :: dmu2,dmu1
      real (PREC), dimension(4) :: dni2,dni1,dmu2t,dmu1t
      real (PREC), dimension(maxmu,12) :: xi
      real (PREC), dimension(maxnu,6) :: eta
      real (PREC), dimension(9,maxnu) :: d1ni,dni
      real (PREC), dimension(9,maxmu) :: d1mu,dmu
      real (PREC), dimension(9,4) :: cint2,cint3l,cint3r,cint4
      real (PREC), dimension(10) :: hmu,hmu_p
      real (PREC), dimension(5) :: exeven,exodd
      real (PREC), dimension(maxmu) :: vmu,vxi,vxisq,vxi1,vxi2,vmu_p,vxi_p,vxisq_p,vxi1_p,vxi2_p,wmu
      real (PREC), dimension(maxnu) :: vni,veta,vetasq,veta1,veta2,vni_p,veta_p,vetasq_p,veta1_p,veta2_p,wni	
      real (PREC) :: r,r2,z1,z2,rinf,cutorb,cutcoul,cutexch

end module discret
