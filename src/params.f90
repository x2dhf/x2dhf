module params
  integer, parameter :: SP = kind(1.0)
  integer, parameter :: DP = kind(1.0d0)
  !  integer, parameter :: DP16 = SELECTED_REAL_KIND (1.d308 )
  integer, parameter :: PREC = DP
  integer, parameter :: PREC16 = DP

! In order to make the 2DHF program the array dimensions in various routines have to be
! set according to the values of the following variables:

! maxnu    = maximum number of grid points in \nu direction 
! maxmu    = maximum number of grid points in \mu direction 

! maxbasis = maximum number of primitive gaussian basis functions used to construct
!            numerical molecular orbitals out of data retrieved from GAUSSIAN outputs

  integer, parameter :: maxnu    = 1500
  integer, parameter :: maxmu    = 2500
  integer, parameter :: maxbasis = 650


! the following values should be left unchanged
  integer, parameter :: maxgrids = 1
  integer, parameter :: maxorb   = 60
  integer, parameter :: maxmpole  = 8

  character*10 :: formint,formfp,formfp64,formfp128
  data formint /'(10i15)  '/
  data formfp64,formfp128/'(5e25.16) ','(5e41.32) '/

  integer :: izero,ione,itwo,ithree,ifour
  data izero,ione,itwo,ithree,ifour /0,1,2,3,4/

  integer :: iinp5,iinp11,iinp12,iinp13,iinp14,iout6,iout21,iout22,iout23,iout24
  data iinp5/5/,iout6/6/,iinp11/11/,iinp12/12/,iinp13/13/,iinp14/14/,iout21/21/,iout22/22/,iout23/23/,iout24/24/

  real (PREC) :: zero,half,tenth,one,two,three,four,five,six,seven,eight,ten,pii,precis
  data zero,half,tenth,one,two,three,four,five,six,seven,eight,ten&
       /0.0_PREC,0.50_PREC,0.10_PREC,1.0_PREC,2.0_PREC,3.0_PREC,4.0_PREC,5.0_PREC,6.0_PREC, &
       7.0_PREC,8.0_PREC,10.0_PREC/

   integer :: incrni,incrmu,maxunit,lengthfp,lengthint,lengthintin,lengthfpin,iprint16,idat,&
       iord_nu_orb,iord_nu_coul,iord_nu_exch,iord_mu_orb,iord_mu_coul,iord_mu_exch,&
       lexchrecl,nsuppl,nsctch,inon

   integer, dimension(1000) :: idbg,iprint

end module params
