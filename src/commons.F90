! SPDX-License-Identifier: GPL-2.0-or-later

! Copyright (C) 1996-2024  Jacek Kobus 

module commons
  use params
#ifdef LIBXC
  use, intrinsic :: iso_c_binding
  ! in case lxclabels is extended make sure that nlxclabels is increased
  ! accordingly
  integer, parameter :: nlxclabels=3
  integer(c_int), dimension(nlxclabels) :: lxcFuncs2use
#endif
  
  real (PREC) :: bohr2ang,au2Debye,ainfcm
  data au2Debye /2.541765_PREC/
  ! 1a.u.=5.29177210671212x10^-11m
  ! data bohr2ang /0.529177249_PREC/
  ! NIST: 0.5291772108(18)
  ! SzS:  0.5291772490085559
  data bohr2ang /0.529177210671212_PREC/
  ! ainfcm: bohr radius in cm
  ! data ainfcm/0.529177249e-08_PREC/
  data ainfcm/0.5291772410671212e-08_PREC/
  
  character*2, dimension(0:118) :: element
  data element/ &
       ' ','H ','He','Li','Be','B ','C ','N ','O ','F ','Ne','Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca',&
       'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ','Zr',&
       'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd',&
       'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg',&
       'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm',&
       'Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og'/

  !   atomic masses of most abundant isotopes extracted from
  !   The 1983 Atomic Mass Evaluation by Wapstra and Audi

  real (PREC), dimension(0:100) :: alphaopt
  real (PREC), dimension(0:118) :: atweight
  real (PREC), dimension(0:103) :: dgsz(0:103)

  ! atomic weights for elements 101-118 taken from https://pubchem.ncbi.nlm.nih.gov/periodic-table/
  data atweight/ &
       0.0_PREC, 1.0078250350_PREC,4.002603240_PREC,7.01600300_PREC,9.01218220_PREC,1.10093054e1_PREC,&
       1.2e1_PREC,1.4003074002e1_PREC,1.599491463e1_PREC,1.899840322e1_PREC,1.99924356e1_PREC,&
       2.29897677e1_PREC,2.39850423e1_PREC,2.69815386e1_PREC,2.79769271e1_PREC,3.09737620e1_PREC, &
       3.19720707e1_PREC,3.4968852721e1_PREC,3.99623837e1_PREC,3.89637074e1_PREC,3.99627906e1_PREC,&
       4.49559100e1_PREC,4.79479473e1_PREC,5.09439617e1_PREC,5.19405098e1_PREC,5.49380471e1_PREC,&
       5.59349393e1_PREC,5.89331976e1_PREC,5.79353462e1_PREC,6.29295989e1_PREC,6.39291448e1_PREC,&
       6.8925580e1_PREC,7.39211774e1_PREC,7.49215942e1_PREC,7.99165196e1_PREC,7.89183361e1_PREC,&
       8.3911507e1_PREC,8.4911794e1_PREC,8.79056188e1_PREC,8.8905849e1_PREC,8.99047026e1_PREC,&
       9.29063772e1_PREC,9.79054073e1_PREC,9.7907215e1_PREC,1.019043485e2_PREC,1.02905500e2_PREC,&
       1.05903478e2_PREC,1.06905092e2_PREC,1.13903357e2_PREC,1.14903882e2_PREC,1.199021991e2_PREC,&
       1.209038212e2_PREC,1.29906229e2_PREC,1.26904473e2_PREC,1.31904144e2_PREC,1.32905429e2_PREC,&
       1.37905232e2_PREC,1.38906347e2_PREC,1.39905433e2_PREC,1.40907647e2_PREC,1.41907719e2_PREC,&
       1.44912743e2_PREC,1.51919728e2_PREC,1.52921225e2_PREC,1.57924019e2_PREC,1.58925342e2_PREC,&
       1.63929171e2_PREC,1.64930319e2_PREC,1.65930290e2_PREC,1.689342120_PREC,1.73938859e2_PREC,&
       1.74940770e2_PREC,1.799465457e2_PREC,1.80947992e2_PREC,1.83950928e2_PREC,1.86955744e2_PREC,&
       1.91961467e2_PREC,1.92962917e2_PREC,1.94964766e2_PREC,1.96966543e2_PREC,2.01970617e2_PREC,&
       2.02972320e2_PREC,2.07976627e2_PREC,2.08980374e2_PREC,2.08982404e2_PREC,2.09987126e2_PREC,&
       2.22017571e2_PREC,2.23019733e2_PREC,2.26025403e2_PREC,2.27027750e2_PREC,2.320380508e2_PREC,&
       2.31035880e2_PREC,2.380507847e2_PREC,2.370481678e2_PREC,2.44064199e2_PREC,2.43061375e2_PREC,&
       2.47070347e2_PREC,2.47070300e2_PREC,2.51079580e2_PREC,2.52082944e2_PREC,2.57095099e2_PREC,&
       2.5809843e2_PREC, 2.5910100e2_PREC, 2.66120e2_PREC,   2.67122e2_PREC,   2.68126e2_PREC,&
       2.69128e2_PREC,   2.70133e2_PREC,   2.691336e2_PREC,  2.77154e2_PREC,   2.82166e2_PREC,&
       2.82169e2_PREC,   2.86179e2_PREC,   2.86182e2_PREC,   2.90192e2_PREC,   2.90196e2_PREC,&
       2.93205e2_PREC,   2.94211e2_PREC,   2.95216e2_PREC/ 
  
  ! recommended alpha values for elements Z=2-41 (see K.Schwarz Phys. Rev. B 5 (1972) 2466-2468)
  ! alpha for hydrogen tuned to produce total energy equat to -0.500000au
  data alphaopt/     0.772980_PREC,0.676400_PREC,0.772980_PREC,0.781470_PREC,0.768230_PREC,0.765310_PREC,&
       0.759280_PREC,0.751970_PREC,0.744470_PREC,0.737320_PREC,0.730810_PREC,0.731150_PREC,0.729130_PREC,&
       0.728530_PREC,0.727510_PREC,0.726200_PREC,0.724750_PREC,0.723250_PREC,0.721770_PREC,0.721170_PREC,&
       0.719840_PREC,0.718410_PREC,0.716950_PREC,0.715560_PREC,0.713520_PREC,0.712790_PREC,0.711510_PREC,&
       0.710180_PREC,0.708960_PREC,0.706970_PREC,0.706730_PREC,0.706900_PREC,0.706840_PREC,0.706650_PREC,&
       0.706380_PREC,0.706060_PREC,0.705740_PREC,0.705530_PREC,0.705040_PREC,0.704650_PREC,0.704240_PREC,&
       0.703830_PREC,0.703830_PREC,0.703830_PREC,0.703830_PREC,0.703830_PREC,0.703830_PREC,0.703830_PREC,&
       0.703830_PREC,0.703830_PREC,0.703830_PREC,0.703830_PREC,0.703830_PREC,0.703830_PREC,0.703830_PREC,&
       0.703830_PREC,0.703830_PREC,0.703830_PREC,0.703830_PREC,0.703830_PREC,0.703830_PREC,0.703830_PREC,&
       0.703830_PREC,0.703830_PREC,0.703830_PREC,0.703830_PREC,0.703830_PREC,0.703830_PREC,0.703830_PREC,&
       0.703830_PREC,0.703830_PREC,0.703830_PREC,0.703830_PREC,0.703830_PREC,0.703830_PREC,0.703830_PREC,&
       0.703830_PREC,0.703830_PREC,0.703830_PREC,0.703830_PREC,0.703830_PREC,0.703830_PREC,0.703830_PREC,&
       0.703830_PREC,0.703830_PREC,0.703830_PREC,0.703830_PREC,0.703830_PREC,0.703830_PREC,0.703830_PREC,&
       0.703830_PREC,0.703830_PREC,0.703830_PREC,0.703830_PREC,0.703830_PREC,0.703830_PREC,0.703830_PREC,&
       0.703830_PREC,0.703830_PREC,0.703830_PREC,0.70383_PREC/

  ! d parameters for Green, Sellin, Zachor model HF potential added values
  ! for Z=1 and 2 to get approximate 1s orbital energies for H and He
  data dgsz/        0.0000_PREC, 0.1000_PREC, 1.9901_PREC, 0.5630_PREC, 0.8580_PREC, 0.9790_PREC, &
       0.8800_PREC, 0.7760_PREC, 0.7080_PREC, 0.5750_PREC, 0.5000_PREC, 0.5610_PREC, 0.6210_PREC, &
       0.7290_PREC, 0.8170_PREC, 0.8680_PREC, 0.8850_PREC, 0.8810_PREC, 0.8620_PREC, 1.0060_PREC, &
       1.1540_PREC, 1.1160_PREC, 1.0600_PREC, 0.9960_PREC, 0.8370_PREC, 0.8660_PREC, 0.8070_PREC, &
       0.7510_PREC, 0.7000_PREC, 0.6060_PREC, 0.6120_PREC, 0.6310_PREC, 0.6490_PREC, 0.6630_PREC, &
       0.6750_PREC, 0.6840_PREC, 0.6890_PREC, 0.7440_PREC, 0.7980_PREC, 0.8550_PREC, 0.8660_PREC, &
       0.8310_PREC, 0.8250_PREC, 0.8550_PREC, 0.8030_PREC, 0.7880_PREC, 0.7370_PREC, 0.7540_PREC, &
       0.7750_PREC, 0.8100_PREC, 0.8410_PREC, 0.8700_PREC, 0.8960_PREC, 0.9190_PREC, 0.9400_PREC, &
       1.0220_PREC, 1.1080_PREC, 1.1500_PREC, 1.0810_PREC, 0.9700_PREC, 0.9380_PREC, 0.9050_PREC, &
       0.8730_PREC, 0.8420_PREC, 0.8620_PREC, 0.8300_PREC, 0.7540_PREC, 0.7280_PREC, 0.7020_PREC, &
       0.6770_PREC, 0.6540_PREC, 0.6650_PREC, 0.6720_PREC, 0.6760_PREC, 0.6790_PREC, 0.6800_PREC, &
       0.6800_PREC, 0.6790_PREC, 0.6610_PREC, 0.6570_PREC, 0.6710_PREC, 0.6900_PREC, 0.7080_PREC, &
       0.7260_PREC, 0.7440_PREC, 0.7610_PREC, 0.7770_PREC, 0.8180_PREC, 0.8590_PREC, 0.8990_PREC, &
       0.9270_PREC, 0.8870_PREC, 0.8800_PREC, 0.8720_PREC, 0.8320_PREC, 0.8220_PREC, 0.8420_PREC, &
       0.8300_PREC, 0.7900_PREC, 0.7780_PREC, 0.7660_PREC, 0.7540_PREC, 0.7420_PREC, 0.755_PREC/

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  character*8, dimension(maxorb)  :: gusym,bond,orbsym
  character*8, dimension(4*maxorb)  :: spin
  character*8, dimension(5*maxorb)  :: spn

  logical lversion,ldatPresent,altSweeps,lomegaz,lbreakci,lhomonucl  
  logical DFT,HF,HFS,LXC,lxcHyb,lxcPolar,lxcSupp,OED,TED,SCMC,HFinput,initAddData
  logical linitFuncsGauss,linitFuncsGaussC,linitFuncsHF,linitFuncsHydrogen,linitFuncsLDA,&
       linitFuncsMolcas,linitFuncsNodat,linitFuncsNoexch,initFuncsOED,linitFuncsOld,linitFuncsQRHF 
  logical lpotCoulomb,lpotCoul2,lpotCoul3,lpotFermi,lpotGauss,&
       lpotGSZ,lpotGSZG,lpotKH,lpotSAP,lspherium,lfefield,lmmoments,&
       lpotHarm2,lpotHarm3,lpotHooke,lintracule,lextracule
  logical loffDiagLM ,hfIncl,homoIncl,lcaoIncl,ldaIncl,ldaSAPIncl,mcsorpt,mmoments,omegaIncl,openShell
  logical lpotmcsor,lorbmcsor,lmcsorpt,lcoulexch,openmp,pthread,tpool
  logical lfixcoul,lfixexch,lfixorb,lfixorb4init,lfastexch,lfastscf,lcheckTotEnergy,&
       lplot,lpthreadpoolq,ltail

  logical inFormatI32,inFormatI64,inFormatR128,outFormatI32,outFormatI64,outFormatR128,&
       inUnFormatted,outUnFormatted  
  logical, dimension(maxorb,maxorb) :: offDiagLM 

  integer (KIND=IPREC) :: nexchp,norb,nel,nexch,nelInput,maxpots
  
  integer (KIND=IPREC) :: iharm2xy,magfield,inclorb,iinterp,imethod,ipot,&
       ienterm,ienterm4init,inoterm,iftail,iplot,istep
  integer (KIND=IPREC) :: iout4pair,iout4dft,iout4kinpot

  integer (KIND=IPREC) :: saveSCFdata,recalcMM
  integer (KIND=IPREC) :: io,ir, maxscf,idftex,idftcorr,ini4,iscmc,&
       ihomon,icanon,nenlast,nnolast,nscf2skip,lsor,mpole,&
       iomega,icase,mulast,ngrids_p,nni_p,mxnmu_p,mxsize_p,npbasis,lmtype,iscf,maxresets,nscfextra
  integer (KIND=IPREC) :: verboseLevel

  integer (KIND=IPREC) :: iform
  integer (KIND=IPREC),dimension(4) ::  iadint2,iadint3l,iadint3r,iadint4
  integer (KIND=IPREC),dimension(5) ::  icolours,i6b
  integer (KIND=IPREC),dimension(0:4,0:maxthreads4mcsor-1) :: nstart,nstop
  integer (KIND=IPREC),dimension(0:maxthreads4mcsor-1) :: nstart7,nstart6a,nstart6b,nstop7,nstop6a,nstop6b      
  integer (KIND=IPREC),dimension(1) :: nmu,ngsize,inpv,ibmu,iemu,ioffs
  integer (KIND=IPREC) :: nmu_p,ibmu_p,iemu_p
  integer (KIND=IPREC) :: nsor4orb,nsor4pot,ireset

  character*8 :: meshOrdering
  
  integer (KIND=IPREC),dimension(30) :: i4b,i5b,i4e,i5e,i4si,i5si,i4ng,i5ng,iadext,iadnor,iadex1,iadex2,iadex3
  integer (KIND=IPREC),dimension(maxorb) :: i1b,i2b,i1e,i2e,i1si,i2si,i1ng,i2ng,i1mu,i2mu,ifixorb,iorn,&
       inhyd,inhydlcao,inDFT,nn,ll,mm,iocc,isymOrb,lock,ige,ihomo,iscforder,itouch,iloop
  integer (KIND=IPREC),dimension(2) :: maxsororb,maxsorpot
  
  integer (KIND=IPREC),dimension(2*maxorb) :: lagraon
  integer (KIND=IPREC),dimension(300) :: iqm
  integer (KIND=IPREC),dimension(maxorb*(maxorb+1)/2) :: i3b,i3e,i3si,i3ng,i3mu,ilc
  integer (KIND=IPREC),dimension(maxorb*(1+maxorb/2)) :: iwexch,i3xrec1,i3xrec2,i3orb1,i3orb2
  integer (KIND=IPREC),dimension(maxbasis) :: lprim,mprim,icgau,ixref
  integer (KIND=IPREC),dimension(maxnu*maxmu) :: indx1nu,indx1mu
  integer (KIND=IPREC),dimension(maxorb,maxorb) :: i3xind,ilc2,k2,nlm
  integer (KIND=IPREC),dimension(9,maxorb) :: mgx
  integer (KIND=IPREC),dimension(240,240) :: imagn

  
  logical :: ldetectNaN,lfixNaN
  integer (KIND=IPREC) :: lxcFuncs
  
  character*4, dimension(10) :: cdftex,cdftcorr
  
  real (PREC), dimension(1) :: rgrid
  real (PREC) :: rgrid_p  
  real (PREC) :: homolevl,homoallign,ffield,fgrad,zcutoff,harm2xy,gammaf,z1atmass,z2atmass

  real (PREC) :: orbEnergyIncThld,recalcMMfactor
  real (PREC) :: alphaf,alphaflxc,trelaxRealOrb,trelaxRealCoul,trelaxRealExch,&
       tmomenReal,trelaxCpu,trelaxCpuOrb,trelaxCpuCoul,trelaxCpuExch,&           
       tortho,trayl,tmomen,tlagra,ttoten,sflagra,dflagra,&
       vkt,vnt,evt,epott,etot,virrat,etotFN,ictot,iready,hallign,&
       wtwoelLXC,dens_threshold,plotThreshold
  
  real (PREC), dimension(10) :: elect,electDA,total,totalDA
  
  real (PREC) :: enkin,ennucel,encoul,enexch,entot,encouldft,enexchdft,edftex,edftcorr,wexhyb,&
       alegk0,epspot,ompot,apot,v0pot,hooke,r_p,z1_p,z2_p,rinf_p,date_p,time_p
  
  integer (KIND=IPREC) :: mpot, nsimp
  
  real (PREC) :: co1lda,co2lda
  real (PREC), dimension(12) :: calp
  real (PREC), dimension(maxorb) :: eza1,eza2,co1,co2,wstorthog,qxm,occ,sign,demax,vk,vn,engt,pnc,div
  real (PREC), dimension(maxorb) :: orbNorm
  real (PREC), dimension(maxorb,maxorb) :: ee
  
  real (PREC), dimension(maxbasis) :: fngau2,shngau,primexp,coeff
  
  real (PREC), dimension(maxorb,maxbasis) :: primcoef
  
  real (PREC), dimension(maxbasis,maxbasis) :: ovl
  
  integer (KIND=IPREC),dimension(maxorb+1) :: nexchpots
  integer (KIND=IPREC),dimension(maxorb,2*maxorb) :: ibexcp,ins1,ins2,ipcs,isyms
  integer (KIND=IPREC),dimension(maxorb,2*maxorb) :: deltam4pot
  integer (KIND=IPREC),dimension(maxorb*(maxorb)/2) :: iorb1mm,iorb2mm,idelmm,ipcmm
  integer (KIND=IPREC) :: nexchmm
end module commons
