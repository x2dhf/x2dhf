module commons8
  use params

  real (PREC) :: bohr2ang,au2Debye,ainfcm
  data au2Debye /2.541765_PREC/
  data bohr2ang /0.529177249_PREC/
  ! ainfcm: bohr radius in cm
  data ainfcm/0.529177249e-08_PREC/

  character*2, dimension(0:100) :: element
  data element/ &
       ' ','H ','He','Li','Be','B ','C ','N ','O ','F ','Ne','Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca',&
       'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ','Zr',&
       'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd',&
       'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg',&
       'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm'/

  !   atomic masses of most abundant isotopes extracted from
  !   The 1983 Atomic Mass Evaluation by Wapstra and Audi

  real (PREC), dimension(0:100) :: atweight,alphaopt
  real (PREC), dimension(0:103) :: dgsz(0:103)

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
       2.47070347e2_PREC,2.47070300e2_PREC,2.51079580e2_PREC,2.52082944e2_PREC,2.57095099e2_PREC/

  !     recommended alpha values for elements Z=2-41 (see K.Schwarz Phys. Rev. B 5 (1972) 2466-2468)
  !     alpha for hydrogen tuned to produce total energy equat to -0.500000au
  ! data alphaopt/0.000000_PREC, 0.000000_PREC,0.772980_PREC,0.781470_PREC,0.768230_PREC,0.765310_PREC,&
    data alphaopt/0.772980_PREC,0.676400_PREC,0.772980_PREC,0.781470_PREC,0.768230_PREC,0.765310_PREC,&
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

  !    d parameters for Green, Sellin, Zachor model HF potential
  !    added values for Z=1 and 2 to get approximate 1s orbital energies
  !    for H and He
  data dgsz/0.0000_PREC, 0.1000_PREC, 1.9901_PREC, 0.5630_PREC, 0.8580_PREC, 0.9790_PREC, &
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
      character*8, dimension(maxorb)  :: gut,bond,orbsym
      character*8, dimension(4*maxorb)  :: spin
      character*8, dimension(5*maxorb)  :: spn

      logical lcaoIncl,ldaIncl,omegaIncl
      integer :: nexchp,norb,nel,nexch,nelInput

      integer :: ifefield,iharm2xy,magfield,inclorb,idump,iinterp,imethod,ipot,ienterm,inoterm,iftail,ifermi,iplot
      integer :: inout32,inout64,inout128,inpform,ioutform

      integer :: io,ir, maxscf,nobckup,islat,idft,idftex,idftcorr,iscmc,ihomon,ibreak,icanon,ilagra, &
           nenlast,nnolast,nscf2skip,maxsor1,maxsor2,maxsor3,lsor,ipoiss,ialtsweeps,isstart,isstop, &
           isstep,mpole,iprtlev,ioo,iomega,icase,mulast, &
           ngrids_p,nni_p,mxnmu_p,mxsize_p,npbasis,lmtype,iscf

      integer ::  iexlorb,iexlcoul,iexlexp,iform

      integer, dimension(4) ::  iadint2,iadint3l,iadint3r,iadint4
      integer, dimension(5) ::  i6b
      integer, dimension(10) :: nmu,ngsize,inpv,ibmu,iemu,ioffs,iorder,nmu_p,ngsize_p,ibmu_p,iemu_p,ioffs_p
      real (PREC), dimension(10) :: rgrid_p
      integer, dimension(20) :: i4b,i5b,i4e,i5e,i4si,i5si,i4ng,i5ng,iadext,iadnor,iadex1,iadex2,iadex3
      integer, dimension(maxorb) :: i1b,i2b,i1e,i2e,i1si,i2si,i1ng,i2ng,i1mu,i2mu,i3nexcp,ifix,iorn,&
           inhyd,inhydlcao,nn,ll,mm,iocc,isymOrb,lock,ige,ihomo,iscforder,maxsororb,maxsorpot,itouch,iloop

      integer, dimension(120) :: lagraon
      integer, dimension(300) :: iqm
      integer, dimension(1830) :: i3b,i3e,i3si,i3ng,i3mu,ilc
      integer, dimension(1860) :: iwexch,i3xrec1,i3xrec2,i3orb1,i3orb2
      integer, dimension(maxbasis) :: lprim,mprim,icgau,ixref
      integer, dimension(4,10) :: ingr1,ingr2
      integer, dimension(60,60) :: i3btv,i3xind,i3xpair,i3brec,i3xk,i3breck,nlmf,nlm,nonorthog
      integer, dimension(9,60) :: mgx
      integer, dimension(240,240) :: imagn

      character*4, dimension(10) :: cdftex,cdftcorr
      
      integer :: nlxclabels,lxcFuncs
      !parameter (nlxclabels=55)
      parameter (nlxclabels=337)
      character*30 lxclabels(nlxclabels)
      integer, dimension(nlxclabels) :: lxcnumbs,lxcFuncs2use

      
      ! in case lxclabels is extended make sure that nlxclabels is increased accordingly (see common8.f90)

      data lxclabels /                   &
           'lda_x',                      & !  Exchange
           'lda_c_wigner',               & !  Wigner parametrization
           'lda_c_rpa',                  & !  Random Phase Approximation
           'lda_c_hl',                   & !  Hedin & Lundqvist
           'lda_c_gl',                   & !  Gunnarson & Lundqvist        
           'lda_c_xalpha',               & !  Slater Xalpha                       
           'lda_c_vwn',                  & !  Vosko, Wilk, & Nusair (5)
           'lda_c_vwn_rpa',              & !  Vosko, Wilk, & Nusair (RPA)
           'lda_c_pz',                   & !  Perdew & Zunger
           'lda_c_pz_mod',               & !  Perdew & Zunger (Modified)   
           'lda_c_ob_pz',                & !  Ortiz & Ballone (PZ)         
           'lda_c_pw',                   & !  Perdew & Wang
           'lda_c_pw_mod',               & !  Perdew & Wang (Modified)     
           'lda_c_ob_pw',                & !  Ortiz & Ballone (PW)         
           'lda_c_2d_amgb',              & !  Attaccalite et al
           'lda_c_2d_prm',               & !  Pittalis, Rasanen & Marques correlation in 2D
           'lda_c_vbh',                  & !  von Barth & Hedin            
           'lda_c_1d_csc',               & !  Casula, Sorella, and Senatore 1D correlation
           'lda_x_2d',                   & !  Exchange in 2D
           'lda_xc_teter93',             & !  Teter 93 parametrization
           'lda_x_1d',                   & !  Exchange in 1D
           'lda_c_ml1',                  & !  Modified LSD (version 1) of Proynov and Salahub
           'lda_c_ml2',                  & !  Modified LSD (version 2) of Proynov and Salahub 
           'lda_c_gombas',               & !  Gombas parametrization
           'lda_c_pw_rpa',               & !  Perdew & Wang fit of the RPA 
           'lda_c_1d_loos',              & !  P-F Loos correlation LDA
           'lda_c_rc04',                 & !  Ragot-Cortona
           'lda_c_vwn_1',                & !  Vosko, Wilk, & Nusair (1)
           'lda_c_vwn_2',                & !  Vosko, Wilk, & Nusair (2)
           'lda_c_vwn_3',                & !  Vosko, Wilk, & Nusair (3)
           'lda_c_vwn_4',                & !  Vosko, Wilk, & Nusair (4)
           'lda_xc_zlp',                 & !  Zhao, Levy & Parr, Eq. (20)
           'lda_k_tf',                   & !  Thomas-Fermi kinetic energy functional
           'lda_k_lp',                   & !  Lee and Parr Gaussian ansatz           
           'lda_xc_ksdt',                & !  Karasiev et al. parametrization
           'lda_c_chachiyo',             & !  Chachiyo simple 2 parameter correlation
           'lda_c_lp96',                 & !  Liu-Parr correlation
           'lda_x_rel',                  & !  Relativistic exchange
           'lda_xc_1d_ehwlrg_1',         & !  LDA constructed from slab-like systems of 1 electron
           'lda_xc_1d_ehwlrg_2',         & !  LDA constructed from slab-like systems of 2 electrons 
           'lda_xc_1d_ehwlrg_3',         & !  LDA constructed from slab-like systems of 3 electrons 
           'lda_x_erf',                  & !  Attenuated exchange LDA (erf)
           'lda_xc_lp_a',                & !  Lee-Parr reparametrization B 
           'lda_xc_lp_b',                & !  Lee-Parr reparametrization B 
           'lda_x_rae',                  & !  Rae self-energy corrected exchange  
           'lda_k_zlp',                  & !  kinetic energy version of ZLP
           'lda_c_mcweeny',              & !  McWeeny 76 
           'lda_c_br78',                 & !  Brual & Rothstein 78 
           'lda_c_pk09',                 & !  Proynov and Kong 2009
           'lda_c_ow_lyp',               & !  Wigner with corresponding LYP parameters 
           'lda_c_ow',                   & !  Optimized Wigner 
           'lda_xc_gdsmfb',              & !  Groth et al. parametrization 
           'lda_c_gk72',                 & !  Gordon and Kim 1972
           'lda_c_karasiev',             & !  Karasiev reparameterization of Chachiyo   
           'lda_k_lp96',                 & !  Liu-Parr kinetic 
           'gga_x_gam',                  & !  GAM functional from Minnesota 
           'gga_c_gam',                  & !  GAM functional from Minnesota            
           'gga_x_hcth_a',               & !  HCTH-A
	   'gga_x_ev93',                 & !  Engel and Vosko
	   'gga_x_bcgp',                 & !  Burke, Cancio, Gould, and Pittalis             
	   'gga_c_bcgp',                 & !  Burke, Cancio, Gould, and Pittalis
	   'gga_x_lambda_oc2_n',         & !  lambda_OC2(N) version of PBE                   
	   'gga_x_b86_r',                & !  Revised Becke 86 Xalpha,beta,gamma (with mod. grad. correction) 
	   'gga_x_lambda_ch_n',          & !  lambda_CH(N) version of PBE                    
	   'gga_x_lambda_lo_n',          & !  lambda_LO(N) version of PBE                    
	   'gga_x_hjs_b88_v2',           & !  HJS screened exchange corrected B88 version
	   'gga_c_q2d',                  & !  Chiodo et al
	   'gga_x_q2d',                  & !  Chiodo et al
	   'gga_x_pbe_mol',              & !  Del Campo, Gazquez, Trickey and Vela (PBE-like) 
	   'gga_k_tfvw',                 & !  Thomas-Fermi plus von Weiszaecker correction
	   'gga_k_revapbeint',           & !  erpolated version of REVAPBE                
	   'gga_k_apbeint',              & !  interpolated version of APBE                   
	   'gga_k_revapbe',              & !  revised APBE                                   
	   'gga_x_ak13',                 & !  Armiento & Kuemmel 2013
	   'gga_k_meyer',                & !  Meyer,  Wang, and Young
	   'gga_x_lv_rpw86',             & !  Berland and Hyldgaard
	   'gga_x_pbe_tca',              & !  PBE revised by Tognetti et al                  
	   'gga_x_pbeint',               & !  PBE for hybrid interfaces
	   'gga_c_zpbeint',              & !  spin-dependent gradient correction to PBEint
	   'gga_c_pbeint',               & !  PBE for hybrid interfaces                          
	   'gga_c_zpbesol',              & !  spin-dependent gradient correction to PBEsol       
	   'gga_xc_opbe_d',              & !  oPBE_D functional of Goerigk and Grimme   
	   'gga_xc_opwlyp_d',            & !  oPWLYP-D functional of Goerigk and Grimme 
	   'gga_xc_oblyp_d',             & !  oBLYP-D functional of Goerigk and Grimme
	   'gga_x_vmt84_ge',             & !  VMT{8,4} with constraint satisfaction with mu = mu_GE  
	   'gga_x_vmt84_pbe',            & !  VMT{8,4} with constraint satisfaction with mu = mu_PBE
	   'gga_x_vmt_ge',               & !  Vela, Medel, and Trickey with mu = mu_GE  
	   'gga_x_vmt_pbe',              & !  Vela, Medel, and Trickey with mu = mu_PBE
	   'gga_c_n12_sx',               & !  N12-SX functional from Minnesota         
	   'gga_c_n12',                  & !  N12 functional from Minnesota
	   'gga_x_n12',                  & !  N12 functional from Minnesota
	   'gga_c_regtpss',              & !  Regularized TPSS correlation (ex-VPBE)
	   'gga_c_op_xalpha',            & !  one-parameter progressive functional (XALPHA version)
	   'gga_c_op_g96',               & !  one-parameter progressive functional (G96 version)
	   'gga_c_op_pbe',               & !  one-parameter progressive functional (PBE version)
	   'gga_c_op_b88',               & !  one-parameter progressive functional (B88 version)
	   'gga_c_ft97',                 & !  Filatov & Thiel correlation
	   'gga_c_spbe',                 & !  PBE correlation to be used with the SSB exchange   
	   'gga_x_ssb_sw',               & !  Swart, Sola and Bickelhaupt correction to PBE
	   'gga_x_ssb',                  & !  Swart, Sola and Bickelhaupt  
	   'gga_x_ssb_d',                & !  Swart, Sola and Bickelhaupt dispersion  
	   'gga_xc_hcth_407p',           & !  HCTH/407+                                
	   'gga_xc_hcth_p76',            & !  HCTH p=7/6                               
	   'gga_xc_hcth_p14',            & !  HCTH p=1/4                               
	   'gga_xc_b97_gga1',            & !  Becke 97 GGA-1                           
	   'gga_c_hcth_a',               & !  HCTH-A
	   'gga_x_bpccac',               & !  BPCCAC (GRAC for the energy)
	   'gga_c_revtca',               & !  Tognetti, Cortona, Adamo (revised)
	   'gga_c_tca',                  & !  Tognetti, Cortona, Adamo
	   'gga_x_pbe',                  & !  Perdew, Burke & Ernzerhof exchange
	   'gga_x_pbe_r',                & !  Perdew, Burke & Ernzerhof exchange (revised)   
	   'gga_x_b86',                  & !  Becke 86 Xalpha,beta,gamma
	   'gga_x_herman',               & !  Herman et al original GGA
	   'gga_x_b86_mgc',              & !  Becke 86 Xalpha,beta,gamma (with mod. grad. correction) 
	   'gga_x_b88',                  & !  Becke 88
	   'gga_x_g96',                  & !  Gill 96
	   'gga_x_pw86',                 & !  Perdew & Wang 86
	   'gga_x_pw91',                 & !  Perdew & Wang 91
	   'gga_x_optx',                 & !  Handy & Cohen OPTX 01
	   'gga_x_dk87_r1',              & !  dePristo & Kress 87 (version R1)
	   'gga_x_dk87_r2',              & !  dePristo & Kress 87 (version R2)               
	   'gga_x_lg93',                 & !  Lacks & Gordon 93
	   'gga_x_ft97_a',               & !  Filatov & Thiel 97 (version A)
	   'gga_x_ft97_b',               & !  Filatov & Thiel 97 (version B) 
	   'gga_x_pbe_sol',              & !  Perdew, Burke & Ernzerhof exchange (solids)    
	   'gga_x_rpbe',                 & !  Hammer, Hansen & Norskov (PBE-like)
	   'gga_x_wc',                   & !  Wu & Cohen
	   'gga_x_mpw91',                & !  Modified form of PW91 by Adamo & Barone 
	   'gga_x_am05',                 & !  Armiento & Mattsson 05 exchange
	   'gga_x_pbea',                 & !  Madsen (PBE-like)
	   'gga_x_mpbe',                 & !  Adamo & Barone modification to PBE
	   'gga_x_xpbe',                 & !  xPBE reparametrization by Xu & Goddard         
	   'gga_x_2d_b86_mgc',           & !  Becke 86 MGC for 2D systems
	   'gga_x_bayesian',             & !  ayesian best fit for the enhancement factor
	   'gga_x_pbe_jsjr',             & !  JSJR reparametrization by Pedroza, Silva & Capelle 
	   'gga_x_2d_b88',               & !  Becke 88 in 2D
	   'gga_x_2d_b86',               & !  Becke 86 Xalpha,beta,gamma
	   'gga_x_2d_pbe',               & !  Perdew, Burke & Ernzerhof exchange in 2D
	   'gga_c_pbe',                  & !  Perdew, Burke & Ernzerhof correlation
	   'gga_c_lyp',                  & !  Lee, Yang & Parr
	   'gga_c_p86',                  & !  Perdew 86
	   'gga_c_pbe_sol',              & !  Perdew, Burke & Ernzerhof correlation SOL          
	   'gga_c_pw91',                 & !  Perdew & Wang 91
	   'gga_c_am05',                 & !  Armiento & Mattsson 05 correlation
	   'gga_c_xpbe',                 & !  xPBE reparametrization by Xu & Goddard             
	   'gga_c_lm',                   & !  Langreth and Mehl correlation
	   'gga_c_pbe_jrgx',             & !  JRGX reparametrization by Pedroza, Silva & Capelle 
	   'gga_x_optb88_vdw',           & !  Becke 88 reoptimized to be used with vdW functional of Dion et al 
	   'gga_x_pbek1_vdw',            & !  PBE reparametrization for vdW                  
	   'gga_x_optpbe_vdw',           & !  PBE reparametrization for vdW 
	   'gga_x_rge2',                 & !  Regularized PBE
	   'gga_c_rge2',                 & !  Regularized PBE                                    
	   'gga_x_rpw86',                & !  refitted Perdew & Wang 86 
	   'gga_x_kt1',                  & !  Exchange part of Keal and Tozer version 1
	   'gga_xc_kt2',                 & !  Keal and Tozer version 2                  
	   'gga_c_wl',                   & !  Wilson & Levy
	   'gga_c_wi',                   & !  Wilson & Ivanov 
	   'gga_x_mb88',                 & !  Modified Becke 88 for proton transfer 
	   'gga_x_sogga',                & !  Second-order generalized gradient approximation 
	   'gga_x_sogga11',              & !  Second-order generalized gradient approximation 2011
	   'gga_c_sogga11',              & !  Second-order generalized gradient approximation 2011
	   'gga_c_wi0',                  & !  Wilson & Ivanov initial version
	   'gga_xc_th1',                 & !  Tozer and Handy v. 1 
	   'gga_xc_th2',                 & !  Tozer and Handy v. 2
	   'gga_xc_th3',                 & !  Tozer and Handy v. 3
	   'gga_xc_th4',                 & !  Tozer and Handy v. 4 
	   'gga_x_c09x',                 & !  C09x to be used with the VdW of Rutgers-Chalmers
	   'gga_c_sogga11_x',            & !  To be used with HYB_GGA_X_SOGGA11_X  
	   'gga_x_lb',                   & !  van Leeuwen & Baerends
	   'gga_xc_hcth_93',             & !  HCTH functional fitted to  93 molecules  
	   'gga_xc_hcth_120',            & !  HCTH functional fitted to 120 molecules  
	   'gga_xc_hcth_147',            & !  HCTH functional fitted to 147 molecules  
	   'gga_xc_hcth_407',            & !  HCTH functional fitted to 407 molecules  
	   'gga_xc_edf1',                & !  Empirical functionals from Adamson, Gill, and Pople
	   'gga_xc_xlyp',                & !  XLYP functional
	   'gga_xc_kt1',                 & !  Keal and Tozer version 1                  
	   'gga_xc_b97_d',               & !  Grimme functional to be used with C6 vdW term
	   'gga_xc_pbe1w',               & !  Functionals fitted for water 
	   'gga_xc_mpwlyp1w',            & !  Functionals fitted for water 
	   'gga_xc_pbelyp1w',            & !  Functionals fitted for water 
	   'gga_x_lbm',                  & !  van Leeuwen & Baerends modified
	   'gga_x_ol2',                  & !  Exchange form based on Ou-Yang and Levy v.2
	   'gga_x_apbe',                 & !  mu fixed from the semiclassical neutral atom   
	   'gga_k_apbe',                 & !  mu fixed from the semiclassical neutral atom   
	   'gga_c_apbe',                 & !  mu fixed from the semiclassical neutral atom       
	   'gga_k_tw1',                  & !  Tran and Wesolowski set 1 (Table II)           
	   'gga_k_tw2',                  & !  Tran and Wesolowski set 2 (Table II)           
	   'gga_k_tw3',                  & !  Tran and Wesolowski set 3 (Table II)           
	   'gga_k_tw4',                  & !  Tran and Wesolowski set 4 (Table II)           
	   'gga_x_htbs',                 & !  Haas, Tran, Blaha, and Schwarz
	   'gga_x_airy',                 & !  Constantin et al based on the Airy gas
	   'gga_x_lag',                  & !  Local Airy Gas
	   'gga_xc_mohlyp',              & !  Functional for organometallic chemistry 
	   'gga_xc_mohlyp2',             & !  Functional for barrier heights 
	   'gga_xc_th_fl',               & !  Tozer and Handy v. FL
	   'gga_xc_th_fc',               & !  Tozer and Handy v. FC  
	   'gga_xc_th_fcfo',             & !  Tozer and Handy v. FCFO 
	   'gga_xc_th_fco',              & !  Tozer and Handy v. FCO 
	   'gga_c_optc',                 & !  Optimized correlation functional of Cohen and Handy
	   'gga_c_pbeloc',               & !  Semilocal dynamical correlation
	   'gga_xc_vv10',                & !  Vydrov and Van Voorhis
	   'gga_c_pbefe',                & !  PBE for formation energies                         
	   'gga_c_optc',                 & !  one-parameter progressive functional (PW91 version)
	   'gga_x_pbefe',                & !  PBE for formation energies                     
	   'gga_x_cap',                  & !  Correct Asymptotic Potential
	   'gga_x_eb88',                 & !  Non-empirical (excogitated) B88 functional of Becke and Elliott 
	   'gga_c_pbe',                  & !  Del Campo, Gazquez, Trickey and Vela (PBE-like)    
	   'gga_k_absp3',                & !  gamma-TFvW form by Acharya et al [g = 1 - 1.513/N^0.35] 
	   'gga_k_absp4',                & !  gamma-TFvW form by Acharya et al [g = l = 1/(1 + 1.332/N^(1/3))] 
	   'gga_c_bmk',                  & !  Boese-Martin for kinetics                
	   'gga_c_tau_hcth',             & !  correlation part of tau-hcth             
	   'gga_c_hyb_tau_hcth',         & !  correlation part of hyb_tau-hcth         
	   'gga_x_beefvdw',              & !  BEEF-vdW exchange
	   'gga_xc_beefvdw',             & !  BEEF-vdW exchange-correlation 
	   'gga_x_pbetrans',             & !  radient-based interpolation between PBE and revPBE
	   'gga_x_chachiyo',             & !  Chachiyo exchange
	   'gga_k_vw',                   & !  von Weiszaecker functional 
	   'gga_k_ge2',                  & !  Second-order gradient expansion (l = 1/9) 
	   'gga_k_golden',               & !  TF-lambda-vW form by Golden (l = 13/45) 
	   'gga_k_yt65',                 & !  TF-lambda-vW form by Yonei and Tomishima (l = 1/5) 
	   'gga_k_baltin',               & !  TF-lambda-vW form by Baltin (l = 5/9) 
	   'gga_k_lieb',                 & !  TF-lambda-vW form by Lieb (l = 0.185909191) 
	   'gga_k_absp1',                & !  gamma-TFvW form by Acharya et al [g = 1 - 1.412/N^(1/3)] 
	   'gga_k_absp2',                & !  gamma-TFvW form by Acharya et al [g = 1 - 1.332/N^(1/3)] 
	   'gga_k_gr',                   & !  gamma-TFvW form by Gazquez and Robles 
	   'gga_k_ludena',               & !  gamma-TFvW form by Ludena 
	   'gga_k_gp85',                 & !  gamma-TFvW form by Ghosh and Parr 
	   'gga_k_pearson',              & !  Pearson
	   'gga_k_ol1',                  & !  Ou-Yang and Levy v.1
	   'gga_k_ol2',                  & !  Ou-Yang and Levy v.2 
	   'gga_k_fr',                   & !  Fuentealba & Reyes (B88 version) 
	   'gga_k_fr',                   & !  Fuentealba & Reyes (PW86 version) 
	   'gga_k_dk',                   & !  DePristo and Kress
	   'gga_k_perdew',               & !  Perdew                                
	   'gga_k_vsk',                  & !  Vitos, Skriver, and Kollar            
	   'gga_k_vjks',                 & !  Vitos, Johansson, Kollar, and Skriver 
	   'gga_k_ernzerhof',            & !  Ernzerhof 
	   'gga_k_lc94',                 & !  Lembarki & Chermette 
	   'gga_k_llp',                  & !  Lee, Lee & Parr 
	   'gga_k_thakkar',              & !  Thakkar 1992
	   'gga_x_wpbeh',                & !  short-range version of the PBE
	   'gga_x_hjs_pbe',              & !  HJS screened exchange PBE version
	   'gga_x_hjs_pbe_sol',          & !  HJS screened exchange PBE_SOL version 
	   'gga_x_hjs_b88',              & !  HJS screened exchange B88 version 
	   'gga_x_hjs_b97x',             & !  HJS screened exchange B97x version 
	   'gga_x_ityh',                 & !  short-range recipe for exchange GGA functionals
	   'gga_x_sfat',                 & !  short-range recipe for exchange GGA functionals
	   'gga_x_sg4',                  & !  Semiclassical GGA at fourth order
	   'gga_c_sg4',                  & !  Semiclassical GGA at fourth order
	   'gga_x_gg99',                 & !  Gilbert and Gill 1999
	   'gga_x_pbepow',               & !  PBE power
	   'gga_x_kgg99',                & !  Gilbert and Gill 1999 (mixed) 
	   'gga_xc_hle16',               & !  high local exchange 2016                 
	   'gga_c_scan_e0',              & !  GGA component of SCAN
	   'gga_c_gapc',                 & !  GapC
	   'gga_c_gaploc',               & !  Gaploc
	   'gga_c_zvpbeint',             & !  nother spin-dependent correction to PBEint
	   'gga_c_zvpbesol',             & !  nother spin-dependent correction to PBEsol       
	   'gga_c_tm_lyp',               & !  Takkar and McCarthy reparametrization 
	   'gga_c_tm_pbe',               & !  Thakkar and McCarthy reparametrization 
	   'gga_c_w94',                  & !  Wilson 94 (Eq. 25)
	   'gga_c_cs1',                  & !  A dynamical correlation functional
	   'gga_x_b88m',                 & !  Becke 88 reoptimized to be used with mgga_c_tau1 
	   'gga_k_pbe3',                 & !  Three parameter PBE-like expansion             
	   'gga_k_pbe4',                 & !  Four  parameter PBE-like expansion             
	   'gga_k_exp4',                 & !  Intermediate form between PBE3 and PBE4
	   'hyb_gga_x_n12_sx',           & !  N12-SX functional from Minnesota
	   'hyb_gga_xc_b97_1p',          & !  version of B97 by Cohen and Handy        
	   'hyb_gga_xc_pbe_mol0',        & !  PBEmol0             
	   'hyb_gga_xc_pbe_sol0',        & !  PBEsol0             
	   'hyb_gga_xc_pbeb0',           & !  PBEbeta0            
	   'hyb_gga_xc_pbe_molb0',       & !  PBEmolbeta0         
	   'hyb_gga_xc_pbe50',	         & !  PBE0 with 50% exx   
	   'hyb_gga_xc_b3pw91',	         & !  The original (ACM) hybrid of Becke
	   'hyb_gga_xc_b3lyp',           & !  The (in)famous B3LYP                  
	   'hyb_gga_xc_b3p86',           & !  Perdew 86 hybrid similar to B3PW91    
	   'hyb_gga_xc_o3lyp',           & !  hybrid using the optx functional
	   'hyb_gga_xc_mpw1k',           & !  mixture of mPW91 and PW91 optimized for kinetics 
	   'hyb_gga_xc_pbeh',            & !  aka PBE0 or PBE1PBE
	   'hyb_gga_xc_b97',             & !  Becke 97
	   'hyb_gga_xc_b97_1',           & !  Becke 97-1                               
	   'hyb_gga_xc_b97_2',           & !  Becke 97-2                               
	   'hyb_gga_xc_x3lyp',           & !  hybrid by Xu and Goddard 
	   'hyb_gga_xc_b1wc',            & !  Becke 1-parameter mixture of WC and PBE
	   'hyb_gga_xc_b97_k',           & !  Boese-Martin for Kinetics                
	   'hyb_gga_xc_b97_3',           & !  Becke 97-3                               
	   'hyb_gga_xc_mpw3pw',          & !  mixture with the mPW functional       
	   'hyb_gga_xc_b1lyp',           & !  Becke 1-parameter mixture of B88 and LYP         
	   'hyb_gga_xc_b1pw91',          & !  Becke 1-parameter mixture of B88 and PW91        
	   'hyb_gga_xc_mpw1pw',          & !  Becke 1-parameter mixture of mPW91 and PW91      
	   'hyb_gga_xc_mpw3lyp',         & !  mixture of mPW and LYP                
	   'hyb_gga_xc_sb98_1a',         & !  Schmider-Becke 98 parameterization 1a    
	   'hyb_gga_xc_sb98_1b',         & !  Schmider-Becke 98 parameterization 1b    
	   'hyb_gga_xc_sb98_1c',         & !  Schmider-Becke 98 parameterization 1c    
	   'hyb_gga_xc_sb98_2a',         & !  Schmider-Becke 98 parameterization 2a    
	   'hyb_gga_xc_sb98_2b',         & !  Schmider-Becke 98 parameterization 2b    
	   'hyb_gga_xc_sb98_2c',         & !  Schmider-Becke 98 parameterization 2c    
	   'hyb_gga_x_sogga11_x',        & !  Hybrid based on SOGGA11 form
	   'hyb_gga_xc_hse03',           & !  the 2003 version of the screened hybrid HSE
	   'hyb_gga_xc_hse06',           & !  the 2006 version of the screened hybrid HSE 
	   'hyb_gga_xc_hjs_pbe',         & !  HJS hybrid screened exchange PBE version 
	   'hyb_gga_xc_hjs_pbe_sol',     & !  HJS hybrid screened exchange PBE_SOL version 
	   'hyb_gga_xc_hjs_b88',         & !  HJS hybrid screened exchange B88 version 
	   'hyb_gga_xc_hjs_b97x',        & !  HJS hybrid screened exchange B97x version 
	   'hyb_gga_xc_cam_b3lyp',       & !  CAM version of B3LYP
	   'hyb_gga_xc_tuned_cam_b3lyp', & !  CAM version of B3LYP tuned for excitations
	   'hyb_gga_xc_bhandh',          & !  Becke half-and-half                              
	   'hyb_gga_xc_bhandhlyp',       & !  Becke half-and-half with B88 exchange            
	   'hyb_gga_xc_mb3lyp_rc04',     & !  B3LYP with RC04 LDA                   
	   'hyb_gga_xc_mpwlyp1m',        & !  MPW with 1 par. for metals/LYP                   
	   'hyb_gga_xc_revb3lyp',        & !  Revised B3LYP                         
	   'hyb_gga_xc_camy_blyp',       & !  BLYP with yukawa screening
	   'hyb_gga_xc_pbe0_13',         & !  PBE0-1/3            
	   'hyb_gga_xc_b3lyps',          & !  B3LYP* functional                     
	   'hyb_gga_xc_wb97',            & !  Chai and Head-Gordon
	   'hyb_gga_xc_wb97x',           & !  Chai and Head-Gordon                     
	   'hyb_gga_xc_lrc_wpbeh',       & !  Long-range corrected functional by Rorhdanz et al 
	   'hyb_gga_xc_wb97x_v',         & !  Mardirossian and Head-Gordon             
	   'hyb_gga_xc_lcy_pbe',         & !  PBE with yukawa screening
	   'hyb_gga_xc_lcy_blyp',        & !  BLYP with yukawa screening
	   'hyb_gga_xc_lc_vv10',         & !  Vydrov and Van Voorhis
	   'hyb_gga_xc_camy_b3lyp',      & !  B3LYP with Yukawa screening
	   'hyb_gga_xc_wb97x_d',         & !  Chai and Head-Gordon                     
	   'hyb_gga_xc_hpbeint',         & !  hPBEint             
	   'hyb_gga_xc_lrc_wpbe',        & !  Long-range corrected functional by Rorhdanz et al 
	   'hyb_gga_xc_b3lyp5',          & !  B3LYP with VWN functional 5 instead of RPA 
	   'hyb_gga_xc_edf2',            & !  Empirical functional from Lin, George and Gill
	   'hyb_gga_xc_cap0',            & !  Correct Asymptotic Potential hybrid
	   'hyb_gga_xc_lc_wpbe',         & !  Long-range corrected functional by Vydrov and Scuseria 
	   'hyb_gga_xc_hse12',           & !  HSE12 by Moussa, Schultz and Chelikowsky 
	   'hyb_gga_xc_hse12s',          & !  Short-range HSE12 by Moussa, Schultz, and Chelikowsky 
	   'hyb_gga_xc_hse_sol',         & !  HSEsol functional by Schimka, Harl, and Kresse 
	   'hyb_gga_xc_cam_qtp_01',      & !  CAM-QTP(01): CAM-B3LYP retuned using ionization potentials of water 
	   'hyb_gga_xc_mpw1lyp',         & !  Becke 1-parameter mixture of mPW91 and LYP       
	   'hyb_gga_xc_mpw1pbe',         & !  Becke 1-parameter mixture of mPW91 and PBE       
	   'hyb_gga_xc_kmlyp',           & !  Kang-Musgrave hybrid                  
	   'hyb_gga_xc_b5050lyp'         & !  Like B3LYP but more exact exchange    
	   /
      
      data lxcnumbs /  &
              1,  2,  3,  4,  5,  6,  7,  8,  9, 10,& 
             11, 12, 13, 14, 15, 16, 17, 18, 19, 20,& 
             21, 22, 23, 24, 25, 26, 27, 28, 29, 30,& 
             31, 43, 50, 51,259,287,289,532,536,537,& 
            538,546,547,548,549,550,551,552,554,573,& 
            574,577,578,579,580, 32, 33, 34, 35, 38,&
             39, 40, 41, 44, 45, 46, 47, 48, 49, 52,&
             53, 54, 55, 56, 57, 58, 59, 60, 61, 62,&
             63, 65, 66, 67, 68, 69, 70, 71, 79, 80,&
             82, 83, 84, 85, 86, 87, 88, 89, 90, 91,&
             92, 93, 94, 95, 96, 97, 98, 99,100,101,&
            102,103,104,105,106,107,108,109,110,111,&
            112,113,114,115,116,117,118,119,120,121,&
            122,123,124,125,126,127,128,129,130,131,&
            132,133,134,135,136,137,138,139,140,141,&
            142,143,144,145,146,147,148,149,150,151,&
            152,153,154,155,156,157,158,159,160,161,&
            162,163,164,165,166,167,170,173,174,175,&
            182,183,184,185,186,187,188,189,190,191,&
            192,193,194,195,196,197,198,199,200,246,&
            255,258,262,265,270,271,272,277,278,280,&
            281,283,285,286,291,298,500,501,502,503,&
            504,505,506,507,508,509,510,511,512,513,&
            514,515,516,517,518,519,520,521,522,523,&
            524,525,526,527,528,529,530,533,534,535,&
            539,544,545,553,555,556,557,558,559,560,&
            561,565,570,595,596,597, 81,266,273,274,&
            275,276,290,401,402,403,404,405,406,407,&
            408,410,411,412,413,414,415,416,417,418,&
            419,420,421,422,423,424,425,426,427,428,&
            429,430,431,432,433,434,435,436,437,453,&
            454,455,456,459,463,464,465,466,467,468,&
            469,470,471,472,473,475,476,477,478,479,&
            480,481,482,483,484,485,572/

      character*80 :: header

      real (PREC) :: hni,hni_p
      real (PREC) :: homolevl,homoallign,ffield,fgrad,zcutoff,harm2xy,gammaf,z1atmass,z2atmass

      real (PREC) :: facmul,exlorb,exlcoul,exlexp,alphaf,diver,trelax,tortho,trayl,tmomen,tlagra,ttoten,sflagra,dflagra,&
           vkt,vnt,evt,epott,etot,virrat,etotFN,ictot,iready,hallign

      real (PREC), dimension(10) :: elect,electDA,total,totalDA

      real (PREC) :: enkin,ennucel,encoul,enexch,entot,encouldft,enexchdft,edftex,edftcorr,alegk0,epspot,ompot,apot,&
           v0pot,hook,r_p,z1_p,z2_p,rinf_p,date_p,time_p

      integer :: mpot, nsimp

      real (PREC), dimension(12) :: calp
      real (PREC), dimension(60) :: eza1,eza2,co1,co2,area,wstorthog,qxm,occ,sign,demax,vk,vn,eng,engi,engt,pnc,div

      real (PREC), dimension(maxbasis) :: fngau2,shngau,primexp,coeff

      real (PREC), dimension(60,maxbasis) :: primcoef

      real (PREC), dimension(maxbasis,maxbasis) :: ovl

end module commons8
