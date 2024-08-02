! SPDX-License-Identifier: GPL-2.0-or-later

! Copyright (C) 1996-2024  Jacek Kobus 

module data4II
  use params

    integer (KIND=IPREC) :: i,icompLAdd,icompLEnc,icompLExp,id1,id2,iform_t,ihit,inpiexit,&
         inzero,iopenshell,iorb,iput,iput1,isum,isum0,isum1,itmp,itmp1,itmp2,itotq,&
         izz1,izz2,j,j1,mt,n,nbsym,next,nmethods,nmutot,nlabels,no,nonortho,maxflags

    real (PREC) :: clo,cloe,co12,fmfield,ftmp,ftmp1,ftmp2,tmp1,tmp2,totchar,totq,z1t,z2t 

    parameter (nmethods=6,nlabels=70,maxflags=40,maxdfts=2)
    character*4 cdftext(maxdfts),cdftcorrt(maxdfts)
    character*8 clabel,clabel1,clabel2,char8
    character*8 labellc(nlabels),cmethod(nmethods)
    character*8 sigma,pi,delta,phi,space,angstrom,endl,guttmp
    character*30 char30
    character*1, dimension(80) :: title
    
    integer (KIND=IPREC),dimension(maxflags) :: id

    data labellc /&
         'break'    , 'config'   , 'conv'     , 'debug'    , 'dft'     ,&
         'fefield'  , 'fermi'    , 'fix'      , 'fixorb'   , 'mmoments',&
         'initial'  , 'inout'    , 'interp'   , 'lcao'     , 'lagra'   ,&
         'mcsor'    , 'method'   , 'multipol' , 'nuclei'   , 'omega'   ,&
         'orbpot'   , 'order'    , 'potgsz'   , 'potgszg'  , 'potharm3',&
         'potkh'    , 'print'    , 'prtevery' , 'scf'      , 'scfexch' ,&
         'sor'      , 'stop'     , 'title'    , 'xalpha'   , 'potharm2',&
         'pothooke' , 'potsap'   , 'lxcpolar' , 'detnan'   , 'fixnan'  ,&
         'mcsor-ce' , 'mcsor-o'  , 'sor4orb'  , 'sor4pot'  , 'out4pair',&
         'slowexch' , 'chktoten' , 'coulexch' , 'fastscf'  , 'ldasap'  ,&
         'omegaopt' , 'tail'     , 'kinpot'   , 'lxcpolar' , 'densthld',&
         'potcoul2' , 'potcoul3' , 'lm'       , 'lm0'      , 'lm1'     ,&
         'lm2'      , 'altsweep' , 'omegaz'   , 'intracul' , 'extracul',&
         'fixpot'   , 'gauss'    , 'grid'     , 'homo'     , 'plot'    /

    data cdftext /'lda','b88'/
    data cdftcorrt /'lyp','vwn'/

    data angstrom/'angstrom'/,cmethod/'hf','oed','hfs','dft','scmc','ted'/,&
         endl/'end'/,space/'    '/
    data sigma/'sigma'/,pi/'pi'/,delta/'delta'/,phi/'phi'/

    data inpiexit/-99999/

end module data4II

