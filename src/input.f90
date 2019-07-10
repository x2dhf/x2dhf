
module input_m
  use params
  use xc_f90_types_m
  use xc_f90_lib_m

    integer :: i,icompLAdd,icompLEnc,icompLExp,id1,id2,iform_t,ig,ihit,inpiexit,&
         inzero,iopenshell,iorb,iput,iput1,isum,isum0,isum1,itmp,itmp1,itmp2,itotq,&
         izz1,izz2,j,j1,mgi,mt,n,nbsym,next,nmethods,nmutot,nlabels,no,nonortho,maxflags

    real (PREC) :: clo,cloe,co12,fmfield,ftmp,ftmp1,ftmp2,tmp1,tmp2,totchar,totq,z1t,z2t 

    parameter (nmethods=5,nlabels=43,maxflags=40)
    character*4 cdftext(10),cdftcorrt(10)
    character*8 clabel,clabel1,clabel2,char8
    character*8 labellc(nlabels),cmethod(nmethods)
    character*8 sigma,pi,delta,phi,space,angstrom,endl,guttmp
    character*30 char30
    
    integer, dimension(maxflags) :: id

    TYPE(xc_f90_pointer_t) :: xc_func
    TYPE(xc_f90_pointer_t) :: xc_info

    
    data labellc /'break','config','conv','debug','dft',&
         'exchio','fefield','fermi','fix','fixorb',&
         'fixcoul','fixexch','gauss','grid','homo',&
         'initial','inout','interp','lcao','lagra',&
         'mcsor','method','multipol','nuclei','omega',&
         'orbpot','order','potgsz','potgszg','poth3',&
         'potkh','print','prtevery','scf','scfexch',&
         'sor','sormcsor','stop','title','xalpha',&
         'potharm', 'pothook','potsap'/


    data cdftext /'lda','b88','pw86','pw91','x','x','x','x','x','x'/
    data cdftcorrt /'lyp','vwn','x','x','x','x','x','x','x','x'/

    data angstrom/'angstrom'/,cmethod/'hf','oed','hfs','dft','scmc'/,&
         endl/'end'/,space/'    '/
    data sigma/'sigma'/,pi/'pi'/,delta/'delta'/,phi/'phi'/

    data inpiexit/-99999/

end module input_m

