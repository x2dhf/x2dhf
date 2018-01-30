module scf
  use params
      real (PREC), dimension(1,1) :: h
      real (PREC), dimension(1200) :: cmulti
      real (PREC), dimension(1830) :: pnx

      real (PREC), dimension(3660) :: excdi,excqu,excoc,exche,exc5,exc6,exc7,exc8,engo
      real (PREC), dimension(7200) :: gec,gca,gxa,gcb,gxb,gc,gx
      real (PREC), dimension(10) :: ssf,ssf_p,diag,rgrid,ovforb,ovfcoul,ovfexch
      real (PREC), dimension(60,60) :: sflagrat
end module scf
