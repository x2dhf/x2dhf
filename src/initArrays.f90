! ### initArrays ###
!
!     Setup various arrays and variables

subroutine initArrays (cw_suppl,cw_sor)
  use params
  use commons8

  integer, dimension(*) ::  cw_sor
  real (PREC), dimension(*) :: cw_suppl

  ! length=nsuppl*nni*mxnmu
  ! allocate(cw_suppl(length),stat=error)
  ! if (stat.ne.0) then
  !    print*,"initArrays: Error! Couldn't allocate memory for array cw_suppl, dim=",length
  ! endif

  !     initialize arrays within common blocks (except for differentiating
  !     arrays which are defined in prepdiff) and check input data

  call initCBlocks

  !     employ input data to determine size and division of cw arrays
  !     initialize address arrays

  call initAddr

  !     Initialises various supplementary arrays of case-dependent lengths

  ! i4b(1)=borb
  ! i4b(2)=bpot
  ! i4b(3)=d
  ! i4b(4)=e
  ! i4b(5)=f0
  ! i4b(6)=f1
  ! i4b(7)=f2
  ! i4b(8)=f3
  ! i4b(9)=f4
  ! i4b(10)=g
  ! i4b(11)=wjac1
  ! i4b(12)=wjac2
  ! i4b(13)=wgt1
  ! i4b(14)=wgt2


  call initSuppl (cw_suppl(i4b( 1)),cw_suppl(i4b( 2)),cw_suppl(i4b( 3)),cw_suppl(i4b( 4)), &
       cw_suppl(i4b( 5)),cw_suppl(i4b( 6)),cw_suppl(i4b( 7)),cw_suppl(i4b( 8)),   &
       cw_suppl(i4b( 9)),cw_suppl(i4b(10)),cw_suppl(i4b(11)),cw_suppl(i4b(12)),   &
       cw_suppl(i4b(13)),cw_suppl(i4b(14)) )

  !     prepare meshes for all subgrids

  call initMesh (cw_sor)

  !     calculate weights of exchange contributions to the total energy
  !     expression for the given open/closed shell scf case


  call initExWeights

end subroutine initArrays
