module card
  use params

  parameter (i40=40, i80=80)

  integer :: jrec,jump


  integer, dimension(i40) :: istrt,inumb
  character*1, dimension(i80) :: ia
  character*80 :: header
end module card

