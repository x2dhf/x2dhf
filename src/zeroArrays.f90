! ### zeroArrays ###
!
!     Zeroise main arrays.
module zeroArrays_m
  implicit none
contains

  subroutine zeroArrays (cw_orb,cw_coul,cw_exch,cw_suppl,cw_sctch,cw_sor)
    use params
    use memory
    use scf
    use zeroArray_m

    !  use commons8

    implicit none

    integer :: n
    integer, dimension(*) :: cw_sor
    real (PREC), dimension(*) ::  cw_orb,cw_coul,cw_exch,cw_suppl,cw_sctch

    call zeroArray(length1,cw_orb)
    call zeroArray(length2,cw_coul)
    call zeroArray(length3,cw_exch)
    stop
    call zeroArray(length4,cw_suppl)
    call zeroArray(length5,cw_sctch)

    do n=1,length6
       cw_sor(n)=0
    enddo

    !     off-diagonal Lagrange multipliers must be zero (see fock)

    n=maxorb*maxorb
    call zeroArray(n,engo)

  end subroutine zeroArrays
end module zeroArrays_m
