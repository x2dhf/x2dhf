! ***************************************************************************
! *                                                                         *
! *   Copyright (C) 1996 Leif Laaksonen, Dage Sundholm                      *
! *   Copyright (C) 1996-2010 Jacek Kobus <jkob@fizyka.umk.pl>              *
! *                                                                         *
! *   This program is free software; you can redistribute it and/or modify  *
! *   it under the terms of the GNU General Public License version 2 as     *
! *   published by the Free Software Foundation.                            *
! *                                                                         *
! ***************************************************************************
! ### initMesh ####
!
!     Establishes meshes for each subgrid and determines coefficients
!     of the interpolation polynomials needed in fill arrays.
!
!     Warning! At present at most 3 grids are allowed. See mesh for
!     further detailes.

!     In version 2.0 only single grid is allowed!

module initMesh_m
  implicit none
contains
  subroutine initMesh (cw_sor)
    use params
    use discret
    use memory
    use commons8

    use mesh_m
    use meshsize_m
    use interpolq_m

    implicit none
    integer :: ib1,ib2,ib3,ib4,ib5,ibeg1,ibeg2,ibeg3,ibeg4,ibeg5, &
         ig,igs,imcase,imut,init,is34,is5,mxnmu16,mxsmu,ngrid16

    integer, dimension(length6) :: cw_sor
    ! FIXME
    !  integer, dimension(*) :: cw_sor

    !     Array cw_sor is used to store index arrays needed to perform
    !     sor relaxation on normal and extended grids. Two cases have to be
    !     considered since all ngrids-1 subgrids can either be interior or
    !     last ones depending on the particular orbital or potential
    !     function relaxed (outer diffused orbitals would need, say, a grid
    !     consisting of the 3 subgrids while for the proper description of the
    !     the core, very much confined orbitals only one subgrid would sufice).

    !     Division of the cw_sor array is as folows:

    !     ind1   - relaxation ordering on the extended grid
    !     ind2   - relaxation ordering on the normal grid
    !     ind1ex - extrapolation ordering on the extended grid
    !     ind2ex - extrapolation ordering on the normal grid
    !     ind3   - number of points for relaxation and extrapolation
    !     mxsize8==(nni+8)*(mxnmu+8)
    !     mxnmu8=mxnmu+8

    !     case 1 - all ngrids subgrids are considered as the last ones;
    !              in effect the last 4 columns of orbitals and potentials
    !	       are not relaxed two arrays of the indx1 type

    !     case 2 - all (ngrids-1) subgrids are considered as the interior ones

    !     cw_sor consists of
    !     - ngrids     arrays of indx1 type (case 1);
    !     - (ngrids-1) arrays of indx1 type (case 2);
    !                  the length of each array is determined by ngsize
    !
    !     - ngrids     arrays of indx2 type (case 1);
    !     - (ngrids-1) arrays of indx2 type (case 2);
    !		   the length of each array is determined by ngsize
    !
    !     - ngrids	   arrays of indx6a (case 1)
    !     - (ngrids-1) arrays of indx6a (case 2)
    !     - ngrids	   arrays of indx6b (case 1)
    !     - (ngrids-1) arrays of indx6b (case 2)
    !     - ngrids	   arrays of indx7 (case 1)
    !     - (ngrids-1) arrays of indx7 (case 2)
    !  		   the length of each array is determined by nmu

    !     ingr1(k,ig) and ingr2(k,ig),ig=1,ngrids, k=1,..,4, contain number of
    !     grids points to be relaxed in each subgrid (k=1) and number of points
    !     for which extrapolation is to be carried out (k=2 - ngrd6a, k=3 - ngrd6b,
    !     k=4 - ngrd7); ingr1 and ingr2 correspond to the extended and normal
    !     case, respectively.

    !     i6b(1) - begining of ind1
    !     i6b(2) - begining of ind2
    !     i6b(3) - begining of index1 (indx6a)
    !     i6b(4) - begining of index2 (indx6b)
    !     i6b(5) - begining of index3 (indx7)

    !     more space is reserved for the ordering arrays to avoid problems when
    !     more than 3 subgrid would be used in the future

    !     determine the largest subgrid length in mu variable

    mxsmu=0
    do ig=1,ngrids
       if (mxsmu.lt.nmu(ig)) mxsmu=nmu(ig)
    enddo

    is34=ngrids*mxsmu
    is5 =ngrids*nni
    mxnmu16=mxnmu+16
    ngrid16=nni*mxnmu16

    i6b(1)=1
    i6b(2)=i6b(1)+2*ngrid16
    i6b(3)=i6b(2)+2*ngrid16
    i6b(4)=i6b(3)+2*is34
    i6b(5)=i6b(4)+2*is5

    !     i6b(4)=i6b(3)+is34
    !     i6b(5)=i6b(4)+is5

    !     determine ordering of grid points assuming the subgrids are the last
    !     ones (grid points in the last 4 columns are excluded from relaxation
    !     process

    imcase=1
    ibeg1=i6b(1)
    ibeg2=i6b(2)
    ibeg3=i6b(3)
    ibeg4=i6b(4)
    ibeg5=i6b(5)
    imut=0
    init=0
    igs=0

    do ig=1,ngrids
       ib1=ibeg1+igs
       ib2=ibeg2+igs
       ib3=ibeg3+imut
       ib4=ibeg4+imut
       ib5=ibeg5+init
       iadext(ig)=ib1
       iadnor(ig)=ib2
       iadex1(ig)=ib3
       iadex2(ig)=ib4
       iadex3(ig)=ib5
       !     call mesh(imcase,ig,ib1,ib2,ib3,ib4,ib5,cw_sor(ib1),cw_sor(ib2),cw_sor(ib3),cw_sor(ib4),cw_sor(ib5))
       call meshsize(imcase,ig)
       call mesh(imcase,ig,cw_sor(ib1),cw_sor(ib2),cw_sor(ib3),cw_sor(ib4),cw_sor(ib5))
       imut=imut+nmu(ig)
       init=init+nni
       igs=igs+ngsize(ig)
    enddo

    !     determine ordering of grid points assuming the subgrids are
    !     the interior ones (only boundary points coressponding to
    !     mu=1, ni=1 and ni=nni are excluded from the relaxation process)

    imcase=2
    ibeg1=i6b(1)+ngrid16
    ibeg2=i6b(2)+ngrid16
    ibeg3=i6b(3)+is34
    ibeg4=i6b(4)+is34
    ibeg5=i6b(5)+is5
    imut=0
    init=0
    igs=0
    do ig=1,ngrids-1
       ib1=ibeg1+igs
       ib2=ibeg2+igs
       ib3=ibeg3+imut
       ib4=ibeg4+imut
       ib5=ibeg5+init
       iadext(ig+ngrids)=ib1
       iadnor(ig+ngrids)=ib2
       iadex1(ig+ngrids)=ib3
       iadex2(ig+ngrids)=ib4
       iadex3(ig+ngrids)=ib5
       call mesh(imcase,ig,cw_sor(ib1),cw_sor(ib2),cw_sor(ib3),cw_sor(ib4),cw_sor(ib5))
       imut=imut+nmu(ig)
       init=init+nni
       igs=igs+ngsize(ig)
    enddo

  end subroutine initMesh
end module initMesh_m
