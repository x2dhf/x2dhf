! ### prepdiff1 ###
!    This routine initializes arrays used for differentiation needed in nuclder

subroutine prepdiff1 
  use params
  use discret
  use commons8

  implicit none
  integer :: i,ib,ie,ig,k
  real (PREC) :: w1,w2
  real (PREC), dimension(9) :: aa1,aa2,a1,a2

  !  derivative coeff. from 8th-order Stirling interpolation formula

  data aa1/3.0_PREC,-32.0_PREC,168.0_PREC,-672.0_PREC,0.0_PREC,672.0_PREC, &
       -168.0_PREC,32.0_PREC,-3.0_PREC/
  data aa2/-9.0_PREC,128.0_PREC,-1008.0_PREC,8064.0_PREC,-14350.0_PREC,   &
       8064.0_PREC,-1008.0_PREC,128.0_PREC,-9.0_PREC /

  !  initialization of the array performing differentiation
  !  based on the 9-point interpolating formula
  
  w1=1.0_PREC/(84.0_PREC)
  w2=1.0_PREC/(504.0_PREC)
  
  !   initialize dni1, dni2, dmu1 and dmu2
  
  do ig=1,ngrids
     do  i=1,4
        dmu2(i,ig)=w2*aa2(i)/(hmu(ig)*hmu(ig))
        dmu1(i,ig)=w1*aa1(i)/hmu(ig)
     enddo
  enddo
  
  do i=1,4
     dni2(i)=w2*aa2(i)/(hni*hni)
     dni1(i)=w1*aa1(i)/hni
  enddo

  !  initialize arrays used in matrix differentiation

  !  initialize dmu

  ib=1
  do ig=1,ngrids
     do k=1,9
        a2(k)=w2*aa2(k)/(hmu(ig)*hmu(ig))
        a1(k)=w1*aa1(k)/hmu(ig)
     enddo

     ie=iemu(ig)
     do i=ib,ie
        do k=1,9
           !          dmu(k,i)=a2(k)+borb(1,i)*a1(k)
           dmu(k,i)=a1(k)
        enddo
     enddo
     ib=ie+1
  enddo

  !  initialize dni 
  do k=1,9
     a2(k)=w2*aa2(k)/(hni*hni)
     a1(k)=w1*aa1(k)/hni
  enddo
  
  do i=1,nni
     do k=1,9
        !        dni(k,i)=a2(k)+d(i,1)*a1(k)
        dni(k,i)=a1(k)
     enddo
  enddo
  
end subroutine prepdiff1
      
