module evalGauss
  use params
  use discret
  use commons8
  implicit none

contains
  subroutine evaluate(bf)

    integer :: ic1, igp, imu, in, inioff, ipb, l1, m1
    real (PREC) :: costh1, d1, expt, fnorm, r1, z, psi1
    real (PREC), dimension(:,:) :: bf
    real (PREC), external :: plegendg

    ! Sanity check
    if(size(bf,1) .ne. nni*mxnmu .or. size(bf,2) .ne. npbasis) then
       write (*,*) 'bf matrix needs to be initialized, got size ',size(bf,1),' x ',size(bf,2)
    end if

    !     loop over basis set functions
    do ipb=1,npbasis
       ! Index of the center of the primitive
       !   ic1=1 basis function at centre Z1
       !   ic1=2 basis function at bond centre
       !   ic1=3 basis function at centre Z2
       ic1 =icgau(ipb)

       ! Primitive exponent
       d1  =primexp(ipb)
       ! (l,m) values
       l1  =lprim(ipb)
       m1  =abs(mprim(ipb))

       ! Normalization: primitive times spherical harmonic
       fnorm=fngau2(ipb)*shngau(ipb)

       ! Loop over grid
       do imu=1,mxnmu
          inioff=(imu-1)*nni
          do in=1,nni
             igp=inioff+in

             ! for each grid point, i.e. for (vmu(imu),vni(ini))
             ! determine its distance |r_i| from the nucleus Z_i and
             ! cosine of the polar angle costh_i between the z axis
             ! and the vector r_i

             z=(r/2.0_PREC)*vxi(imu)*veta(in)

             if (ic1.eq.1) then
                ! Center on Z1
                r1=(r/2.0_PREC)*(vxi(imu)+veta(in))
                z=z+r/2.0_PREC

             elseif (ic1.eq.3) then
                ! Center on Z2
                r1=(r/2.0_PREC)*(vxi(imu)-veta(in))
                z=z-r/2.0_PREC

             elseif (ic1.eq.2) then
                ! Bond center
                r1=(r/2.0_PREC)*sqrt(vxisq(imu)+vetasq(in)-1.0_PREC)

             else
                write (*,*) 'invalid center ',ic1
                stop
             endif

             if (r1.lt.precis) then
                costh1=0.0_PREC
             else
                costh1=z/r1
                ! Sanity check
                if(costh1 < -1.0_PREC) costh1 = -1.0_PREC
                if(costh1 > 1.0_PREC) costh1 = 1.0_PREC
             endif

             ! Calculate exponential
             expt=exp(-d1*r1*r1)

             ! Value of basis function is
             if (r1.lt.precis) then
                if(l1.gt.0) then
                   psi1=0.0_PREC
                else
                   psi1=fnorm*expt*plegendg(l1,m1,costh1)
                end if
             else
                psi1=fnorm*r1**l1*expt*plegendg(l1,m1,costh1)
             endif

             ! Store value
             bf(igp,ipb)=psi1
          enddo
       enddo
    enddo
  end subroutine evaluate
end module evalGauss
