SUBROUTINE rpn2_qad(mx,meqn,maux,mbc, idir, iface, &
     qf,qc,auxf,auxc,amdq,apdq)

  IMPLICIT NONE
  INTEGER mx, meqn, maux, mbc, idir, iface
  DOUBLE PRECISION qf(meqn,mx), qc(meqn,mx)
  DOUBLE PRECISION auxf(maux,mx), auxc(maux,mx)

  DOUBLE PRECISION amdq(meqn,mx)  
  DOUBLE PRECISION apdq(meqn,mx)

!!  # Must use edge velocities
  integer color_equation
  common /eqn_comm/ color_equation

  !! automatically allocated
  DOUBLE PRECISION qvc(meqn), qvf(meqn)
  DOUBLE PRECISION auxvc_center(maux), auxvf_center(maux)
  DOUBLE PRECISION fluxf(meqn), fluxc(meqn)

  DOUBLE PRECISION fd
  INTEGER m,i, iface_cell, sgn

  !! idir refers to direction of the Riemann solver
  !! iface refers to the face of the Cartesian grid

  DO i = 1,mx
    do m = 1,maux
        auxvf_center(m) = auxf(m,i)
        auxvc_center(m) = auxc(m,i)
    end do

    do m = 1,meqn
        qvf(m) = qf(m,i)
        qvc(m) = qc(m,i)
    end do


    !! Compute flux in fine grid ghost cell at left edge of grid
    if (idir .eq. 0) then
        iface_cell = 1-iface    !! Left edge of Cartesian grid is right edge of ghost cell
    else
        iface_cell = 5-iface    !! swap bottom and top edges
    endif

    if (color_equation .eq. 0) then
        call rpn2_cons_update_manifold(meqn,maux,idir,iface_cell,qvf,auxvf_center,fluxf)
        call rpn2_cons_update_manifold(meqn,maux,idir,iface_cell,qvc,auxvc_center,fluxc)
    else
        call rpn2_cons_update_zero(meqn,maux,idir,iface_cell,qvf,auxvf_center,fluxf)
        call rpn2_cons_update_zero(meqn,maux,idir,iface_cell,qvc,auxvc_center,fluxc)
    endif        

    do m = 1,meqn
        fd = fluxf(m) - fluxc(m)
        apdq(m,i) = 0.5*fd  
        amdq(m,i) = 0.5*fd
    end do

  ENDDO

END SUBROUTINE rpn2_qad




SUBROUTINE  rpn2_cons_update_manifold(meqn,maux,   &
                                      idir, iface, q,  &
                                      auxvec_center,  flux)

  IMPLICIT NONE

  INTEGER meqn,maux,idir, m, iface
  DOUBLE PRECISION q(meqn), flux(meqn)
  DOUBLE PRECISION auxvec_center(maux)
  DOUBLE PRECISION urot, g
  INTEGER k

  !! # Get cell-centered velocity projected to face
  !! # 'iface' (in 0,1,2,3)

  urot = auxvec_center(2+iface)
  if (iface .eq. 3) then
    if (abs(urot) .gt. 1e-8) then
!!        write(6,*) 'urot = ', urot
    endif
  endif

  !! Scaling for fluxes (since this is supposed to look like a scaled apdq, amdq)
  g = auxvec_center(6+iface)  !! Edge length


  !! #  f(q) = (n dot u)*q
  DO m = 1,meqn
     !! # Don't multiply by edgelength (scaling is done elsewhere)
     flux(m) = g*urot*q(m)
  ENDDO

END SUBROUTINE rpn2_cons_update_manifold

SUBROUTINE  rpn2_cons_update_zero(meqn,maux,   &
                                  idir, iface, q,  &
                                  auxvec_center,  flux)

  IMPLICIT NONE
  INTEGER meqn,maux,idir, iface
  DOUBLE PRECISION q(meqn), flux(meqn)
  DOUBLE PRECISION auxvec_center(maux)

  INTEGER m

  !! #  f(q) = (n dot u)*q
  DO m = 1,meqn
     !! # No flux function available for equations in non-conservative form
     flux(m) = 0
  ENDDO

END SUBROUTINE rpn2_cons_update_zero
