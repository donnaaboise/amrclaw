!! For each coarse-fine interface, a Riemann problem between an inner
!! ghost cell value on the fine grid and cell value in the adjacent coarse
!! cell must be solved and added to corresponding location in
!! **node(ffluxptr, mptr)** for conservative fix later
!!
!! -------------------------------------------------------------
!!

SUBROUTINE qad(valbig,mitot,mjtot,nvar, &
     svdflx,qc1d,lenbc,lratiox,lratioy,hx,hy, &
     maux,auxbig,auxc1d,delt,mptr)

  USE amr_module, only : timemult, nghost, max1d, maxaux, mwaves,nestlevel,rnode, &
                          auxtype, node, method
  IMPLICIT NONE

  INTEGER mitot, mjtot, nvar, lenbc, lratiox, lratioy,maux
  INTEGER mptr
  DOUBLE PRECISION hx,hy,delt

  DOUBLE PRECISION, target :: valbig(nvar,mitot,mjtot)
  DOUBLE PRECISION, target :: auxbig(maux,mitot,mjtot)
  DOUBLE PRECISION qc1d(nvar,lenbc)
  DOUBLE PRECISION svdflx(nvar,lenbc)
  DOUBLE PRECISION auxc1d(maux,lenbc)

  double precision, pointer :: q(:,:,:)
  double precision, pointer :: aux(:,:,:)

  !!
  !! ::::::::::::::::::::::::::: QAD ::::::::::::::::::::::::::::::::::
  !!  are added in to coarse grid value, as a conservation fixup.
  !!  Done each fine grid time step. If source terms are present, the
  !!  coarse grid value is advanced by source terms each fine time step too.

  !!  No change needed in this sub. for spherical mapping: correctly
  !!  mapped vals already in bcs on this fine grid and coarse saved
  !!  vals also properly prepared
  !!
  !! Side 1 is the left side of the fine grid patch.  Then go around clockwise.
  !! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  !!
  !!      # local storage
  !!      # note that dimension here are bigger than dimensions used
  !!      # in rp2, but shouldn't matter since wave is not used in qad
  !!      # and for other arrays it is only the last parameter that is wrong
  !!      #  ok as long as meqn, mwaves < maxvar

  INTEGER max1dp1
  PARAMETER (max1dp1 = max1d+1)
  DOUBLE PRECISION ql(nvar,max1dp1)
  DOUBLE PRECISION qr(nvar,max1dp1)
  DOUBLE PRECISION wave(nvar,mwaves,max1dp1)
  DOUBLE PRECISION s(mwaves,max1dp1)
  DOUBLE PRECISION amdq(nvar,max1dp1)
  DOUBLE PRECISION apdq(nvar,max1dp1)
  DOUBLE PRECISION :: auxlbig(maxaux*max1dp1)
  DOUBLE PRECISION :: auxrbig(maxaux*max1dp1)

!!  DOUBLE PRECISION, POINTER :: auxl(:,:)
!!  DOUBLE PRECISION, POINTER :: auxr(:,:)

  !!
  !!  WARNING: auxl,auxr dimensioned at max possible, but used as if
  !!  they were dimensioned as the real maux by max1dp1. Would be better
  !!  of course to dimension by maux by max1dp1 but this wont work if maux=0
  !!  So need to access using your own indexing into auxl,auxr.

  INTEGER iaddaux
  iaddaux(iaux,i) = (i-1)*maux + iaux

  !!
  !!      aux is auxiliary array with user parameters needed in Riemann solvers
  !!          on fine grid corresponding to valbig
  !!      auxc1d is coarse grid stuff from around boundary, same format as qc1d
  !!      auxl, auxr are work arrays needed to pass stuff to rpn2
  !!      maux is the number of aux variables, which may be zero.
  !!

  DOUBLE PRECISION tgrid
  INTEGER nc, nr, level, index, l
  INTEGER i,j,ma,ivar, lind, ncrse, ic, jc, jfine, influx
  INTEGER iaux, ifine

  INTEGER mx,my,mbc,meqn
  DOUBLE PRECISION dt, dx, dy

  tgrid = rnode(timemult, mptr)
  nr = mitot-2*nghost
  nc = mjtot-2*nghost
  level = node(nestlevel, mptr)
  index = 0

  mbc = nghost
  mx = nr
  my = nc
  meqn = nvar
  dt = delt
  dx = hx
  dy = hy

  q(1:meqn,1-mbc:mx+mbc,1-mbc:my+mbc) => valbig

  if (maux .gt. 0) then
    aux(1:maux,1-mbc:mx+mbc,1-mbc:my+mbc) => auxbig
!!    auxl(1:maux,1:max1dp1) => auxlbig
!!    auxr(1:maux,1:max1dp1) => auxrbig
  endif


  !! --------
  !!  side 1
  !! --------
  !!
  DO j = 1, my
     IF (maux .GT. 0) THEN
        DO ma = 1,maux
           IF (auxtype(ma) .EQ. "xleft") THEN
              !! # Assuming velocity at left-face, this fix
              !! # preserves conservation in incompressible flow:
                auxlbig(iaddaux(ma,j+1)) = aux(ma,1,j)
           ELSE
              !! # Normal case -- we set the aux arrays
              !! # from the cell corresponding  to q
                auxlbig(iaddaux(ma,j+1)) = aux(ma,0,j)
           ENDIF
        ENDDO
     ENDIF
     DO ivar = 1, meqn
        ql(ivar,j+1) = q(ivar,0,j)
     ENDDO
  ENDDO

  lind = 0
  ncrse = my/lratioy
  DO jc = 1, ncrse
     index = index + 1
     DO l = 1, lratioy
        lind = lind + 1
        IF (maux .GT. 0) THEN
           DO ma = 1,maux
                auxrbig(iaddaux(ma,lind)) = auxc1d(ma,index)
           ENDDO
        ENDIF
        DO ivar = 1, meqn
           qr(ivar,lind) = qc1d(ivar,index)
        ENDDO
     ENDDO
  ENDDO


  CALL rpn2(1,max1dp1-2*mbc,nvar,mwaves,maux,mbc, &
       my+1-2*mbc,ql,qr,auxlbig,auxrbig,wave,s,amdq,apdq)
  !!
  !! we have the wave. for side 1 add into sdflxm
  !!

  influx = 0
  DO j = 1, nc/lratioy
     influx  = influx + 1
     jfine = (j-1)*lratioy
     DO  ivar = 1, meqn
        DO  l = 1, lratioy
           svdflx(ivar,influx) = svdflx(ivar,influx) &
                + amdq(ivar,jfine+l+1) * hy * delt &
                + apdq(ivar,jfine+l+1) * hy * delt
        ENDDO
     ENDDO
  ENDDO

  !! --------
  !!  side 2
  !! --------
  !!
  IF (mjtot .EQ. 2*nghost+1) THEN
     !! # a single row of interior cells only happens when using the
     !! # 2d amrclaw code to do a 1d problem with refinement.
     !! # (feature added in Version 4.3)
     !! # skip over sides 2 and 4 in this case
     go to 299
  ENDIF

       do 210 i = nghost+1, mitot-nghost
        if (maux .gt. 0) then
          do 205 ma = 1,maux
             auxrbig(iaddaux(ma,i-nghost)) = auxbig(ma,i,mjtot-nghost+1)
 205         continue
          endif
        do 210 ivar = 1, nvar
            qr(ivar,i-nghost) = valbig(ivar,i,mjtot-nghost+1)
 210    continue

       lind = 0
       ncrse = (mitot-2*nghost)/lratiox
       do 220 ic = 1, ncrse
         index = index + 1
         do 225 l = 1, lratiox
         lind = lind + 1
         if (maux.gt.0) then
            do 224 ma=1,maux
             if (auxtype(ma).eq."yleft") then
!!                # Assuming velocity at bottom-face, this fix
!!                # preserves conservation in incompressible flow:
                 ifine = (ic-1)*lratiox + nghost + l
                 auxlbig(iaddaux(ma,lind+1)) = auxbig(ma,ifine,mjtot-nghost+1)
               else
                 auxlbig(iaddaux(ma,lind+1)) = auxc1d(ma,index)
               endif
  224          continue
            endif
         do 225 ivar = 1, nvar
 225         ql(ivar,lind+1) = qc1d(ivar,index)
 220    continue

       call rpn2(2,max1dp1-2*nghost,nvar,mwaves,maux,nghost, &
                   nr+1-2*nghost,ql,qr,auxlbig,auxrbig,wave,s,amdq,apdq)
!!
!! we have the wave. for side 2. add into sdflxp
!!
       do 230 i = 1, nr/lratiox
          influx  = influx + 1
          ifine = (i-1)*lratiox
          do 240 ivar = 1, nvar
            do 250 l = 1, lratiox
              svdflx(ivar,influx) = svdflx(ivar,influx) &
                          - amdq(ivar,ifine+l+1) * hx * delt &
                          - apdq(ivar,ifine+l+1) * hx * delt
 250         continue
 240       continue
 230    continue

 299  continue



!!  DO i = 1, mx
!!     IF (maux .GT. 0) THEN
!!        DO ma = 1,maux
!!           auxr(ma,i) = aux(ma,i,my+1)
!!        ENDDO
!!     ENDIF
!!     DO  ivar = 1, meqn
!!        qr(ivar,i) = q(ivar,i,my+1)
!!     ENDDO
!!  ENDDO

!!  lind = 0
!!  ncrse = mx/lratiox
!!  DO ic = 1, ncrse
!!     index = index + 1
!!     DO l = 1, lratiox
!!        lind = lind + 1
!!        IF (maux .GT. 0) THEN
!!           DO  ma = 1,maux
!!              IF (auxtype(ma) .EQ. "yleft") THEN
!!                 !! # Assuming velocity at bottom-face, this fix
!!                 !! # preserves conservation in incompressible flow:
!!                 ifine = (ic-1)*lratiox + l
!!                 auxl(ma,lind+1) = aux(ma,ifine,my+1)
!!              ELSE
!!                 auxl(ma,lind+1) = auxc1d(ma,index)
!!              ENDIF
!!           ENDDO
!!        ENDIF
!!        DO  ivar = 1, meqn
!!           ql(ivar,lind+1) = qc1d(ivar,index)
!!        ENDDO
!!     ENDDO
!!  ENDDO

!!  CALL rpn2(2,max1dp1-2*mbc,meqn,mwaves,maux,mbc, &
!!       mx+1-2*mbc,ql,qr,auxl,auxr,wave,s,amdq,apdq)
!!  !!
!!  !! we have the wave. for side 2. add into sdflxp
!!  !!
!!  DO i = 1, nr/lratiox
!!     influx  = influx + 1
!!     ifine = (i-1)*lratiox
!!     DO ivar = 1, nvar
!!        DO l = 1, lratiox
!!           svdflx(ivar,influx) = svdflx(ivar,influx) &
!!                - amdq(ivar,ifine+l+1) * hx * delt &
!!                - apdq(ivar,ifine+l+1) * hx * delt
!!        ENDDO
!!     ENDDO
!!  ENDDO



  !! --------
  !!  side 3
  !! --------
  !!
  DO j = nghost+1, mjtot-nghost
     IF (maux.GT.0) THEN
        DO ma = 1,maux
           auxrbig(iaddaux(ma,j-nghost)) = auxbig(ma,mitot-nghost+1,j)
        ENDDO
     ENDIF
     DO ivar = 1, nvar
        qr(ivar,j-nghost) = valbig(ivar,mitot-nghost+1,j)
     ENDDO
  ENDDO

  lind = 0
  ncrse = (mjtot-2*nghost)/lratioy
  DO jc = 1, ncrse
     index = index + 1
     DO l = 1, lratioy
        lind = lind + 1
        IF (maux .GT. 0) THEN
           DO ma = 1,maux
              IF (auxtype(ma).EQ."xleft") THEN
                 !! # Assuming velocity at left-face, this fix
                 !! # preserves conservation in incompressible flow:
                 jfine = (jc-1)*lratioy + nghost + l
                 auxlbig(iaddaux(ma,lind+1)) =  &
                      auxbig(ma,mitot-nghost+1,jfine)
              ELSE
                 auxlbig(iaddaux(ma,lind+1)) = auxc1d(ma,index)
              ENDIF
           ENDDO
        ENDIF
        DO ivar = 1, nvar
           ql(ivar,lind+1) = qc1d(ivar,index)
        ENDDO
     ENDDO
  ENDDO

  CALL rpn2(1,max1dp1-2*nghost,nvar,mwaves,maux,nghost, &
       nc+1-2*nghost,ql,qr,auxlbig,auxrbig,wave,s,amdq,apdq)
  !!
  !! we have the wave. for side 3 add into sdflxp
  !!
  DO j = 1, nc/lratioy
     influx = influx + 1
     jfine = (j-1)*lratioy
     DO ivar = 1, nvar
        DO l = 1, lratioy
           svdflx(ivar,influx) = svdflx(ivar,influx) &
                - amdq(ivar,jfine+l+1) * hy * delt &
                - apdq(ivar,jfine+l+1) * hy * delt
        ENDDO
     ENDDO
  ENDDO

  !! --------
  !!  side 4
  !! --------
  !!
  IF (mjtot .EQ. 2*nghost+1) THEN
     !! # a single row of interior cells only happens when using the
     !! # 2d amrclaw code to do a 1d problem with refinement.
     !! # (feature added in Version 4.3)
     !! # skip over sides 2 and 4 in this case
     go to 499
  ENDIF

  DO  i = nghost+1, mitot-nghost
     IF (maux.GT.0) THEN
        DO ma = 1,maux
           IF (auxtype(ma).EQ."yleft") THEN
              !! # Assuming velocity at bottom-face, this fix
              !! # preserves conservation in incompressible flow:
              auxlbig(iaddaux(ma,i-nghost+1)) = auxbig(ma,i,nghost+1)
           ELSE
              auxlbig(iaddaux(ma,i-nghost+1)) = auxbig(ma,i,nghost)
           ENDIF
        ENDDO
     ENDIF
     DO ivar = 1, nvar
        ql(ivar,i-nghost+1) = valbig(ivar,i,nghost)
     ENDDO
  ENDDO

  lind = 0
  ncrse = (mitot-2*nghost)/lratiox
  DO ic = 1, ncrse
     index = index + 1
     DO l = 1, lratiox
        lind = lind + 1
        IF (maux .GT. 0) THEN
           DO ma=1,maux
              auxrbig(iaddaux(ma,lind)) = auxc1d(ma,index)
           ENDDO
        ENDIF
        DO  ivar = 1, nvar
           qr(ivar,lind) = qc1d(ivar,index)
        ENDDO
     ENDDO
  ENDDO

  CALL rpn2(2,max1dp1-2*nghost,nvar,mwaves,maux,nghost, &
       nr+1-2*nghost,ql,qr,auxlbig,auxrbig,wave,s,amdq,apdq)
  !!
  !! we have the wave. for side 4. add into sdflxm
  !!
  DO i = 1, nr/lratiox
     influx  = influx + 1
     ifine = (i-1)*lratiox
     DO ivar = 1, nvar
        DO l = 1, lratiox
           svdflx(ivar,influx) = svdflx(ivar,influx) &
                + amdq(ivar,ifine+l+1) * hx * delt &
                + apdq(ivar,ifine+l+1) * hx * delt
        ENDDO
     ENDDO
  ENDDO

499 continue

  !!      # for source terms:
  IF (method(5) .NE. 0) THEN   ! should I test here if index=0 and all skipped?
     !!   call src1d(nvar,nghost,lenbc,qc1d,maux,auxc1d,tgrid,delt)
     !! # how can this be right - where is the integrated src term used?
  ENDIF

  RETURN
END SUBROUTINE qad
