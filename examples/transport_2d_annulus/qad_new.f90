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

  INTEGER mitot, mjtot, nvar, lenbc, lratiox, lratioy, maux
  INTEGER mptr
  DOUBLE PRECISION hx,hy,delt

  DOUBLE PRECISION, TARGET :: valbig(nvar,mitot,mjtot)
  DOUBLE PRECISION, TARGET :: auxbig(maux,mitot,mjtot)
  DOUBLE PRECISION qc1d(nvar,lenbc)
  DOUBLE PRECISION svdflx(nvar,lenbc)
  DOUBLE PRECISION auxc1d(maux,lenbc)

  DOUBLE PRECISION, POINTER :: q(:,:,:)
  DOUBLE PRECISION, POINTER :: aux(:,:,:)

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

  DOUBLE PRECISION, TARGET :: auxlbig(maxaux*max1dp1)
  DOUBLE PRECISION, TARGET :: auxrbig(maxaux*max1dp1)

  DOUBLE PRECISION, POINTER :: auxl(:,:)
  DOUBLE PRECISION, POINTER :: auxr(:,:)

  !!
  !!  WARNING: auxl,auxr dimensioned at max possible, but used as if
  !!  they were dimensioned as the real maux by max1dp1. Would be better
  !!  of course to dimension by maux by max1dp1 but this wont work if maux=0
  !!  So need to access using your own indexing into auxl,auxr.

  !!  INTEGER iaddaux
  !!  iaddaux(iaux,i) = (i-1)*maux + iaux

  !!
  !!      aux is auxiliary array with user parameters needed in Riemann solvers
  !!          on fine grid corresponding to valbig
  !!      auxc1d is coarse grid stuff from around boundary, same format as qc1d
  !!      auxl, auxr are work arrays needed to pass stuff to rpn2
  !!      maux is the number of aux variables, which may be zero.
  !!

  DOUBLE PRECISION tgrid
  INTEGER nc, nr, level, index, l
  INTEGER i,j,ma,lind, ncrse, ic, jc, ifine, jfine
  INTEGER iaux

  INTEGER mx,my,mbc,meqn, mxc, myc, mq
  DOUBLE PRECISION dt, dx, dy, delta_fix
  integer iface, idir
  logical prt

!!  return 

  tgrid = rnode(timemult, mptr)
  nr = mitot-2*nghost
  nc = mjtot-2*nghost
  level = node(nestlevel, mptr)
  index = 0

  !! Rename variables to use Clawpack convention
  mbc = nghost
  mx = nr
  my = nc
  meqn = nvar
  dt = delt
  dx = hx
  dy = hy

  mxc = mx/lratiox
  myc = my/lratioy

  !! Redimension arrays to use indexing that starts at 1-mbc, etc
  q(1:meqn,1-mbc:mx+mbc,1-mbc:my+mbc) => valbig

  if (maux .gt. 0) then
     aux(1:maux,1-mbc:mx+mbc,1-mbc:my+mbc) => auxbig
     auxl(1:maux,1:max1dp1) => auxlbig
     auxr(1:maux,1:max1dp1) => auxrbig
  endif


  !! Counter for indexing into 1d arrays of coarse grid values 
  index = 0

  !! --------
  !!  side 1
  !! --------

  !! All looping is over fine grid
  DO j = 1,my
     IF (maux .GT. 0) THEN
        DO ma = 1,maux
           IF (auxtype(ma) .EQ. "xleft") THEN
              !! # Assuming velocity at left-face, this fix
              !! # preserves conservation in incompressible flow:
              auxl(ma,j) = aux(ma,1,j)
           ELSE
              !! # Normal case -- we set the aux arrays
              !! # from the cell corresponding  to q
              auxl(ma,j) = aux(ma,0,j)
           ENDIF
        ENDDO
     ENDIF
     DO mq = 1,meqn
        ql(mq,j) = q(mq,0,j)
     ENDDO
  ENDDO

  DO jc = 1,myc
     DO l = 1,lratioy
        jfine = (jc-1)*lratioy + l
        IF (maux .GT. 0) THEN
           DO ma = 1,maux
              auxr(ma,jfine) = auxc1d(ma,index + jc)
           ENDDO
        ENDIF
        DO mq = 1,meqn
           qr(mq,jfine) = qc1d(mq,index + jc)
        ENDDO
     ENDDO
  ENDDO



  !! side 1 (iface = 0)
  !! qf = ql
  !! qc = qr

  idir = 0
  iface = 0
  CALL rpn2_qad(my,meqn,maux,mbc, idir, iface, &
                ql,qr,auxl,auxr,amdq,apdq)


  !!
  !! we have the wave. for side 1 add into sdflxm
  !!

  DO jc = 1,myc
     DO l = 1,lratioy
        jfine = (jc-1)*lratioy + l
        DO  mq = 1,meqn
           delta_fix = amdq(mq,jfine) + apdq(mq,jfine)
           svdflx(mq,index+jc) = svdflx(mq,index+jc) + dy*dt*delta_fix
        ENDDO
     ENDDO
  ENDDO
  index  = myc

  !! --------
  !!  side 2
  !! --------
  !!
  IF (my .EQ. 1) THEN
     !! # a single row of interior cells only happens when using the
     !! # 2d amrclaw code to do a 1d problem with refinement.
     !! # (feature added in Version 4.3)
     !! # skip over sides 2 and 4 in this case
     go to 299
  ENDIF


  DO i = 1,mx
     IF (maux .GT. 0) THEN
        DO ma = 1,maux
           auxr(ma,i) = aux(ma,i,my+1)
        ENDDO
     ENDIF
     DO mq = 1,meqn
        qr(mq,i) = q(mq,i,my+1)
     ENDDO
  ENDDO

  DO ic = 1, mxc
     DO l = 1,lratiox
        ifine = (ic-1)*lratiox + l
        IF (maux .GT. 0) THEN
           DO  ma = 1,maux
              IF (auxtype(ma) .EQ. "yleft") THEN
                 auxl(ma,ifine) = aux(ma,ifine,my+1)
              ELSE
                 auxl(ma,ifine) = auxc1d(ma,index + ic)
              ENDIF
           ENDDO
        ENDIF
        DO  mq = 1,meqn
           ql(mq,ifine) = qc1d(mq,index + ic)
        ENDDO
     ENDDO
  ENDDO

  !! side 2 (iface = 4  (top edge))
  !! qf = qr
  !! qc = ql

  idir = 1
  iface = 3
  CALL rpn2_qad(mx,meqn,maux,mbc, idir, iface, &
                qr,ql,auxr,auxl,amdq,apdq)

  !! we have the wave. for side 2. add into sdflxp

  DO ic = 1,mxc
     DO l = 1,lratiox
        ifine = (ic-1)*lratiox + l
        DO mq = 1,meqn
           delta_fix = amdq(mq,ifine) + apdq(mq,ifine)
           svdflx(mq,index+ic) = svdflx(mq,index+ic) + dx*dt*delta_fix
        ENDDO
     ENDDO
  ENDDO
  index  = index + mxc

299 continue



  !! --------
  !!  side 3
  !! --------
  !!
  DO j = 1, my
     IF (maux .GT. 0) THEN
        DO ma = 1,maux
           auxr(ma,j) = aux(ma,mx+1,j)
        ENDDO
     ENDIF
     DO mq = 1, meqn
        qr(mq,j) = q(mq,mx+1,j)
     ENDDO
  ENDDO

  DO jc = 1,myc
     DO l = 1,lratioy
        jfine = (jc-1)*lratioy + l
        IF (maux .GT. 0) THEN
           DO ma = 1,maux
              IF (auxtype(ma).EQ."xleft") THEN
                 auxl(ma,jfine) = aux(ma,mx+1,jfine)
              ELSE
                 auxl(ma,jfine) = auxc1d(ma,index + jc)
              ENDIF
           ENDDO
        ENDIF
        DO mq = 1, meqn
           ql(mq,jfine) = qc1d(mq,index + jc)
        ENDDO
     ENDDO
  ENDDO


  !! side 3 (iface = 1  (right edge))
  !! qf = qr
  !! qc = ql

  idir = 0
  iface = 1
  CALL rpn2_qad(my,meqn,maux,mbc, idir, iface, &
                qr,ql,auxr,auxl,amdq,apdq)


  !!
  !! we have the wave. for side 3 add into sdflxp
  !!
  DO jc = 1, myc
     DO l = 1, lratioy
        jfine = (jc-1)*lratioy + l
        DO mq = 1, meqn
           delta_fix = amdq(mq,jfine) + apdq(mq,jfine)
           svdflx(mq,index + jc) = svdflx(mq,index + jc) - dy*dt*delta_fix
        ENDDO
     ENDDO
  ENDDO
  index  = index + myc

  !! --------
  !!  side 4
  !! --------
  !!
  IF (my .EQ. 1) THEN
     !! # a single row of interior cells only happens when using the
     !! # 2d amrclaw code to do a 1d problem with refinement.
     !! # (feature added in Version 4.3)
     !! # skip over sides 2 and 4 in this case
     go to 499
  ENDIF

  DO  i = 1, mx
     IF (maux .GT. 0) THEN
        !! Is this conditional needed?  Loop won't do anything if maux == 0
        DO ma = 1,maux
           IF (auxtype(ma) .EQ. "yleft") THEN
              !! But this left edge may not coincide with the coarse grid left edge              
              auxl(ma,i) = aux(ma,i,1)
           ELSE
              auxl(ma,i) = aux(ma,i,0)
           ENDIF
        ENDDO
     ENDIF
     DO mq = 1, meqn
        ql(mq,i) = q(mq,i,0)
     ENDDO
  ENDDO

  DO ic = 1,mxc
     DO l = 1,lratiox
        ifine = (ic-1)*lratiox + l
        IF (maux .GT. 0) THEN
           DO ma = 1,maux
              auxr(ma,ifine) = auxc1d(ma,index + ic)
           ENDDO
        ENDIF
        DO  mq = 1, meqn
           qr(mq,ifine) = qc1d(mq,index + ic)
        ENDDO
     ENDDO
  ENDDO

  !! side 4 (iface = 2  (right edge))
  !! qf = ql
  !! qc = qr  

  idir = 1
  iface = 2
  CALL rpn2_qad(mx,meqn,maux,mbc, idir, iface, &                
                ql,qr,auxl,auxr,amdq,apdq)


  !!
  !! we have the wave. for side 4. add into sdflxm
  !!
  DO ic = 1,mxc
     DO l = 1,lratiox
        ifine = (ic-1)*lratiox + l
        DO mq = 1,meqn
           delta_fix = amdq(mq,ifine) + apdq(mq,ifine)
           svdflx(mq,index + ic) = svdflx(mq,index + ic) + dx*dt*delta_fix
        ENDDO
     ENDDO
  ENDDO

499 continue

  !!      # for source terms:
  IF (method(5) .NE. 0) THEN   ! should I test here if index=0 and all skipped?
     !!   call src1d(meqn,nghost,lenbc,qc1d,maux,auxc1d,tgrid,delt)
     !! # how can this be right - where is the integrated src term used?
  ENDIF

  RETURN
END SUBROUTINE qad
