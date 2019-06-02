      subroutine annulus46_rpt2adv_manifold(ixy,maxm,meqn,
     &      mwaves,mbc,mx,ql,qr,aux1,aux2,aux3,
     &      imp,asdq,bmasdq,bpasdq,maux)
      implicit none

      integer ixy,maxm, meqn, mwaves, mbc, mx, imp, m, maux
      double precision     ql(meqn,1-mbc:maxm+mbc)
      double precision     qr(meqn,1-mbc:maxm+mbc)
      double precision   asdq(meqn,1-mbc:maxm+mbc)
      double precision bmasdq(meqn,1-mbc:maxm+mbc)
      double precision bpasdq(meqn,1-mbc:maxm+mbc)
      double precision   aux1(maux,1-mbc:maxm+mbc)
      double precision   aux2(maux,1-mbc:maxm+mbc)
      double precision   aux3(maux,1-mbc:maxm+mbc)

      integer i, i1, idir, ii, jj, idx
      double precision vm, vp, gm,gp
      double precision face_data(-1:1,-1:1,0:3)

c     # Idir is direction of transverse solve
      idir = 2-ixy  
      do i = 2-mbc,mx+mbc
         i1 = i-2+imp    !#  =  i-1 for amdq,  i for apdq


         do m = 0,3
             if (ixy .eq. 1) then
                 do ii = -1,1
                     face_data(ii,-1,m) = aux1(i+ii,2+m)
                     face_data(ii, 0,m) = aux2(i+ii,2+m) 
                     face_data(ii, 1,m) = aux3(i+ii,2+m)
                 enddo              
             else
                 do jj = -1,1
                     face_data(-1,jj,m)  = aux1(i+jj,2+m)
                     face_data( 0,jj,m)  = aux2(i+jj,2+m) 
                     face_data( 1,jj,m)  = aux3(i+jj,2+m)
                 enddo
             endif
         end do



         idx = imp - 2
         if (ixy .eq. 1) then
            gm = aux2(8,i1)   
            gp = aux3(8,i1)
            vm = gm*(face_data(idx,-1,3) + face_data(idx,0,2))/2.d0
            vp = gp*(face_data(idx, 0,3) + face_data(idx,1,2))/2.d0
         else
            gm = aux2(6,i1)
            gp = aux3(6,i1)
            vm = gm*(face_data(-1,idx,1) + face_data(0,idx,0))/2.d0
            vp = gp*(face_data( 0,idx,1) + face_data(1,idx,0))/2.d0
         endif


c        # Velocities are stored in [2,3,4,5] (corresponding to faces [0,1,2,3])         
c         if (ixy .eq. 1) then
c            vm = gm*(aux1(i1,4) + aux2(i1,3))/2.d0
c            vp = gp*(aux2(i1,4) + aux3(i1,3))/2.d0
c         else
c            vm = gm*(aux1(i1,2) + aux2(i1,1))/2.d0
c            vp = gp*(aux2(i1,2) + aux3(i1,1))/2.d0
c         endif

          if (.false.) then
              write(6,200) 'imp = ', imp
              m = 2
              do jj = 1,-1,-1
                  write(6,100) face_data(-1,jj,m), face_data(-1,jj,m), 
     &               face_data(-1,jj,m)              
              enddo
              write(6,*) ' '
          endif

100       format(3F16.8)
200       format(A,I5)

         do m = 1,meqn
            bmasdq(m,i) = min(vm, 0.d0) * asdq(m,i)
            bpasdq(m,i) = max(vp, 0.d0) * asdq(m,i)
         enddo
      enddo

      return
      end
