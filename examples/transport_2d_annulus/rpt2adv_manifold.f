c      subroutine rpt2adv(ixy,maxm,meqn,
c     &      mwaves,mbc,mx,ql,qr,aux1,aux2,aux3,
c     &      imp,asdq,bmasdq,bpasdq)
      subroutine rpt2adv_manifold(ixy,imp,maxm,meqn,mwaves,maux,
     &               mbc,mx,ql,qr,
     &               aux1,aux2,aux3,asdq,bmasdq,bpasdq)
      implicit none

      integer ixy,maxm, maux, meqn, mwaves, mbc, mx, imp, m
      double precision     ql(meqn,1-mbc:maxm+mbc)
      double precision     qr(meqn,1-mbc:maxm+mbc)
      double precision   asdq(meqn,1-mbc:maxm+mbc)
      double precision bmasdq(meqn,1-mbc:maxm+mbc)
      double precision bpasdq(meqn,1-mbc:maxm+mbc)
      double precision   aux1(maux,1-mbc:maxm+mbc)
      double precision   aux2(maux,1-mbc:maxm+mbc)
      double precision   aux3(maux,1-mbc:maxm+mbc)

      integer iface, i, i1

c     # ixy = 1 : Normal solve is in x direction; Transverse is in y dir.   
c     #           so we need y face values for transverse solve   
c     # ixy = 2 : Normal solve is in y direction; Transverse is in x dir.      
c     #           so we need x face values for transverse solve   

c     # (1,2), (2,0) --> iface = (2/-1)*(ixy-2) = 2*(2-ixy)

c     # iface = 2*(2-ixy)
      if (ixy .eq. 1) then
         iface = 2  !! Bottom edge
      else
         iface = 0  !! left edge
      endif

      do i = 2-mbc,mx+mbc
         i1 = i-2+imp    !#  =  i-1 for amdq,  i for apdq
         do m = 1,meqn
            bmasdq(m,i) = min(aux2(1+iface,i1), 0.d0) * asdq(m,i)
            bpasdq(m,i) = max(aux3(1+iface,i1), 0.d0) * asdq(m,i)
         enddo
      enddo

      return
      end
