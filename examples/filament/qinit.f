c     =====================================================
       subroutine qinit(meqn,mbc,mx,my,
     &      xlower,ylower,dx,dy,q,maux,aux)
c     =====================================================

c     # Set initial conditions for q.
c     # Sample scalar equation with data that is piecewise constant with
c     # q = 1.0  if  0.1 < x < 0.6   and   0.1 < y < 0.6
c     #     0.1  otherwise

       implicit none

       integer meqn, mbc, mx, my, maux
       double precision xlower, ylower, dx, dy
       double precision q(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
       double precision aux(maux,1-mbc:mx+mbc, 1-mbc:my+mbc)

       integer i, j, mq
       double precision xlow, ylow, w

       do mq = 1,meqn
          do i = 1-mbc,mx+mbc
             do j = 1-mbc,my+mbc
                xlow = xlower + (i-1)*dx
                ylow = ylower + (j-1)*dy
                call cellave(xlow,ylow,dx,dy,w)
                q(mq,i,j) = w
             enddo
          enddo
       enddo

       end

      double precision function  fdisc(xc,yc)
      implicit none

      double precision xc,yc, xp, yp, zp

      double precision r

      r = sqrt((xc-0.5d0)**2 + (yc-1.d0)**2)

      fdisc = r-0.25d0
      end
