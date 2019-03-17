      subroutine mapc2m_annulus(xc,yc,xp,yp,zp)
      implicit none

      double precision xc,yc,xp,yp,zp, r

      double precision pi, pi2
      common /compi/ pi, pi2

      double precision beta
      common /annulus_comm/ beta


      r = beta + (1-beta)*yc
      xp = r*cos(pi2*xc)
      yp = r*sin(pi2*xc)
      zp = 0


      end
