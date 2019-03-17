      double precision function  fdisc(blockno,xc,yc)
      implicit none

      double precision xc,yc, xp, yp, zp
      integer blockno
      double precision th, tp

      double precision pi, pi2
      common /compi/ pi, pi2

      double precision beta
      common /annulus_comm/ beta

      double precision init_radius
      common /initradius_comm/ init_radius

      double precision x0, y0, r0, r, theta, ravg

      call mapc2m_annulus(xc,yc,xp,yp,zp)

c     # Torus or annulus
c      th = atan2(yp,xp)
c      tp = abs(th - pi/2.d0)
c      fdisc = tp - pi/8.d0

      r0 = init_radius
      ravg = (1 + beta)/2.d0
      x0 = -ravg
      y0 = 0

c      theta = pi2*xc
c      r = beta + (1-beta)*yc
c      if (abs(theta-pi) .lt. pi/12.d0 .and. abs(r-ravg) .lt. r0) then
c            fdisc = -1.d0
c      else
c            fdisc = 1
c      endif

      r = sqrt((xp - x0)**2 + (yp-y0)**2)
      fdisc = r - r0

      end
