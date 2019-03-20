      subroutine mapc2m_annulus(ac,bc,xp,yp,zp)
      implicit none

      double precision ac,bc, xc,yc,xp,yp,zp, r
      double precision l1(4)

      double precision pi, pi2
      common /compi/ pi, pi2

      integer mapping
      common /mapping_comm/ mapping

      double precision beta
      common /annulus_comm/ beta


      call annulus_transform_coordinates(ac,bc,xc,yc,mapping)

      r = beta + (1-beta)*yc
      xp = r*cos(pi2*xc)
      yp = r*sin(pi2*xc)
      zp = 0


      end
