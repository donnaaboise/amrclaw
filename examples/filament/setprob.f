      subroutine setprob
      implicit double precision (a-h,o-z)
      character(len=25) fname

      double precision figsize(2)

      common /compi/ pi
      common /comtikz/ figsize

c
c     # compute pi, used in psi.f
      pi = 4.d0 * datan(1.d0)
c
c     # save 2*pi and tperiod in common block for use in b4step2:
c
      pi2 = 2.d0*pi
c

      iunit = 7
      fname = 'setprob.data'
      call opendatafile(iunit, fname)
      read(7,*) figsize(1)
      read(7,*) figsize(2)
      close(7)

      return
      end
