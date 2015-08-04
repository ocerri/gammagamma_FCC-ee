ccc   gg --> SM higgs subprocess amplitude
      subroutine higgs(mx,pp,mm,pm,mp)
      implicit double precision (a-z)
      complex*16 pp,mm,pm,mp

      include 'pi.f'
      include 'ewpars.f'

      xh=(mt/mx)**2
      fh=-2d0*(dasin(1d0/(2d0*dsqrt(xh))))**2
      ih=3d0*xh*(2d0+(4d0*xh-1d0)*fh)

      lambdacaph=0.25d0
      nfh=5d0

      b1h=(33d0-2d0*nfh)/(12d0*pi)
      b11h=(153d0-19d0*nfh)/(2d0*pi*(33d0-2d0*nfh))

      hmpp=ih*mx**2*alphas(mx**2)/4d0/pi/v*2d0/3d0
      kf=(1d0+alphas(mx**2)/pi*(pi**2+11d0/2d0))

cccccccccccc

      pp=hmpp*dsqrt(kf)
      mm=pp
      mp=0d0
      pm=0d0
      
      return
      end







