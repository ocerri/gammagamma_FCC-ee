ccc   spin-summed gamma gamma --> W+W- cross section
      subroutine ww(mx,u,t,out)
      implicit double precision(a-y)

      include 'ewpars.f'
      include 'pi.f'
      include 'scorr.f'
      include 'mom.f'
      include 'norm.f'

      beta=dsqrt(1d0-4d0*mw**2/mx**2)
      gam=mx**2/mw**2
      costt=(t-u)/beta/mx**2

      if(scorr)then
         out=(4d0*pi*alpha/(1d0-beta**2*costt**2))**2
         out=out/4d0*9d0
         out=out*beta/(32d0*pi*mx**2)
      else
         out=(1d0-8d0*(2d0/3d0+1d0/gam)/(1d0-beta**2*costt**2)
     &        +32d0*(1d0/3d0+1d0/gam**2)/(1d0-beta**2*costt**2)**2)
         out=out*3d0*pi*alpha**2*beta/mx**2
      endif

      out=out*2d0
      out=out*conv

      return
      end
