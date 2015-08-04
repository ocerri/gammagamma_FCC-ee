ccc   spin-summed gamma gamma --> l+l- cross section
      subroutine ll(mx,u,t,out)
      implicit double precision(a-y)

      include 'mq.f'
      include 'ewpars.f'
      include 'pi.f'
      include 'mom.f'
      include 'norm.f'

      beta=dsqrt(1d0-4d0*mq**2/mx**2)
      gam=mx**2/mq**2
      costt=(t-u)/beta/mx**2

      out=(1d0+2d0*beta**2*(1d0-beta**2)*(1d0-costt**2)
     &     -beta**4*costt**4)/(1d0-beta**2*costt**2)**2
      out=out*2d0*pi*alpha**2*beta/mx**2

      out=out*2d0
      out=out*conv

      return
      end
