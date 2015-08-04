ccc   gg --> gamma gamma subprocess amplitude
      subroutine gamgam(p,mu,u,t,pp,mm,pm,mp)
      implicit double precision (a-z)
      complex*16 pp,mm,pm,mp
      integer p

      include 'pi.f'
      include 'zi.f'
      include 'vars.f'
      include 'ewpars.f'
      include 'partonmom2.f'

      cphi=p1(1)/dsqrt(p1(1)**2+p1(2)**2)
      sphi=p1(2)/dsqrt(p1(1)**2+p1(2)**2)

      if(mu.gt.mb)then
         qf=11d0/9d0
      else
         qf=10d0/9d0
      endif

      alphaem=1d0/137d0
      norm=4d0*qf*alphaem*alphas(mu**2/4d0)
      sh=mx**2

      if(p.eq.1)then   ! ++
         mm=1d0
         pp=-0.5d0*(t**2+u**2)/sh**2*((dlog(t/u))**2+pi**2)
     &        -(t-u)/sh*dlog(t/u)-1d0
         mp=1d0
         pm=1d0
      elseif(p.eq.2)then  ! --
         mm=-0.5d0*(t**2+u**2)/sh**2*((dlog(t/u))**2+pi**2)
     &        -(t-u)/sh*dlog(t/u)-1d0
         pp=1d0
         mp=1d0
         pm=1d0
      elseif(p.eq.3)then  ! -+
         mm=1d0
         pp=1d0
         mp=-0.5d0*(t**2+sh**2)/u**2*((dlog(-t/sh))**2+
     &        2d0*zi*pi*dlog(-t/sh))-(t-sh)/u*(dlog(-t/sh)+zi*pi)-1d0
         pm=-0.5d0*(u**2+sh**2)/t**2*((dlog(-sh/u))**2+
     &        2d0*zi*pi*dlog(-sh/u))-(sh-u)/t*(dlog(-sh/u)+zi*pi)-1d0
      elseif(p.eq.4)then  ! +-
         mm=1d0
         pp=1d0
         mp=-0.5d0*(u**2+sh**2)/t**2*((dlog(-sh/u))**2+
     &        2d0*zi*pi*dlog(-sh/u))-(sh-u)/t*(dlog(-sh/u)+zi*pi)-1d0
         pm=-0.5d0*(t**2+sh**2)/u**2*((dlog(-t/sh))**2+
     &        2d0*zi*pi*dlog(-t/sh))-(t-sh)/u*(dlog(-t/sh)+zi*pi)-1d0
      endif
      
      pp=pp*norm
      mm=mm*norm
      pm=pm*norm
      mp=mp*norm

cccccccc

      pm=pm*(cphi-zi*sphi)**2
      mp=mp*(cphi+zi*sphi)**2
      
      return
      end
