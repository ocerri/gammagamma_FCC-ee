      subroutine lightlightpol(p,mu,u,t,pp,mm,pm,mp)
      implicit double precision (a-z)
      complex*16 pp,mm,pm,mp
      integer p

      include 'pi.f'
      include 'zi.f'
      include 'vars.f'
      include 'ewpars.f'
      include 'norm.f'

ccccc   include phase for pm/mp????

      if(mu.gt.mb)then
         qf=35d0/81d0
      else
         qf=34d0/81d0
      endif

      qf=qf*3d0
      ql=3d0
  
      norm=(8d0*(qf+ql)*alpha**2)**2
      norm=norm/16d0/pi/mx**2
      norm=norm/4d0
      norm=norm*conv
      
      norm=dsqrt(norm)

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
      
      return
      end
