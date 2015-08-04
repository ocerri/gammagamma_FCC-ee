ccc   gg --> eta'eta' subprocess amplitude
      subroutine etapetap(p,mx,u,t,pp,mm,pm,mp)
      implicit double precision (a-z)
      complex*16 pp,mm,pm,mp
      integer p
 
      include 'mixing.f'
      include 'partonmom2.f'
      include 'zi.f'

cccccc

      phi=datan(p1(2)/p1(1))
      cphi=p1(1)/dsqrt(p1(1)**2+p1(2)**2)
      sphi=p1(2)/dsqrt(p1(1)**2+p1(2)**2)

ccccccc

      call eta(p,mx,u,t,ppqqs,ppqg,ppgg,pmqqs,pmqqns)

      pp=dcos(thetap1)**2*(ppqqs+ppqg+ppgg)*fpi1**2
      mm=pp
      pm=fpi1**2*dcos(thetap1)**2*pmqqs
     &+(fpi8**2*dsin(thetap8)**2+fpi1**2*dcos(thetap1)**2)*pmqqns   
      mp=pm

      pm=pm*(cphi+zi*sphi)**2
      mp=mp*(cphi-zi*sphi)**2

      return
      end
