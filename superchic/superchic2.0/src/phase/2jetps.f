ccc   generates two-body phase space for particles of equal mass mq
      subroutine twojetps(mx,mq,ps,u,t)
      implicit double precision(a-y)
      double precision px(4),pcm(4)
      double precision pboo(4)
      integer i

      include 'mt.f'
      include 'partonmom2.f'
      include 'pi.f'
      include 'mom.f'

      call r2455(rphi)
      call r2455(rtheta)

      phi=2d0*pi*rphi
      ctheta=-1d0+2d0*rtheta
      stheta=dsqrt(1d0-ctheta**2)
      beta=dsqrt(1d0-4d0*mq**2/mx**2)

      pcm(4)=mx/2d0
      pcm(1)=stheta*dcos(phi)*pcm(4)*beta
      pcm(2)=stheta*dsin(phi)*pcm(4)*beta
      pcm(3)=ctheta*pcm(4)*beta
      
      u=mq**2-mx*(pcm(4)+pcm(3))
      t=mq**2-mx*(pcm(4)-pcm(3))

      do i=1,4
         px(i)=q(i,5)
         p1(i)=pcm(i)
      enddo

      p2(4)=pcm(4)
      do i=1,3
         p2(i)=-pcm(i)
      enddo

      call boost(mx,px,pcm,pboo)

      do i=1,4
         q(i,6)=pboo(i)
         q(i,7)=q(i,5)-q(i,6)
      enddo

      ps=4d0*pi*beta

      mt=dsqrt(mq**2+p1(1)**2+p1(2)**2)

      return
      end
