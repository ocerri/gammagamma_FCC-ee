ccc   generates two-body phase space for particles of mass m1,m2
      subroutine twobody(j,in,i1,i2,m1,m2,wt)
      implicit double precision(a-y)
      integer in,i1,i2,npart
      integer i,j
      double precision pcm(4),pboo(4),px(4)
      double precision am(100),pout(4,100)
  
      include 'mom.f'
      include 'wtinit.f'
      include 'partonmom4.f'

      npart=2
      ein=dsqrt(q(4,in)**2-q(3,in)**2-q(2,in)**2-q(1,in)**2)
      am(1)=m1
      am(2)=m2

      call rambo(npart,ein,am,pout,wt)

      do i=1,4
         px(i)=q(i,in)
         pcm(i)=pout(i,1)
      enddo 

      if(in.eq.5.or.in.eq.6)then
         do i=1,4
            paa(i)=pcm(i)
         enddo
      elseif(in.eq.7)then
         do i=1,4
            pbb(i)=pcm(i)
         enddo
      endif

      call boost(ein,px,pcm,pboo)
      
      do i=1,4
         q(i,i1)=pboo(i)
         q(i,i2)=q(i,in)-q(i,i1)
      enddo

      wt=wt/wt2i(j)

      return
      end
