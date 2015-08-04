ccc   calculates normalization factor when using RAMBO 
ccc   for three-body phase space 
      subroutine threebodyinit(ein,m1,m2,m3)
      implicit double precision(a-y)
      integer nrun,npart,n
      double precision am(100),pout(4,100)

      include 'wtinit.f'
  
      nrun=100000

      npart=3
      am(1)=m1
      am(2)=m2
      am(3)=m3

      sum=0d0

      do n=1,nrun
         call rambo(npart,ein,am,pout,wt)
         sum=sum+wt
      enddo

      wt3i=sum/dfloat(nrun) 

      return
      end
