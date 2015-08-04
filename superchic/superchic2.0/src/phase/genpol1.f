ccc   generates polarization vectors for chi_1
      subroutine genpol1(in,echi1)
      implicit double precision(a-y)
      double precision n1(4),n2(4)
      double precision nchi(4)
      integer i,in
      complex*16 echi1(3,4)

      include 'zi.f'
      include 'mom.f'

      m=dsqrt(q(4,in)**2-q(3,in)**2-q(2,in)**2-q(1,in)**2)

      pnorm=dsqrt(q(1,in)**2+q(2,in)**2+q(3,in)**2)
      do i=1,3
         nchi(i)=q(i,in)/pnorm
      enddo

      n1(1)=nchi(2)/dsqrt(nchi(1)**2+nchi(2)**2)
      n1(2)=-nchi(1)/dsqrt(nchi(1)**2+nchi(2)**2)
      n1(3)=0d0

      n2(1)=n1(2)*nchi(3)-n1(3)*nchi(2)
      n2(2)=n1(3)*nchi(1)-n1(1)*nchi(3)
      n2(3)=n1(1)*nchi(2)-n1(2)*nchi(1)

      echi1(1,4)=0d0
      echi1(2,4)=0d0

      do i=1,3
         echi1(1,i)=(n1(i)+zi*n2(i))/dsqrt(2d0)
         echi1(2,i)=-(n1(i)-zi*n2(i))/dsqrt(2d0)
      enddo

      echi1(3,4)=dsqrt(q(4,in)**2-m**2)/m
      do i=1,3
         echi1(3,i)=q(4,in)*nchi(i)/m
      enddo

      return
      end
