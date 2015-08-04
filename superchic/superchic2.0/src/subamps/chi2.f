ccc   gg --> chi_0 subprocess amplitude
      subroutine chi2(p,mqq,mxx,q1,q2,echi2,out)
      implicit double precision (a-z)
      double precision q1(2),q2(2)
      double precision qt1(4),qt2(4)
      complex*16 out,echi2(5,4,4),cpp
      integer p,i

      include 'pi.f'
      include 'zi.f'
      include 'vars.f'
      include 'mom.f'
      include 'quarkonia.f'

      qt1(4)=0d0
      qt2(4)=0d0
      qt1(3)=0d0
      qt2(3)=0d0

      do i=1,2
         qt1(i)=q1(i)
         qt2(i)=q2(i)
      enddo

      qt1sq=qt1(1)**2+qt1(2)**2
      qt2sq=qt2(1)**2+qt2(2)**2
      q1q2=(mx**2+qt1sq+qt2sq)/2d0

      cchi=dsqrt(pi*mx**3*gamchi0/3d0)

      cpp=s*(qt1(1)*qt2(1)*echi2(p,1,1)+qt1(1)*qt2(2)*echi2(p,1,2)
     &+qt1(2)*qt2(1)*echi2(p,2,1)+qt1(2)*qt2(2)*echi2(p,2,2))
      cpp=cpp-2d0*(qt1(1)*qt2(1)+qt1(2)*qt2(2))*
     &(q(3,1)*q(3,2)*echi2(p,3,3)+q(4,1)*q(4,2)*echi2(p,4,4)
     &-q(3,1)*q(4,2)*echi2(p,3,4)-q(4,1)*q(3,2)*echi2(p,4,3))
 
      cpp=cpp*cchi*dsqrt(2d0)*mx/s
      cpp=cpp/(2d0*mqq*mx+qt1sq+qt2sq)**2*4d0
      cpp=cpp*dsqrt(mx/mxx)
      cpp=cpp*mx**2/2d0

      out=cpp

      return
      end
