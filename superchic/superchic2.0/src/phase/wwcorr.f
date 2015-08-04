ccc   spin correlations in W+W- production (leptonic decays)
      subroutine wwcorr(mx,u,t,out)
      implicit double precision(a-y)
      double precision rhoW1(3),rhoW2(3),rhoggp(3,3)
      integer mm,jj,nn

      include 'ewpars.f'
      include 'pi.f'
      include 'partonmom2.f'
      include 'partonmom4.f'

      cost1=(p1(1)*paa(1)+p1(2)*paa(2)+p1(3)*paa(3))
     &/dsqrt((paa(1)**2+paa(2)**2+paa(3)**2)*
     &(p1(1)**2+p1(2)**2+p1(3)**2))
      sint1=dsqrt(1d0-cost1**2)

      cost2=(p2(1)*pbb(1)+p2(2)*pbb(2)+p2(3)*pbb(3))
     &/dsqrt((pbb(1)**2+pbb(2)**2+pbb(3)**2)*
     &(p2(1)**2+p2(2)**2+p2(3)**2))
      sint2=dsqrt(1d0-cost2**2)

      beta=dsqrt(1d0-4d0*mw**2/mx**2)
      gam=mx**2/mw**2
      costt=(t-u)/beta/mx**2
      sintt=dsqrt(1d0-costt**2)

      rhoggp(1,1)=8d0*(16d0+(4d0+gam)**2*sintt**4)/gam**2  
      do mm=2,3
         lam=2d0*(-1.5d0+dble(mm-1))
         rhoggp(mm,1)=64d0*sintt**2*(1d0+lam**2*costt**2)/gam
         rhoggp(1,mm)=64d0*sintt**2*(1d0+lam**2*costt**2)/gam
      enddo
      do mm=2,3
         do jj=2,3
            lam3=2d0*(-1.5d0+dble(mm-1))
            lam4=2d0*(-1.5d0+dble(jj-1))
            rhoggp(mm,jj)=8d0*beta**2*(lam3+lam4)**2+
     &           (-8d0*(1d0+lam3*lam4)+4d0*gam*(1d0+lam3*lam4))**2
     &           /(2d0*gam**2)+8d0*(lam3-lam4)**2*costt**2+
     &           (8d0*(1d0+lam3*lam4)+2d0*gam*(1d0-lam3*lam4)
     &           -8d0*(1d0+lam3*lam4)*costt**2
     &           +2d0*gam*(1d0-lam3*lam4)*costt**2)**2/(2d0*gam**2)
         enddo
      enddo
      
      rhoW1(1)=-sint1/dsqrt(2d0)
      rhoW1(2)=(1d0+cost1)/2d0
      rhoW1(3)=(1d0-cost1)/2d0
      rhoW2(1)=-sint2/dsqrt(2d0)
      rhoW2(2)=(1d0+cost2)/2d0
      rhoW2(3)=(1d0-cost2)/2d0
      
      wt1=0d0

      do mm=1,3
         do nn=1,3     
            wt1=wt1+rhoggp(mm,nn)*rhoW1(mm)**2*rhoW2(nn)**2
         enddo
      enddo

      out=wt1

      return
      end
