ccccc Photoproduction form factors (assume Pomeron coupling universal
ccccc between eigenstates
      subroutine formfacphot(io,t1,t2,out)
      implicit double precision(a-y)
      integer i1,i2,io

      include 'bpsi.f'
      include 'photo.f'
      include 'mp.f'
      include 'pi.f'
      
      if(prot.eq.1)then
         wt=dexp(-bpsi*t2)
         qsq=(xgam**2*mp**2+t1)/(1d0-xgam)
      else
         wt=dexp(-bpsi*t1) 
         qsq=(xgam**2*mp**2+t2)/(1d0-xgam)
      endif

      qsqmin=mp**2*xgam**2/(1d0-xgam)
      ge=1d0/(1d0+qsq/0.71d0)**4
      mum=7.78d0
      gm=ge*mum
      fe=(4d0*mp**2*ge+qsq*gm)/(4d0*mp**2+qsq)
      fm=gm

      ww1=fe/qsq**2
      ww2=xgam**2/2d0*fm/qsq

      ww1=ww1/pi/(1d0-xgam)/137d0*2d0
      out1=dsqrt(wt*ww1)

      ww2=ww2/pi/(1d0-xgam)/137d0*2d0
      out2=dsqrt(wt*ww2)

      if(io.eq.1)then
         out=out1
      else
         out=out2
      endif

      return
      end
