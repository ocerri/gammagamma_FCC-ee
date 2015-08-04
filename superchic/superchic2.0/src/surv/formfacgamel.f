ccccc EPA form factors (electron)
      subroutine formfacgamel(io,t1,t2,out)
      implicit double precision(a-y)
      integer i1,i2,io

      include 'photo.f'
      include 'ewpars.f'
      include 'pi.f'
      include 'x.f'
      
      qsq=(x1**2*me**2+t1)/(1d0-x1)
      qsqp=(x2**2*me**2+t2)/(1d0-x2)      

      fe=1d0
      fm=1d0

cccccccccc

      qsqmin=me**2*x1**2/(1d0-x1)
      ww1p=((1d0-x1)*(1d0-qsqmin/qsq)*fe)/qsq
      ww1pa=x1**2/2d0*fm/qsq

      ww1=((1d0-x1)*fe)/qsq**2/(1d0-x1)
      ww1=ww1/pi/(1d0-x1)/137d0

      ww1p=ww1p/pi/(1d0-x1)/137d0
      ww1pa=ww1pa/pi/(1d0-x1)/137d0


ccccccccccc

      qsqmin=me**2*x2**2/(1d0-x2)

      ww2p=((1d0-x2)*(1d0-qsqmin/qsqp)*fe)
     &/qsqp
      ww2pa=x2**2/2d0*fm/qsqp

      ww2=((1d0-x2)*fe)/qsqp**2/(1d0-x2)
      ww2=ww2/pi/(1d0-x2)/137d0

      ww2p=ww2p/pi/(1d0-x2)/137d0
      ww2pa=ww2pa/pi/(1d0-x2)/137d0

      if(io.eq.1)then
         out=dsqrt(ww1*ww2)
      elseif(io.eq.2)then
         out=dsqrt(ww1p*ww2pa+ww2p*ww1pa+ww1pa*ww2pa)
      endif

      return
      end
