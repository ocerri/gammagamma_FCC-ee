      subroutine cut(icut)
      implicit double precision(a-y)
      double precision p1(4),p2(4),p3(4)
      integer icut,jflag,i
      logical accut

      include 'gencuts.f'
      include 'vars.f'
      include 'mom.f'
      include 'proc.f'
      include 'jetalg.f'
      include 'decay.f'
      include 'pi.f'

      icut=0

ccccccccccccccccccccccccccccccccccccccccccccc
ccc
cc    Place user-defined cuts here if needed
cc
cc    if(..)return - return for failed cut
cc
cccccccccccccccccccccccccccccccccccccccccccccc

      pt1=dsqrt(q(1,3)**2+q(2,3)**2)
      pt2=dsqrt(q(1,4)**2+q(2,4)**2)

      ptx=dsqrt(q(1,5)**2+q(2,5)**2)

c      if(ptx.gt.1.5d0)return

c      print*,q1sq

c      if(pt2**2.gt.0.4d0)return

c      if(pt1.gt.100d-3)return
c      if(pt2.gt.100d-3)return

      if(decay4)then

         if(proc.eq.41)then
            
            et1=dsqrt(q(1,7)**2+q(2,7)**2)
            et2=dsqrt(q(1,8)**2+q(2,8)**2)
            et3=dsqrt(q(1,9)**2+q(2,9)**2)
            et4=dsqrt(q(1,10)**2+q(2,10)**2)
            
            pmod1=dsqrt(q(1,7)**2+q(2,7)**2+q(3,7)**2)
            pmod2=dsqrt(q(1,8)**2+q(2,8)**2+q(3,8)**2)
            pmod3=dsqrt(q(1,9)**2+q(2,9)**2+q(3,9)**2)
            pmod4=dsqrt(q(1,10)**2+q(2,10)**2+q(3,10)**2)
            
            eta1=0.5d0*dlog((pmod1+q(3,7))/(pmod1-q(3,7)))
            eta2=0.5d0*dlog((pmod2+q(3,8))/(pmod2-q(3,8)))
            eta3=0.5d0*dlog((pmod3+q(3,9))/(pmod3-q(3,9)))
            eta4=0.5d0*dlog((pmod4+q(3,10))/(pmod4-q(3,10)))

            elseif(proc.eq.17.or.proc.eq.18.or.proc.eq.19
     &           .or.proc.eq.42.or.proc.eq.43)then

            et1=dsqrt(q(1,8)**2+q(2,8)**2)
            et2=dsqrt(q(1,9)**2+q(2,9)**2)
            et3=dsqrt(q(1,10)**2+q(2,10)**2)
            et4=dsqrt(q(1,11)**2+q(2,11)**2)
            
            pmod1=dsqrt(q(1,8)**2+q(2,8)**2+q(3,8)**2)
            pmod2=dsqrt(q(1,9)**2+q(2,9)**2+q(3,9)**2)
            pmod3=dsqrt(q(1,10)**2+q(2,10)**2+q(3,10)**2)
            pmod4=dsqrt(q(1,11)**2+q(2,11)**2+q(3,11)**2)
            
            eta1=0.5d0*dlog((pmod1+q(3,8))/(pmod1-q(3,8)))
            eta2=0.5d0*dlog((pmod2+q(3,9))/(pmod2-q(3,9)))
            eta3=0.5d0*dlog((pmod3+q(3,10))/(pmod3-q(3,10)))
            eta4=0.5d0*dlog((pmod4+q(3,11))/(pmod4-q(3,11)))

            
            endif

      if(et1.lt.ptamin4)return
      if(et2.lt.ptbmin4)return
      if(et3.lt.ptcmin4)return
      if(et4.lt.ptdmin4)return
      if(eta1.gt.etaamax4)return
      if(eta2.gt.etabmax4)return
      if(eta3.gt.etacmax4)return
      if(eta4.gt.etadmax4)return
      if(eta1.lt.etaamin4)return
      if(eta2.lt.etabmin4)return
      if(eta3.lt.etacmin4)return
      if(eta4.lt.etadmin4)return

      elseif(dps.eq.2.or.dps.eq.12.or.decay2)then


      pmod6=dsqrt(q(1,6)**2+q(2,6)**2+q(3,6)**2)
      pmod7=dsqrt(q(1,7)**2+q(2,7)**2+q(3,7)**2)
      eta6=0.5d0*dlog((pmod6+q(3,6))/(pmod6-q(3,6)))
      eta7=0.5d0*dlog((pmod7+q(3,7))/(pmod7-q(3,7)))
      y6=0.5d0*dlog((q(4,6)+q(3,6))/(q(4,6)-q(3,6)))
      y7=0.5d0*dlog((q(4,7)+q(3,7))/(q(4,7)-q(3,7)))
      pt6=dsqrt(q(1,6)**2+q(2,6)**2)
      pt7=dsqrt(q(1,7)**2+q(2,7)**2)


      if(pt6.lt.ptamin)return
      if(pt7.lt.ptbmin)return
      if(eta6.lt.etaamin)return
      if(eta7.lt.etabmin)return
      if(eta6.gt.etaamax)return
      if(eta7.gt.etabmax)return


      

      pmod8=dsqrt(q(1,8)**2+q(2,8)**2+q(3,8)**2)
      pmod9=dsqrt(q(1,9)**2+q(2,9)**2+q(3,9)**2)
      pmod10=dsqrt(q(1,10)**2+q(2,10)**2+q(3,10)**2)
      pmod11=dsqrt(q(1,11)**2+q(2,11)**2+q(3,11)**2)

c      y8=0.5d0*dlog((q(4,8)+q(3,8))/(q(4,8)-q(3,8)))
c      y9=0.5d0*dlog((q(4,9)+q(3,9))/(q(4,9)-q(3,9)))
c      y10=0.5d0*dlog((q(4,10)+q(3,10))/(q(4,10)-q(3,10)))
c      y11=0.5d0*dlog((q(4,11)+q(3,11))/(q(4,11)-q(3,11)))

      y8=0.5d0*dlog((pmod8+q(3,8))/(pmod8-q(3,8)))
      y9=0.5d0*dlog((pmod9+q(3,9))/(pmod9-q(3,9)))
      y10=0.5d0*dlog((pmod10+q(3,10))/(pmod10-q(3,10)))
      y11=0.5d0*dlog((pmod11+q(3,11))/(pmod11-q(3,11)))

c$$$      if(y8.lt.2d0)return
c$$$      if(y9.lt.2d0)return
c$$$      if(y8.gt.4.5d0)return
c$$$      if(y9.gt.4.5d0)return
c$$$      if(y10.lt.2d0)return
c$$$      if(y11.lt.2d0)return
c$$$      if(y10.gt.4.5d0)return
c$$$      if(y11.gt.4.5d0)return

      elseif(dps.eq.3)then

      et1=dsqrt(q(1,6)**2+q(2,6)**2)
      et2=dsqrt(q(1,7)**2+q(2,7)**2)
      et3=dsqrt(q(1,8)**2+q(2,8)**2)

      pmod1=dsqrt(q(1,6)**2+q(2,6)**2+q(3,6)**2)
      pmod2=dsqrt(q(1,7)**2+q(2,7)**2+q(3,7)**2)
      pmod3=dsqrt(q(1,8)**2+q(2,8)**2+q(3,8)**2)

      eta1=0.5d0*dlog((pmod1+q(3,6))/(pmod1-q(3,6)))
      eta2=0.5d0*dlog((pmod2+q(3,7))/(pmod2-q(3,7)))
      eta3=0.5d0*dlog((pmod3+q(3,8))/(pmod3-q(3,8)))

      if(et1.lt.ptamin3)return
      if(et2.lt.ptbmin3)return
      if(et3.lt.ptcmin3)return
      if(eta1.gt.etaamax3)return
      if(eta2.gt.etabmax3)return
      if(eta3.gt.etacmax3)return
      if(eta1.lt.etaamin3)return
      if(eta2.lt.etabmin3)return
      if(eta3.lt.etacmin3)return
     
cccccccccc kt/anti-kt/durham alg

      if(jalg.eq.'kt')then

      d1=q(1,6)**2+q(2,6)**2
      dmin1=d1
      d2=q(1,7)**2+q(2,7)**2
      dmin1=min(d2,dmin1)
      d3=q(1,8)**2+q(2,8)**2
      dmin1=min(d3,dmin1)

      elseif(jalg.eq.'antikt')then

      d1=1d0/(q(1,6)**2+q(2,6)**2)
      dmin1=d1
      d2=1d0/(q(1,7)**2+q(2,7)**2)
      dmin1=min(d2,dmin1)
      d3=1d0/(q(1,8)**2+q(2,8)**2)
      dmin1=min(d3,dmin1)

      elseif(jalg.eq.'Durham')then

      d1=1d0
      d2=1d0
      d3=1d0
      dmin1=1d0

      endif

ccccccccc  

      d12=min(d1,d2)*((eta1-eta2)**2+dphi(6,7)**2)/rjet**2
      dmin2=d12
      d13=min(d1,d3)*((eta1-eta3)**2+dphi(6,8)**2)/rjet**2
      dmin2=min(d13,dmin2)
      d23=min(d2,d3)*((eta2-eta3)**2+dphi(7,8)**2)/rjet**2
      dmin2=min(d23,dmin2)

      if(dmin1.lt.dmin2)then
         jflag=1
      else
         jflag=2
      endif

      if(jflag.eq.2)return


cccccccccccccccccccccc

      d34=eta1-eta2
      d35=eta1-eta3
      d45=eta2-eta3

      a34=sinh(d34/2d0)**2/(cosh(d35/2d0)**2+cosh(d45/2d0)**2)
      a35=sinh(d35/2d0)**2/(cosh(d34/2d0)**2+cosh(d45/2d0)**2)
      a45=sinh(d45/2d0)**2/(cosh(d35/2d0)**2+cosh(d34/2d0)**2)      

      phicut=60d0*pi/180d0
      acut=0.7d0
      ccut=4d0

      amin=0.7d0
      amax=1.3d0

c      if(dphi(6,7).lt.phicut.and.a34.lt.acut)return
c      if(dphi(6,8).lt.phicut.and.a35.lt.acut)return
c      if(dphi(7,8).lt.phicut.and.a45.lt.acut)return

c$$$      if(dphi(6,7).lt.phicut.and.cosh(d34).lt.ccut)return
c$$$      if(dphi(6,8).lt.phicut.and.cosh(d35).lt.ccut)return
c$$$      if(dphi(7,8).lt.phicut.and.cosh(d45).lt.ccut)return
c$$$c$$$
c$$$      if(dphi(6,7).gt.phicut.and.dphi(6,8).gt.phicut.and.
c$$$     &     dphi(7,8).gt.phicut)return

c      if(cosh(d34).lt.ccut.and.cosh(d35).lt.ccut.and.
c     &     cosh(d45).lt.ccut)return

c      if(dabs(d34).lt.2.6d0.and.dabs(d45).lt.2.6d0.and.dabs(d35).
c     &     lt.2.6d0)return

c      if(a34.lt.amin.and.a35.lt.amin.and.a45.lt.amin)return
c      if(a34.gt.amax.and.a35.gt.amax.and.a45.gt.amax)return

c$$$      accut=.true.
c$$$
c      if(a34.gt.amin.and.a34.lt.amax)accut=.false.
c      if(a35.gt.amin.and.a35.lt.amax)accut=.false.
c      if(a45.gt.amin.and.a45.lt.amax)accut=.false.
c$$$
c$$$      if(a34.gt.amin.and.a34.lt.amax)then
c$$$         if(dphi(6,7).lt.phicut)return
c$$$      endif
c$$$      if(a35.gt.amin.and.a35.lt.amax)then
c$$$         if(dphi(6,8).lt.phicut)return
c$$$      endif
c$$$      if(a45.gt.amin.and.a45.lt.amax)then
c$$$         if(dphi(7,8).lt.phicut)return
c$$$      endif
c$$$
c      if(accut)return

      elseif(decay3)then

      et1=dsqrt(q(1,6)**2+q(2,6)**2)
      et2=dsqrt(q(1,8)**2+q(2,8)**2)
      et3=dsqrt(q(1,9)**2+q(2,9)**2)

      pmod1=dsqrt(q(1,6)**2+q(2,6)**2+q(3,6)**2)
      pmod2=dsqrt(q(1,8)**2+q(2,8)**2+q(3,8)**2)
      pmod3=dsqrt(q(1,9)**2+q(2,9)**2+q(3,9)**2)

      eta1=0.5d0*dlog((pmod1+q(3,6))/(pmod1-q(3,6)))
      eta2=0.5d0*dlog((pmod2+q(3,8))/(pmod2-q(3,8)))
      eta3=0.5d0*dlog((pmod3+q(3,9))/(pmod3-q(3,9)))

  
      if(et1.lt.ptamin3)return
      if(et2.lt.ptbmin3)return
      if(et3.lt.ptcmin3)return
      if(eta1.gt.etaamax3)return
      if(eta2.gt.etabmax3)return
      if(eta3.gt.etacmax3)return
      if(eta1.lt.etaamin3)return
      if(eta2.lt.etabmin3)return
      if(eta3.lt.etacmin3)return
           
      
  
      elseif(dps.eq.1)then


      endif

       icut=1

      return
      end
 
