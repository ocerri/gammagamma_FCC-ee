ccc   calculates CEP cross section
      function cs(rarr,wgt)
      implicit double precision (a-z)
      complex*16 wt(10)
      double precision rarr(10),ran(5)
      integer i,p

      include 'polvecs.f'
      include 'gencuts.f'
      include 'pi.f'
      include 'unweighted.f'
      include 'bin.f'
      include 'x.f'
      include 'zi.f'
      include 'vars.f'
      include 'survpars.f'
      include 'mom.f'
      include 'mandelstam.f'
      include 'range.f'
      include 'polarization.f'
      include 'proc.f'
      include 'mq.f'
      include 'norm.f'
      include 'meta.f'
      include 'mres.f'
      include 'decay.f'
      include 'mp.f'
      include 'quarkonia.f'
      include 'scorr.f'
      include 'photo.f'
      include 'bpsi.f'
      include 'mpip.f'
      include 'gamma.f'
      include 'ewpars.f'
      include 'eff.f'
      include 'brs.f'
      include 'wmax.f'
      include 'unwsurv.f'
      include 'survin.f'
      include 'inparticle.f'

      wtt=0d0

      if(dps.eq.1)then
         mx=mres
      else
         rm=rarr(5)
         mx=mmin+(mmax-mmin)*rm
         if(mmin.lt.1d-3)mmin=1d-3
         msub=1d0/1d0*(1d0/mmax**1d0+
     &           (1d0/mmin**1d0-1d0/mmax**1d0)*rm)
         mx=(1d0/msub/1d0)**(1d0/1d0)
         if(mx.lt.2d0*mq)goto 777
      endif

      if(inparticle.eq.'prot')mpp=mp
      if(inparticle.eq.'el')mpp=me

      if(photo)then
         
         r1=rarr(2)
         r3=rarr(3)
         r4=rarr(4)
         r5=rarr(5)
         
         call r2455(r2)
         
         ptmax=dsqrt(5d0)
         ptmin=0d0
         
         pt2sq=(ptmax-ptmin)*r1+ptmin
         pt2sq=pt2sq**2   
            
         phi1=2d0*pi*r2
         phi2=2d0*pi*r3+phi1
         
         rmx=dsqrt(pt2sq+mx**2)
         
         xgmin=(mx/rts)**2
         ypmax=dlog(xgmin**2*mpp**2+ptmax**2)
         ypmin=dlog(xgmin**2*mpp**2)
         yp=(ypmax-ypmin)*r4+ypmin
         
         wtpt=2d0*dsqrt(pt2sq)*(ptmax-ptmin)
         
         if(r5.gt.0.5d0)then  
            pt1sq=dexp(yp)-xgmin**2*mpp**2
            wt1=xgmin**2*mpp**2+pt1sq
         else                    
            pt1sq=pt2sq         
            pt2sq=dexp(yp)-xgmin**2*mpp**2       
            wt1=xgmin**2*mpp**2+pt2sq
         endif
         
         pt2x=dsqrt(pt2sq)*dcos(phi2)
         pt2y=dsqrt(pt2sq)*dsin(phi2)
         pt1x=dsqrt(pt1sq)*dcos(phi1)
         pt1y=dsqrt(pt1sq)*dsin(phi1)
         
         ptxx=(pt1x+pt2x)**2+(pt1y+pt2y)**2
         rmx=dsqrt(ptxx+mx**2)
            
         ymax1=ymax
         ymin1=ymin
         ycut=dlog(rts/rmx)
            
         if(ymax.gt.ycut) ymax=ycut
         if(-ymin.gt.ycut) ymin=-ycut
            
         ry=rarr(1)
         yx=ymin+(ymax-ymin)*ry
            
         wty=ymax-ymin

         ymin=ymin1
         ymax=ymax1

         if(r5.gt.0.5d0)then    ! photon emitted from q(1,k) 
            xgam=rmx*dexp(yx)/rts ! photon mom. fraction 
            wpsi=dsqrt(xgam*s)  ! proton-photon cms energy
            xglu=(rmx/wpsi)**2  ! gluon mom. fraction
            x1=xgam
            x2=xglu
            prot=1
         else                   ! photon emitted from q(1,k)
            xgam=rmx*dexp(-yx)/rts
            wpsi=dsqrt(xgam*s)  
            xglu=(rmx/wpsi)**2
            x2=xgam
            x1=xglu
            prot=2
         endif
         
         bpsi=bpsi0+4d0*alphapb*dlog(wpsi/w0b)
         
      elseif(gamma)then

         r2=rarr(2)
         r3=rarr(3)
         r4=rarr(4)
           
         call r2455(r1)
         
         phi1=2d0*pi*r1
         phi2=2d0*pi*r2+phi1
         
         if(inparticle.eq.'prot')then
            ptmax=dsqrt(3d0)       
         elseif(inparticle.eq.'el')then
            ptmax=dsqrt(50d0)       
         endif

         ptmin=0d0
            
         xgmin=(mx/rts)**2

         ypmax=dlog(xgmin**2*mpp**2+ptmax**2)
         ypmin=dlog(xgmin**2*mpp**2)
            
         yp=(ypmax-ypmin)*r3+ypmin
         ypp=(ypmax-ypmin)*r4+ypmin
         
         pt1sq=dexp(yp)-xgmin**2*mpp**2
         pt2sq=dexp(ypp)-xgmin**2*mpp**2

         if(pt1sq.lt.0d0)then
            pt1sq=0d0
         endif
         
         if(pt2sq.lt.0d0)then
            pt2sq=0d0
         endif
         
         pt2x=dsqrt(pt2sq)*dcos(phi2)
         pt2y=dsqrt(pt2sq)*dsin(phi2)
         pt1x=dsqrt(pt1sq)*dcos(phi1)
         pt1y=dsqrt(pt1sq)*dsin(phi1)
         
         ptxx=(pt1x+pt2x)**2+(pt1y+pt2y)**2
         rmx=dsqrt(ptxx+mx**2)
         
         ymax1=ymax
         ymin1=ymin
         ycut=dlog(rts/rmx)
         
         if(ymin.gt.0d0.and.ycut.lt.ymin)goto 777
         if(ymax.lt.0d0.and.ycut.gt.ymax)goto 777
            
         if(ymax.gt.ycut) ymax=ycut
         if(-ymin.gt.ycut) ymin=-ycut
            
         ry=rarr(1)
         yx=ymin+(ymax-ymin)*ry
         
         wty=ymax-ymin
         
         ymin=ymin1
         ymax=ymax1
         
         x1=rmx*dexp(yx)/rts    ! photon 1 mom. fraction 
         x2=rmx*dexp(-yx)/rts   ! photon 2 mom. fraction 
             
         qsq=(x1**2*mpp**2+pt1sq)/(1d0-x1)
         qsqp=(x2**2*mpp**2+pt2sq)/(1d0-x2)
         wpsi=dsqrt(x1*x2*s)
         
      else
         
         r1=rarr(2)
         r2=rarr(3)
         r4=rarr(4)
      
         call r2455(r3)
         
         ptmax=dsqrt(2d0)
         
         pt1sq=r1*ptmax
         pt2sq=r2*ptmax
         
         pt1sq=pt1sq**2
         pt2sq=pt2sq**2
         
         phi1=2d0*pi*r3
         phi2=2d0*pi*r4+phi1
            
         pt1x=dsqrt(pt1sq)*dcos(phi1)
         pt1y=dsqrt(pt1sq)*dsin(phi1)
         pt2x=dsqrt(pt2sq)*dcos(phi2)
         pt2y=dsqrt(pt2sq)*dsin(phi2)
         
         ptxsq=(pt1x+pt2x)**2+(pt1y+pt2y)**2
         rmx=dsqrt(mx**2+ptxsq)

cccccccccccccccccccccc
      
         ymax1=ymax
         ymin1=ymin
         ycut=dlog(rts/rmx)
         
         if(ymin.gt.0d0.and.ycut.lt.ymin)goto 777
         if(ymax.lt.0d0.and.ycut.gt.ymax)goto 777
         
         if(ymax.gt.ycut) ymax=ycut
         if(-ymin.gt.ycut) ymin=-ycut
         
         ry=rarr(1)
         yx=ymin+(ymax-ymin)*ry
         
         ymin=ymin1
         ymax=ymax1
         
         x1=rmx/rts*dexp(yx)
         x2=rmx/rts*dexp(-yx)
 
      endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccc

      aa1=(1d0-x1)*rts/dsqrt(2d0)
      aa2=(1d0-x2)*rts/dsqrt(2d0)
      cc1=0.5d0*(pt2sq+mpp**2)
      cc2=0.5d0*(pt1sq+mpp**2)

c     impose massive on-shell condition by solving
c                   p1+ + cc1/p2- = aa1
c                   p2- + cc2/p1+ = aa2 

      root1sq=(cc1-cc2-aa1*aa2)**2-4d0*cc2*aa1*aa2
      root2sq=(cc2-cc1-aa1*aa2)**2-4d0*cc1*aa1*aa2
      if(root1sq.le.0d0.or.root2sq.le.0d0)then
         weight=0d0
         goto 777
      endif
      p1p=(cc2-cc1+aa1*aa2+dsqrt(root1sq))/(2d0*aa2)
      p2m=(cc1-cc2+aa1*aa2+dsqrt(root2sq))/(2d0*aa1)
      p1m=(pt1sq+mpp**2)/(2d0*p1p)
      p2p=(pt2sq+mpp**2)/(2d0*p2m)

      q(1,3)=pt1x
      q(2,3)=pt1y
      q(3,3)=(p1p-p1m)/dsqrt(2d0)
      q(4,3)=(p1p+p1m)/dsqrt(2d0)

      q(1,4)=pt2x
      q(2,4)=pt2y
      q(3,4)=(p2p-p2m)/dsqrt(2d0)
      q(4,4)=(p2p+p2m)/dsqrt(2d0)

      do i=1,4
         q(i,5)=q(i,1)+q(i,2)-q(i,3)-q(i,4)
      enddo

ccccccccccccccccccccccccccccccccc

         if(dps.eq.2)then
            call twojetps(mx,mq,ps,uh,th)
         elseif(dps.eq.3)then
            call threejetps(mx,mq,ps)
         elseif(dps.eq.12)then
            if(proc.eq.18)then
               call twojetpsm(mx,mpsi,mpsip,ps,uh,th)
            else
               call twojetpsm(mx,meta,metap,ps,uh,th)
            endif
         endif

cccccccccccccccccc

cccc decays

ccccccccccccccccccc

         wt2=1d0
         wt3=1d0

         if(proc.eq.1)then
            call twobody(1,5,6,7,mb,mb,wt2)
         elseif(proc.eq.17.or.proc.eq.18.or.proc.eq.19)then
            call twobody(1,6,8,9,mmu,mmu,wt2a)
            call twobody(1,7,10,11,mmu,mmu,wt2b)
            wt2=wt2a*wt2b
         elseif(proc.eq.20.or.proc.eq.21.or.proc.eq.22)then
            call twobody(1,5,6,7,0d0,mpsi,wt2a)
            call twobody(2,7,8,9,mmu,mmu,wt2b)
            wt2=wt2a*wt2b
         elseif(proc.gt.22.and.proc.lt.28)then
            call twobody(1,5,6,7,m2b,m2b,wt2)
         elseif(proc.eq.29.or.proc.eq.30.or.proc.eq.31)then
            call twobody(1,5,6,7,0d0,mups,wt2a)
            call twobody(2,7,8,9,mmu,mmu,wt2b)
            wt2=wt2a*wt2b
         elseif(proc.gt.32.and.proc.lt.37)then
            call twobody(1,5,6,7,m2b,m2b,wt2)
         elseif(proc.eq.38.or.proc.eq.39.or.proc.eq.40)then
            call twobody(1,5,6,7,mmu,mmu,wt2)
         elseif(proc.eq.41)then
            call threebody(5,6,7,8,mpsi,mpip,mpip,wt3)
            call twobody(1,6,9,10,mmu,mmu,wt2a)
            wt3=wt3*wt2a
         elseif(proc.eq.42)then
            call twobodyw(6,8,9,0d0,mmu)
            call twobodyw(7,10,11,0d0,mmu)
         elseif(proc.eq.43)then
            call twobodyw(6,8,9,0d0,me)
            call twobodyw(7,10,11,0d0,me)
         endif

ccccccccccccccccccc  cuts ccccccccccccccccc

         neff0=neff0+1

         if(gencuts)then
            call cut(icut)
            if(icut.eq.0)goto 777
         endif

         neff=neff+1

ccccccccccccccccccccccccccccccccccccccccccc

          if(proc.eq.21.or.proc.eq.24.or.proc.eq.26)then
             call genpol1(5,echi1)
          elseif(proc.eq.22.or.proc.eq.25.or.proc.eq.27)then
             call genpol2
          endif

          if(proc.eq.30.or.proc.eq.33.or.proc.eq.35)then
             call genpol1(5,echi1)
          elseif(proc.eq.31.or.proc.eq.34.or.proc.eq.36)then
             call genpol2
          endif

ccccccccc
          
 456      if(photo)then
             call schimcphot(pt1x,pt1y,pt2x,pt2y,wt)
          elseif(gamma)then
             call schimcgam(pt1x,pt1y,pt2x,pt2y,wt)
          else
             call wtgen
             call schimc(pt1x,pt1y,pt2x,pt2y,wt)
          endif

         wtt=0d0
         do p=1,pol
            wtt=wtt+cdabs(wt(p))**2
         enddo   

         wtpol=1d0


         if(scorr)then
            if(proc.eq.17.or.proc.eq.18.or.proc.eq.19)then
               call jpsidecay(wt,wtt)
            endif
            if(proc.eq.20.or.proc.eq.29)then
               call chic0decay3(wtc0)
               wtt=wtt*wtc0
            endif
            if(proc.eq.21.or.proc.eq.30)then
               do i=4,6
                  wt(i)=conjg(wt(i-3))
               enddo
               call chic1decay3(wt,wtt)
            endif
           if(proc.eq.22.or.proc.eq.31)then
              do i=6,10
                 wt(i)=conjg(wt(i-5))
              enddo
              call chic2decay3(wt,wtt)
           endif
           if(proc.eq.24.or.proc.eq.33)then
              do i=4,6
                 wt(i)=conjg(wt(i-3))
              enddo
              call chic1decay2s(wt,wtt)
           endif
           if(proc.eq.25.or.proc.eq.34)then
              do i=6,10
                 wt(i)=conjg(wt(i-5))
              enddo
              call chic2decay2s(wt,wtt)
           endif
           if(proc.eq.26.or.proc.eq.35)then
              do i=4,6
                 wt(i)=conjg(wt(i-3))
              enddo
              call chic1decay2f(wt,wtt)
           endif
           if(proc.eq.27.or.proc.eq.36)then
              do i=6,10
                 wt(i)=conjg(wt(i-5))
              enddo
              call chic2decay2f(wt,wtt)
           endif
           if(proc.eq.38.or.proc.eq.39.or.proc.eq.40)then
              call jpsidecayphot(wtpol)
           endif
           if(proc.eq.42.or.proc.eq.43)then
              call wwcorr(mx,uh,th,wtpol)
           endif
         endif
    
         wtt=wtt*wt3*wt2

         if(decays)then
            do i=1,nbr
               wtt=wtt*br(i)
            enddo
         endif
          
         if(photo)then
            wtt=wtt*wty
            wtt=wtt*wtpt
            wtt=wtt*(ypmax-ypmin)*wt1
            wtt=wtt*(wpsi/w0)**delta*normp*bpsi
            if(scorr)wtt=wtt*wtpol
         elseif(gamma)then
            wtt=wtt*(ypmax-ypmin)**2
            wtt=wtt*(xgmin**2*mpp**2+pt1sq)*(xgmin**2*mpp**2+pt2sq)
            wtt=wtt*wty
            wtt=wtt*2d0/mx
            wtt=wtt*mx**2*(1d0/mmin**1d0-1d0/mmax**1d0)/1d0
            if(scorr)wtt=wtt*wtpol
         else
            wtt=wtt*(ymax-ymin)
            wtt=wtt*4d0*ptmax**2*dsqrt(pt1sq*pt2sq)*pi**2
         endif

         wtt=wtt/sym

         if(photo)goto 888
         if(gamma)goto 888

ccccccccccccc 1 body phase space

         if(dps.eq.1)then
            wtt=wtt/(16d0**2*pi**5)
         endif

         if(dps.gt.1)then
            wtt=wtt/(64d0*pi**2)*ps
            wtt=wtt*2d0*mx
            wtt=wtt/(16d0**2*pi**6)
            wtt=wtt*mx**2*(1d0/mmin**1d0-1d0/mmax**1d0)/1d0
         endif

         wtt=wtt*conv*surv

 888     val=wtt*wgt
         if(bin)then
         call binit(val)
         endif

         if(calcmax)then
            if(wmax.lt.wtt*wgt*ren)then
               wmax=wtt*wgt*ren
            endif
         endif

         if(unw)then
            call r2455(runw)
            call unweight(wtt*wgt*ren,runw)
         endif

 777     cs=wtt

      return
      end

   
