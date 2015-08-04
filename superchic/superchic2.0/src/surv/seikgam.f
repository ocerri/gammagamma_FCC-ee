ccc   integrates bare + screened amplitude over k_t 
ccc   (two-photon induced processes)
      subroutine schimcgam(p1x,p1y,p2x,p2y,outg)
      implicit double precision(a-y)
      implicit complex*16(z)
      integer jx,jy,i1,i2,p,i,nphi,nqt,jqt,jphi
      complex*16 pp,mm,pm,mp
      complex*16 out(4,10),x0(10),x00(10),outg(10)
      double precision outgd(10)
      complex*16 ppa(10),mma(10),pma(10),mpa(10)
      complex*16 screen(2,2)
      integer nni
      complex*16 outt(0:200,0:10)

      include 'nchan.f'
      include 'surv.f'
      include 'vars.f'
      include 'survpars.f'
      include 'polarization.f'
      include 'photo.f'
      include 'bpsi.f'
      include 'proc.f'
      include 'zi.f'
      include 'mandelstam.f'
      include 'pi.f'
      include 'nsurv.f'
      include 'inparticle.f'

      nphi=s2int
      nqt=s2int*4

      hphi=2d0*pi/dble(nphi)
      hqt=2d0/dble(nqt)
      
      zoutg=0d0
      
      do i1=0,nqt
         do p=0,10
            outt(i1,p)=0d0
         enddo
      enddo

      do p=1,pol
         outg(p)=0d0
         do i=1,4
            out(i,p)=0d0
         enddo
      enddo

      if(proc.eq.42.or.proc.eq.43.or.proc.eq.48.or.proc.eq.49)then
         call ww(mx,uh,th,pinc)
         do p=1,pol                       
            call wwpol(p,uh,th,pp,mm,pm,mp)
            ppa(p)=pp
            mma(p)=mm
            pma(p)=pm
            mpa(p)=mp     
         enddo
      elseif(proc.eq.44.or.proc.eq.45.or.proc.eq.46.or.
     &        proc.eq.50.or.proc.eq.51.or.proc.eq.52)then
         call ll(mx,uh,th,pinc)
         do p=1,pol                       
            call llpol(p,mx,uh,th,pp,mm,pm,mp)
            ppa(p)=pp
            mma(p)=mm
            pma(p)=pm
            mpa(p)=mp       
         enddo
      elseif(proc.eq.47.or.proc.eq.53)then
         call lightlight(mx,uh,th,pinc)
         do p=1,pol                       
            call lightlightpol(p,mx,uh,th,pp,mm,pm,mp)
            ppa(p)=pp
            mma(p)=mm
            pma(p)=pm
            mpa(p)=mp
         enddo
      endif

      if(sfac)then

         do jqt=0,nqt

            qt=(dble(jqt))*hqt

            tp2=qt**2
            
            do i1=1,nch
               do i2=1,nch
                  call screeningint(i1,i2,tp2,sc,sc1)  
                  screen(i1,i2)=sc
               enddo
            enddo

            do jphi=1,nphi
               
               phiq=(dble(jphi)-0.5d0)*hphi
               
               tpx=qt*dcos(phiq)
               tpy=qt*dsin(phiq)
               wt=hphi*qt*hqt
               
               p1xp=p1x-tpx
               p1yp=p1y-tpy
               t12=p1xp**2+p1yp**2
               p2xp=tpx+p2x
               p2yp=tpy+p2y
               t22=p2xp**2+p2yp**2
               
               call formfacgam(1,t12,t22,x00p)
               call formfacgam(2,t12,t22,x00p2)

           do p=1,pol

              do i1=1,nch
                 do i2=1,nch
                    
                    x0(p)=x00p*pp0(i1)*pp0(i2)/dble(nch)**2/gaa(i1)
     &                   /gaa(i2)

                       zout=-0.5d0*(ppa(p)+mma(p))*(p1xp*p2xp+p1yp*p2yp)
     &                   -0.5d0*zi*(ppa(p)-mma(p))*(p1xp*p2yp-p2xp*p1yp)
     &                      +0.5d0*(p1xp*p2xp-p1yp*p2yp
     &                      +zi*(p1xp*p2yp+p1yp*p2xp))*mpa(p)
     &                      +0.5d0*(p1xp*p2xp-p1yp*p2yp
     &                      -zi*(p1xp*p2yp+p1yp*p2xp))*pma(p)
                     
                       zout=zout*2d0
                    
                       outt(jqt,p)=outt(jqt,p)+x0(p)*wt
     &                      *screen(i1,i2)*zout

                       if(p.eq.1)then
                          x00p2=x00p2*pp0(i1)*pp0(i2)/dble(nch)**2
     &                         /gaa(i1)/gaa(i2)
                          outt(jqt,0)=outt(jqt,0)+x00p2*wt*screen(i1,i2)
                       endif

                 enddo
              enddo
           enddo
      enddo
     
      enddo

      do jqt=0,nqt-2,2
         do p=1,pol
            outg(p)=outg(p)+(outt(jqt,p)+outt(jqt+1,p)*4d0
     &           +outt(jqt+2,p))/3d0
         enddo
         zoutg=zoutg+(outt(jqt,0)+outt(jqt+1,0)*4d0
     &           +outt(jqt+2,0))/3d0
      enddo

      endif

      t11=p1x**2+p1y**2
      t22=p2x**2+p2y**2
  
      if(inparticle.eq.'prot')then
         call formfacgam(1,t11,t22,x00p)
         call formfacgam(2,t11,t22,x00p2)
      elseif(inparticle.eq.'el')then
         call formfacgamel(1,t11,t22,x00p)
         call formfacgamel(2,t11,t22,x00p2)
      endif

      zoutg=zoutg+x00p2*dsqrt(pinc)  
      
      outgg1=cdabs(zoutg)**2
      
      do p=1,pol
         zout=-0.5d0*(ppa(p)+mma(p))*(p1x*p2x+p1y*p2y)
     &        -0.5d0*zi*(ppa(p)-mma(p))*(p1x*p2y-p2x*p1y)
     &        +0.5d0*(p1x*p2x-p1y*p2y+zi*(p1x*p2y+p1y*p2x))*mpa(p)
     &        +0.5d0*(p1x*p2x-p1y*p2y-zi*(p1x*p2y+p1y*p2x))*pma(p)
         
         outg(p)=outg(p)+zout*x00p*2d0
         dbl=cdabs(outg(p))
         dbl=dsqrt(dbl**2+outgg1/dble(pol))
         outg(p)=dbl
         enddo
     
      return
      end
