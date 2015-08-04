ccc   integrates bare + screened amplitude over k_t 
ccc   (QCD induced processes)
      subroutine schimc(p1x,p1y,p2x,p2y,out)
      implicit double precision(a-y)
      integer jx,jy,i1,i2,p,nphi,nqt,jqt,jphi,nnphi
      complex*16 out(10),x0(10),x00(10),x0p(10)
      complex*16 screen(2,2)
      complex*16 outt(0:200,10)

      include 'nchan.f'
      include 'mt.f'
      include 'surv.f'
      include 'vars.f'
      include 'survpars.f'
      include 'polarization.f'
      include 'nsurv.f'
      include 'pi.f'

      call setmu(mu)

      hx=2d0/dble(nx)

      do p=1,pol
         out(p)=0d0
      enddo

      nphi=s2int
      nqt=s2int*4

      if(nphi.lt.10)nphi=10

      hphi=2d0*pi/dble(nphi)
      hqt=1.5d0/dble(nqt)

      do i1=0,nqt
         do p=1,10
            outt(i1,p)=0d0
         enddo
      enddo

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

            if(jqt.eq.0)then
               nnphi=1
            else
               nnphi=nphi
            endif

           do jphi=1,nnphi
             
               phiq=(dble(jphi)-0.5d0)*hphi
               
               tpx=qt*dcos(phiq)
               tpy=qt*dsin(phiq)
               wt=hphi*qt*hqt

               if(jqt.eq.0)wt=qt*hqt

            p1xp=p1x-tpx
            p1yp=p1y-tpy
            t12=p1xp**2+p1yp**2
            p2xp=tpx+p2x
            p2yp=tpy+p2y
            t22=p2xp**2+p2yp**2

          call bare(mu,p1xp,p1yp,p2xp,p2yp,x0p)

           do p=1,pol
              
              do i1=1,nch
                 do i2=1,nch
                    
                    x0(p)=x0p(p)*pp0(i1)*pp0(i2)/dble(nch)**2
                    x0(p)=x0(p)*dexp(-((t12+0.08d0+bb0(i1))*bex(i1))
     &                   **cc0(i1)+(bex(i1)*(bb0(i1)+0.08d0))**cc0(i1))
                    x0(p)=x0(p)*dexp(-((t22+0.08d0+bb0(i2))*bex(i2))
     &                   **cc0(i2)+(bex(i2)*(bb0(i2)+0.08d0))**cc0(i2)) 
                    
                    outt(jqt,p)=outt(jqt,p)+x0(p)*wt*screen(i1,i2)
                    
                 enddo
              enddo

           enddo

      enddo
      enddo

      do jqt=0,nqt-2,2
         do p=1,pol
            out(p)=out(p)+(outt(jqt,p)+outt(jqt+1,p)*4d0
     &           +outt(jqt+2,p))/3d0
         enddo
      enddo

      endif

      t11=p1x**2+p1y**2
      t22=p2x**2+p2y**2
  
      call formfac(t11,t22,x00p)
      call bare(mu,p1x,p1y,p2x,p2y,x00)

      do p=1,pol
         out(p)=out(p)+x00(p)*x00p
      enddo

      return
      end
