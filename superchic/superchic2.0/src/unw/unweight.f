ccc   writes event information to array for unweighted generation
      subroutine unweight(wt,r)
      implicit double precision(a-y)
      integer i,j

      include 'record.f'
      include 'unweighted.f'
      include 'mom.f'
      include 'proc.f'
      include 'decay.f'
      include 'hepevt.f'
      include 'leshouches.f'
      include 'wmax.f'

         if(wt/wmax.gt.r)then

           evnum=evnum+1
          
            if(proc.eq.17.or.proc.eq.18.or.proc.eq.19)then
               nup=11
            elseif(decay4)then
               if(proc.eq.41)then
                  nup=10
               elseif(proc.eq.42.or.proc.eq.43)then
                  nup=11
               endif
            elseif(dps.eq.2.or.decay2)then
               nup=7
            elseif(dps.eq.3)then
               nup=8
            elseif(decay3)then
               nup=9
            elseif(dps.eq.1)then
               nup=5
            endif

            do j=3,nup
               do i=1,4
                 evrec(evnum,j,i)=q(i,j)
               enddo
            enddo

         endif

         nhep=nup

      return
      end
