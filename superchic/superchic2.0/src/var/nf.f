ccc   sets number of active flavours
      function nf(qsq)
      implicit double precision (a-z)

      if(qsq.lt.1.96d0)then
      nf=3d0
      elseif(qsq.lt.22.56d0)then
      nf=4d0
      else
      nf=5d0
      endif

      return
      end
