ccc   gg --> rhorho subprocess amplitude (x,y dependent)
      subroutine rhorhoxy(p,x,y,cost,out)
      implicit double precision (a-z)
      integer p

      a=(1d0-x)*(1d0-y)+x*y
      b=(1d0-x)*(1d0-y)-x*y

      if(p.eq.1)then
      out=-cost*(1d0+cost)/(a**2-b**2*cost**2)
      out=out*(4d0/3d0*b**2-3d0/2d0*a)
      else
      out=cost*(1d0-cost)/(a**2-b**2*cost**2)
      out=out*(4d0/3d0*b**2-3d0/2d0*a)
      endif

      pp=0d0
      mm=0d0
     
      return
      end

