c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
c     This file is part of the Test Set for IVP solvers
c     http://www.dm.uniba.it/~testset/
c
c        Car Axis problem (in index 3 formulation)
c        index 3 DAE of dimension 10
c
c     DISCLAIMER: see
c     http://www.dm.uniba.it/~testset/disclaimer.php
c
c     The most recent version of this source file can be found at
c     http://www.dm.uniba.it/~testset/src/problems/caraxis.f
c
c     This is revision
c     $Id: caraxis.F,v 1.2 2006/10/02 10:29:13 testset Exp $
c
c-----------------------------------------------------------------------
c----------------------------------------------------------------------

c-----------------------------------------------------------------------
c     residual function
c-----------------------------------------------------------------------

      SUBROUTINE carres(X,Y,YPRIME,cj,DELTA,ier,RPAR,IPAR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	    parameter(N=10)

      DIMENSION Y(N),DELTA(N),IPAR(*),RPAR(*),YPRIME(N)
      double precision m,eps, k
C
      M = 10d0
      eps = 1d-2
      k = m*eps*eps/2d0
      
      call carfunc(n,x,y,delta,rpar,ipar)
      do i=1,4
         delta(i) =yprime(i)-delta(i)
      enddo   
      do i=5,8
         delta(i) = k*yprime(i)- delta(i)
      enddo
      do i=9,10
         delta(i) = -delta(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine carfunc(neqn,t,y,f,rpar,ipar)
      integer neqn,ierr,ipar(*), i
      double precision t,y(neqn),f(neqn),rpar(*)      
      double precision M,L,L0,w,r,xb,yb,Ll,Lr,eps,g,
     +                 xl,yl,xr,yr,lam1,lam2

      eps = 1d-2
      M   = 10d0
      L   = 1d0
      L0  = 0.5d0
      r   = 0.1d0
      w   = 10d0
      g   = 1d0
      yb  = r*sin(w*t)
      xb  = sqrt(L*L-yb*yb)

      do 10 i=1,4
         f(i) = y(i+4)
   10 continue

      xl   = y(1)
      yl   = y(2)
      xr   = y(3)
      yr   = y(4)
      lam1 = y(9)
      lam2 = y(10)

      Ll = sqrt(xl**2+yl**2)
      Lr = sqrt((xr-xb)**2+(yr-yb)**2)

      f(5)  =(L0-Ll)*xl/Ll +lam1*xb+2d0*lam2*(xl-xr)
      f(6)  =(L0-Ll)*yl/Ll +lam1*yb+2d0*lam2*(yl-yr)-M*eps*eps*g/2d0
      f(7)  =(L0-Lr)*(xr-xb)/Lr -2d0*lam2*(xl-xr)
      f(8)  =(L0-Lr)*(yr-yb)/Lr -2d0*lam2*(yl-yr)-M*eps*eps*g/2d0

      f(9)  = xb*xl+yb*yl
      f(10) = (xl-xr)**2+(yl-yr)**2-L*L

      return
      end

c-----------------------------------------------------------------------
      subroutine carsoln(neqn, y)
      integer neqn
      double precision  y(neqn)

      y(  1) =  0.4934557842755629d-001
      y(  2) =  0.4969894602303324d+000
      y(  3) =  0.1041742524885400d+001
      y(  4) =  0.3739110272652214d+000
      y(  5) = -0.7705836840321485d-001
      y(  6) =  0.7446866596327776d-002
      y(  7) =  0.1755681574942899d-001
      y(  8) =  0.7703410437794031d+000
      y(  9) = -0.4736886750784630d-002
      y( 10) = -0.1104680411345730d-002

      return
      end
c----------------------------------------------------------------------------
 