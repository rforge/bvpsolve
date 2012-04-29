
c ===================================================================================
c print R-messages
c ===================================================================================
C just a string
      subroutine rprint(msg)
      character (len=*) msg
           call dblepr(msg, -1, 0, 0)
      end subroutine 

C printing with one integer and a double
      subroutine rprintid(msg, i1, d1)
      character (len=*) msg
      double precision d1
      integer i1
        call dblepr(msg, -1, d1, 1)
        call intpr(" ", -1, i1, 1)
      end subroutine 

C printing with one double
      subroutine rprintd1(msg, d1)
      character (len=*) msg
      double precision d1
        call dblepr(msg, -1, d1, 1)
      end subroutine 

C printing with two doubles
      subroutine rprintd2(msg, d1, d2)
      character (len=*) msg
      double precision DBL(2), d1, d2
        DBL(1) = d1
        DBL(2) = d2
        call dblepr(msg, -1, DBL, 2)
      end subroutine 

C printing with one integer
      subroutine rprinti1(msg, i1)
      character (len=*) msg
      integer i1
        call intpr(msg, -1, i1, 1)
      end subroutine 

      subroutine rprinti2(msg, i1, i2)
      character (len=*) msg
      INTEGER IN(2), i1, i2
        IN(1) = i1
        IN(2) = i2
        call intpr(msg, -1, IN, 2)
      end subroutine 

      subroutine rprinti3(msg, i1, i2, i3)
      character (len=*) msg
      INTEGER IN(3), i1, i2, i3
        IN(1) = i1
        IN(2) = i2
        IN(3) = i3
        call intpr(msg, -1, IN, 3)
      end subroutine 
