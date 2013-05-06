
      subroutine dblepr_k(label, nchar, data, ndata)
      integer nchar, ndata
      character*(*) label
      double precision data(ndata)
      integer nc
      nc = nchar
      if(nc .lt. 0) nc = len(label)
      call dblepr(label, nc, data, ndata)
      end

      subroutine intpr_k(label, nchar, data, ndata)
      integer nchar, ndata
      character*(*) label
      integer data(ndata)
      integer nc
      nc = nchar
      if(nc .lt. 0) nc = len(label)
      call intpr(label, nc, data, ndata)
      end

      subroutine logpr(msg,len,ldum)
      integer len
      character (len=*) msg
      logical ldum
      character (len=8) logi

      if (ldum) then
         logi ='  TRUE '
      else
         logi ='  FALSE'
      endif



      call dblepr(MSG//logi,-1,0,0)


      end subroutine
c ===================================================================================
c print R-messages
c ===================================================================================
C just a string
      subroutine rprint(msg)
      character (len=*) msg
           call dblepr(msg, -1, 0, 0)
      end subroutine 

C printing with one logical
      subroutine rprintl1(msg, d1)
      character (len=*) msg
      logical d1

      call logpr(msg, -1, d1)

      end subroutine
C printing with one integer and a double
      subroutine rprintid(msg, i1, d1)
      character (len=*) msg
      double precision d1
      integer i1
        call dblepr(msg, -1, d1, 1)
        call intpr(" ", -1, i1, 1)
      end subroutine 

      subroutine rprintid3(msg, i1, d1, d2, d3)
      character (len=*) msg
      double precision dbl(3), d1, d2, d3
      integer i1
        DBL(1) = d1
        DBL(2) = d2
        DBL(3) = d3

        call dblepr(msg, -1, dbl, 3)
        call intpr(" ", -1, i1, 1)
      end subroutine 

      subroutine rprintd3(msg, d1, d2, d3)
      character (len=*) msg
      double precision dbl(3), d1, d2, d3
        DBL(1) = d1
        DBL(2) = d2
        DBL(3) = d3

        call dblepr(msg, -1, dbl, 3)
      end subroutine 

C printing with one double
      subroutine rprintd1(msg, d1)
      character (len=*) msg
      double precision d1, DBL(1)
        DBL(1) = d1
        call dblepr(msg, -1, DBL, 1)
      end subroutine 

C printing with two doubles
      subroutine rprintd2(msg, d1, d2)
      character (len=*) msg
      double precision DBL(2), d1, d2
        DBL(1) = d1
        DBL(2) = d2
        call dblepr_k(msg, -1, DBL, 2)
      end subroutine 

C printing with one integer
      subroutine rprinti1(msg, i1)
      character (len=*) msg
      integer i1, IN(1)
        IN(1) = i1
        call intpr_k(msg, -1, IN, 1)
      end subroutine 

      subroutine rprinti2(msg, i1, i2)
      character (len=*) msg
      INTEGER IN(2), i1, i2
        IN(1) = i1
        IN(2) = i2
        call intpr_k(msg, -1, IN, 2)
      end subroutine 

      subroutine rprinti3(msg, i1, i2, i3)
      character (len=*) msg
      INTEGER IN(3), i1, i2, i3
        IN(1) = i1
        IN(2) = i2
        IN(3) = i3
        call intpr_k(msg, -1, IN, 3)
      end subroutine 
