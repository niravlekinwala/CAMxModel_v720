      subroutine rassgn2d(ncol,nrow,io,jo,nmesh,ncolf,nrowf,cval,fval)
c
c----CAMx v7.20 220430
c
c     RASSGN2D assigns real fine grid values from coarse grid
c
c      Copyright 1996 - 2022
c     Ramboll
c          
c     Modifications:
c        none
c
c     Input arguments:
c        ncol              number of columns in the parent grid
c        nrow              number of rows in the parent grid
c        io                starting i index for the fine grid
c        jo                starting j index for the fine grid
c        nmesh             mesh number
c        ncolf             number of columns in fine grid
c        nrowf             number of rows in fine grid
c        cval              cell centered value on coarse grid
c
c     Output arguments:
c        fval              cell centered value on coarse grid
c
c     Subroutine called:
c        none
c
c     Called by:
c        CAMx
c        STARTUP
c
      real cval(ncol,nrow)
      real fval(ncolf,nrowf)
c
c-----Entry point
c
        do 40 jfin = 2,nrowf-1
          j = (jfin - 2)/nmesh + jo
          do 30 ifin = 2,ncolf-1
            i = (ifin - 2)/nmesh + io
            fval(ifin,jfin) = cval(i,j)
  30      continue
  40    continue
c
      return
      end
