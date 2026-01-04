      subroutine init_Lslice(iproc_id,igrd,ioff,joff)
c
      use pigsty
      use node_mod
c
      implicit none
c
c----CAMx v7.20 220430
c
c     This routine initializes the LSlice array for the first time step
c     of the day. This is only needed to identify the slice that owns 
c     the PiG when LVISIG is turned on. The AVEPIG routine needs to know
c     which slice owns each PiG.
c
c      Copyright 1996 - 2022
c     Ramboll
c
c     Modifications:
c
c     Input arguments:
c        iproc_id            process ID
c        igrd                grid index 
c        ioff                X offset of this slice in complete domain
c        ioff                Y offset of this slice in complete domain
c
c     Output arguments:
c
c     Subroutines Called:
c        PIGCOORD
c
c     Called by:
c        AVEPIG
c
      include "camx.prm"
      include "flags.inc"
c
c-----Argunment delarations
c
      integer iproc_id
      integer igrd
      integer ioff
      integer joff
c
c------Local variables
c
      integer n, i, j, idum
      real    xpig, ypig
c
c-----Entry point
c
      do 10 n = 1,npig
        if( lmpi ) Lslice(igrd,iproc_id,n) = 0
c
c-----Skip puffs that are not in the current grid
c
        if(ingrd(n) .NE. igrd) goto 10
c
        xpig = (xpigf(n) + xpigb(n))/2.
        ypig = (ypigf(n) + ypigb(n))/2.
        call pigcoord(xpig,ypig,i,j,idum)
c
c--- adjust cell index and skip if PiG not in this slice ---
c
        i = i-ioff
        j = j-joff
        if( i .LT. ia .OR. i .GT. iz ) goto 10
        if( j .LT. ja .OR. j .GT. jz ) goto 10
        if( lmpi ) Lslice(igrd,iproc_id,n) = 1
   10 continue
c
      return
      end
