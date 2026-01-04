      subroutine findtrop(nz,temp,hght,ktrop)
c
c----CAMx v7.20 220430
c
c     FINDTROP finds the tropopause in each vertical column based in input T/P
c     and layer heights. The method is based on the WMO definition of
c     the tropopause where temperature gradient changes from <-2K/km to
c     >-2K/km over at least a 2 km depth.
c
c      Copyright 1996 - 2022
c     Ramboll
c
c     Modifications:
c        none
c
c     Input arguments:
c        nz                  Number of layers
c        temp                Temperature (K)
c        hght                Layer interface heights (m)
c
c     Output arguments:
c        ktrop               Layer index containing tropopause
c
c     Routines called:
c        none
c
c     Called by:
c        STRATO3
c
      implicit none
      integer nz
      integer ktrop
      real temp(nz),hght(nz)

      integer k,kp,km
      real dtdz1,dtdz2,zmid,zmidkp,zmidkm,zp2k,zm2k
c
c-----Entry point
c
c-----Scan vertical T-profile above 5 km
c
      do k = 1,nz-1
        if (hght(k).lt.5000.) cycle 
        zmid = (hght(k) + hght(k-1))/2.
        zp2k = zmid + 2000.          ! Alt 2km above this level
        do kp = k+1,nz
          zmidkp = (hght(kp) + hght(kp-1))/2.
          if (zmidkp.ge.zp2k .or. kp.eq.nz) exit
        enddo
        dtdz1 = (temp(kp) - temp(k))/(zmidkp - zmid)
c
c-----Check T-grad above and below layer and set layer index if trop found
c
        if (dtdz1.gt.-0.002) then
          ktrop = k
          goto 200
        endif
      enddo
c
c-----If we went through whole profile, set trop to top layer
c
      ktrop = nz
 200  continue

      return
      end
