      subroutine chkpartial()
      use grid
      use tracer
      use filunit
      implicit none
c
c----CAMx v7.20 220430
c
c     This routione checks for the integrity of the data read in the 
c     partial source area maps. The sum of the fractions in each grid
c     cell must be with an acceptable range of 100%. If it within the 
c     range then the values are normalized to ensure that the total 
c     equals 100%, to machine precision.
c
c      Copyright 1996 - 2022
c     Ramboll
c          
c     Modifications:
c
c     Input arguments:
c
c     Output arguments:
c
c     Routines called:
c
c     Called by:
c
      include "camx.prm"
c
c-----External Functions
c
      integer istrln
c
c-----Local parameters ----
c
      real FUZZ
      parameter( FUZZ = 0.02 )
c
c-----Local variables ----
c
      integer icl, jcl, igrp, igrd, ipart
      real    sum_frac
c
c-----Entry point
c
c  --- loop over grids ---
c
      do igrd=1,ngrid
c
c  --- loop over groups, initializing the count ---
c
         do igrp=1,ngroup
c
c  --- loop over cells ---
c
            do icl=1,nxcell(igrd) 
              do jcl=1,nycell(igrd) 
                sum_frac = 0.
                do ipart=1,npartial(igrp,igrd)
                   sum_frac = sum_frac + frcmap(igrp,ipart,igrd,icl,jcl)
                enddo
c
c  --- if sum is zero then use regular model source map ---
c
                if( sum_frac .EQ. 0. ) then
                   frcmap(igrp,1,igrd,icl,jcl) = 1.0
                   igrmap(igrp,1,igrd,icl,jcl) = igrmap(0,1,igrd,icl,jcl)
                   sum_frac = 1.0
                endif
c
c  --- if sum is outside range, write an error ---
c
                if( ABS(1.0 - sum_frac) .GT. FUZZ) then
                   write(iout,'(//,a)') 'ERROR in CHKPARTIAL:'
                   write(iout,'(/,1X,2A)') 'The sum of the fractional ',
     &                                      'grid values do not add to 100%'
                   write(iout,'(1X,A,I5,A,I5)') 'Grid: ',igrd,
     &                                          '   Emissions group: ',igrp
                   write(iout,'(1X,A,I3,A,I3,A)') 'Cell: (',icl,',',jcl,')'
                   write(iout,'(1X,A,F8.2,A)') 'The total for this cell is: ',
     &                                                       sum_frac*100.,'%'
                   call camxerr()
                 endif
c
c  --- normalize the values to make sure it is exactly 100% ---
c
                 
                do ipart=1,npartial(igrp,igrd)
                   frcmap(igrp,ipart,igrd,icl,jcl) = 
     &                       frcmap(igrp,ipart,igrd,icl,jcl) / sum_frac
                enddo
c
c  --- next grid cell ---
c
              enddo
            enddo
c
c ---- next group ---
c
         enddo
c
c ---- next grid ---
c
      enddo
c
      return
      end
