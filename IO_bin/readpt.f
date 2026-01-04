      subroutine readpt( )
      use filunit
      use grid
      use chmstry
      use camxcom
      use ptemiss
c
c----CAMx v7.20 220430
c
c     READPT reads the time-variant records of the point source file and
c     cycles through to current time/date to load point emission rates
c
c      Copyright 1996 - 2022
c     Ramboll
c          
c     Modifications:
c     05/01/03  Time span of emissions must now match emiss update interval
c     04/21/04  Now reads effective plume height which may overide plume
c               height calculation
c     04/14/04  Removed uncessary and errant check for plume rise distribution
c               flags
c     12/07/14  Revised for VBS emissions
c     01/08/16  Updated for Revised VBS emissions
c
c     Input arguments:
c        none
c
c     Output arguments:
c        none
c
c     Routines Called:
c        none
c
c     Called by:
c        CAMx
c
      include 'camx.prm'
      include 'flags.inc'
      include 'vbs.inc'
c
      character*4 ptspec(10)
      character*10 ptspcname
c
      real emispts(MXPTSRC)
c
      integer nptsrc_in
c
      data pi /3.1415927/
c
c-----Entry point
c
c----Initialize point emissions to zero
c
      ptemis = 0.
      do ifile=1,npoint_files
         if( is_netcdf_iptem(ifile) ) cycle
         kount = 1
 100     continue
         read(iptem(ifile),end=900) idat1,tim1,idat2,tim2
         if(NINT(1000*tim2) .EQ. 0) then
           tim2 = 24.
           idat2 = idat1
         endif
         ichktm1 = NINT( 1000*(tim1) )
         if( le1day ) then
            ichktm2 = NINT( 1000*(tim2) )
         else
            ichktm2 = NINT( 1000*(tim2)+24000*(idat2-idat1) )
         endif
         if( ichktm2 .EQ. 0 ) ichktm2 = 24000
         ichkems = NINT( 1000*(dtems/60.) )
         if( (ichktm2 - ichktm1) .NE. ichkems ) then
             write(iout,'(//,a)')'ERROR in READPT:'
             write(iout,*) 'Time interval in point emissions file does'
             write(iout,*)  ' not match emissions update time interval.'
             write(iout,*) '   Beginning Date/Time (Hour): ',idat1,tim1
             write(iout,*) '   Ending Date/Time    (Hour): ',idat2,tim2
             write(iout,*) '   Emiss Input interval(Hour): ',dtems/60.
             call camxerr()
         endif
         tim1 = 100.*aint(tim1) + 60.*amod(tim1,1.)
         tim2 = 100.*aint(tim2) + 60.*amod(tim2,1.)
         read(iptem(ifile)) idum,nptsrc_in
c
c  --- calculatate the total points and make sure it matches ----
c
         if( nptsrc_in .NE. nptsrc_files(ifile) ) then 
            write(iout,'(//,a)') 'ERROR in READPT:'
            write(iout,*)'Number of points read from point source files',
     &               ' differs from header value'
            write(iout,*)'Record: ',nptsrc_in,' Header: ',nptsrc_files(ifile) 
            call camxerr()
         endif 
c
         idx_start = idx_start_pts(ifile)
         read (iptem(ifile)) (idum,idum,idum,flowrat(n),effph(n),
     &                    n=idx_start+1,idx_start+nptsrc_in)
c
c----Verify that plume distribution override is valid
c
         do n=idx_start+1,idx_start+nptsrc_in
           if(flowrat(n) .LT. 0. .AND. effph(n) .LT. 0. ) then
              if( ABS(effph(n)) .LT. ABS(flowrat(n)) ) then
                 write(iout,'(//,a)')'ERROR in READPT:'
                 write(iout,*) 'Invalid values found for flow rate and plume height.'
                 write(iout,*) 'The plume distribution override is triggered. But the '
                 write(iout,*) 'value for flow rate (base) is larger than value for ',
     &                      'plume height (top).'
                 call camxerr()
              endif
           endif
         enddo
c
         do 50 ll = 1,nptspc(ifile) 
           read(iptem(ifile)) idum,(ptspec(i),i=1,10),(emispts(n),
     &                                 n=idx_start+1,idx_start+nptsrc_in)
           write(ptspcname,'(10A1)') ptspec
c   --- skip if not being modelled ---
c
           if(idxpnt_files(ifile,ll) .LE. 0 ) cycle
c
c----Load into the global array ---
c
           if(  lvbs .AND. LVBSPREPROC) then
             if( ptspcname .EQ. 'POA_OP    ' ) then
                 do l = 0, NVOLBIN
                    do ispc=1,nemspc
                       if( emspcname(ispc) .EQ. spname(kpap_c(l)) ) then
                          do n=idx_start+1,idx_start+nptsrc_in
                             ptemis(n,ispc) = emispts(n) * poa_op_ef(l)
                          enddo
                       endif 
                    enddo
                 enddo
                 goto 50
             endif
             if( ptspcname .EQ. 'POA_BB    ' ) then
                 do l = 0, NVOLBIN
                    do ispc=1,nemspc
                       if( emspcname(ispc) .EQ. spname(kpfp_c(l)) ) then
                          do n=idx_start+1,idx_start+nptsrc_in
                             ptemis(n,ispc) = emispts(n) * poa_bb_ef(l)
                          enddo
                       endif
                    enddo
                 enddo
                goto 50
             endif
           endif
           do n=idx_start+1,idx_start+nptsrc_in
             ptemis(n,idxpnt_files(ifile,ll)) = emispts(n)
           enddo
 50      continue
         write(iout,'(a40,2(f7.0,i8.5))')
     &       'Read point source file at ',tim1,idat1,tim2,idat2 
         call flush(iout)
c
c-----Check times only if LE1DAY = T, otherwise check both time and date
c
         if (le1day) then
           if (abs(tim1-time).lt.0.01 .and. tim2.gt.time) goto 200
           if (tim1-time.ge.0.01) goto 900
         else
           if ((idat1.lt.date .or.
     &      (idat1.eq.date .and. abs(tim1-time).lt.0.01)) .and.
     &      (idat2.gt.date .or.
     &      (idat2.eq.date .and. tim2.gt.time)))
     &    goto 200
         endif
         goto 100
c 
c-----Convert emission rates from moles/(dtems-hours) to moles/s for gas 
c     or g/(dtems-hours) to g/s for aero species 
c 
 200     do 10 l = 1,nemspc
           do n=idx_start+1,idx_start+nptsrc_in
             ptemis(n,l) = ptemis(n,l)/(60.*dtems) 
           enddo 
 10      continue 
c 
         goto 999
c
c-----End of file reached; if 1-day emissions requested, rewind and read 
c     through header once more.  Otherwise, report error and exit
c
 900     continue
         if (le1day) then
           if (kount.ge.2) then
             write(iout,'(//,a)')'ERROR in READPT:'
             write(iout,*)'Cannot match model time with point source time.'
             call camxerr()
           endif
           rewind(iptem(ifile))
           read(iptem(ifile)) idum 
           read(iptem(ifile)) dum  
           read(iptem(ifile)) idum  
           read(iptem(ifile)) idum 
           read(iptem(ifile)) idum 
           read(iptem(ifile)) dum
           kount = kount + 1
           goto 100
        else
           write(iout,'(//,a)')'ERROR in READPT:'
           write(iout,*)'End of point source file reached.  Make sure the '
           write(iout,*)
     &            'file is for the correct day and contains all hours.'
           call camxerr()
         endif
c
 999     continue
      enddo
c
      return
      end
