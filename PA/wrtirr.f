      subroutine wrtirr(iendat,endtim)
      use filunit
      use grid
      use camxfld
      use chmstry
      use procan
      use tracer
c
c----CAMx v7.20 220430
c
c     This routine writes to the output file for the Integrated Reaction
c     Rates (IRR) data for the Process Analysis algorithm.  Each record 
c     contains the all of the data for once cell and species for the 
c     specified time period.
c     This routine also calls the subroutines that will write the 
c     gridded Chemical Process Analysis data.
c
c      Copyright 1996 - 2022
c     Ramboll
c
c     Modifications:
c        8/25/06    CPA output files now all UAM format, one file per grid
c
c     Input arguments:
c        endtim     ending time for this period
c        iendat     ending date for this period
c
c     Subroutines Called:
c       WRTCGCPA
c       WRTFGCPA
c
c
c     Called by:
c        CAMx
c
      include 'camx.prm'
      include 'flags.inc'
      include 'ncf_iodat.inc'
c
c-----Argument declarations
c
      integer iendat
      real    endtim 
c
c-----Local variables
c
      character*200 action
      integer       num_dims
      logical       loutspec(MXSPEC+MXTRSP)
c
c-----Entry point
c
c  --- loop over the number of sub-domain cells ---
c
      do icel=1,npa_cels
c
c  --- write out the record for this cell ---
c
           write(irr_unit,ERR=7000) iendat, endtim, 
     &              ipadom(icel), ipanst(icel),
     &                  ipax(icel), ipay(icel), ipaz(icel),
     &                                  (cirr(icel,i),i=1,nirrrxn)
c
c  --- next subdomain cell ---
c
      enddo
c
c  --- call routine to write the gridded CPA arrays ---
c
      if( lsfcfl ) then
         do igrd = 1,ngrid        
            if( .NOT. lcdfout ) then
                call wrtcpa(igrd,iendat,endtim,igrd,ncol(igrd),nrow(igrd),
     &                           nlay(igrd),ntotsp,ptconc(ipsa3d(igrd)))
            else
                nlayav = nlay(igrd)
                num_dims = 4
                if( .NOT. l3davg(igrd) ) then
                    nlayav = 1
                    num_dims = 3
                endif
                write(action,'(A,I2)') 'Writing CPA file for grid: ',igrd
                call ncf_wrt_data_tstep(action,iowsfc(igrd),ntotsp)
                loutspec = .TRUE.
                call ncf_wrt_data_species(action,iowsfc(igrd),igrd,num_dims,
     &               ncol(igrd),nrow(igrd),nlay(igrd),nlay(igrd),nlayav,
     &                 ntotsp,ptname,loutspec,ptconc(ipsa3d(igrd)),
     &                            height(iptr3d(igrd)))
            endif
         enddo
      endif
      goto 9999
c
c----Error messages
c
 7000 continue
      write(iout,'(//,a)') 'ERROR in WRTIRR:'
      write(iout,'(1X,2A)',ERR=9999) 'Writing data to the output ',
     &                         'Integrated Reaction Rates (.irr) file.'
      write(iout,'(10X,A,I8.5,5X,A,F8.1)') 
     &      'Date: ',iendat,'Time: ',endtim
      call camxerr()
c
 9999 continue
      return
      end
