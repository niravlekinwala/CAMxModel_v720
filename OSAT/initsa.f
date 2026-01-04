C**** INITSA
c
      subroutine initsa(version,nxx,nxy,idate,begtim,jdate,endtim)
      use grid
      use tracer
      use rtracchm
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routine calls all of the routines that set up the run for the
c   passive tracer algorithm.  The first thing done on routine OPENPT
c   is to set the flag for using passive tracer alogithm.  If this
c   flag is false coming out of the subroutine then all other routines
c   are skipped.
c
c      Copyright 1996 - 2022
c     Ramboll 
c
c     Argument description:
c      Inputs:
c         version C  version number of the model
c         nxx     I  number of cells in the X-direction
c         nxy     I  number of cells in the Y-direction
c         idate   I  beginning date of the simulation (YYJJJ)
c         begtim  R  beginning time of the simulation
c         jdate   I  ending date of the simulation (YYJJJ)
c         endtim  R  ending time of the simulation
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     05/27/96   --gwilson--    Original development
c     11/06/01   --cemery--     Input dates are now Julian
c     01/30/02   --gwilson--    Added code for RTRAC probing tool
c     07/11/07   --cemery--     Added code for RTCMC probing tool
c     07/16/07   --bkoo--       Added check for HDDM
c     03/01/16   --gwilson--    Added partial source area map
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
c
c-----------------------------------------------------------------------
c    Arguement declarations:
c-----------------------------------------------------------------------
c
      character*20 version
      integer      nxx(*)
      integer      nxy(*)
      integer      idate
      integer      jdate
      real         begtim
      real         endtim
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      integer ip, ic, igrd, igrp, jfin, ifin, ipart
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c  ---- set the global variables to the passed arguments ---
c
      do igrd=1,ngrid
        nxcell(igrd) = nxx(igrd)
        nycell(igrd) = nxy(igrd)
      enddo 
c
c  --- if doing RTRAC, call routine to read chemistry definitions file ---
c 
      if( tectyp .EQ. RTRAC ) then
         call rdchmrt(version)
         ntrtim = 0 
         nsaspc = ntotsp
         ntrcls = 1
         iptcls(1) = 1
         nptcls(1) = ntotsp
         ipttim = ntotsp + 1
c
c  --- if doing RTCMC, call routine to read chemistry definitions file
c      and generate the chemistry solver ---
c
      elseif( tectyp .EQ. RTCMC ) then
         call rdschm(ldbg)
         call gencmc(ldbg)
         ntrtim = 0
         nsaspc = ntotsp
         ntrcls = 1
         iptcls(1) = 1
         nptcls(1) = ntotsp
         ipttim = ntotsp + 1
         call alloc_tracer_vdep(ngrid,ncol,nrow,ntotsp)
c
c  ---- call routine to read the source mapping file ---
c
      else
          if( nregin .GT. 0 ) then
              call resmap()
              if( lpartial ) then
                 call rdpartial_map()
                 call chkpartial()
              endif
c
c  --- now assign any values that were not read in ---
c
              do ip = 1,ngrid
                do 20 ic = 1,nchdrn(ip)
                  igrd = idchdrn(ic,ip)
                  do igrp=0,ngroup
                     if( .NOT. lmapfl(0,igrd) ) then
                        do jfin = 1,nrow(igrd)
                          j = (jfin - 2)/nmesh(igrd) + j1(igrd)
                          do ifin = 1,ncol(igrd)
                             i = (ifin - 2)/nmesh(igrd) + i1(igrd)
                             do ipart=1,npartial(igrp,igrd)
                               igrmap(igrp,ipart,igrd,ifin,jfin) =
     &                                          igrmap(igrp,ipart,ip,i,j)
                               frcmap(igrp,ipart,igrd,ifin,jfin) =
     &                                           frcmap(igrp,ipart,ip,i,j)
                             enddo
                          enddo
                        enddo
                     endif
                  enddo
   20           continue
              enddo
          endif
c
c   --- read the header of the initial conditions file to get init tracer species ---
c
          if( lsa_ioric ) then
             if( is_netcdf_ioric ) then
                call ncf_sa_icprep(begtim,idate,begtim,idate)
             else
                call sa_icprep(begtim,idate,endtim,jdate)
             endif
          endif
c
c   --- read the header of the initial conditions file to get init tracer species ---
c
          if( lsa_iorbc ) then
             if( is_netcdf_iorbc ) then
                call ncf_sa_bcprep(begtim,idate,endtim,jdate)
             else
                call sa_bcprep(begtim,idate,endtim,jdate)
             endif
             if( lsa_iortc ) then
                 if( is_netcdf_iorbc ) then
                     call sa_tcprep(begtim,idate,endtim,jdate)
                 else
                     call ncf_sa_tcprep(begtim,idate,endtim,jdate)
                 endif
             endif
          endif
c
c   ---- call routine to set the species flags ----
c
          call stabsa(idate,begtim,jdate,endtim)
c
c   ---- call routine to set up the tracer species lists ----
c
          if( ltrace ) then
            call specsa(idate,begtim,jdate,endtim)
          else if( lddm .OR. lhddm ) then
            call specddm( )
          endif
      endif
c
c  ---- return to the calling routine ---
c
      goto 9999
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
