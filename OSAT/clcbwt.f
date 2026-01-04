c**** CLCBWT
c
      subroutine clcbwt(idate,btim,jdate,etim,ncolx,nrowy,nlays)
      use filunit
      use chmstry
      use bndary
      use camxcom
      use tracer
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routine calculates the weighted reactivity factor for VOC
c   species for the boundary conditions.  The mass is weighted by layer
c   thickness giving the weighted average for the cell.  The average
c   over all cells is then calculated.  The averages are calculated for
c   the entire bounday and for each boundary seperately.
c
c      Copyright 1996 - 2022
c     Ramboll
c
c     Argument declarations:
c        idate   I   beginning date of simulation (YYJJJ)
c        btim    R   beginning time of simulation
c        jdate   I   ending date of simulation (YYJJJ)
c        etim    R   ending time of simulation
c        ncolx    I   number of columns in coarse grid
c        nrowy    I   number of rows in coarse grid
c        nlays    I   number of layers in coarse grid
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     05/29/96   --gwilson--    Original development
c     11/06/01   --cemery--     Input dates are now Julian
c     11/4/09    -cemery-       Removed input top concentrations
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
      integer   idate
      real      btim
      integer   jdate
      real      etim
      integer   ncolx
      integer   nrowy
      integer   nlays
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      integer      ndate, jdlast, i, j, numbc_tracers, numic_tracers
      integer      idtnow, ioff
      real         ttime, ttlast, timnow
c
      real sumvoc(MXTRCLS)
      real sumkoh(MXTRCLS)
      real summir(MXTRCLS)
      real sumyld(MXTRCLS)
      real yldhvoc(MXTRCLS), yldlvoc(MXTRCLS)
      real consum(MXSPEC,0:IDXBTP)
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c   --- call routine to get the file pointers to the proper place ---
c
      if( is_netcdf_ibc ) then
          call ncf_bndprep(btim,idate,etim,jdate)
      else
          call bndprep(btim,idate,etim,jdate)
      endif
c
c   --- set ending date and time to be consistent with how time is 
c       counted here ----
c      
  222 continue
      ndate = idate
      ttime = btim/100.0
      jdlast = jdate
      ttlast = etim/100.0
      if( ttlast .EQ. 0. ) then
         jdlast = jdlast - 1
         ttlast = 24.0
      endif
      if( MOD(jdlast,1000) .EQ. 0 ) then
        jdlast = (INT(jdlast/1000)-1)*1000 + 365
      endif
c
c  --- initialize the dates and times ---
c
      idtnow = ndate
      timnow = ttime
c
      if( .NOT. is_netcdf_ibc ) then
         call rdsumbc(idtnow,timnow,jdlast,ttlast,ncolx,nrowy,nlays,consum)
      else
         call ncf_rdsumbc(idtnow,timnow,jdlast,ttlast,ncolx,nrowy,nlays,consum)
         bnddate_time_tflag = 0
         bnddate_time_etflag = 0
      endif
c
c   --- calculate the fractions ---
c
      if( lbndry .AND. .NOT. lsa_iorbc ) then
         do j=1,IDXBTP
            do icls=1,ntrcls
               sumvoc(icls) = 0.
               sumkoh(icls) = 0.
               summir(icls) = 0.
               sumyld(icls) = 0.
               yldhvoc(icls) = 0.
               yldlvoc(icls) = 0.
            enddo
            do i=1,nspec
               if( consum(i,j) .GT. 0. ) then
                   do icls=1,ntrcls
                      sumvoc(icls) = sumvoc(icls) + 
     &                                     consum(i,j) * trspmap(i,icls)
                      sumkoh(icls) = sumkoh(icls) + 
     &                        consum(i,j) * rkohrt(i)  * trspmap(i,icls)
                      summir(icls) = summir(icls) + 
     &                        consum(i,j) * rmirrt(i)  * trspmap(i,icls)
                      if( yhratmap(i,icls) .GT. 0. .OR.
     &                    ylratmap(i,icls) .GT. 0. ) then
                          sumyld(icls) = sumyld(icls) + consum(i,j)
                          yldhvoc(icls) = yldhvoc(icls) + 
     &                                     consum(i,j) * yhratmap(i,icls)
                          yldlvoc(icls) = yldlvoc(icls) + 
     &                                     consum(i,j) * ylratmap(i,icls)
                      endif
                   enddo
               endif
            enddo
            do icls=1,ntrcls
               if( sumvoc(icls) .GT. 0. ) then
                   wtkoh(iptcls(icls)+j) = sumkoh(icls) / sumvoc(icls)
                   wtmir(iptcls(icls)+j) = summir(icls) / sumvoc(icls)
               else
                   wtkoh(iptcls(icls)+j) = 0.
                   wtmir(iptcls(icls)+j) = 0.
               endif
               if( sumyld(icls) .GT. 0. ) then
                   yhrates(iptcls(icls)+j) = yldhvoc(icls) / sumyld(icls)
                   ylrates(iptcls(icls)+j) = yldlvoc(icls) / sumyld(icls)
               else
                   yhrates(iptcls(icls)+j) = 0.
                   ylrates(iptcls(icls)+j) = 0.
               endif
            enddo
         enddo
      else
         do icls=1,ntrcls
            sumvoc(icls) = 0.
            sumkoh(icls) = 0.
            summir(icls) = 0.
            sumyld(icls) = 0.
            yldhvoc(icls) = 0.
            yldlvoc(icls) = 0.
         enddo
         do i=1,nspec
            if( consum(i,0) .GT. 0. ) then
                do icls=1,ntrcls
                    sumvoc(icls) = sumvoc(icls) + 
     &                                     consum(i,0) * trspmap(i,icls)
                    sumkoh(icls) = sumkoh(icls) + 
     &                        consum(i,0) * rkohrt(i)  * trspmap(i,icls)
                    summir(icls) = summir(icls) + 
     &                        consum(i,0) * rmirrt(i)  * trspmap(i,icls)
                    if( yhratmap(i,icls) .GT. 0. .OR.
     &                  ylratmap(i,icls) .GT. 0. ) then
                       sumyld(icls) = sumyld(icls) + consum(i,0)
                       yldhvoc(icls) = yldhvoc(icls) + 
     &                                     consum(i,0) * yhratmap(i,icls)
                       yldlvoc(icls) = yldlvoc(icls) + 
     &                                     consum(i,0) * ylratmap(i,icls)
                    endif
                enddo
            endif
         enddo
         numbc_tracers = 1
         if( lsa_iorbc ) numbc_tracers = ncls_iorbc
         numic_tracers = 1
         if( lsa_ioric ) numic_tracers = ncls_ioric
         do icls=1,ntrcls
            do ioff=1,numbc_tracers
               if( sumvoc(icls) .GT. 0. ) then
                  wtkoh(iptcls(icls)+numic_tracers+ioff-1) = 
     &                                    sumkoh(icls) / sumvoc(icls)
                  wtmir(iptcls(icls)+numic_tracers+ioff-1) = 
     &                                    summir(icls) / sumvoc(icls)
               else
                  wtkoh(iptcls(icls)+numic_tracers+ioff-1) = 0.
                  wtmir(iptcls(icls)+numic_tracers+ioff-1) = 0.
               endif
               if( sumyld(icls) .GT. 0. ) then
                  yhrates(iptcls(icls)+numic_tracers+ioff-1) = 
     &                                    yldhvoc(icls) / sumyld(icls)
                  ylrates(iptcls(icls)+numic_tracers+ioff-1) = 
     &                                    yldlvoc(icls) / sumyld(icls)
               else
                  yhrates(iptcls(icls)+numic_tracers+ioff-1) = 0.
                  ylrates(iptcls(icls)+numic_tracers+ioff-1) = 0.
               endif
            enddo
            do ioff=1,ncls_iortc
               if( sumvoc(icls) .GT. 0. ) then
                  wtkoh(iptcls(icls)+numic_tracers+numbc_tracers+ioff-1) =
     &                                    sumkoh(icls) / sumvoc(icls)
                  wtmir(iptcls(icls)+numic_tracers+numbc_tracers+ioff-1) =
     &                                    summir(icls) / sumvoc(icls)
               else
                  wtkoh(iptcls(icls)+numic_tracers+numbc_tracers+ioff-1) = 0.
                  wtmir(iptcls(icls)+numic_tracers+numbc_tracers+ioff-1) = 0.
               endif
               if( sumyld(icls) .GT. 0. ) then
                  yhrates(iptcls(icls)+numic_tracers+numbc_tracers+ioff-1) =
     &                                    yldhvoc(icls) / sumyld(icls)
                  ylrates(iptcls(icls)+numic_tracers+numbc_tracers+ioff-1) =
     &                                    yldlvoc(icls) / sumyld(icls)
               else
                  yhrates(iptcls(icls)+numic_tracers+numbc_tracers+ioff-1) = 0.
                  ylrates(iptcls(icls)+numic_tracers+numbc_tracers+ioff-1) = 0.
               endif
            enddo
         enddo
      endif
c
c  --- return to the calling routine ---
c
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c    Format statements:
c-----------------------------------------------------------------------
c
 8000 format(A10,F10.0)
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
