***** EMPREPSA
c
      subroutine emprepsa(idate,btim,igrid,numpts_lst)
      use filunit
      use grid
      use chmstry
      use ptemiss
      use tracer
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routine reads the source group emissions files and moves the
c   file pointer to the first hour to be used in this run.  Each point
c   source and surace emissions files respresenting a source group is
c   read.
c
c      Copyright 1996 - 2022
c     Ramboll
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'flags.inc'
c
c-----------------------------------------------------------------------
c    Arguement declarations:
c-----------------------------------------------------------------------
c
      integer   idate
      real      btim
      integer   igrid
      integer   numpts_lst
c
c-----------------------------------------------------------------------
c    Local parameters:
c-----------------------------------------------------------------------
c
      integer MXOUT
      parameter (MXOUT = MAX(MXSPEC,MXTRSP))
c
c-----------------------------------------------------------------------
c    External functions:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      character*200 fname
      character*10  namein, cfile
      character*4   ifname(10), ifnote(60)
      integer       ifile, nseg, nspect, igroup
      integer       ispc, ibgdat, iendat
      integer       idateb, idatee, nxcelt, nycelt, nzceltm npts_in_list
      integer       idum, izonet, iounit, isegm, npoint, idxpt_base
      integer       i, j, ndate, seg4(4), idxfile, num_emsfiles, num_ptsfiles
      integer       ii, jj, igrd, buffer_offset
      real          tutmx, tutmy, tdxcel, tdycel
      real          tbegti, tendti, begtim, endtim, dum
      real          txorig, tyorig, ttime, tdummy
c
      character*4 iname(10,MXOUT)
      real        txloc (MXPTSRC)
      real        tyloc (MXPTSRC)
      real        tstkhgt(MXPTSRC)
      real        tstkdiam(MXPTSRC)
      real        tstktemp(MXPTSRC)
      real        tstkvelo(MXPTSRC)
c
      common /comemprepsa/ txloc, tyloc, tstkhgt, tstkdiam, 
     &                     tstktemp, tstkvelo
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c   --- allocate array for checking list ---
c
      if( lptsrc ) then
         xlocpt = xloc
         ylocpt = yloc
         idx_last_found = 1
      endif
c
c   --- set the date and time ---
c
      ndate = idate
      ttime = btim/100
      if( ttime .EQ. 24.0 ) then
         ndate = ndate + 1
         ttime = 0.
         if( MOD(ndate,1000) .GT. 365 ) then
            if( MOD(INT(ndate/1000),4) .EQ. 0 ) then
               if( MOD(ndate,1000) .EQ. 367 )
     &                        ndate = (INT(ndate/1000)+1)*1000 + 1
            else
               ndate = (INT(ndate/1000)+1)*1000 + 1
            endif
         endif
      endif
c
c   --- initialize postion zero, which is for the model emissions ----
c
      if( ltrace .AND. (tectyp .NE. RTRAC .AND. tectyp .NE. RTCMC) ) then
         do idxfile=1,nemiss_files(igrid)
            if( larsrc .AND. ltrace ) then
                ltemfl(igrid,0,idxfile) = .TRUE.
            else
                ltemfl(igrid,0,idxfile) = .FALSE.
            endif
         enddo
      endif
c
c  --- calculate number of sources in base inventory ---
c
      if( tectyp .NE. SA ) then
         numpts_lst = 0
         do i=1,npoint_files
            numpts_lst = numpts_lst + nptsrc_files(i)
         enddo
      endif
c
c   --- loop over the source grouping files ----
c
      igroup_start = 1
      if( ltrace .AND. (tectyp .NE. RTRAC .AND. tectyp .NE. RTCMC) ) 
     &                                               igroup_start = 0
      do 10 igroup=igroup_start,ngroup
           num_emsfiles = num_iortem(igrid,igroup)
           if( igroup .EQ. 0 ) num_emsfiles = nemiss_files(igrid)
           do 11 idxfile=1,num_emsfiles
c
c   --- only read if surface emissions file is supplied ---
c
             if( larsrc ) then
                if( igroup .EQ. 0 ) then
                   if( is_netcdf_iarem(igrid,idxfile)) cycle
                   iounit = iarem(igrid,idxfile)
                   write(fname,'(A,I3)') 'EMISSIONS -- UNIT ',
     &                                            iarem(igrid,idxfile)
                   buffer_offset = buffer_offset_iarem(igrid,idxfile)
                else
                   if( is_netcdf_iortem(igrid,igroup,idxfile)) cycle
                   iounit = iortem(igrid,igroup,idxfile)
                   fname = temfil(igrid,igroup,idxfile)
                endif
c
c   --- rewind the file ----
c
                rewind(iounit,ERR=7005)
c
c   --- read the header of the emissions file ---
c
                read(iounit,ERR=7000) ifname, ifnote, nseg, nspect, 
     &                                 idateb, tbegti, idatee, tendti
                if( nspect .GT. MXTRSP ) then
                  write(iout,'(//,A)') 'ERROR in EMPREPSA:'
                  write(iout,*) 'A parameter in the camx.prm is not ',
     &                                        'sufficiently large.'
                  write(iout,*) 
     &                    'Please change the value for parameter: MXTRSP'
                  write(iout,*) 'It should be set to a value of at least: ',
     &                                                             nspect
                  call flush(iout)
                  call camxerr()
                endif
c
                read(iounit,ERR=7000) tutmx, tutmy, izonet, txorig, 
     &                 tyorig, tdxcel, tdycel, nxcelt, nycelt, nzcelt
                read(iounit,ERR=7000) seg4
                read(iounit,ERR=7000) 
     &                   ((iname(i,ispc),i=1,10),ispc=1,nspect)
                if( igroup .GT. 0 ) then
                   buffer_offset_iortem(igrid,igroup,idxfile) = 0
                   if( nxcelt .EQ. ncol(igrid)-2 .AND. nycelt .EQ. nrow(igrid)-2 )
     &                                buffer_offset_iortem(igrid,igroup,idxfile) = 1
                   buffer_offset = buffer_offset_iortem(igrid,igroup,idxfile)
                endif
c
c   --- set up index for species names ---
c
                nspcem(igrid,igroup,idxfile) = nspect
                do 20 ispc=1,nspcem(igrid,igroup,idxfile)
                   write(namein,'(10A1)') (iname(i,ispc),i=1,10)
                   if(namein .EQ.'HNO2      ') namein = 'HONO      '
                   if(namein .EQ.'HCHO      ' .AND. kHCHO .eq. nspec+1 )
     &                                            namein = 'FORM      '
                   do j=1,nspec
                     if(namein .EQ. spname(j) ) then
                         idxems(igrid,igroup,idxfile,ispc) = j
                         goto 20
                      endif
                   enddo
  20            continue
c
c   --- read the data until the correct position in the file is reached ---
c
  111           continue
                read(iounit,ERR=7002,END=7003) ibgdat, begtim, 
     &                                              iendat, endtim
                if(NINT(endtim) .EQ. 0) then
                  endtim = 24.
                  iendat = ibgdat
                endif
                if( le1day ) then
                   if( INT(begtim*100) .EQ. INT(ttime*100) ) goto 222
                else
                   if( ibgdat .EQ. ndate .AND. 
     &                   INT(begtim*100) .EQ. INT(ttime*100) ) goto 222
                endif
                do ispc=1,nspcem(igrid,igroup,idxfile)
                    read(iounit,ERR=7002,END=7003) idum, 
     &                 (dum,i=1,10),((dum,i=1+buffer_offset,ncol(igrid)-buffer_offset),
     &                                     j=1+buffer_offset,nrow(igrid)-buffer_offset)
                enddo
                goto 111
  222           continue
                backspace(iounit)
             endif
   11     continue
c
c   --- point sources are only for the coarse grid -- grid #1 ---
c
          if( igrid .GT. 1 ) cycle
c
c   --- don't have group 0 with point sources ---
c
          if( igroup .EQ. 0 ) cycle
c
c   --- only read if elevated point source emissions file is supplied ---
c
          do idxfile=1,num_iortpt(igroup)
             if( ltptfl(igroup,idxfile) .AND. lptsrc ) then
                if( is_netcdf_iortpt(igroup,idxfile) ) cycle
                iounit = iortpt(igroup,idxfile)
                fname = tptfil(igroup,idxfile)
c
c   --- rewind the file ----
c
                rewind(iounit,ERR=7005)
c
c   --- read the header of the emissions file ---
c
                read(iounit,ERR=7000) ifname, ifnote, nseg, nspect, 
     &                                 idateb, tbegti, idatee, tendti
                if( nspect .GT. MXTRSP ) goto 7001
c
c   --- check to see if this is a compact file ---
c
                write(cfile,'(10A1)') ifname
                lcompactpt(igroup,idxfile) = .FALSE.
                if( cfile .EQ. 'PTSOURCECP' ) lcompactpt(igroup,idxfile) = .TRUE.
c
                read(iounit,ERR=7000) tutmx, tutmy, izonet, txorig, 
     &                 tyorig, tdxcel, tdycel, nxcelt, nycelt, nzcelt
                read(iounit,ERR=7000) seg4
                read(iounit,ERR=7000) 
     &                   ((iname(i,ispc),i=1,10),ispc=1,nspect)
c
c   --- set up index for species names ---
c
                nspcpt(igroup,idxfile) = nspect
                do 30 ispc=1,nspcpt(igroup,idxfile)
                   write(namein,'(10A1)') (iname(i,ispc),i=1,10)
                   if(namein .EQ.'HNO2      ') namein = 'HONO      '
                   if(namein .EQ.'HCHO      ' .AND. kHCHO .eq. nspec+1 )
     &                                            namein = 'FORM      '
                   do j=1,nspec
                     if(namein .EQ. spname(j) ) then
                         idxpts(igroup,idxfile,ispc) = j
                         goto 30
                      endif
                    enddo
  30            continue
                read(iounit,ERR=7000) isegm, npoint
                if( npoint .GT. MXPTSRC ) goto 7006
c
                read(iounit,ERR=7002) (txloc(i),tyloc(i),tstkhgt(i),tstkdiam(i),
     &                                     tstktemp(i),tstkvelo(i),i=1,npoint)
                if( tectyp .NE. SA ) then
                   do i=1,npoint
                      tstkvelo(i) = tstkvelo(i)/3600.
c
c  --- make sure this source is in the base inventory and flag it ---
c
                      numpts_lst = numpts_lst + 1
                      if( llatlon ) then
                          xstk(numpts_lst,1) = txloc(i) - xorg
                          ystk(numpts_lst,1) = tyloc(i) - yorg
                      else
                          xstk(numpts_lst,1) = txloc(i)/1000. - xorg
                          ystk(numpts_lst,1) = tyloc(i)/1000. - yorg
                      endif
                      do igrd = 2,ngrid
                         xstk(numpts_lst,igrd) = 
     &                      xstk(numpts_lst,1) - (inst1(igrd)-1)*delx
                         ystk(numpts_lst,igrd) = 
     &                      ystk(numpts_lst,1) - (jnst1(igrd)-1)*dely
                         ii = 2 + FLOOR(xstk(numpts_lst,igrd)/
     &                                 delx*FLOAT( meshold(igrd) ) )
                         jj = 2 + FLOOR(ystk(numpts_lst,igrd)/
     &                                dely*FLOAT( meshold(igrd) ) )
                         if (ii .GT. 1 .AND. ii .LT. ncol(igrd) .AND.
     &                        jj .GT. 1 .AND. jj .LT. nrow(igrd)) then
                            nosrc(igrd) = nosrc(igrd) + 1
                            idsrc(nosrc(igrd),igrd) = numpts_lst
                            isrc(nosrc(igrd),igrd) = ii
                            jsrc(nosrc(igrd),igrd) = jj
                         endif
                     enddo
                     xlocpt(numpts_lst) = txloc(i)
                     ylocpt(numpts_lst) = tyloc(i)
                     xloc(numpts_lst) = txloc(i)
                     yloc(numpts_lst) = tyloc(i)
                     hstk(numpts_lst) = tstkhgt(i)
                     dstk(numpts_lst) = tstkdiam(i)
                     tstk(numpts_lst) = tstktemp(i)
                     vstk(numpts_lst) = tstkvelo(i)
                     lpiglet(numpts_lst) = .FALSE.
                     if( dstk(numpts_lst) .LT. 0. ) lpiglet(numpts_lst) = .TRUE.
                     idx_point_in_list(igroup,idxfile,i) = numpts_lst
                   enddo
                endif
c
c   --- read the data until the correct position in the file is reached ---
c
  333           continue
                read(iounit,ERR=7002,END=7003) ibgdat, begtim, 
     &                                              iendat, endtim
                if(NINT(endtim) .EQ. 0) then
                  endtim = 24.
                  iendat = ibgdat
                endif
                if( le1day ) then
                   if( INT(begtim*100) .EQ. INT(ttime*100) ) goto 444
                else
                   if( ibgdat .EQ. ndate .AND. 
     &                INT(begtim*100) .EQ. INT(ttime*100) ) goto 444
                endif
                read(iounit,ERR=7002) isegm, npoint
                read(iounit,ERR=7002) (idum,idum,idum,dum,dum,i=1,npoint)
                do ispc=1,nspcpt(igroup,idxfile)
                    read(iounit,ERR=7002,END=7003) idum,(dum,i=1,10),
     &                                                (dum,i=1,npoint)
                enddo
                goto 333
  444           continue
                backspace(iounit)
c
             endif
          enddo
   10 continue
c
c  ---- return to calling routine ---
c
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,A)') 'ERROR in EMPREPSA:'
      write(iout,'(1X,3A)') 'Reading header of ',
     &                     'emissions file: ',fname(:istrln(fname))
      call camxerr()
c
 7001 continue
      write(iout,'(//,A)') 'ERROR in EMPREPSA:'
      write(iout,*) 'A parameter in the camx.prm is not ',
     &                                        'sufficiently large.'
      write(iout,*)
     &                    'Please change the value for parameter: MXTRSP'
      write(iout,*) 'It should be set to a value of at least: ',
     &                                                             nspect
      call flush(iout)
      call camxerr()
c
 7002 continue
      write(iout,'(//,A)') 'ERROR in EMPREPSA:'
      write(iout,'(1X,2A)') 'Reading emissions grid in file: ',
     &                                        fname(:istrln(fname))
      call camxerr()
c
 7003 continue
      write(iout,'(//,A)') 'ERROR in EMPREPSA:'
      write(iout,'(1X,3A)') 'Premature end-of-file reading ',
     &                     'emissions file: ',fname(:istrln(fname))
      write(iout,'(1X,A)') 'Make sure date/time matches the simulation '
      write(iout,'(1X,A)') 'or turn on single emissions flag.'
      call camxerr()
c
 7004 continue
      write(iout,'(//,A)') 'ERROR in EMPREPSA:'
      if( lddm .OR. lhddm ) then
         write(iout,'(1X,2A,/,A)') 'Point source found in DDM ',
     &               'point source file is not found ',' in the base inventory.'
      else
         write(iout,'(1X,2A,/,A)') 'Point source found in Source Apportionment ',
     &               'point source file is not found ',' in the base inventory.'
      endif
      write(iout,'(1X,2A)') '   File: ',fname(:istrln(fname))
      write(iout,'(1X,A,I10)') '   Point ID: ',i
      call camxerr()
c
 7005 continue
      if( igroup .EQ. 0 ) then
         if( iarem(igrid,1) .GT. 0 ) then
            write(iout,'(//,A)') 'ERROR in EMPREPSA:'
            write(iout,'(1X,2A)') 'Cannot access emissions file',
     &                             ' provided for OSAT processing. '
            write(iout,'(1X,2A)') 'If it is the same as the file ',
     &            'used for the regular model, you may have to make a '
            write(iout,'(1X,A)')'copy and specify the name of the copy.'
         else
            write(iout,'(//,A)') 'ERROR in EMPREPSA:'
            write(iout,'(1X,2A)') 'Flexi-nesting the emissions inputs',
     &              ' is not supported with a '
            write(iout,'(1X,2A)') 'Probing Tools application.'
            write(iout,'(1X,2A)') 'You must supply an emissions file ',
     &                                                 'for each grid.'
         endif
      else
         write(iout,'(//,A)') 'ERROR in EMPREPSA:'
         write(iout,'(1X,3A)') 'Cannot access emissions file',
     &     ' provided for OSAT processing: ',fname(:istrln(fname))
      endif
      call camxerr()
c
 7006 continue
      write(iout,'(//,A)') 'ERROR in EMPREPSA:'
      write(iout,*) 'A parameter in the camx.prm is not ',
     &                                        'sufficiently large.'
      write(iout,*) 'Please change the value for parameter: MXPTSRC'
      write(iout,*) 'It should be set to a value of at least: ',npoint
      call flush(iout)
      call camxerr()
c
c-----------------------------------------------------------------------
c    Format statements:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
