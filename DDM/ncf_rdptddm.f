c***** NCF_RDPTDDM
c
      subroutine ncf_rdptddm(ndate,ttime)
      use chmstry
      use grid
      use filunit
      use ptemiss
      use tracer
      implicit none
c
c----CAMx v7.20 220430
c
c
c      Copyright 1996 - 2022
c     Ramboll
c
c-----------------------------------------------------------------------
c
c   This routine reads one hour of point source emissions for the DDM 
c   process and fills the approproate arrays.  This routine is for NETCDF 
c   files. 
c    Argument descriptions:
c     Outputs:
c     Inputs:
c       ndate     I    julian day of current hour
c       ttime     R    current hour
c       
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     09/26/19  Original development
c
c-----------------------------------------------------------------------
c   Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'flags.inc'
      include 'camx.inc'
      include 'netcdf.inc'
c       
c-----------------------------------------------------------------------
c   Argument declarations:
c-----------------------------------------------------------------------
c
      integer ndate
      real    ttime
c
c-----------------------------------------------------------------------
c   External functions:
c-----------------------------------------------------------------------
c
      integer ncf_get_tstep
      integer istrln
c
c-----------------------------------------------------------------------
c   Local variables:
c-----------------------------------------------------------------------
c
      character*200 fname, action
      character*20  this_var
      character*10  cname
      integer       nedge, igroup, idxfile, iounit, this_tstep
      integer       ierr, idxpt, this_time_tflag, this_time_etflag
      integer       imod, idxbase, icel, jcel, igrd, ig, ii, jj
      integer       imap, iptr, this_varid,  data_start(2), data_count(2)
      integer       iddm, ispc
      real          xloccrs, yloccrs, xloctmp, yloctmp, emstmp, pi
      logical       luse
c
      real,    allocatable, dimension(:) ::  emispts
      integer, allocatable, dimension(:) ::  sa_region
      real,    allocatable, dimension(:) ::  effph_in
      real,    allocatable, dimension(:) ::  plume_bot
      real,    allocatable, dimension(:) ::  plume_top
c
      data pi /3.1415927/
c
c-----------------------------------------------------------------------
c   Entry point:
c-----------------------------------------------------------------------
c
c  --- set the number of BC edges --
c
      if( lbndry ) then
        nedge = 5
      else
        nedge = 1
      endif
c
c   --- loop over all of the groups ----
c
      do 10 igroup=1,ngroup
        do 11 idxfile=1,num_iortpt(igroup)
c
c   --- skip if filename not supplied ---
c
          if( .NOT. is_netcdf_iortpt(igroup,idxfile) .OR. .NOT. 
     &            ltptfl(igroup,idxfile) .OR. .NOT. lptsrc ) goto 11
          write(action,'(2A,I3,A,I3)') 'Reading the DDM NetCDF point source',
     &            ' emissions file for Group: ',igroup,' File: ',idxfile
c
c   --- set the unit number for surface emissions file ---
c
          iounit = iortpt(igroup,idxfile)
          fname = tptfil(igroup,idxfile)
c
c  ---- allocate the temporary array ---
c
          allocate( emispts(nptsrc_safile(igroup,idxfile)) )
          allocate( sa_region(nptsrc_safile(igroup,idxfile)) )
          allocate( effph_in(nptsrc_safile(igroup,idxfile)) )
          allocate( plume_bot(nptsrc_safile(igroup,idxfile)) )
          allocate( plume_top(nptsrc_safile(igroup,idxfile)) )
          plume_bot = 0.
          plume_top = 0.
c
c  ---- get the region override flag ---
c
          sa_region = 0
          this_var = 'saoverride'
          ierr = nf_inq_varid(iounit, this_var, this_varid)
          if( ierr .NE. NF_NOERR ) cycle
          ierr = nf_get_var_int(iounit,this_varid,sa_region)
          if( ierr .NE. NF_NOERR ) goto 7001
c
c  ---- get the index for timestep containing this time ----
c
          this_tstep = ncf_get_tstep(iounit,action,ndate,ttime,
     &               this_time_tflag,this_time_etflag,le1day,.TRUE.)
c
c   --- set the indexes for what to read ---
c
          data_start(1) = 1
          data_count(1) = nptsrc_safile(igroup,idxfile)
          data_start(2) = this_tstep
          data_count(2) = 1
c
c  ---- plumerise ---
c
          this_var = 'plumerise'
          ierr = nf_inq_varid(iounit, this_var, this_varid)
          if( ierr .EQ. NF_NOERR) then
             ierr = nf_get_vara_real(iounit,this_varid,data_start,
     &                                  data_count,effph_in)
             if( ierr .NE. NF_NOERR) goto 7001
          endif
c
c  ---- load plume bottom, if available ---
c
          this_var = 'plume_bottom'
          ierr = nf_inq_varid(iounit, this_var, this_varid)
          if( ierr .EQ. NF_NOERR) then
             ierr = nf_get_vara_real(iounit,this_varid,data_start,
     &                                  data_count,plume_bot)
             if( ierr .NE. NF_NOERR) goto 7001
          endif
c
c  ---- load plume top, if available ---
c
          this_var = 'plume_top'
          ierr = nf_inq_varid(iounit, this_var, this_varid)
          if( ierr .EQ. NF_NOERR) then
             ierr = nf_get_vara_real(iounit,this_varid,data_start,
     &                                  data_count,plume_top)
             if( ierr .NE. NF_NOERR) goto 7001
          endif
c
c   --- check integrity of data ---
c
          do idxpt=1,nptsrc_safile(igroup,idxfile)
             if( plume_bot(idxpt) .GT. 0 .AND. 
     &                    plume_top(idxpt) .LE. 0. ) goto 7002
             if( plume_bot(idxpt) .LE. 0 .AND. 
     &                    plume_top(idxpt) .GT. 0. ) goto 7002
             if( plume_bot(idxpt) .GT. 0. .AND.
     &           plume_bot(idxpt) .GT. plume_top(idxpt) ) goto 7002
          enddo
c          
c   --- read the emissions for this hour ---
c
          do ispc=1,nspec
c
c  --- see if this species is used by DDM ---
c
             luse = .FALSE.
             cname = spname(ispc)
             do iddm=1,nemddm
                 luse = .FALSE.
                 if( emddmsp(iddm) .EQ. cname ) luse = .TRUE.
                 if( emddmsp(iddm) .EQ. NAMALL  ) luse = .TRUE.
                 if( emddmsp(iddm) .EQ. NAMVOC
     &                          .AND. lvocsp(ispc) ) luse = .TRUE.
                 if( emddmsp(iddm) .EQ. NAMNOX
     &                          .AND. lnoxsp(ispc) ) luse = .TRUE.
                 if( emddmsp(iddm) .EQ. NAMHRV
     &                                 .AND. lhrvoc(ispc) ) luse = .TRUE.
c
c  --- if this DDM species matches this modeled species, read and load it ---
c
                 if( .NOT. luse ) cycle
                 this_var = cname                 
                 ierr = nf_inq_varid(iounit, this_var, this_varid)
                 if( ierr .NE. NF_NOERR) cycle
                 ierr = nf_get_vara_real(iounit,this_varid,data_start,
     &                                               data_count,emispts)
                 if( ierr .NE. NF_NOERR) goto 7000
                 do 30 idxpt=1,nptsrc_safile(igroup,idxfile)
                    idxbase = idx_point_in_list(igroup,idxfile,idxpt)
                    if( effph_in(idxpt) .GT. 0. )
     &                           effph(idxbase) = -effph_in(idxpt)
                    if( plume_bot(idxpt) .GT. 0. )
     &                           effph(idxbase) = -plume_bot(idxpt)
                    if( plume_top(idxpt) .GT. 0. )
     &                           flowrat(idxbase) = -plume_top(idxpt)
                    if( llatlon ) then
                       xloccrs = xlocpt(idxbase) - xorg
                       yloccrs = ylocpt(idxbase) - yorg
                    else
                       xloccrs = xlocpt(idxbase)/1000. - xorg
                       yloccrs = ylocpt(idxbase)/1000. - yorg
                    endif
                    icel = 1 + FLOOR( xloccrs/delx )
                    jcel = 1 + FLOOR( yloccrs/dely )
                    if( icel .LT. 1 .OR. icel .GT. ncol(1) ) goto 30
                    if( jcel .LT. 1 .OR. jcel .GT. nrow(1) ) goto 30
c
c   --- find out if a nest contains this source  ---
c
                    igrd = 1
                    do ig = 2,ngrid
                       xloctmp = xloccrs - (inst1(ig)-1)*delx
                       yloctmp = yloccrs - (jnst1(ig)-1)*dely
                       if( xloctmp .GT. 0. .AND. 
     &                                    yloctmp .GT. 0. ) then
                            ii = 2 + FLOOR( xloctmp/delx *
     &                                       FLOAT( meshold(ig) ) )
                            jj = 2 + FLOOR( yloctmp/dely *
     &                                              FLOAT( meshold(ig) ) )
                            if( ii .GT. 1 .AND. jj .GT. 1 .AND. 
     &                                     ii .LT. ncol(ig) .AND. 
     &                                             jj .LT. nrow(ig) ) then
                               igrd = ig
                               icel = ii
                               jcel = jj
                            endif
                        endif
                    enddo
c
c  --- get the region for this cell from mapping array ----
c
                    imap = igrmap(0,1,igrd,icel,jcel)
                    if( sa_region(idxpt) .GT. 0 ) then
                       if( .NOT. lptoverride ) goto 7001
                       imap = ABS( sa_region(idxpt) )
                    endif
                    if( imap .LE. 0 .OR. imap .GT. nregin ) goto 30
c
c  --- convert to emissions time and put into array ---
c
                    emstmp = emispts(idxpt)/(60.*dtems)
                    iptr = iptddm(ispc) + nicddm + 
     &                    nbcddm*nedge + (iddm-1)*ngroup*nregin +
     &                                        (imap-1)*ngroup + igroup - 1
                       if( emstmp .GT. bdnl(ispc) )
     &                                sapnts(idxbase,iptr) = 
     &                                     sapnts(idxbase,iptr) + emstmp
c
c  --- next point ---
c
   30            continue
c
c  --- next DDM species ---
c
              enddo
c
c  --- next model species ---
c
          enddo
c
c  --- deallocate local arrays --
c
          deallocate( emispts )
          deallocate( sa_region )
          deallocate( effph_in )
          deallocate( plume_bot )
          deallocate( plume_top )
c
c  --- get the next file ---
c
  11    continue
  10  continue
c
      goto 9999
c
c----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,a)') 'ERROR in NCF_RDPTDDM: '
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot find variable in file: ',
     &                                      this_var(:istrln(this_var))
      write(iout,'(2A)') 'File: ',fname(:istrln(fname))
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7001 continue
      write(iout,'(//,a)') 'ERROR in NCF_RDPTDDM: '
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot read data for variable: ',
     &                                      this_var(:istrln(this_var))
      write(iout,'(2A)') 'File: ',fname(:istrln(fname))
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7002 continue
      write(iout,'(//,a)') 'ERROR in NCF_RDPTDDM::'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Point source found with invalid plume ',
     &                                        'distribution values.'
      write(iout,'(A,I10)') 'Point source number in file: ',idxpt
      write(iout,'(A,F10.4)') 'Plume bottom: ',plume_bot(idxpt)
      write(iout,'(A,F10.4)') 'Plume top   : ',plume_top(idxpt)
      call camxerr()
c
c----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
