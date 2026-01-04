C**** NCF_READPT
c
      subroutine ncf_readpt(ifile,iunit,numpts,numspcs,pntems)
      use chmstry
      use camxcom
      use filunit
      use ptemiss
      implicit none
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c     NCF_READPT reads the NetCDF pint source emissions file 
c     loads the data into global arrays for emissions injection.
c
c      Copyright 1996 - 2022
c     Ramboll
c
c     Input arguments:
c        ifile               index of file in list for this grid
c        iunit               file unit for area emissions file
c        numpnts             total number of point sources
c        numspcs             number of emissions species
c
c     Output arguments:
c        pntems              point source emissions rate (mole/s)
c
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     02/20/17   --gwilson--    Original development
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'flags.inc'
      include 'vbs.inc'
      include 'netcdf.inc'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      integer ifile 
      integer iunit
      integer numpts
      integer numspcs
      real    pntems(numpts,numspcs)
c
c-----------------------------------------------------------------------
c    External functions:
c-----------------------------------------------------------------------
c
      integer istrln
      integer ncf_get_tstep
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      character*200 action
      character*16  this_var
      integer       this_tstep, ierr, ispc, ipt, ibin
      integer       idx_start, this_time_tflag, this_time_etflag
      integer       this_varid, data_start(2), data_count(2)
      real          pi
c
      real, allocatable, dimension(:) ::  emispts
      real, allocatable, dimension(:) ::  plume_bot
      real, allocatable, dimension(:) ::  plume_top
c
c-----------------------------------------------------------------------
c    Data statements:
c-----------------------------------------------------------------------
c
      data pi /3.1415927/
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
      if( .NOT. is_netcdf_iptem(ifile) ) goto 9999
      write(action,'(2A,I3)') 'Reading the NetCDF point source',
     &                               ' emissions file. File: ',ifile
      idx_start = idx_start_pts(ifile)
c
c  ---- allocate the temporary array ---
c
      allocate( emispts(nptsrc_files(ifile)) )
      allocate( plume_bot(nptsrc_files(ifile)) )
      allocate( plume_top(nptsrc_files(ifile)) )
      emispts = 0.
      plume_bot = 0.
      plume_top = 0.
c
c  ---- get the index for timestep containing this time ----
c
      this_tstep = ncf_get_tstep(iunit,action,date,time,
     &               this_time_tflag,this_time_etflag,le1day,.TRUE.)
c
c   --- set the indexes for what to read ---
c
      data_start(1) = 1
      data_count(1) = nptsrc_files(ifile)
      data_start(2) = this_tstep
      data_count(2) = 1

      do ipt=1,nptsrc_files(ifile)
         effph(ipt+idx_start) = 0.
         flowrat(ipt+idx_start) = 0.
      enddo
c
c  ---- plumerise ---
c
      this_var = 'plumerise'
      ierr = nf_inq_varid(iunit, this_var, this_varid)
      if( ierr .EQ. NF_NOERR) then
          ierr = nf_get_vara_real(iunit,this_varid,data_start,
     &                                  data_count,emispts)
          if( ierr .NE. NF_NOERR) goto 7001
c
c  --- Put the data into the global array ----
c
          do ipt=1,nptsrc_files(ifile)
             if( emispts(ipt) .GT. 0. ) effph(ipt+idx_start) = -emispts(ipt)
          enddo
      endif
c
c  ---- loop over the list of emissions spcies ---
c
      do ispc=1,nspec
          if( lemmap(ispc) .LE. 0 ) cycle
c
c  --- load data into temporary arrray ---
c
          this_var = emspcname(lemmap(ispc))
          ierr = nf_inq_varid(iunit, this_var, this_varid)
          if( ierr .NE. NF_NOERR) cycle
          ierr = nf_get_vara_real(iunit,this_varid,data_start,
     &                                               data_count,emispts)
          if( ierr .NE. NF_NOERR) goto 7001
c
c-----Put the data into the global array ----
c
          do ipt=1,nptsrc_files(ifile)
             pntems(idx_start+ipt,lemmap(ispc)) = emispts(ipt)/(60.*dtems)
          enddo
c
c ---- next species ---
c
      enddo
c
c ---- handle VBS as special case ---
c
      if( lvbs .AND. LVBSPREPROC ) then
          this_var = 'POA_OP    '
          ierr = nf_inq_varid(iunit, this_var, this_varid)
          if( ierr .EQ. NF_NOERR) then
             ierr = nf_get_vara_real(iunit,this_varid,data_start,
     &                                               data_count,emispts)
             if( ierr .NE. NF_NOERR) goto 7001
             do ibin = 0, NVOLBIN
                do ispc=1,nemspc
                   if( emspcname(ispc) .EQ. spname(kpap_c(ibin)) ) then
                     do ipt=1,nptsrc_files(ifile)
                        pntems(idx_start+ipt,ispc) = 
     &                     ( emispts(ipt) * poa_op_ef(ibin) ) / (60.*dtems)
                     enddo
                   endif
                enddo
             enddo
          endif
          this_var = 'POA_BB    '
          ierr = nf_inq_varid(iunit, this_var, this_varid)
          if( ierr .EQ. NF_NOERR) then
             ierr = nf_get_vara_real(iunit,this_varid,data_start,
     &                                               data_count,emispts)
             if( ierr .NE. NF_NOERR) goto 7001
             do ibin = 0, NVOLBIN
                do ispc=1,nemspc
                   if( emspcname(ispc) .EQ. spname(kpfp_c(ibin)) ) then
                     do ipt=1,nptsrc_files(ifile)
                        pntems(idx_start+ipt,ispc) = 
     &                ( emispts(ipt) * poa_bb_ef(ibin) ) / (60.*dtems)
                     enddo
                   endif
                enddo
             enddo
          endif
      endif
c
c   --- get plume bottom, if available ---
c
      this_var = 'plume_bottom'
      ierr = nf_inq_varid(iunit, this_var, this_varid)
      if( ierr .EQ. NF_NOERR) then
          ierr = nf_get_vara_real(iunit,this_varid,data_start,
     &                                  data_count,plume_bot)
          if( ierr .NE. NF_NOERR) goto 7001
      endif
c
c   --- get plume top, if available ---
c
      this_var = 'plume_top'
      ierr = nf_inq_varid(iunit, this_var, this_varid)
      if( ierr .EQ. NF_NOERR) then
          ierr = nf_get_vara_real(iunit,this_varid,data_start,
     &                                  data_count,plume_top)
          if( ierr .NE. NF_NOERR) goto 7001
      endif
c
c  --- Put the data into the global array ----
c
      do ipt=1,nptsrc_files(ifile)
          if( plume_top(ipt) .GT. 0. ) then
             if( plume_bot(ipt) .GE. plume_top(ipt) ) goto 7002
             if( plume_bot(ipt) .LT. 0. ) goto 7002
             effph(ipt+idx_start) = -plume_bot(ipt)
             flowrat(ipt+idx_start) = -plume_top(ipt)
          endif
      enddo
c
c  ---- deallocate the temporary array ---
c
      deallocate( emispts )
      deallocate( plume_bot )
      deallocate( plume_top )
c
c  --- write message to out file ---
c
      write(iout,'(a35,i5,a,f7.0,i8.5)') 'Read point source file ',
     &                                         ifile,' at ',time,date
c
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,a)') 'ERROR in NCF_READPT: '
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot find variable id for: ',
     &                                      this_var(:istrln(this_var))
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7001 continue
      write(iout,'(//,a)') 'ERROR in NCF_READPT:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot read data for variable: ',
     &                                      this_var(:istrln(this_var))
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7002 continue
      write(iout,'(//,a)') 'ERROR in NCF_READPT:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Point source found with invalid plume ',
     &                                        'distribution values.'
      write(iout,'(A,I10)') 'Point source number in file: ',ipt
      write(iout,'(A,F10.4)') 'Plume bottom: ',plume_bot(ipt)
      write(iout,'(A,F10.4)') 'Plume top   : ',plume_top(ipt)
      call camxerr()
c
 9999 continue
      return
      end
