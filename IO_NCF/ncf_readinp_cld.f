C**** NCF_READINP_CLD
c
      subroutine ncf_readinp_cld(igrd,ncol,nrow,nlay,cldwtr,ranwtr,
     &                      snowtr,gplwtr,cldod,cldph,cigfrc,cigtim,
     &                                    cigwtr,cigpcp,cigent,cigdet)
      use filunit
      implicit none
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routines reads the cloud data in the 3D met file.
c 
c      Copyright 1996 - 2022
c     Ramboll
c           
c     Input arguments:
c        igrd                grid index
c        ncol                number of columns
c        nrow                number of rows
c        nlay                number of layers
c
c     Output arguments:
c        cldwtr              cloud water content (g/m3)
c        ranwtr              rain water content (g/m3)
c        snowtr              snow water content (g/m3)
c        gplwtr              graupel water content (g/m3)
c        cldod               cloud optical depth
c        cldph               cloud water pH
c        cigfrc              convective cloud fraction (unitless)
c        cigtim              convective cloud timescale (s)
c        cigwtr              convective cloud water (g/m3)
c        cigpcp              convective cloud precip (g/m3)
c        cigent              convective cloud entrainment (kg/m2/s)
c        cigdet              convective cloud detrainment (kg/m2/s)
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
      include 'camx.inc'
      include 'flags.inc'
      include 'netcdf.inc'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      integer igrd
      integer ncol
      integer nrow
      integer nlay
      real    cldwtr(ncol,nrow,nlay)
      real    ranwtr(ncol,nrow,nlay)
      real    snowtr(ncol,nrow,nlay)
      real    gplwtr(ncol,nrow,nlay)
      real    cldod(ncol,nrow,nlay)
      real    cldph(ncol,nrow,nlay)
      real    cigfrc(ncol,nrow)
      real    cigtim(ncol,nrow)
      real    cigwtr(ncol,nrow,nlay)
      real    cigpcp(ncol,nrow,nlay)
      real    cigent(ncol,nrow,nlay)
      real    cigdet(ncol,nrow,nlay)
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
      character*20  this_var
      integer       this_varid, ierr, data_start_3d(4), data_count_3d(4)
      integer       data_start_2d(3), data_count_2d(3)
      integer       this_tstep, nfound, i, j, k
      integer       this_time_tflag, this_time_etflag
      real          whr, wmn, atim
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
      call zeros(cldwtr,ncol*nrow*nlay)
      call zeros(ranwtr,ncol*nrow*nlay)
      call zeros(snowtr,ncol*nrow*nlay)
      call zeros(gplwtr,ncol*nrow*nlay)
      call zeros(cldod,ncol*nrow*nlay)
c
      if( .NOT. is_netcdf_i3dmet(igrd) ) goto 9999
c
      write(action,'(A,I3)') 'Reading Cloud data in 3D met file for grid ',igrd
c
c --- Find the date/time at the next update interval ---
c
      whr = aint(time/100.)
      wmn = anint(amod(time,100.))
      atim = 100.*whr + wmn
c
c   --- get the timestep that encompasses this date/time ---
c
      this_tstep = ncf_get_tstep(i3dmet(igrd),action,date,atim,
     &                this_time_tflag,this_time_etflag,.FALSE.,.FALSE.)
c
c   --- set the indexes for what to read ---
c
      data_start_3d(1) = 1
      data_count_3d(1) = ncol
      data_start_3d(2) = 1
      data_count_3d(2) = nrow
      data_start_3d(3) = 1
      data_count_3d(3) = nlay
      data_start_3d(4) = this_tstep
      data_count_3d(4) = 1
c
      data_start_2d(1) = 1
      data_count_2d(1) = ncol
      data_start_2d(2) = 1
      data_count_2d(2) = nrow
      data_start_2d(3) = this_tstep
      data_count_2d(3) = 1
c
c ---- cloud water content ----
c
      this_var = "cloudwater"
      ierr = nf_inq_varid(i3dmet(igrd), this_var, this_varid)
      if( ierr .EQ. NF_NOERR) then
         ierr = nf_get_vara_real(i3dmet(igrd),this_varid,data_start_3d,
     &                                             data_count_3d,cldwtr)
         if( ierr .NE. NF_NOERR) goto 7001
      endif
c
c ---- rain water content ----
c
      this_var = "rainwater"
      ierr = nf_inq_varid(i3dmet(igrd), this_var, this_varid)
      if( ierr .EQ. NF_NOERR) then
         ierr = nf_get_vara_real(i3dmet(igrd),this_varid,data_start_3d,
     &                                             data_count_3d,ranwtr)
         if( ierr .NE. NF_NOERR) goto 7001
      endif
c
c ---- snow water content ----
c
      this_var = "snowater"
      ierr = nf_inq_varid(i3dmet(igrd), this_var, this_varid)
      if( ierr .EQ. NF_NOERR) then
         ierr = nf_get_vara_real(i3dmet(igrd),this_varid,data_start_3d,
     &                                             data_count_3d,snowtr)
         if( ierr .NE. NF_NOERR) goto 7001
      endif
c
c ---- graupel water content ----
c
      this_var = "grplwater"
      ierr = nf_inq_varid(i3dmet(igrd), this_var, this_varid)
      if( ierr .EQ. NF_NOERR) then
         ierr = nf_get_vara_real(i3dmet(igrd),this_varid,data_start_3d,
     &                                             data_count_3d,gplwtr)
         if( ierr .NE. NF_NOERR) goto 7001
      endif
c
c ---- cloud optical depth ---
c
      this_var = "cloudod"
      ierr = nf_inq_varid(i3dmet(igrd), this_var, this_varid)
      if( ierr .EQ. NF_NOERR) then
         ierr = nf_get_vara_real(i3dmet(igrd),this_varid,data_start_3d,
     &                                             data_count_3d,cldod)
         if( ierr .NE. NF_NOERR) goto 7001
      endif
c
c-----Initialize cloud pH
c
      cldph = 5.
c
c --- if not doing convection in cloud just exit ---
c
      if( .NOT. lcig(igrd) ) goto 9999
c
c ---- convective cloud fraction ----
c
      nfound = 0
      this_var = "kf_cldfrac"
      ierr = nf_inq_varid(i3dmet(igrd), this_var, this_varid)
      if( ierr .EQ. NF_NOERR) then
         nfound = nfound + 1
         ierr = nf_get_vara_real(i3dmet(igrd),this_varid,data_start_2d,
     &                                             data_count_2d,cigfrc)
         if( ierr .NE. NF_NOERR) goto 7001
      endif
c
c ----  convective cloud timescale ----
c
      this_var = "kf_tscale"
      ierr = nf_inq_varid(i3dmet(igrd), this_var, this_varid)
      if( ierr .EQ. NF_NOERR) then
         nfound = nfound + 1
         ierr = nf_get_vara_real(i3dmet(igrd),this_varid,data_start_2d,
     &                                             data_count_2d,cigtim)
         if( ierr .NE. NF_NOERR) goto 7001
      endif
c
c ---- convective cloud water ---
c
      this_var = "kf_cldwater"
      ierr = nf_inq_varid(i3dmet(igrd), this_var, this_varid)
      if( ierr .EQ. NF_NOERR) then
         nfound = nfound + 1
         ierr = nf_get_vara_real(i3dmet(igrd),this_varid,data_start_3d,
     &                                             data_count_3d,cigwtr)
         if( ierr .NE. NF_NOERR) goto 7001
      endif
c
c ---- convective cloud precip ---
c
      this_var = "kf_pcpwater"
      ierr = nf_inq_varid(i3dmet(igrd), this_var, this_varid)
      if( ierr .EQ. NF_NOERR) then
         nfound = nfound + 1
         ierr = nf_get_vara_real(i3dmet(igrd),this_varid,data_start_3d,
     &                                             data_count_3d,cigpcp)
         if( ierr .NE. NF_NOERR) goto 7001
      endif
c
c ---- convective cloud entrainment ---
c
      this_var = "kf_entrain"
      ierr = nf_inq_varid(i3dmet(igrd), this_var, this_varid)
      if( ierr .EQ. NF_NOERR) then
         nfound = nfound + 1
         ierr = nf_get_vara_real(i3dmet(igrd),this_varid,data_start_3d,
     &                                             data_count_3d,cigent)
         if( ierr .NE. NF_NOERR) goto 7001
      endif
c
c ---- convective cloud detrainment ---
c
      this_var = "kf_detrain"
      ierr = nf_inq_varid(i3dmet(igrd), this_var, this_varid)
      if( ierr .EQ. NF_NOERR) then
         nfound = nfound + 1
         ierr = nf_get_vara_real(i3dmet(igrd),this_varid,data_start_3d,
     &                                             data_count_3d,cigdet)
         if( ierr .NE. NF_NOERR) goto 7001
      endif
c
c  --- make sure everything is here for convection in cloud ---
c
      if( nfound .NE. 6 ) then
         lcig(igrd) = .false.
         write(iout,'(//,a)')'WARNING in NCF_READINP_CLD:'
         write(iout,'(a)') action
         write(iout,'(2a)')'Subgrid convective cloud data were not',
     &                       ' found on file'
         write(iout,'(2a)')'Convective treatment will not be',
     &                       ' applied for this grid'
         write(iout,*)
      endif
c
c  --- write message to out file ---
c
      write(iout,'(a40,f7.0,i8.5,a,i3)') 
     &    'Read cloud data from 3D met file at ',atim,date,' grid',igrd
      call flush(iout)
c
c  --- successful completion ---
c
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,a)') 'ERROR in NCF_READINP_CLD:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot find variable id for: ',
     &                                      this_var(:istrln(this_var))
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7001 continue
      write(iout,'(//,a)') 'ERROR in NCF_READINP_CLD:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot read data for variable: ',
     &                                      this_var(:istrln(this_var))
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 9999 continue
      return
      end
