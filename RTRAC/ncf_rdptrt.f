c***** NCF_RDPTRT.F
c
      subroutine ncf_rdptrt(ndate,ttime)
      use tracer
      use ptemiss
      use filunit
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
c   This routine reads one hour of emissions for the RTRAC process
c   and fills the approproate arrays.  The emissions file for one grid
c   but each emissions groups is read. This version is for NetCDF files.
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
c     06/13/19   --gwilson--   Originial development
c
c-----------------------------------------------------------------------
c   Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'camx.inc'
      include 'flags.inc'
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
      integer istrln
      integer ncf_get_tstep
c
c-----------------------------------------------------------------------
c   Local variables:
c-----------------------------------------------------------------------
c
      character*200 fname, action
      character*10  this_var, cname
      integer       ierr, iounit, this_tstep, idxpt, i, j
      integer       this_time_tflag, this_time_etflag, idx_start
      integer       data_start(2), data_count(2), this_varid, ispc
c
      real,    allocatable, dimension(:) ::  emispts
      integer, allocatable, dimension(:) ::  sa_region
      real,    allocatable, dimension(:) ::  effph_in
c
c-----------------------------------------------------------------------
c   Edtry point:
c-----------------------------------------------------------------------
c
      action = 'Reading RTRAC point source file.'
c
c   --- skip if filename not supplied ---
c
      if( .NOT. is_netcdf_iortpt(1,1) .OR. .NOT. ltptfl(1,1) 
     &                             .OR. .NOT. lptsrc ) goto 9999
c
c  --- initialize emissions to zero ---
c
      do i = 1,nptsrc
        do j = 1,ntotsp
         sapnts(i,j) = 0.
        enddo
      enddo
c
c   --- set the unit number for surface emissions file ---
c
      iounit = iortpt(1,1)
      fname = tptfil(1,1)
c
c   --- allocate the local array ---
c
      allocate( emispts(nptsrc_safile(1,1)) )
      allocate( sa_region(nptsrc_safile(1,1)) )
      allocate( effph_in(nptsrc_safile(1,1)) )
c
c  ---- get the index for timestep containing this time ----
c
      this_tstep = ncf_get_tstep(iounit,action,ndate,ttime,
     &               this_time_tflag,this_time_etflag,le1day,.TRUE.)
c
c   --- set the indexes for what to read ---
c
      data_start(1) = 1
      data_count(1) = nptsrc_safile(1,1)
      data_start(2) = this_tstep
      data_count(2) = 1
      idx_start = nptsrc - nptsrc_safile(1,1)
c
c  ---- plumerise ---
c
      this_var = 'plumerise'
      ierr = nf_inq_varid(iounit, this_var, this_varid)
      if( ierr .NE. NF_NOERR) goto 7000
      ierr = nf_get_vara_real(iounit,this_varid,data_start,
     &                                  data_count,effph_in)
      if( ierr .NE. NF_NOERR) goto 7001
c
c  ---- put into global array ---
c
      do idxpt = 1,nptsrc_safile(1,1)
          effph(idxpt+idx_start) = effph_in(idxpt)
      enddo
c
c   --- read the emissions for this hour ---
c
      do ispc=1,ntotsp
c
c   --- if species not in file just skip it ---
c
         this_var = ptname(ispc)
         ierr = nf_inq_varid(iounit, this_var, this_varid)
         if( ierr .NE. NF_NOERR) cycle
c
c  --- read data and load it into the global array ---
c
         ierr = nf_get_vara_real(iounit,this_varid,data_start,
     &                                               data_count,emispts)
         if( ierr .NE. NF_NOERR) goto 7000
         do idxpt = 1,nptsrc_safile(1,1)
             sapnts(idxpt+idx_start,ispc) = emispts(idxpt)/(60.*dtems)
         enddo
c
c  ---- next species ---
c
      enddo
c
c  --- if not found in model species list, write a message ----     
c   
cae   write(idiag,'(1X,4A)') 'Species in RTRAC ',
cae  &            'point source file: ',cname(:istrln(cname)),
cae  &            ' not found in species list ... Skipping.'
c
c  --- get the next file ---
c
      goto 9999
c
c----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,a)') 'ERROR in NCF_RDPTRT: '
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot find variable in file: ',
     &                                      this_var(:istrln(this_var))
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7001 continue
      write(iout,'(//,a)') 'ERROR in NCF_RDPTRT: '
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot read data for variable: ',
     &                                      this_var(:istrln(this_var))
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
c----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
c
      return
      end
