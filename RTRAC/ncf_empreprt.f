*** NCF_EMPREPRT
c
      subroutine ncf_empreprt(igrid,begtim,begdate,endtim,enddate,buffer_offset)
      use ptemiss
      use grid
      use tracer
      use filunit
      implicit none
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routine reads the emissions files for the RTRAC and moves the 
c   file pointer to the first hour to be used in this run.
c
c      Copyright 1996 - 2022
c     Ramboll
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c      09/09/19   --gwilson-- original development
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'flags.inc'
      include 'netcdf.inc'
c
c-----------------------------------------------------------------------
c    Arguement declarations:
c-----------------------------------------------------------------------
c
      integer igrid
      real    begtim
      integer begdate
      real    endtim
      integer enddate
      integer buffer_offset
c
c-----------------------------------------------------------------------
c    External functions:
c-----------------------------------------------------------------------
c
      integer istrln
c
c-----------------------------------------------------------------------
c    Local parameters:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      character*200 action, fname
      character*10  emfil, pntfil, name_in, this_var
      integer       ierr, iounit, this_dimid, this_varid, ipt
      integer       idxpt_base, igrd, ii, jj, numpts
      logical       lno_buffer_cells
c
      real,    allocatable, dimension(:) :: xlocin
      real,    allocatable, dimension(:) :: ylocin
      real,    allocatable, dimension(:) :: hstkin
      real,    allocatable, dimension(:) :: dstkin
      real,    allocatable, dimension(:) :: tstkin
      real,    allocatable, dimension(:) :: vstkin
      real,    allocatable, dimension(:) :: iarray
c
c-----------------------------------------------------------------------
c    Data statments:
c-----------------------------------------------------------------------
c
      data emfil /'EMISSIONS '/
      data pntfil /'PTSOURCE  '/
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c
c   --- only read if surface emissions file is supplied ---
c
      if( .NOT. is_netcdf_iortem(igrid,1,1) ) goto 111
      if( .NOT. ltemfl(igrid,1,1) ) goto 111
      if( .NOT. larsrc ) goto 111

      iounit = iortem(igrid,1,1)
      fname = temfil(igrid,1,1)
      write(action,'(2A,I3)') 'Reading the NetCDF gridded RTRAC emissions ',
     &                                              'file for grid: ',igrid
c
c --- get the type of file to make sure it is a emissions file ---
c
      this_var = 'CAMx_NAME'
      ierr = nf_get_att_text(iounit, NF_GLOBAL, this_var, name_in)
      if( ierr .NE. NF_NOERR ) then
         this_var = 'NAME'
         ierr = nf_get_att_text(iounit, NF_GLOBAL, this_var, name_in)
         if( ierr .NE. NF_NOERR ) goto 7006
      endif
      if( name_in(:9) .NE. emfil(:9) ) goto 7001
c
c --- call routine to make sure grid defintion is consistent ---
c
      lno_buffer_cells = .TRUE.
      if( igrd .EQ. 1 ) lno_buffer_cells = .FALSE.
      call ncf_chk_griddef(iounit,action,igrid,.FALSE.,.TRUE.,
     &                         lno_buffer_cells,.FALSE.,buffer_offset)
c
c --- call routine to make sure file spans the episode ---
c
      call ncf_chk_tstep(iounit,action,begdate,begtim,enddate,endtim,le1day)
c
c   --- only read if elevated point source emissions file is supplied ---
c
  111 continue
      if( .NOT. is_netcdf_iortpt(1,1) ) goto 9999
      if( .NOT. ltptfl(1,1) ) goto 9999
      if( .NOT. lptsrc ) goto 9999
c
c   --- process the file ----
c
      action = 'Reading NetCDF point source file for RTRAC.'
      iounit = iortpt(1,1)
      fname = tptfil(1,1)
c
c --- get the type of file to make sure it is a point source file ---
c
      this_var = 'CAMx_NAME'
      ierr = nf_get_att_text(iounit, NF_GLOBAL, this_var, name_in)
      if( ierr .NE. NF_NOERR ) then
         this_var = 'NAME'
         ierr = nf_get_att_text(iounit, NF_GLOBAL, this_var, name_in)
         if( ierr .NE. NF_NOERR ) goto 7000
      endif
      if( name_in(:istrln(name_in)) .NE. pntfil(:istrln(pntfil)) ) goto 7003
c
c --- call routine to make sure file spans the episode ---
c
      call ncf_chk_tstep(iounit,action,begdate,begtim,
     &                                        enddate,endtim,le1day)
c
c --- get the number of points (used as dimension for columns) ---
c
      ierr = nf_inq_dimid(iounit, "COL", this_dimid )
      if( ierr .NE. NF_NOERR ) goto 7002
      ierr = nf_inq_dimlen(iounit,this_dimid,
     &                                    nptsrc_safile(1,1))
      if( ierr .NE. NF_NOERR ) goto 7002
      numpts = nptsrc - nptsrc_safile(1,1)
c
c  --- allocate the local arrays ---
c
      allocate( xlocin(nptsrc_safile(1,1)) )
      allocate( ylocin(nptsrc_safile(1,1)) )
      allocate( hstkin(nptsrc_safile(1,1)) )
      allocate( dstkin(nptsrc_safile(1,1)) )
      allocate( tstkin(nptsrc_safile(1,1)) )
      allocate( vstkin(nptsrc_safile(1,1)) )
      allocate( iarray(nptsrc_safile(1,1)) )
c
      this_var = "xcoord"
      ierr = nf_inq_varid(iounit,this_var,this_varid)
      if( ierr .NE. NF_NOERR ) goto 7004
      ierr = nf_get_var_real(iounit,this_varid,xlocin)
      if( ierr .NE. NF_NOERR ) goto 7005
c
      this_var = "ycoord"
      ierr = nf_inq_varid(iounit,this_var,this_varid)
      if( ierr .NE. NF_NOERR ) goto 7004
      ierr = nf_get_var_real(iounit,this_varid,ylocin)
      if( ierr .NE. NF_NOERR ) goto 7005
c
c --- get pig flag ---
c
      this_var = "pigflag"
      ierr = nf_inq_varid(iounit,this_var,this_varid)
      if( ierr .NE. NF_NOERR ) goto 7004
      ierr = nf_get_var_int(iounit,this_varid,iarray)
      if( ierr .NE. NF_NOERR ) goto 7005
c
      this_var = "stkheight"
      ierr = nf_inq_varid(iounit,this_var,this_varid)
      if( ierr .NE. NF_NOERR ) goto 7004
      ierr = nf_get_var_real(iounit,this_varid,hstkin)
      if( ierr .NE. NF_NOERR ) goto 7005
c
      this_var = "stkdiam"
      ierr = nf_inq_varid(iounit,this_var,this_varid)
      if( ierr .NE. NF_NOERR ) goto 7004
      ierr = nf_get_var_real(iounit,this_varid,dstkin)
      if( ierr .NE. NF_NOERR ) goto 7005
c
      this_var = "stktemp"
      ierr = nf_inq_varid(iounit,this_var,this_varid)
      if( ierr .NE. NF_NOERR ) goto 7004
      ierr = nf_get_var_real(iounit,this_varid,tstkin)
      if( ierr .NE. NF_NOERR ) goto 7005
c
      this_var = "stkspeed"
      ierr = nf_inq_varid(iounit,this_var,this_varid)
      if( ierr .NE. NF_NOERR ) goto 7004
      ierr = nf_get_var_real(iounit,this_varid,vstkin)
      if( ierr .NE. NF_NOERR ) goto 7005
c
c  --- make sure this source is in the base inventory and flag it ---
c
      do ipt=1,nptsrc_safile(1,1)
          vstkin(ipt) = vstkin(ipt)/3600.
          numpts = numpts + 1 
          if (llatlon) then
            xstk(numpts,1) = xlocin(ipt) - xorg
            ystk(numpts,1) = ylocin(ipt) - yorg
          else
            xstk(numpts,1) = xlocin(ipt)/1000. - xorg
            ystk(numpts,1) = ylocin(ipt)/1000. - yorg
          endif
          do igrd = 2,ngrid
             xstk(numpts,igrd) = xstk(numpts,1) - (inst1(igrd) - 1)*delx
             ystk(numpts,igrd) = ystk(numpts,1) - (jnst1(igrd) - 1)*dely
             ii = 2 + FLOOR(xstk(numpts,igrd)/delx*FLOAT( meshold(igrd) ) )
             jj = 2 + FLOOR(ystk(numpts,igrd)/dely*FLOAT( meshold(igrd) ) )
             if (ii.gt.1 .and. ii.lt.ncol(igrd) .and.
     &                    jj.gt.1 .and. jj.lt.nrow(igrd)) then
                nosrc(igrd) = nosrc(igrd) + 1
                idsrc(nosrc(igrd),igrd) = numpts
                isrc(nosrc(igrd),igrd) = ii
                jsrc(nosrc(igrd),igrd) = jj
             endif
          enddo
          xloc(numpts) = xlocin(ipt)
          yloc(numpts) = ylocin(ipt)
          xlocpt(numpts) = xlocin(ipt)
          ylocpt(numpts) = ylocin(ipt)
          hstk(numpts) = hstkin(ipt)
          dstk(numpts) = ABS(dstkin(ipt))
          tstk(numpts) = tstkin(ipt)
          vstk(numpts) = vstkin(ipt)
          lpiglet(numpts) = .FALSE.
          idx_point_in_list(1,1,ipt) = numpts
          if( iarray(ipt) .GT. 0 ) lpiglet(numpts) = .TRUE.
      enddo
c
c  --- deallocate the local arrays ---
c
      deallocate( xlocin )
      deallocate( ylocin )
      deallocate( hstkin )
      deallocate( dstkin )
      deallocate( tstkin )
      deallocate( vstkin )
      deallocate( iarray )
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
      write(iout,'(//,a)') 'ERROR in NCF_EMPREPRT:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot get necessary global attribute: ',
     &                                      this_var(:istrln(this_var))
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7001 continue
      write(iout,'(//,a)') 'ERROR in NCF_EMPREPRT:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(A)') 'Input file is not the correct type.'
      write(iout,'(2A)') 'Looking for type: ',emfil(:istrln(emfil))
      write(iout,'(2A)') 'Found in file   : ',name_in(:istrln(name_in))
      call camxerr()
c
 7002 continue
      write(iout,'(//,a)') 'ERROR in NCF_EMPREPRT:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(A)') 'Cannot find the dimension id for number ',
     &                                      'of point sources (NCOL)'
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7003 continue
      write(iout,'(//,a)') 'ERROR in NCF_EMPREPRT:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(A)') 'Input file is not the correct type.'
      write(iout,'(2A)') 'Looking for type: ',pntfil(:istrln(emfil))
      write(iout,'(2A)') 'Found in file   : ',name_in(:istrln(name_in))
      call camxerr()
c
 7004 continue
      write(iout,'(//,a)') 'ERROR in NCF_EMPREPRT:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot find variable id for: ',
     &                                      this_var(:istrln(this_var))
      write(iout,'(1X,2A)') 'File: ',fname(:istrln(fname))
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7005 continue
      write(iout,'(//,a)') 'ERROR in NCF_EMPREPRT:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot read data for variable: ',
     &                                      this_var(:istrln(this_var))
      write(iout,'(1X,2A)') 'File: ',fname(:istrln(fname))
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7006 continue
      write(iout,'(//,A)') 'ERROR in NCF_EMPREPRT:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)',ERR=9999) 'Cannot find global variable for ',
     &                                       'type of file: CAMx_NAME'
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
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
c
      return
      end
