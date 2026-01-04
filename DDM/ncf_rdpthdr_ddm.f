      subroutine ncf_rdpthdr_ddm(igroup,idxfile,begtim,begdate,
     &                              endtim,enddate,numpts_tot)
      use filunit
      use ptemiss
      use tracer
      use grid
      implicit none
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c     NCF_RDPTHDR_DDM reads the header information from the global section
c     of the NetCDF point source file. This version is for DDM.
c
c      Copyright 1996 - 2022
c     Ramboll
c
c      Argument description:
c       Inputs:
c         igroup     I group number for this file
c         idxfile    I index for this file in arrays
c         begtim     R model begin time
c         begdat     I model begin date (YYJJJ)
c         endtim     R model end time
c         enddate    I model end date (YYJJJ)
c         numpts_tot I number of points read so far
c       Outputs:
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
      include 'netcdf.inc'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      integer igroup
      integer idxfile
      real    begtim
      integer begdate
      real    endtim
      integer enddate
      integer numpts_tot
c
c-----------------------------------------------------------------------
c    External functions:
c-----------------------------------------------------------------------
c
      integer istrln
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      character*200 action, fname
      character*10  pntfil, name_in, this_var
      integer       this_dimid, this_varid, ierr, tmpmap(MXSPEC), numpts
      integer       ipt, iounit, idxpt_base, ii, jj, igrd, idx_start, i
      integer       numpts_in, itmp
c
      real,    allocatable, dimension(:) :: txloc
      real,    allocatable, dimension(:) :: tyloc
      real,    allocatable, dimension(:) :: tstkhght
      real,    allocatable, dimension(:) :: tstkdiam
      real,    allocatable, dimension(:) :: tstktemp
      real,    allocatable, dimension(:) :: tstkvelo
      real,    allocatable, dimension(:) :: iarray
      logical, allocatable, dimension(:) :: already_found
c
c-----------------------------------------------------------------------
c    Data statments:
c-----------------------------------------------------------------------
c
      data pntfil /'PTSOURCE  '/
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c   --- assume no sources have been found ---
c
      allocate( already_found(MXPTSRC) )
      already_found = .FALSE.
c
c   --- loop over number of files ---
c
      iounit = iortpt(igroup,idxfile)
      if( .NOT. is_netcdf_iortpt(igroup,idxfile) ) goto 9999
      write(action,'(2A,I3)')  'Reading NetCDF DDM point source file.',
     &                                               ' Group: ',igroup
      fname = tptfil(igroup,idxfile)
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
      if( name_in(:8) .NE. pntfil(:8) ) goto 7001
c
c --- call routine to make sure file spans the episode ---
c
      call ncf_chk_tstep(iounit,action,begdate,begtim,
     &                                        enddate,endtim,le1day)
c
c --- call routine to setup species mappping array ---
c
      call ncf_set_emiss_mapping(iounit,action)
c
c --- get the number of points (used as dimension for columns) ---
c
      ierr = nf_inq_dimid(iounit, "COL", this_dimid )
      if( ierr .NE. NF_NOERR ) goto 7002
      ierr = nf_inq_dimlen(iounit,this_dimid,numpts_in)
      if( ierr .NE. NF_NOERR ) goto 7002
      idx_start_sapts(igroup,idxfile) = numpts_tot
      numpts_tot = numpts_tot + numpts_in
c
c  --- check that parameters are the same ---
c
      nptsrc_safile(igroup,idxfile) = numpts_in
      allocate( txloc(numpts_in) )
      allocate( tyloc(numpts_in) )
      allocate( tstkhght(numpts_in) )
      allocate( tstkdiam(numpts_in) )
      allocate( tstktemp(numpts_in) )
      allocate( tstkvelo(numpts_in) )
      allocate( iarray(numpts_in) )
      this_var = "xcoord"
      ierr = nf_inq_varid(iounit,this_var,this_varid)
      if( ierr .NE. NF_NOERR ) goto 7005
      ierr = nf_get_var_real(iounit,this_varid,txloc)
      if( ierr .NE. NF_NOERR ) goto 7006
c
      this_var = "ycoord"
      ierr = nf_inq_varid(iounit,this_var,this_varid)
      if( ierr .NE. NF_NOERR ) goto 7005
      ierr = nf_get_var_real(iounit,this_varid,tyloc)
      if( ierr .NE. NF_NOERR ) goto 7006
c
      this_var = "stkheight"
      ierr = nf_inq_varid(iounit,this_var,this_varid)
      if( ierr .NE. NF_NOERR ) goto 7005
      ierr = nf_get_var_real(iounit,this_varid,tstkhght)
      if( ierr .NE. NF_NOERR ) goto 7006
c
      this_var = "stkdiam"
      ierr = nf_inq_varid(iounit,this_var,this_varid)
      if( ierr .NE. NF_NOERR ) goto 7005
      ierr = nf_get_var_real(iounit,this_varid,tstkdiam)
      if( ierr .NE. NF_NOERR ) goto 7006
c
      this_var = "stktemp"
      ierr = nf_inq_varid(iounit,this_var,this_varid)
      if( ierr .NE. NF_NOERR ) goto 7005
      ierr = nf_get_var_real(iounit,this_varid,tstktemp)
      if( ierr .NE. NF_NOERR ) goto 7006
c
      this_var = "stkspeed"
      ierr = nf_inq_varid(iounit,this_var,this_varid)
      if( ierr .NE. NF_NOERR ) goto 7005
      ierr = nf_get_var_real(iounit,this_varid,tstkvelo)
      if( ierr .NE. NF_NOERR ) goto 7006
      do ipt=1,numpts_in
        tstkvelo(ipt) = tstkvelo(ipt)/3600.
      enddo
c
c --- get pig flag ---
c
      this_var = "pigflag"
      ierr = nf_inq_varid(iounit,this_var,this_varid)
      if( ierr .NE. NF_NOERR ) goto 7005
      ierr = nf_get_var_int(iounit,this_varid,iarray)
      if( ierr .NE. NF_NOERR ) goto 7006
      idx_start = idx_start_sapts(igroup,idxfile)
      do ipt=1,numpts_in
c
c  --- make sure this source is in the base inventory and flag it ---
c
          if( llatlon ) then
              xstk(idx_start+ipt,1) = txloc(ipt) - xorg
              ystk(idx_start+ipt,1) = tyloc(ipt) - yorg
          else
              xstk(idx_start+ipt,1) = txloc(ipt)/1000. - xorg
              ystk(idx_start+ipt,1) = tyloc(ipt)/1000. - yorg
          endif
          xloc(idx_start+ipt) = txloc(ipt)
          yloc(idx_start+ipt) = tyloc(ipt)
          do igrd = 2,ngrid
             xstk(idx_start+ipt,igrd) =
     &                    xstk(idx_start+ipt,1) - (inst1(igrd)-1)*delx
             ystk(idx_start+ipt,igrd) =
     &                    ystk(idx_start+ipt,1) - (jnst1(igrd)-1)*dely
             ii = 2 + FLOOR(xstk(idx_start+ipt,igrd)/
     &                        delx*FLOAT( meshold(igrd) ) )
             jj = 2 + FLOOR(ystk(idx_start+ipt,igrd)/
     &                        dely*FLOAT( meshold(igrd) ) )
             if (ii .GT. 1 .AND. ii .LT. ncol(igrd) .AND.
     &                    jj .GT. 1 .AND. jj .LT. nrow(igrd)) then
                nosrc(igrd) = nosrc(igrd) + 1
                idsrc(nosrc(igrd),igrd) = idx_start+ipt
                isrc(nosrc(igrd),igrd) = ii
                jsrc(nosrc(igrd),igrd) = jj
             endif
         enddo
         xlocpt(idx_start+ipt) = txloc(ipt)
         ylocpt(idx_start+ipt) = tyloc(ipt)
         hstk(idx_start+ipt) = tstkhght(ipt)
         dstk(idx_start+ipt) = ABS(tstkdiam(ipt))
         tstk(idx_start+ipt) = tstktemp(ipt)
         vstk(idx_start+ipt) = tstkvelo(ipt)
         lpiglet(idx_start+ipt) = .FALSE.
         if( iarray(ipt) .GT. 0 ) lpiglet(idx_start+ipt) = .TRUE.
         idx_point_in_list(igroup,idxfile,ipt) = idx_start+ipt
      enddo
      deallocate( txloc )
      deallocate( tyloc )
      deallocate( tstkhght )
      deallocate( tstkdiam )
      deallocate( tstktemp )
      deallocate( tstkvelo )
      deallocate( iarray )
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,a)') 'ERROR in NCF_RDPTHDR_DDM:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot get necessary global attribute: ',
     &                                      this_var(:istrln(this_var))
      write(iout,'(1X,2A)') 'File: ',fname(:istrln(fname))
      call camxerr()
c
 7001 continue
      write(iout,'(//,a)') 'ERROR in NCF_RDPTHDR_DDM:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(A)') 'Input file is not the correct type.'
      write(iout,'(2A)') 'Looking for type: ',pntfil(:istrln(pntfil))
      write(iout,'(2A)') 'Found in file   : ',name_in(:istrln(name_in))
      write(iout,'(1X,2A)') 'File: ',fname(:istrln(fname))
      call camxerr()
c
 7002 continue
      write(iout,'(//,a)') 'ERROR in NCF_RDPTHDR_DDM:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(A)') 'Cannot find the dimension id for number ',
     &                                      'of point sources (NCOL)'
      write(iout,'(1X,2A)') 'File: ',fname(:istrln(fname))
      call camxerr()
c
 7003 continue
      write(iout,'(//,a)') 'ERROR in NCF_RDPTHDR_DDM:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(A)') 'Cannot find the dimension value for number ',
     &                                      'of point sources (NCOL)'
      write(iout,'(1X,2A)') 'File: ',fname(:istrln(fname))
      call camxerr()
c
 7004 continue
      write(iout,'(//,a)') 'ERROR in NCF_RDPTHDR_DDM:'
      write(iout,'(1X,2A,/,A)') 'Point source found in Source Apportionment ',
     &               'point source file is not found ',' in the base inventory.'
      write(iout,'(1X,2A)') '   File: ',fname(:istrln(fname))
      write(iout,'(1X,A,I10)') '   Point ID: ',ipt
      call camxerr()
c
 7005 continue
      write(iout,'(//,a)') 'ERROR in NCF_RDPTHDR_DDM:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot find variable id for: ',
     &                                      this_var(:istrln(this_var))
      write(iout,'(1X,2A)') 'File: ',fname(:istrln(fname))
      call camxerr()
c
 7006 continue
      write(iout,'(//,a)') 'ERROR in NCF_RDPTHDR_DDM:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot read data for variable: ',
     &                                      this_var(:istrln(this_var))
      write(iout,'(1X,2A)') 'File: ',fname(:istrln(fname))
      call camxerr()
c
 7007 continue
      write(iout,'(//,a)') 'ERROR in NCF_RDPTHDR_DDM:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'List of point sources does not match the ',
     &                                          'regular model file.'
      write(iout,'(A,I10,A)') 'Location of point source: ',ipt,
     &                                             ' does not match.'
      write(iout,'(1X,2A)') 'File: ',fname(:istrln(fname))
      call camxerr()
c
 7008 continue
      write(iout,'(//,a)') 'ERROR in NCF_RDPTHDR_DDM:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'PiG flag for sources in file is ',
     &                       'inconsistent with regular model file.'
      write(iout,'(A,I10,A)') 'Flag for point source: ',ipt,
     &                                             ' does not match.'
      write(iout,'(1X,2A)') 'File: ',fname(:istrln(fname))
      call camxerr()
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      deallocate( already_found )
      return
      end
