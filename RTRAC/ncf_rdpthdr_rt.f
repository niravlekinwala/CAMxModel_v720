      subroutine ncf_rdpthdr_rt(begtim,begdate,endtim,enddate)
      use filunit
      use ptemiss
      use tracer
      implicit none
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c     NCF_RDPTHDR_RT reads the header information from the global section
c     of the NetCDF point source file. This version is for RTRAC.
c
c      Copyright 1996 - 2022
c     Ramboll
c
c      Argument description:
c       Inputs:
c         begtim  R model begin time
c         begdat  I model begin date (YYJJJ)
c         endtim  R model end time
c         enddate I model end date (YYJJJ)
c       Outputs:
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     06/13/19   --gwilson--    Original development
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
      real    begtim
      integer begdate
      real    endtim
      integer enddate
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
      integer       ipt, iounit, idxpt_base, itmp
c
      real, allocatable, dimension(:) :: txloc
      real, allocatable, dimension(:) :: tyloc
      real, allocatable, dimension(:) :: tstkhght
      real, allocatable, dimension(:) :: tstkdiam
      real, allocatable, dimension(:) :: tstktemp
      real, allocatable, dimension(:) :: vstktemp
      real, allocatable, dimension(:) :: iarray
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
c   --- loop over number of files ---
c
      iounit = iortpt(1,1)
      if( .NOT. is_netcdf_iortpt(1,1) ) goto 9999
      write(action,'(2A,I3)')  'Reading NetCDF RTRAC point source file.'
      fname = tptfil(1,1)
c
c --- get the type of file to make sure it is a point source file ---
c
      this_var = 'CAMx_NAME'
      ierr = nf_get_att_text(iounit, NF_GLOBAL, this_var, name_in)
      if( ierr .NE. NF_NOERR ) then
         this_var = 'NAME'
         ierr = nf_get_att_text(iounit, NF_GLOBAL, this_var, name_in)
         if( ierr .NE. NF_NOERR ) goto 7009
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
      call ncf_set_emsrt_mapping(iounit,action)
c
c --- get the number of points (used as dimension for columns) ---
c
      ierr = nf_inq_dimid(iounit, "COL", this_dimid )
      if( ierr .NE. NF_NOERR ) goto 7002
      ierr = nf_inq_dimlen(iounit,this_dimid,numpts)
      if( ierr .NE. NF_NOERR ) goto 7002
c
c  --- check that parameters are the same ---
c
      allocate( txloc(numpts) )
      allocate( tyloc(numpts) )
      allocate( tstkhght(numpts) )
      allocate( tstkdiam(numpts) )
      allocate( tstktemp(numpts) )
      allocate( vstktemp(numpts) )
      allocate( iarray(numpts) )
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
      ierr = nf_get_var_real(iounit,this_varid,vstktemp)
      if( ierr .NE. NF_NOERR ) goto 7006
      do ipt=1,numpts
        vstktemp(ipt) = vstktemp(ipt)/3600.
      enddo
c
c --- get pig flag ---
c
      this_var = "pigflag"
      ierr = nf_inq_varid(iounit,this_var,this_varid)
      if( ierr .NE. NF_NOERR ) goto 7005
      ierr = nf_get_var_int(iounit,this_varid,iarray)
      if( ierr .NE. NF_NOERR ) goto 7006
      do ipt=1,numpts
c
c  --- make sure this source is in the base inventory and flag it ---
c
            idx_point_in_list(1,1,ipt) = idxpt_base
            if( iarray(ipt) .GT. 0 ) lpigsa(idxpt_base) = .TRUE.
         endif
      enddo
      deallocate( txloc )
      deallocate( tyloc )
      deallocate( tstkhght )
      deallocate( tstkdiam )
      deallocate( tstktemp )
      deallocate( vstktemp )
      deallocate( iarray )
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,a)') 'ERROR in NCF_RDPTHDR_RT:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot get necessary global attribute: ',
     &                                      this_var(:istrln(this_var))
      write(iout,'(1X,2A)') 'File: ',fname(:istrln(fname))
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7001 continue
      write(iout,'(//,a)') 'ERROR in NCF_RDPTHDR_RT:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(A)') 'Input file is not the correct type.'
      write(iout,'(2A)') 'Looking for type: ',pntfil(:istrln(pntfil))
      write(iout,'(2A)') 'Found in file   : ',name_in(:istrln(name_in))
      write(iout,'(1X,2A)') 'File: ',fname(:istrln(fname))
      call camxerr()
c
 7002 continue
      write(iout,'(//,a)') 'ERROR in NCF_RDPTHDR_RT:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(A)') 'Cannot find the dimension id for number ',
     &                                      'of point sources (NCOL)'
      write(iout,'(1X,2A)') 'File: ',fname(:istrln(fname))
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7003 continue
      write(iout,'(//,a)') 'ERROR in NCF_RDPTHDR_RT:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(A)') 'Cannot find the dimension value for number ',
     &                                      'of point sources (NCOL)'
      write(iout,'(1X,2A)') 'File: ',fname(:istrln(fname))
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7004 continue
      write(iout,'(//,a)') 'ERROR in NCF_RDPTHDR_RT:'
      write(iout,'(1X,2A,/,A)') 'Point source found in Source Apportionment ',
     &               'point source file is not found ',' in the base inventory.'
      write(iout,'(1X,2A)') '   File: ',fname(:istrln(fname))
      write(iout,'(1X,A,I10)') '   Point ID: ',ipt
      call camxerr()
c
 7005 continue
      write(iout,'(//,a)') 'ERROR in NCF_RDPTHDR_RT:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot find variable id for: ',
     &                                      this_var(:istrln(this_var))
      write(iout,'(1X,2A)') 'File: ',fname(:istrln(fname))
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7006 continue
      write(iout,'(//,a)') 'ERROR in NCF_RDPTHDR_RT:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot read data for variable: ',
     &                                      this_var(:istrln(this_var))
      write(iout,'(1X,2A)') 'File: ',fname(:istrln(fname))
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7007 continue
      write(iout,'(//,a)') 'ERROR in NCF_RDPTHDR_RT:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'List of point sources does not match the ',
     &                                          'regular model file.'
      write(iout,'(A,I10,A)') 'Location of point source: ',ipt,
     &                                             ' does not match.'
      write(iout,'(1X,2A)') 'File: ',fname(:istrln(fname))
      call camxerr()
c
 7008 continue
      write(iout,'(//,a)') 'ERROR in NCF_RDPTHDR_RT:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'PiG flag for sources in file is ',
     &                       'inconsistent with regular model file.'
      write(iout,'(A,I10,A)') 'Flag for point source: ',ipt,
     &                                             ' does not match.'
      write(iout,'(1X,2A)') 'File: ',fname(:istrln(fname))
      call camxerr()
c
 7009 continue
      write(iout,'(//,A)') 'ERROR in NCF_RDPTHDR_RT:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)',ERR=9999) 'Cannot find global variable for ',
     &                                       'type of file: CAMx_NAME'
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
