c**** NCF_WRT_VARS_BASE
c
      subroutine ncf_wrt_vars_base(action,iounit)
      use filunit
      use grid
      implicit none
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routine writes the variable definitions and descriptions to 
c    the NetCDF file
c
c      Copyright 1996 - 2022
c     Ramboll
c      Argument description:
c       Inputs:
c           action  C  name of file to open
c           iounit  I  NetCDF file ID of file
c       Outputs:
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     02/20/17   --gwilson--    Original development
c     03/15/21   --gwilson--    Changed TOPO to float from double
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'netcdf.inc'
      include 'ncf_iodat.inc'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      character*(*) action
      integer       iounit
      integer       nspcs
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
      integer ierr, grid_dimid(2), z_dimid(4), time_dimid(3), this_varid
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
      grid_dimid(1) = ncf_col_dimid
      grid_dimid(2) = ncf_row_dimid
c
      z_dimid(1) = ncf_col_dimid
      z_dimid(2) = ncf_row_dimid
      z_dimid(3) = ncf_lay_dimid
      z_dimid(4) = ncf_tstep_dimid
c
      time_dimid(1) = ncf_date_time_dimid
      time_dimid(2) = ncf_var_dimid
      time_dimid(3) = ncf_tstep_dimid
c
      ierr = nf_def_var(iounit, "X", NF_DOUBLE, 1, 
     &                                      ncf_col_dimid, this_varid)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_text(iounit,this_varid,'units',
     &                                    istrln(ncf_x_units),ncf_x_units)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_att_text(iounit,this_varid,'long_name',
     &                              istrln(ncf_x_long_name),ncf_x_long_name)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_att_text(iounit,this_varid,'var_desc',
     &                                istrln(ncf_x_var_desc),ncf_x_var_desc)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_def_var(iounit, "Y", NF_DOUBLE, 1, 
     &                                      ncf_row_dimid, this_varid)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_att_text(iounit,this_varid,'units',
     &                                    istrln(ncf_y_units),ncf_y_units)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_att_text(iounit,this_varid,'long_name',
     &                            istrln(ncf_y_long_name),ncf_y_long_name)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_att_text(iounit,this_varid,'var_desc',
     &                             istrln(ncf_y_var_desc),ncf_y_var_desc)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_def_var(iounit, "layer", NF_INT, 1, 
     &                                ncf_lay_dimid, this_varid)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_att_text(iounit,this_varid,'units',
     &                           istrln(ncf_layer_units),ncf_layer_units)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_att_text(iounit,this_varid,'long_name',
     &                   istrln(ncf_layer_long_name),ncf_layer_long_name)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_att_text(iounit,this_varid,'var_desc',
     &                     istrln(ncf_layer_var_desc),ncf_layer_var_desc)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_def_var(iounit, "TFLAG", NF_INT, 3, 
     &                                  time_dimid, this_varid)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_att_text(iounit,this_varid,'units',
     &                          istrln(ncf_tflag_units),ncf_tflag_units)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_att_text(iounit,this_varid,'long_name',
     &                  istrln(ncf_tflag_long_name),ncf_tflag_long_name)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_att_text(iounit,this_varid,'var_desc',
     &                     istrln(ncf_tflag_var_desc),ncf_tflag_var_desc)
      if( ierr .NE. NF_NOERR ) goto 7000

      ierr = nf_def_var(iounit, "ETFLAG", NF_INT, 3, 
     &                                 time_dimid, this_varid)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_att_text(iounit,this_varid,'units',
     &                        istrln(ncf_etflag_units),ncf_etflag_units)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_att_text(iounit,this_varid,'long_name',
     &                istrln(ncf_etflag_long_name),ncf_etflag_long_name)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_att_text(iounit,this_varid,'var_desc',
     &                  istrln(ncf_etflag_var_desc),ncf_etflag_var_desc)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_def_var(iounit, "longitude", NF_DOUBLE, 2, 
     &                              grid_dimid, this_varid)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_att_text(iounit,this_varid,'units',
     &                  istrln(ncf_longitude_units),ncf_longitude_units)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_att_text(iounit,this_varid,'long_name',
     &          istrln(ncf_longitude_long_name),ncf_longitude_long_name)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_att_text(iounit,this_varid,'var_desc',
     &             istrln(ncf_longitude_var_desc),ncf_longitude_var_desc)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_att_text(iounit,this_varid,'coordinates',
     &        istrln(ncf_longitude_coordinates),ncf_longitude_coordinates)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_def_var(iounit, "latitude", NF_DOUBLE, 2, 
     &                                grid_dimid, this_varid)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_att_text(iounit,this_varid,'units',
     &                   istrln(ncf_latitude_units),ncf_latitude_units)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_att_text(iounit,this_varid,'long_name',
     &          istrln(ncf_latitude_long_name),ncf_latitude_long_name)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_att_text(iounit,this_varid,'var_desc',
     &             istrln(ncf_latitude_var_desc),ncf_latitude_var_desc)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_att_text(iounit,this_varid,'coordinates',
     &       istrln(ncf_latitude_coordinates),ncf_latitude_coordinates)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_def_var(iounit, "topo", NF_FLOAT, 2, 
     &                                grid_dimid, this_varid)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_att_text(iounit,this_varid,'units',
     &                   istrln(ncf_topo_units),ncf_topo_units)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_att_text(iounit,this_varid,'long_name',
     &          istrln(ncf_topo_long_name),ncf_topo_long_name)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_att_text(iounit,this_varid,'var_desc',
     &             istrln(ncf_topo_var_desc),ncf_topo_var_desc)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_att_text(iounit,this_varid,'coordinates',
     &       istrln(ncf_topo_coordinates),ncf_topo_coordinates)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_def_var(iounit, "z", NF_FLOAT, 4, 
     &                                         z_dimid, this_varid)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_att_text(iounit,this_varid,'units',
     &                                 istrln(ncf_z_units),ncf_z_units)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_att_text(iounit,this_varid,'long_name',
     &                        istrln(ncf_z_long_name),ncf_z_long_name)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_att_text(iounit,this_varid,'var_desc',
     &                           istrln(ncf_z_var_desc),ncf_z_var_desc)
      if( ierr .NE. NF_NOERR ) goto 7000
      ierr = nf_put_att_text(iounit,this_varid,'coordinates',
     &                     istrln(ncf_z_coordinates),ncf_z_coordinates)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,a)') 'ERROR in NCF_WRT_VARS_BASE:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(A)') 'Cannot create file variable for grid variables.'
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
 
