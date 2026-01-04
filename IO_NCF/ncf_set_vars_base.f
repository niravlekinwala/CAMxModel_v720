c**** NCF_SET_VARS_BASE
c
      subroutine ncf_set_vars_base()
      implicit none
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routine sets the variable definitions and descriptions to 
c   the NetCDF file
c
c      Copyright 1996 - 2022
c     Ramboll
c      Argument description:
c       Inputs:
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
      include 'ncf_iodat.inc'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c    External functions:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
      ncf_x_units = "km"
      ncf_x_long_name = "X coordinate"
      ncf_x_var_desc = "X cartesian distance from projection origin"
      if( llatlon ) then
           ncf_x_units = "degrees"
           ncf_x_var_desc = "Longitude degrees east"
      endif
      ncf_y_units = "km"
      ncf_y_long_name = "Y coordinate"
      ncf_y_var_desc = "Y cartesian distance from projection origin"
      if( llatlon ) then
           ncf_y_units = "degrees"
           ncf_y_var_desc = "Longitude degrees north"
      endif
      ncf_layer_units = "Layer index"
      ncf_layer_long_name = "Model layer"
      ncf_layer_var_desc = "Model layer"
      ncf_longitude_units = "Degrees east"
      ncf_longitude_long_name = "Longitude"
      ncf_longitude_var_desc = "Longitude degrees east"
      ncf_longitude_coordinates = "latitude longitude"
      ncf_latitude_units = "Degrees north"
      ncf_latitude_long_name = "Latitude"
      ncf_latitude_var_desc = "Latitude degrees north"
      ncf_latitude_coordinates = "latitude longitude"
      ncf_topo_units = "m MSL"
      ncf_topo_long_name = "topographic elevation"
      ncf_topo_var_desc = "topographic elevation m above sea level"
      ncf_topo_coordinates = "latitude longitude"
      ncf_z_units = "m"
      ncf_z_long_name = "Layer height"
      ncf_z_var_desc = "Layer interface heights AGL"
      ncf_z_coordinates = "latitude longitude"
c
      ncf_tflag_units = "YYYYDDD,HHMMSS"
      ncf_tflag_long_name = "Start time flag"
      ncf_tflag_var_desc = "Timestep start date and time"
      ncf_etflag_units = "YYYYDDD,HHMMSS"
      ncf_etflag_long_name = "End time flag"
      ncf_etflag_var_desc = "Timestep end date and time"
c
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
 
