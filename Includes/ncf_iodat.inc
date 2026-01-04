c----CAMx v7.20 220430
c  
c     NC_IODAT.INC contains NCF file global attributes and metadata
c                            
c      Copyright 1996 - 2022
c     Ramboll
c            
c     Modifications:  
c        02/01/17   Origininal Development
c        07/23/18   Added ncf_bidi_nh3_drydep
c        04/28/21   Added strat_o3_profile
c
c-----------------------------------------------------------------------
c   Parameters for number of "base" variables: (appear in each file):
c-----------------------------------------------------------------------
c
c   NCF_BASE_VARS  -- number of the base variables
c
      integer NCF_BASE_VARS
c
      parameter (NCF_BASE_VARS = 9 )
c
c-----------------------------------------------------------------------
c   Dimension definitions:
c-----------------------------------------------------------------------
c
c   ncf_date_time       -- number of variables in time-stamp (always 2 - 
c                          YYYYDDD and HHMMSS)
c   ncf_date_time_dimid -- NetCDF dimension ID for this dimension
c   ncf_lay             -- number of vertical layers in modeling grid
c   ncf_lay_dimid       -- NetCDF dimension ID for this dimension
c   ncf_col             -- number of columns in modeling grid
c   ncf_col_dimid       -- NetCDF dimension ID for this dimension
c   ncf_row             -- number of rows in modeling grid
c   ncf_row_dimid       -- NetCDF dimension ID for this dimension
c   ncf_var             -- number of variables in data section
c   ncf_var_dimid       -- NetCDF dimension ID for this dimension
c   ncf_tstep_dimid     -- NetCDF dimension ID for TSTEP dimension
c
c-----------------------------------------------------------------------
c
      integer ncf_date_time
      integer ncf_date_time_dimid
      integer ncf_lay
      integer ncf_lay_dimid
      integer ncf_col
      integer ncf_col_dimid
      integer ncf_row
      integer ncf_row_dimid
      integer ncf_var
      integer ncf_var_dimid
      integer ncf_tstep_dimid

      common /ncf_dom_dat/ ncf_date_time, ncf_date_time_dimid, ncf_lay, 
     &                     ncf_lay_dimid, ncf_col, ncf_col_dimid,
     &                     ncf_row, ncf_row_dimid, ncf_var, ncf_var_dimid,
     &                     ncf_tstep_dimid
c
c-----------------------------------------------------------------------
c  Domain definition variables:
c-----------------------------------------------------------------------
c
c   ncf_x_units               -- units of coordinate values
c   ncf_x_long_name           -- long name of this variable
c   ncf_x_var_desc            -- description of this variable
c   ncf_y_units               -- units of coordinate values
c   ncf_y_long_name           -- long name of this variable
c   ncf_y_var_desc            -- description of this variable
c   ncf_layer_units           -- units of layer values
c   ncf_layer_long_name       -- long name of this variable
c   ncf_layer_var_desc        -- description of this variable
c   ncf_longitude_units       -- units of longitude values
c   ncf_longitude_long_name   -- long name of this variable
c   ncf_longitude_var_desc    -- description of this variable
c   ncf_longitude_coordnates  -- coordinates of this variable
c   ncf_latitude_units        -- units of longitude values
c   ncf_latitude_long_name    -- long name of this variable
c   ncf_latitude_var_desc     -- description of this variable
c   ncf_latitude_coordnates   -- coordinates of this variable
c   ncf_topo_units            -- units of longitude values
c   ncf_topo_long_name        -- long name of this variable
c   ncf_topo_var_desc         -- description of this variable
c   ncf_topo_coordnates       -- coordinates of this variable
c   ncf_z_units               -- units of layer heights
c   ncf_z_long_name           -- long name of this variable
c   ncf_z_var_desc            -- description of this variable
c
      character*20 ncf_x_units
      character*20 ncf_x_long_name
      character*60 ncf_x_var_desc
      character*20 ncf_y_units
      character*20 ncf_y_long_name
      character*60 ncf_y_var_desc
      character*20 ncf_layer_units
      character*20 ncf_layer_long_name
      character*60 ncf_layer_var_desc
      character*20 ncf_longitude_units
      character*20 ncf_longitude_long_name
      character*60 ncf_longitude_var_desc
      character*60 ncf_longitude_coordinates
      character*20 ncf_latitude_units
      character*20 ncf_latitude_long_name
      character*60 ncf_latitude_var_desc
      character*60 ncf_latitude_coordinates
      character*20 ncf_topo_units
      character*20 ncf_topo_long_name
      character*60 ncf_topo_var_desc
      character*60 ncf_topo_coordinates
      character*20 ncf_z_units
      character*20 ncf_z_long_name
      character*60 ncf_z_var_desc
      character*60 ncf_z_coordinates
c
      common /ncf_dom_chr/ ncf_x_units, ncf_x_long_name, ncf_x_var_desc,
     &                     ncf_y_units, ncf_y_long_name, ncf_y_var_desc,
     &                     ncf_layer_units, ncf_layer_long_name, 
     &                     ncf_layer_var_desc, ncf_longitude_units, 
     &                     ncf_longitude_long_name, ncf_longitude_var_desc, 
     &                     ncf_longitude_coordinates, ncf_latitude_units, 
     &                     ncf_latitude_long_name, ncf_latitude_var_desc, 
     &                     ncf_latitude_coordinates, ncf_topo_units, 
     &                     ncf_topo_long_name, ncf_topo_var_desc, 
     &                     ncf_topo_coordinates, ncf_z_units, ncf_z_long_name, 
     &                     ncf_z_var_desc, ncf_z_coordinates
c
c-----------------------------------------------------------------------
c  Time stamp variables:
c-----------------------------------------------------------------------
c
c   ncf_tflag             -- array of date/times for each output interval
c   ncf_tflag_units       -- units of the time stamp variable
c   ncf_tflag_long_name   -- long name of the time stamp variable
c   ncf_tflag_var_desc    -- description of the time stamp variable
c   ncf_tflag_varid       -- variable ID for the time stamp variable
c   ncf_etflag            -- array of date/times for each output interval
c   ncf_etflag_units      -- units of the time stamp variable
c   ncf_etflag_long_name  -- long name of the time stamp variable
c   ncf_etflag_var_desc   -- description of the time stamp variable
c   ncf_etflag_varid      -- variable ID for the time stamp variable
c
      integer,     allocatable, dimension(:,:) :: ncf_tflag
      character*20 ncf_tflag_units 
      character*20 ncf_tflag_long_name
      character*60 ncf_tflag_var_desc
      integer      ncf_tflag_varid
      integer,     allocatable, dimension(:,:) :: ncf_etflag
      character*20 ncf_etflag_units 
      character*20 ncf_etflag_long_name
      character*60 ncf_etflag_var_desc
      integer      ncf_etflag_varid

      common /ncf_time_chr/ ncf_tflag_units, ncf_tflag_long_name, 
     &                      ncf_tflag_var_desc, ncf_etflag_units, 
     &                      ncf_etflag_long_name, ncf_etflag_var_desc
c
      common /ncf_time_dat/ ncf_tflag_varid, ncf_etflag_varid
c
c-----------------------------------------------------------------------
c  Species name variables:
c-----------------------------------------------------------------------
c
c   ncf_species_units       -- units of species
c   ncf_species_long_name   -- long_name of species
c   ncf_species_var_desc    -- desciptions of species
c   ncf_species_coordinates -- desciptions of species
c   ncf_species_varid       -- array of variable ID for each species
c   spavnam                 --  array of names that will be written to the output files
c
      character*20 ncf_species_units 
      character*20 ncf_species_long_name
      character*60 ncf_species_var_desc
      character*60 ncf_species_coordinates
      character*10 spavnam(MXSPEC)
c
      common /ncf_spec_chr/ ncf_species_units, ncf_species_long_name,
     &                      ncf_species_var_desc, ncf_species_coordinates,
     &                      spavnam
c
c-----------------------------------------------------------------------
c  Global attributes for IO/API (domain definition):
c-----------------------------------------------------------------------
c
c   ncf_xorig     -- X-coordinate of domain origin
c   ncf_yorig     -- Y-coordinate of domain origin
c   ncf_xcell     -- cell width in X direction
c   ncf_ycell     -- cell width in Y direction
c   ncf_iutm      -- UTM zone of grid
c   ncf_istag     -- 
c   ncf_cproj     -- type of projection (CAMx codes)
c   ncf_gdtyp     -- tyoe of projection (IOAPI codes)
c   ncf_xcent     -- longitude pole
c   ncf_ycent     -- latitude pole
c   ncf_p_alp     -- tree latitude1
c   ncf_p_bet     -- tree latitude2
c   ncf_p_gam     -- same as longitude pole
c   ncf_nlays     -- number of layers
c   ncf_nrows     -- number of rows
c   ncf_ncols     -- number of columns
c   ncf_nthik     -- number of boundary layer cells
c   ncf_nvars     -- number of variables
c   ncf_nsteps    -- number of variables
c   ncf_cur_tstep -- current time step
c
      real    ncf_xorig
      real    ncf_yorig
      real    ncf_xcell
      real    ncf_ycell
      integer ncf_iutm
      integer ncf_istag
      integer ncf_cproj
      integer ncf_gdtyp
      real    ncf_xcent
      real    ncf_ycent
      real    ncf_p_alp
      real    ncf_p_bet
      real    ncf_p_gam
      integer ncf_nlays
      integer ncf_nrows
      integer ncf_ncols
      integer ncf_nthik
      integer ncf_nvars
      integer ncf_nsteps
      integer ncf_cur_tstep
c
      common /ncf_glob_dom_dat/ ncf_xorig, ncf_yorig, ncf_xcell, ncf_ycell,
     &                          ncf_iutm, ncf_istag, ncf_cproj, ncf_gdtyp,
     &                          ncf_xcent, ncf_ycent, ncf_p_alp, ncf_p_bet,
     &                          ncf_p_gam, ncf_nlays, ncf_nrows, ncf_ncols,
     &                          ncf_nthik, ncf_nvars, ncf_nsteps, ncf_cur_tstep
c
c-----------------------------------------------------------------------
c  Global attributes for IO/API (file description):
c-----------------------------------------------------------------------
c
c   ncf_cdate              -- date of file creation (YYYYJJJ)
c   ncf_ctime              -- time of file creation
c   ncf_wdate              -- date of file (YYYYJJJ)
c   ncf_wtime              -- time of file 
c   ncf_name               -- CAMx file type
c   ncf_note               -- simulation note
c   ncf_itzon              -- time zone of simulation
c   ncf_ftype              -- IO-API file type (always 1 = CUSTOM3)
c   ncf_vgtyp              -- IO-API grid type (always 6 = H: m above ground) 
c   ncf_vgtop              -- IO-API grid tio for sigma coordinates (always -9999 = missing)
c   ncf_vglvls             -- IO-API levels 
c   ncf_gdnam              -- IO-API grid name (always CAMx vX.YY)
c   ncf_upnam              -- IO-API last program writing this file (always CAMx vX.YY)
c   ncf_ioapi_ver          -- IO-API version - needs to exit for M3TOOLS and must start with IOAPI-
c   ncf_execid             -- name of program that created the file (needed for M3TOOLS)
c   ncf_filedesc           -- CAMx file tyoe (same as ncf_name)
c   ncf_sdate              -- date of simulation start (YYYYJJJ)
c   ncf_stime              -- time of simulation start
c   ncf_tstep              -- time step (hours times HHMMSS)
c   ncf_grid_id            -- CAMx grid ID
c   ncf_i_grid_start       -- I-cell index of start of grid in master 
c   ncf_i_grid_end         -- I-cell index of dns of grid in master
c   ncf_j_grid_start       -- J-cell index of start of grid in master
c   ncf_j_grid_end         -- J-cell index of dns of grid in master
c   ncf_grid_mesh_factor   -- meshing factor for grid in master
c   ncf_flexi_nest         -- flexi-nesting flag (0 = false, 1 = true) 
c   ncf_advection          -- CAMx advection scheme
c   ncf_chem_solver        -- name of CAMx chemical mechanism
c   ncf_pig                -- CAMx PiG flag
c   ncf_probing_tool       --  CAMx probing tools flag
c   ncf_chemistry          -- CAMx chemistry flag (0 = false, 1 = true)
c   ncf_total_species      -- total species in mechanism
c   ncf_radical_species    -- number of radical species in mechanism
c   ncf_gas_species        -- number of gas species in mechanism
c   ncf_pm_species         -- number of PM species in mechanism
c   ncf_reactions          -- number of reactions in mechanism
c   ncf_drydep             -- CAMx dryp deposition flag
c   ncf_wetdep             -- CAMx wet deposition flag (0 = false, 1 = true)
c   ncf_acm2               -- CAMx ACM2 option (0 = false, 1 = true)
c   ncf_cig_model          -- CAMx cloud in grid option (0 = false, 1 = true)
c   ncf_surface_model      -- chemical surface model flag (0 = false, 1 = true)
c   ncf_inline_ix_emiss    -- inline emissions for halogens flag (0 = false, 1 = true)
c   ncf_bidi_nh3_drydep    -- bi-di nh3 drydep flag (0 = false, 1 = true)
c   ncf_strat_o3_profile   -- stratospheric ozone profile flag (0 = false, 1 = true)
c   ncf_super_stepping     -- super stepping flag (0 = false, 1 = true)
c   ncf_gridded_emiss      -- flag to indicate if gridded emissions are used (0 = false, 1 = true)
c   ncf_point_emiss        --  flag to indicate if point source emissions are used (0 = false, 1 = true)
c   ncf_ignore_emiss_dates -- flag to ignore dates on emissions files (0 = false, 1 = true)
c   ncf_output_3d          -- flag to output 3D concentrations (0 = false, 1 = true)
c   ncf_pig_sample_grid    -- flag to use PiG sampling grid (0 = false, 1 = true)
c   ncf_pig_sample_bckgnd  -- flag to use grid concs in PiG sampleing grid (0 = false, 1 = true)
c   ncf_pig_sample_grid_id -- ID of CAMx grid containing PiG sampling grid
c   ncf_i_sample_start     -- I-cell index of start of PiG sampling grid
c   ncf_i_sample_end       -- I-cell index of end of PiG sampling grid
c   ncf_j_sample_start     -- J-cell index of start of PiG sampling grid
c   ncf_j_sample_end       -- J-cell index of end of PiG sampling grid
c   ncf_sample_mesh_factor -- meshing factor of PiG sampling grid
c   ncf_conventions        -- NetCDF convnetions in this file
c   ncf_history            -- history of file access
c
      integer       ncf_cdate
      integer       ncf_ctime
      integer       ncf_wdate
      integer       ncf_wtime
      character*10  ncf_name
      character*60  ncf_note
      integer       ncf_itzon
      integer       ncf_ftype
      integer       ncf_vgtyp
      real          ncf_vgtop
      real          ncf_vglvls(MXLAYER+1)
      character*16  ncf_gdnam
      character*16  ncf_upnam
      character*16  ncf_ioapi_ver
      character*10  ncf_execid
      character*10  ncf_filedesc
      integer       ncf_sdate
      integer       ncf_stime
      integer       ncf_tstep
      integer       ncf_grid_id
      integer       ncf_i_grid_start
      integer       ncf_i_grid_end
      integer       ncf_j_grid_start
      integer       ncf_j_grid_end
      integer       ncf_grid_mesh_factor
      integer       ncf_flexi_nest
      character*10  ncf_advection
      character*10  ncf_chem_solver
      character*10  ncf_pig
      character*10  ncf_probing_tool
      character*140 ncf_chemistry
      integer       ncf_total_species
      integer       ncf_radical_species
      integer       ncf_gas_species
      integer       ncf_pm_species
      integer       ncf_reactions
      character*10  ncf_drydep
      integer       ncf_wetdep
      integer       ncf_acm2
      integer       ncf_cig_model
      integer       ncf_surface_model
      integer       ncf_inline_ix_emiss
      integer       ncf_bidi_nh3_drydep
      integer       ncf_strat_o3_profile
      integer       ncf_super_stepping
      integer       ncf_gridded_emiss
      integer       ncf_point_emiss
      integer       ncf_ignore_emiss_dates
      integer       ncf_output_3d
      integer       ncf_pig_sample_grid
      integer       ncf_pig_sample_bckgnd
      integer       ncf_pig_sample_grid_id
      integer       ncf_i_sample_start
      integer       ncf_i_sample_end
      integer       ncf_j_sample_start
      integer       ncf_j_sample_end
      integer       ncf_sample_mesh_factor
      character*10  ncf_conventions
      character*60  ncf_history

      common /ncf_glob_chr/ ncf_name, ncf_note, ncf_gdnam, ncf_upnam, ncf_ioapi_ver, 
     &                      ncf_execid, ncf_filedesc, ncf_advection, ncf_chem_solver, 
     &                      ncf_pig, ncf_probing_tool, ncf_chemistry, ncf_drydep, 
     &                      ncf_conventions, ncf_history

      common /ncf_glob_dat/ ncf_cdate, ncf_ctime, ncf_wdate, ncf_wtime, 
     &                      ncf_itzon, ncf_ftype, ncf_vgtyp, ncf_vgtop, ncf_vglvls,
     &                      ncf_sdate, ncf_stime, ncf_tstep, ncf_grid_id, 
     &                      ncf_i_grid_start, ncf_i_grid_end, ncf_j_grid_start, 
     &                      ncf_j_grid_end, ncf_grid_mesh_factor, ncf_flexi_nest, 
     &                      ncf_total_species, ncf_radical_species, ncf_gas_species,
     &                      ncf_pm_species, ncf_reactions, ncf_wetdep, ncf_acm2, 
     &                      ncf_cig_model, ncf_surface_model, ncf_inline_ix_emiss, 
     &                      ncf_bidi_nh3_drydep,ncf_strat_o3_profile,
     &                      ncf_super_stepping, ncf_gridded_emiss, ncf_point_emiss, 
     &                      ncf_ignore_emiss_dates, ncf_output_3d, 
     &                      ncf_pig_sample_grid, ncf_pig_sample_bckgnd,
     &                      ncf_pig_sample_grid_id, ncf_i_sample_start, ncf_i_sample_end,
     &                      ncf_j_sample_start, ncf_j_sample_end, ncf_sample_mesh_factor
c
c-----------------------------------------------------------------------
c  Variables for chunking parameters:
c-----------------------------------------------------------------------
c
c   NCF_CACHESIZE            -- total size of the raw data chunk cache in MegaBytes
c   NCF_NELEMS               -- number slots in the per-variable chunk cache 
c                              (should be a prime number larger than the number of chunks in the cache)
c   NCF_PREEMPTION           -- indicates how much chunks that have been fully read are favored for
c                               preemption (0 means fully read chunks are treated no differently 
c                               than other chunks, while a value of 100 means fully read chunks 
c                               are always preempted before other chunks
c   NCF_SHUFFLE              -- If non-zero, turn on the shuffle filter
c   NCF_DEFLATE              -- If non-zero, turn on the deflate filter at the level specified by the 
c                               deflate_level parameter
c   NCF_DEFLATE_LEVEL        -- 0 = no deflate and 9 = slowest, but best deflate
c   NCF_CHUNK_SIZE_VAR_X  -- number of chunks in X direction
c   NCF_CHUNK_SIZE_VAR_Y  -- number of chunks in Y direction
c  
      integer NCF_CACHESIZE
      integer NCF_NELEMS
      integer NCF_PREEMPTION
      integer NCF_SHUFFLE
      integer NCF_DEFLATE
      integer NCF_DEFLATE_LEVEL
      integer NCF_CHUNK_SIZE_VAR_X
      integer NCF_CHUNK_SIZE_VAR_Y
c
      parameter( NCF_CACHESIZE        = 32 )
      parameter( NCF_NELEMS           = 37 )
      parameter( NCF_PREEMPTION       = 100 )
      parameter( NCF_SHUFFLE          = 0 )
      parameter( NCF_DEFLATE          = 1 )
      parameter( NCF_DEFLATE_LEVEL    = 2 )
      parameter( NCF_CHUNK_SIZE_VAR_X = 2 )
      parameter( NCF_CHUNK_SIZE_VAR_Y = 2 )
