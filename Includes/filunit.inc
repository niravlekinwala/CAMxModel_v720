c----CAMx v7.20 220430
c 
c     FILUNIT.COM contains all model I/O unit numbers
c                           
c      Copyright 1996 - 2022
c     Ramboll Environ
c           
c     Modifications: 
c        8/30/02    Cloud file now contains rain record, rain unit removed.
c                   Cloud/rain and water vapor files can be read in for all
c                   nests
c        1/10/03    Added deposition file units
c        10/12/04   Added file root name as common variable
c        8/23/06    Instantaneous restart files reduced to 1 per grid type
c        8/25/06    Average and deposition output files now one file per grid
c        11/4/09    Removed input top concentrations
c        1/04/11    Revised for new met input format
c        04/30/13   Added surface model
c        4/07/14    Added top con file
c        09/02/14   Added subgrid convective model
c
c-----------------------------------------------------------------------
c     File units for ASCII and Fortran binary I/O files:
c
c     icur_unit       -- unit number of the current file being opened
c     iout            -- message output file
c     idiag           -- diagnostic output file
c     imass           -- mass summary output file
c     iconc           -- coarse grid instantaneous concentration output file
c     ifconc          -- fine grid instantaneous concentration output file
c     ipig            -- PiG output file
c     ichem           -- chemistry parameters input file
c     iphot           -- photolysis lookup input file
c     iic             -- initial conditions input file
c     is_netcdf_iic   -- flag for NetCDF format for 3D Met file
c     ibc             -- boundary conditions input file
c     is_netcdf_ibc   -- flag for NetCDF format for 3D Met file
c     itc             -- top concentrations input file
c     is_netcdf_itc   -- flag for NetCDF format for 3D Met file
c     io3col          -- ozone column input file
c     irstc           -- coarse grid restart input file
c     irstf           -- fine grid restart input file
c     irstp           -- PiG restart file
c     filroot         -- Root output file name
c     npoint_files    -- number of point source files to read
c     nlayers_ems     -- number of layers in emissions files (1 or nlay)
c-----------------------------------------------------------------------
c
      character*200 filroot
      integer       icur_unit
      integer       iout
      integer       idiag
      integer       imass
      integer       iconc
      integer       ifconc
      integer       ipig
      integer       ichem
      integer       iphot
      integer       iic
      integer       ibc
      integer       itc
      integer       io3col
      integer       irstc
      integer       irstf
      integer       irstp
      logical       is_netcdf_iic
      logical       is_netcdf_ibc
      logical       did_update_ibc
      logical       is_netcdf_itc
      logical       did_update_itc
      integer       npoint_files
      integer       nlayers_ems
c
      common /funit/ icur_unit, iout, idiag, imass, iconc, ifconc,
     &               ipig, ichem, iphot, iic, ibc, itc,
     &               io3col, irstc, irstf, irstp,
     &               is_netcdf_iic, is_netcdf_ibc, is_netcdf_itc,
     &               did_update_ibc, did_update_itc, 
     &               npoint_files, nlayers_ems
      common /fuchr/ filroot
c
c-----------------------------------------------------------------------
c     Variables for trackingcurrent  timestep for NetCDF files
c
c     bnddate_time_tflag -- current tflag value for BC file
c     enddate_time_tflag -- current etflag value for BC file
c     topdate_time_tflag -- current tflag value for TC file
c     topdate_time_tflag -- current etflag value for TC file
c
      real bnddate_time_tflag    
      real bnddate_time_etflag    
      real topdate_time_tflag    
      real topdate_time_etflag    
c
      common /tflag_dat/ bnddate_time_tflag, bnddate_time_etflag,
     &                   topdate_time_tflag, topdate_time_etflag
c
c========================= Process Analysis Begin ==============================
c
c     ipr_unit  -- unit number of Integrated Process Rates output file
c     irr_unit  -- unit number of Integrated Reaction Rates output file
c
      integer ipr_unit
      integer irr_unit
c
      common /paunit/ ipr_unit, irr_unit
c
c========================= Process Analysis End ==============================
c
c
c-----------------------------------------------------------------------
c     File units for Fortran binary I/O files:
c     These are dynamically allocated.
c
c     iavg                 -- average concentration output file
c     idep                 -- deposition output file
c     ismin                -- surface model input file
c     ismout               -- surface model output file
c     iarem                -- area emission input file
c     is_netcdf_iarem      -- flag for NetCDF format for 3D Met file
c     iptem                -- point source input file
c     buffer_offset_iarem  -- beginning I cell of grid to read in file 
c     isurf                -- 2D Surface input file
c     is_netcdf_isurf      -- flag for NetCDF format for 3D Met file
c     i3dmet               -- 3D Met input file
c     is_netcdf_i3dmet     -- flag for NetCDF format for 3D Met file
c     i2dmet               -- 2D Met input file
c     is_netcdf_i2dmet     -- flag for NetCDF format for 3D Met file
c     ikv                  -- 3D vertical diffusivity input file
c     is_netcdf_ikv        -- flag for NetCDF format for 3D Met file
c     icld                 -- 3D cloud/rain input file
c     isample              -- 2D PiG sampling grid output file
c     lcig                 -- Subgrid convection flag by grid
c     nemiss_files         -- number of gridded emissions files for each grid
c-----------------------------------------------------------------------
c
      integer, allocatable, dimension(:)   :: iavg
      integer, allocatable, dimension(:)   :: idep
      integer, allocatable, dimension(:)   :: ismin
      integer, allocatable, dimension(:)   :: ismout
      integer, allocatable, dimension(:,:) :: iarem
      logical, allocatable, dimension(:,:) :: is_netcdf_iarem
      integer, allocatable, dimension(:,:) :: buffer_offset_iarem
      integer, allocatable, dimension(:)   :: iptem
      logical, allocatable, dimension(:)   :: is_netcdf_iptem
      integer, allocatable, dimension(:)   :: isurf
      logical, allocatable, dimension(:)   :: is_netcdf_isurf
      integer, allocatable, dimension(:)   :: i3dmet
      logical, allocatable, dimension(:)   :: is_netcdf_i3dmet
      integer, allocatable, dimension(:)   :: i2dmet
      logical, allocatable, dimension(:)   :: is_netcdf_i2dmet
      integer, allocatable, dimension(:)   :: ikv
      logical, allocatable, dimension(:)   :: is_netcdf_ikv
      integer, allocatable, dimension(:)   :: icld
      integer, allocatable, dimension(:)   :: isample
      logical, allocatable, dimension(:)   :: lcig
      integer, allocatable, dimension(:)   :: nemiss_files
c
c-----Number of variables on the met input files
c
c     n3dmet   -- 3D Met file
c     n2dmet   -- 2D Met file
c     nkvmet   -- 3D Kv file
c     ncldmet  -- 3D Cloud file
c     nsrfvar  -- 2D Surface file
c
      integer, allocatable, dimension(:) :: n3dmet
      integer, allocatable, dimension(:) :: n2dmet
      integer, allocatable, dimension(:) :: nkvmet
      integer, allocatable, dimension(:) :: ncldmet
      integer, allocatable, dimension(:) :: nsrfvar
c
      integer itag
      common /tagdat/ itag
