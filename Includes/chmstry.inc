c----CAMx v7.20 220430
c  
c     CHMSTRY.COM contains all chemistry variables 
c                            
c      Copyright 1996 - 2022
c     Ramboll
c            
c     Modifications:  
c       4/4/00     Added aerosol deposition variables to /aerochm/
c       1/9/02     Aerosol size cut points and density now a function of
c                  species
c       8/20/02    Added minimum CWC to define presence of clouds
c      12/12/02    Expanded species list for Mechanism 4
c       1/10/03    Added array for deposition output species names
c      10/14/04    Modifications for mechanism 10
c      11/30/06    Split into two include files to accamodate block-data
c      12/29/06    Expanded species list for the updated SOA scheme
c      01/08/06    Expanded species list for Mechanism 6 (CB05)
c                  Now ipigrxn is a scalar pointer to NO self-reaction
c      07/04/07    Added pointer to hydrolysis of N2O5
c      12/15/08    Added code to handle averaging of radicals
c      03/29/11    Support in-line TUV with aerosol optical depth
c      04/30/13    Added surface model
c       4/07/14    Added top con file
c       04/15/21   Added data structure to keep molecular weight
c 
c-----------------------------------------------------------------------
c    Variables for the number of species in input files:
c
c    narspc        --  number of species in the surface emissions file
c    nptspc        --  number of species in the point emissions file
c    nptsrc_files  --  number of points in each of the point emissions file
c    emspcname     --  names of emssions species 
c    idxems_files  --  index of species into master list for each surface emissions file
c    idxpnt_files  --  index of species into master list for each point emissions file
c    lgas          --  TRUE if species is a gas species
c    idx_start_pts --  index into start of global array for each file
c
c-----------------------------------------------------------------------
c
       integer,      allocatable, dimension(:,:)   :: narspc
       integer,      allocatable, dimension(:)     :: nptspc
       integer,      allocatable, dimension(:)     :: nptsrc_files
       character*10, allocatable, dimension(:)     :: emspcname
       integer,      allocatable, dimension(:,:,:) :: idxems_files
       integer,      allocatable, dimension(:,:)   :: idxpnt_files
       logical,      allocatable, dimension(:)     :: lgas
       integer,      allocatable, dimension(:)     :: idx_start_pts
c
c-----------------------------------------------------------------------
c     Variables for storing chemical reaction data:
c
c     rk     -- reaction rate constant (ppm/hr)
c     ltdep  -- flag to determine if rate constant is temperature dependent
c     lpdep  -- flag to determine if rate constant is pressure dependent
c     bdnl   -- lower vound value for each modeled species (ppm)
c
c-----------------------------------------------------------------------
c
      logical, allocatable, dimension(:) :: ltdep
      logical, allocatable, dimension(:) :: lpdep
      real,    allocatable, dimension(:) :: bdnl
c
      real    rk(10000)
      common /chmratep/ rk
c$omp threadprivate(/chmratep/)
c
c-----------------------------------------------------------------------
c     Variables for photolysis data:
c
c     idphot1  -- ID of primary photolysis reactions
c     idphot2  -- ID of secondary (scaled) photolysis reactions 
c     idphot3  -- ID of primary photolysis reaction to scale to obtain
c                 the secondary photolysis reaction
c     phtscl   -- photolysis reaction scaling factor
c
c-----------------------------------------------------------------------
c
      integer, allocatable, dimension(:) :: idphot1
      integer, allocatable, dimension(:) :: idphot2
      integer, allocatable, dimension(:) :: idphot3
      real,    allocatable, dimension(:) :: phtscl
c 
c-----------------------------------------------------------------------
c     Variables for species names:
c
c     spname  --  name of each modeled species
c     spavg   --  name of output species listed in control file
c     depsp   --  name of each deposition species output to file
c
c-----------------------------------------------------------------------
c
      character*10, allocatable, dimension(:) :: spname
      character*10, allocatable, dimension(:) :: spavg
      character*10, allocatable, dimension(:) :: depsp
c 
c-----------------------------------------------------------------------
c     Variables for mapping input species to internal model order:
c
c     lbcmap   -- mapping of species in the boundary condition file
c     ltcmap   -- mapping of species in the top concentrations file
c     lavmap   -- mapping of species written to average file
c     lemmap   -- mapping of species in the point source emissions file
c     licmap   -- mapping of species in the initial conditions file
c     ldepmap  -- mapping of species written to deposition file
c
c-----------------------------------------------------------------------
c
      integer, allocatable, dimension(:) :: lbcmap
      integer, allocatable, dimension(:) :: ltcmap
      integer, allocatable, dimension(:) :: lavmap
      integer, allocatable, dimension(:) :: lemmap
      integer, allocatable, dimension(:) :: licmap
      integer, allocatable, dimension(:) :: ldepmap
c
c-----------------------------------------------------------------------
c     Variables for chemistry lookup tables:
c
c     tempr  -- temperature table
c     presr  -- pressure table
c     rktbl  -- temperature/pressure-dependent rate constant table
c     htint  -- height AGL table
c     zenint -- zenith angle table
c     prkn   -- reaction rate table
c
c-----------------------------------------------------------------------
c      
      real, allocatable, dimension(:)           :: tempr
      real, allocatable, dimension(:)           :: presr
      real, allocatable, dimension(:)           :: htint
      real, allocatable, dimension(:)           :: zenint
      real, allocatable, dimension(:,:,:)       :: rktbl
      real, allocatable, dimension(:,:,:,:,:,:) :: prkn
c
c-----------------------------------------------------------------------
c     Variables to define parameters for each chemical species:
c
c     henry0       -- Henry's Law constant at STP (molar/atm)
c     tfact        -- Temperature dependence of Henry's Law constant (1/K)
c     diffrat      -- Species diffusivity
c     f0           -- Species reactivity parameter
c     rscale       -- Species scaling factor for surface resistance
c     roprt        -- Aerosol density (g/m3)
c     bext         -- Aerosol dry extinction efficiency (m2/ug)
c     ssa          -- Aerosol single scatter albedo
c     rhadj        -- Aerosol hygroscopic growth flag (1=apply RH adjustment)
c     dcut         -- Aerosol size bin cut points (um)
c     mole_weight  -- Molecular weight of each species
c
c-----------------------------------------------------------------------
c
      real, allocatable, dimension(:) :: henry0
      real, allocatable, dimension(:) :: tfact
      real, allocatable, dimension(:) :: diffrat
      real, allocatable, dimension(:) :: f0
      real, allocatable, dimension(:) :: rscale
      real, allocatable, dimension(:) :: mole_weight
c
      real,    allocatable, dimension(:)   :: roprt
      real,    allocatable, dimension(:,:) :: dcut
      real,    allocatable, dimension(:)   :: bext
      real,    allocatable, dimension(:)   :: ssa
      integer, allocatable, dimension(:)   :: rhadj
c
c-----------------------------------------------------------------------
c     Variables to define parameters for each surface model species:
c
c     smspc    -- Surface model species name
c     smpre    -- Precursor species name
c     smprd    -- Product species name
c     idsmsp   -- Pointer from surface model species to core species
c     idsmpre  -- Pointer from precursor species to surface model species
c     idsmprd  -- Pointer from product species to surface model species
c     smskrat  -- Reaction rate constant on soil (1/min)
c     smsjrat  -- Pholysis rate constant on soil (1/min)
c     smvkrat  -- Reaction rate constant on veg (1/min)
c     smvjrat  -- Pholysis rate constant on veg (1/min)
c     smikrat  -- Reaction rate constant on snow (1/min)
c     smijrat  -- Pholysis rate constant on snow (1/min)
c     smssrb   -- Soil equilib sorption ratio (unitless: soil/air)
c     smvsrb   -- Vegegation equilib sorption ratio (unitless: veg/air)
c     smisrb   -- Snow/ice equilib sorption ratio (unitless: ice/air)
c     smlch    -- Soil leaching rate (1/min)
c     smpen    -- Vegetation penetration rate (1/min)
c     smmlt    -- Snow/ice melt/penetration rate (1/min)
c
c-----------------------------------------------------------------------
c
      character*10, allocatable, dimension(:) :: smspc
      character*10, allocatable, dimension(:) :: smpre
      character*10, allocatable, dimension(:) :: smprd
      integer, allocatable, dimension(:) :: idsmsp
      integer, allocatable, dimension(:) :: idsmpre
      integer, allocatable, dimension(:) :: idsmprd
      real, allocatable, dimension(:) :: smskrat
      real, allocatable, dimension(:) :: smsjrat
      real, allocatable, dimension(:) :: smvkrat
      real, allocatable, dimension(:) :: smvjrat
      real, allocatable, dimension(:) :: smikrat
      real, allocatable, dimension(:) :: smijrat
      real, allocatable, dimension(:) :: smssrb
      real, allocatable, dimension(:) :: smvsrb
      real, allocatable, dimension(:) :: smisrb
      real, allocatable, dimension(:) :: smlch
      real, allocatable, dimension(:) :: smpen
      real, allocatable, dimension(:) :: smmlt
c
c----------------------------------------------------------------------
c    Variables for controlling calls to aerosol routines
c
c     time_aero  -- next time to call aerosol routines for each grid (HHMM)
c     date_aero  -- next date to call aerosol routines for each grid (YYJJJ)
c     aero_dt    -- incremented (actual) time (hr) between calls to aerosol 
c                   routines for each grid
c
      integer, allocatable, dimension(:) :: date_aero
      real,    allocatable, dimension(:) :: time_aero
      real,    allocatable, dimension(:) :: aero_dt
