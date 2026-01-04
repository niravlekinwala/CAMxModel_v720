**** TRACER
c
c----CAMx v7.20 220430
c
c
c----------------------------------------------------------------------
c
c    Include file for tracer parameters and data structures used
c    in the source apportionment version of the CAMx
c
c      Copyright 1996 - 2022
c     Ramboll Environ
c
c-----------------------------------------------------------------------
c   LOG:
c-----------------------------------------------------------------------
c
c     01/03/96   -gwilson-  original development
c     02/10/97   -cemery-   added SA output file root name
c     04/28/97   -gwilson-  added OPPAT and GOAT
c     05/22/97   -gwilson-  added APCA
c     11/17/99   -gwilson-  added DDM
c     01/17/02   -gwilson-  added RTRAC
c     07/31/02   -gwilson-  removed OPPAT
c     09/20/03   -gwilson-  changed the way pointers are stored to 
c                           support many more classes in PSAT
c     07/20/05   -cemery-   Moved most PiG sampling grid variables
c                           to PIGSTY.COM, moved sampling grid parameters
c                           to CAMX.PRM, moved sampling grid file
c                           variables to FILUNIT.COM
c     08/25/05   -cemery-   Revamped PiG/OSAT region/group pointers
c     8/23/06    -cemery-   Instantaneous restart files reduced to 1 per
c                           grid type
c     03/15/09   -gwilson-  Added code for deposition output for tracers
c     8/25/06    -cemery-   Surface output files now one file per grid
c     12/29/06   -bkoo-     Revised for the updated SOA scheme
c     07/11/07   -cemery-   added RTCMC
c     07/16/07   -bkoo-     Revised for HDDM
c                           Added HRVOC
c     06/11/08   -bkoo-     Added rate constant sensitivity
c     07/16/08   -bkoo-     Added DDM turn-off flag
c     03/15/09   -gwilson-  Added code for deposition output for tracers
c     10/29/09   -cemery-   added RTRAC surface model
c     11/4/09    -cemery-   Removed input top concentrations
c     11/12/09   -gwilson-  Added factor for applying new type of top
c                           boundary conditions.
c     11/06/12   -gwilson-  fixed Wall of Cells recrptors for MPI
c     08/23/13   -bkoo-     Added DDM-PM parameters
c     12/20/13   -gwilson-  Added the compact point source file
c     03/18/14   -bkoo-     Added parameters for benzene SOA
c     04/20/15   -gwilson-  Added partial source area map
c     09/11/15   -bkoo-     Revised for SA v3
c     08/25/16   -bkoo-     Updated PSAT & DDM-PM parameters for new SOAP
c     11/09/16   -cemery-   Added Baker APCA point source override option
c     07/23/18   -bkoo-     Added parameter for Bi-Di NH3 source group
c     08/09/18   -bkoo-     Added rate term sensitivity
c     10/26/18   -bkoo-     Updated JDP indexes for new PM-DDM species
c     01/06/19   -cemery-   Added PFE/PMN/PMG/PK/PCA/PAL/PSI/PTI
c     01/09/19   -cemery-   Added DMS
c     04/19/21   -gwilson-  Addded molecular weight for probing tools
c
c-----------------------------------------------------------------------
c    Parameters for Source Apportionment tech type:
c-----------------------------------------------------------------------
c
c
c   VERSA    C  character string for the version of the model
c   VERDDM   C  character string for the version of the model
c   VERHDDM  C  character string for the version of the model
c   VERPA    C  character string for the version of the model
c   VERRTRAC C  character string for the version of the model
c   VERRTCMC C  character string for the version of the model
c
      character*80 VERSA
      character*80 VERDDM
      character*80 VERHDDM
      character*80 VERPA
      character*80 VERRTRAC
      character*80 VERRTCMC
c
      parameter(VERSA=                'Source Apportionment, SA 220530')
      parameter(VERDDM=          'De-coupled Direct Method, DDM 220530')
      parameter(VERHDDM=
     &                'High-order Decoupled Direct Method, HDDM 220530')
      parameter(VERPA=                    'Process Anaylsis, PA 220530')
      parameter(VERRTRAC=              'Reactive Tracers, RTRAC 220530')
      parameter(VERRTCMC=              'Reactive Tracers, RTCMC 220530')
c
c-----------------------------------------------------------------------
c    Parameters and variables for I/O units:
c-----------------------------------------------------------------------
c
c    MXTEMF   I  maximum number of source grouping emission files
c    IDXCRS   I  index of course grid in arrays
c    IDXFIN   I  index of fine grid in arrays
c
      integer   MXTEMF
      integer   IDXCRS
      integer   IDXFIN
c
      parameter( MXTEMF = 100 )
      parameter( IDXCRS = 1 )
      parameter( IDXFIN = 2 )
c
c    iormap                 I unit number of the source area mapping file
c    iorrcp                 I unit number of the file with the receptor definition
c    ioric                  I unit number of the IC file for DDM
c    lsa_ioric              L .TRUE. if using an SA IC file
c    num_ioric	            I number of total species in SA IC file
c    spc_ioric              C names of species in SA IC file
c    ncls_ioric	            I number of species per class in SA IC file
c    is_netcdf_ioric        L flag for NetCDF format for DDM IC file
c    iorbc                  I unit number of the BC file for DDM
c    lsa_iorbc              L .TRUE. if using an SA BC file
c    num_iorbc	            I number of total species in SA tC file
c    spc_iorbc              C names of species in SA IC file
c    ncls_iorbc	            I number of species per class in SA IC file
c    is_netcdf_iorbc        L flag for NetCDF format for DDM BC file
c    sa_bndtim              R current time of SA boundary condition file timestep
c    sa_bnddate             I current date of SA boundary condition file timestep
c    lrd_sa_bc_hdr          L .TRUE. if BC header must be read
c    iortc                  I unit number of the BC file for DDM
c    lsa_iortc              L .TRUE. if using an SA TC file
c    ncls_iortc	            I number of species per class in SA TC file
c    lrd_sa_tc_hdr          L .TRUE. if TC header must be read
c    is_netcdf_iortc        L flag for NetCDF format for DDM TC file
c    iorchm                 I unit mumber of the chemistry parameters for RTRAC/RTCMC
c    iorini                 I unit number of instantaneous file used to initalize
c    iortem                 I unit number of the source grouping files (surface)
c    is_netcdf_iortem       L flag for NetCDF format for surface emissions
c    buffer_offset_iortem   I beginning I cell for grid to read in file
c    iortpt                 I unit number of the source grouping files (elevated)
c    is_netcdf_iortpt       L flag for NetCDF format for point emissions
c    lcompactpt             L flag to determine if PTSOURCE file is compact
c    iowcon                 I unit number of the first output instantaneous file
c    iowsfc                 I unit number of the surface tracer concentrations file
c    iowptdep               I unit number of the tracer deposition outputs
c    iowrcp                 I unit number of output receptor concentrations
c    iowsmp                 I unit number of the RTRAC/RTCMC sampling grid concentration file
c    iorrtsrf               I unit number of the input RTRAC surface model mass file
c    is_netcdf_iorrtsrf     L flag for NetCDF format for point emissions
c    iowrtsrf               I unit number of the output RTRAC surface model mass file
c    num_iortem             I number of emissions files for each group/grid
c    num_iortpt             I number of points source files for each group
c
      integer, allocatable, dimension(:,:)   :: iormap
      integer, allocatable, dimension(:,:)   :: iortpt
      logical, allocatable, dimension(:,:)   :: lcompactpt
      integer, allocatable, dimension(:,:,:) :: iortem
      integer, allocatable, dimension(:)     :: iowsfc
      integer, allocatable, dimension(:)     :: iowptdep
      integer, allocatable, dimension(:,:)   :: num_iortem
      integer, allocatable, dimension(:)     :: num_iortpt
c
      integer, allocatable, dimension(:)     :: iowsmp
      integer, allocatable, dimension(:)     :: iorrtsrf
      integer, allocatable, dimension(:)     :: iowrtsrf
c
      character*10, allocatable, dimension(:) :: spc_ioric
      character*10, allocatable, dimension(:) :: spc_iorbc
c
      logical, allocatable, dimension(:,:)   :: is_netcdf_iortpt
      logical, allocatable, dimension(:,:,:) :: is_netcdf_iortem
      logical, allocatable, dimension(:)     :: is_netcdf_iorrtsrf
c
      integer, allocatable, dimension(:,:,:) :: buffer_offset_iortem
c
      integer iorini(IDXCRS:IDXFIN)
      integer iowcon(IDXCRS:IDXFIN)
      integer iowrcp
      integer ioric
      logical lsa_ioric
      integer num_ioric
      integer ncls_ioric
      integer iorbc
      logical lsa_iorbc
      integer num_iorbc
      integer ncls_iorbc
      real    sa_bndtim
      integer sa_bnddate
      integer iortc
      logical lsa_iortc
      integer ncls_iortc
      integer iorchm
      logical is_netcdf_ioric
      logical is_netcdf_iorbc
      logical is_netcdf_iortc
      logical lrd_sa_bc_hdr
      logical lrd_sa_tc_hdr
c
      integer iorrcp
c
      common /ptiount/ iorrcp, iorini, ioric, lsa_ioric, num_ioric, 
     &                 ncls_ioric, iorbc, lsa_iorbc, num_iorbc, 
     &                 ncls_iorbc, sa_bndtim, sa_bnddate, iortc, 
     &                 lsa_iortc, ncls_iortx, iorchm, iowcon, iowrcp, 
     &                 is_netcdf_ioric, is_netcdf_iorbc, is_netcdf_iortc, 
     &                 lrd_sa_bc_hdr, lrd_sa_tc_hdr
c
c-----------------------------------------------------------------------
c    Variables for filenames:
c-----------------------------------------------------------------------
c
c   flrtsa   C   output file root name
c   mapfil   C   filename for the source region map
c   lmapfl   L   flag to dertimine if source region map has been provided
c   rcpfil   C   filename for the receptor definition file
c   lrcpfil  L   flag to determine if receptor file is to be used
c   icfil    C   filename for the IC file for DDM or RTRAC/RTCMC
c   licfil   L   flag to determine if IC file was supplied for RTRAC/RTCMC 
c   bcfil    C   filename for the BC file for DDM or RTRAC/RTCMC
c   lbcfil   L   flag to determine if BC file was supplied for RTRAC/RTCMC 
c   ltcfil   L   flag to determine if TC file was supplied for RTRAC/RTCMC 
c   chmfil   C   filename for the chemistry parameters file for RTRAC/RTCMC
c   inifil   C   filename for the instantaneous file used for initializing
c                (one for each grid, zero element is for course grid)
c   temfil   C   array of filenames for source group emissions (surface)
c   ltemfl   L   flag to determine if emissions file was supplied for 
c                the source grouping (position 0 is for the regular
c                emisions file)
c   tptfil   C   array of filenames for source group emissions (elevated)
c   ltptfl   L   flag to determine if emissions file was supplied for 
c                the source grouping (position 0 is for the regular
c                emisions file)
c   confil   C   filename for output instantaneous concentrations
c   sfcfil   C   filename for the surface tracer concentration file
c   ptdepfil C   filename for the tracer deposition output file
c   lsfcfl   L   flag to generate a surface tracer output file
c   avgfil   C   filename for the average tracer concentrations
c   smpfil   C   filename for the RTRAC/RTCMC sampling grid concentration file
c   rtsrfin  C   filename for the input RTRAC surface mass file
c   rtsrfout C   filename for the output RTRAC surface mass file
c   lfirst   L   flag for firs time interval of average concs
c   verson   C   string containing the version of the model being run
c
      character*200, allocatable, dimension(:,:)   :: mapfil
      logical,       allocatable, dimension(:,:)   :: lmapfl
      character*200, allocatable, dimension(:,:,:) :: temfil
      logical,       allocatable, dimension(:,:,:) :: ltemfl
      character*200, allocatable, dimension(:,:)   :: tptfil
      logical,       allocatable, dimension(:,:)   :: ltptfl
      character*200, allocatable, dimension(:)     :: smpfil
      character*200, allocatable, dimension(:)     :: rtsrfin
      character*200, allocatable, dimension(:)     :: rtsrfout
      character*200,  allocatable, dimension(:)    :: sfcfil
      character*200,  allocatable, dimension(:)    :: ptdepfil
c
      character*200 flrtsa
      character*200 rcpfil
      logical       lrcpfil
      character*200 icfil
      logical       licfil
      character*200 bcfil
      logical       lbcfil
      character*200 tcfil
      logical       ltcfil
      character*200 chmfil
      character*200 inifil(IDXCRS:IDXFIN)
      character*200 confil(IDXCRS:IDXFIN)
      logical       lsfcfl
      character*200 avgfil
      logical       lfirst
      character*80  verson
c     
      common /filchr/ flrtsa, rcpfil, icfil, bcfil, tcfil, chmfil, 
     &                inifil, confil, avgfil, verson 
      common /fildat/ lrcpfil, licfil, lbcfil, ltcfil, lsfcfl, lfirst
c
c-----------------------------------------------------------------------
c    Parameters for releasing timing tracers:
c-----------------------------------------------------------------------
c
c   TIMEMS   R   concentration (PPM) of timing tracer release emissions
c
      real      TIMEMS
c
      parameter( TIMEMS = 10.0E-07 )
c
c-----------------------------------------------------------------------
c    Parameters for array indices:
c-----------------------------------------------------------------------
c
c    MXCLTYP  I   maximum number class types for SA
c    MXTRCLS  I   maximum number of tracer classes actually used
c    MXALCLS  I   number of tracer classes supported
c    ITRVOC   I   index of VOC emissions
c    ITRO3N   I   index of NOx ozone reactivity
c    ITRO3V   I   index of VOC ozone reactivity
c    ITROON   I   index of odd-oxygen in RGN from O3N
c    ITROOV   I   index of odd-oxygen in RGN from O3V
c    ITRSO2   I   index of SO2 emissions
c    ITRPS4   I   index of particulate sulfate
c    ITRDMS   I   index of DMS
c    ITRNIT   I   index of reactive gas nitrogen w/o odd-oxygen
c    ITRRGN   I   index of reactive gas nitrogen w/ odd-oxygen
c    ITRTPN   I   index of gaseous peroxyl acetyl nitrate
c    ITRNTR   I   index of organic nitrates (RNO3)
c    ITRHN3   I   index of gaseous nitric acid
c    ITRPN3   I   index of particulate nitrate from primary nitrate + RNO3
c    ITRNH3   I   index of gaseous amonnia
c    ITRPN4   I   index of particulate ammonium
c    ITRARO   I   index of Aromatic Secondary Organic Aerosol Precursors
c    ITRISP   I   index of Isoprene Secondary Organic Aerosol Precursors
c    ITRTRP   I   index of Monoterpene Secondary Organic Aerosol Precursors
c    ITRSQT   I   index of Sesquiterpene Secondary Organic Aerosol Precursors
c    ITRCG1   I   index of condensable gas associated with the CG1/SOA1 pair
c    ITRCG2   I   index of condensable gas associated with the CG2/SOA2 pair
c    ITRCG3   I   index of condensable gas associated with the CG3/SOA3 pair
c    ITRCG4   I   index of condensable gas associated with the CG4/SOA4 pair
c    ITRPO1   I   index of SOA associated with the CG1/SOA1 pair
c    ITRPO2   I   index of SOA associated with the CG2/SOA2 pair
c    ITRPO3   I   index of SOA associated with the CG3/SOA3 pair
c    ITRPO4   I   index of SOA associated with the CG4/SOA4 pair
c    ITRPPA   I   index of polymerized SOA (anthropogenic)
c    ITRPPB   I   index of polymerized SOA (biogenic)
c    ITRPEC   I   index of primary elemenatal carbon
c    ITRPOA   I   index of primary organic carbon
c    ITRPFC   I   index of primal fine crustal PM
c    ITRPFN   I   index of other primary fine particulate
c    ITRPCC   I   index of primary course crustal particulate
c    ITRPCS   I   index of primary course particulate
c    ITRPFE   I   index of primary Fe
c    ITRPMN   I   index of primary Mn
c    ITRPMG   I   index of primary Mg
c    ITRPK    I   index of primary K
c    ITRPCA   I   index of primary Ca
c    ITRPAL   I   index of primary Al
c    ITRPSI   I   index of primary Si
c    ITRPTI   I   index of primary Ti
c    ITRHG0   I   index of primary elemental mercury
c    ITRHG2   I   index of primary and secondary oxygenated mercury
c    ITRPHG   I   index of primary particulate mercury
c
      integer MXTRCLS
      integer MXALCLS
      integer ITRVOC
      integer ITRO3N
      integer ITRO3V
      integer ITROON
      integer ITROOV
      integer ITRSO2
      integer ITRPS4
      integer ITRDMS
      integer ITRNIT
      integer ITRRGN
      integer ITRTPN
      integer ITRNTR
      integer ITRHN3
      integer ITRPN3
      integer ITRNH3
      integer ITRPN4
      integer ITRARO
      integer ITRISP
      integer ITRTRP
      integer ITRSQT
      integer ITRCG1
      integer ITRCG2
      integer ITRCG3
      integer ITRCG4
      integer ITRPO1
      integer ITRPO2
      integer ITRPO3
      integer ITRPO4
      integer ITRPPA
      integer ITRPPB
      integer ITRPEC
      integer ITRPOA
      integer ITRPFC
      integer ITRPFN
      integer ITRPCC
      integer ITRPCS
      integer ITRPFE
      integer ITRPMN
      integer ITRPMG
      integer ITRPK
      integer ITRPCA
      integer ITRPAL
      integer ITRPSI
      integer ITRPTI
      integer ITRHG0
      integer ITRHG2
      integer ITRPHG
c
      parameter( MXTRCLS  = 47 )
      parameter( MXALCLS  = 47 )
      parameter( ITRVOC   =  1 )
      parameter( ITRO3N   =  2 )
      parameter( ITRO3V   =  3 )
      parameter( ITROON   =  4 )
      parameter( ITROOV   =  5 )
      parameter( ITRSO2   =  6 )
      parameter( ITRPS4   =  7 )
      parameter( ITRDMS   =  8 )
      parameter( ITRNIT   =  9 )
      parameter( ITRRGN   = 10 )
      parameter( ITRTPN   = 11 )
      parameter( ITRNTR   = 12 )
      parameter( ITRHN3   = 13 )
      parameter( ITRPN3   = 14 )
      parameter( ITRNH3   = 15 )
      parameter( ITRPN4   = 16 )
      parameter( ITRARO   = 17 )
      parameter( ITRISP   = 18 )
      parameter( ITRTRP   = 19 )
      parameter( ITRSQT   = 20 )
      parameter( ITRCG1   = 21 )
      parameter( ITRCG2   = 22 )
      parameter( ITRCG3   = 23 )
      parameter( ITRCG4   = 24 )
      parameter( ITRPO1   = 25 )
      parameter( ITRPO2   = 26 )
      parameter( ITRPO3   = 27 )
      parameter( ITRPO4   = 28 )
      parameter( ITRPPA   = 29 )
      parameter( ITRPPB   = 30 )
      parameter( ITRPEC   = 31 )
      parameter( ITRPOA   = 32 )
      parameter( ITRPFC   = 33 )
      parameter( ITRPFN   = 34 )
      parameter( ITRPCC   = 35 )
      parameter( ITRPCS   = 36 )
      parameter( ITRPFE   = 37 )
      parameter( ITRPMN   = 38 )
      parameter( ITRPMG   = 39 )
      parameter( ITRPK    = 40 )
      parameter( ITRPCA   = 41 )
      parameter( ITRPAL   = 42 )
      parameter( ITRPSI   = 43 )
      parameter( ITRPTI   = 44 )
      parameter( ITRHG0   = 45 )
      parameter( ITRHG2   = 46 )
      parameter( ITRPHG   = 47 )
c
c-----------------------------------------------------------------------
c    Variables for gridded tracer emissions data:
c-----------------------------------------------------------------------
c
c   saemis    R   gridded array of tracer emissions from the surface
c                 emissions files
c   pttop     R   gridded array of tracer top concentrations
c   xlocpt    R   location of the point source in X direction
c   ylocpt    R   location of the point source in Y direction
c   ntrtim    I   number of times per day that a timing release is done
c   lreles    L   flag for determing if a new timing tracer should
c                 be released
c   nreles    I   number of the current timing release
c   ipigmap   I   region in which each PiG source is located
c   ipiggrp   I   group in which each PiG source is located
c   lpigsa    L   flag to determine if source is a PiG source
c                   
      integer, allocatable, dimension(:)   :: ipigmap
      integer, allocatable, dimension(:)   :: ipiggrp
      real,    allocatable, dimension(:)   :: saemis
      real,    allocatable, dimension(:)   :: pttop
      real,    allocatable, dimension(:)   :: xlocpt
      real,    allocatable, dimension(:)   :: ylocpt
      logical, allocatable, dimension(:)   :: lpigsa
c
      integer   nreles
      integer   ntrtim
      logical   lreles
      logical   lptoverride
c
      common /empdat/  lreles, nreles, ntrtim, lptoverride
c
c-----------------------------------------------------------------------
c    Parameters for gridded tracer concentration data:
c-----------------------------------------------------------------------
c
c   BNDLPT    R   lower bound for tracer concentrations
c
      real      BNDLPT
c
      parameter( BNDLPT = 1.0E-16 )
c
      real       EPSYLD ! very small non-zero number to ensure all the SOA precursors counted
      parameter( EPSYLD = 1.0E-25 )
c
c-----------------------------------------------------------------------
c    Variables for gridded tracer concentration data:
c-----------------------------------------------------------------------
c
c   ptconc    R   gridded array of tracer concentrations by layer(SA)
c                 or gridded array of sensitivities by layer (DDM)
c   ptavrg    R   average tracer concentrations at surface (used for
c                 hourly peak receptor)
c   ptvdep    R   diffusion velocities for tracer species
c   senrad    R   gridded array of radical sensitivities by layer
c   sapnts    R   array of tracer emissions for each point source
c
      real, allocatable, dimension(:)   :: ptconc
      real, allocatable, dimension(:)   :: ptavrg
      real, allocatable, dimension(:)   :: ptvdep
      real, allocatable, dimension(:)   :: senrad
      real, allocatable, dimension(:,:) :: sapnts
c
c-----------------------------------------------------------------------
c    Variables for tracer names:
c-----------------------------------------------------------------------
c
c   ptname    C    name of passive tracer species (imbedded into this name
c                 will be the type of tracer, the source region from 
c                 where it originated and the source grouping emissions
c                 file from which it originated)
c   ptnameout C    name of tracer species that will be written to file
c
c   ptdepname C    name of SA deposition species
c   nsaspc    I    number of tracer species currently in the simulation
c   npttim    I    number of timing tracer species (including future timings)
c   ntotsp    I    number of total tracer species or number of species 
c                 times number of parameters (DDM)
c   ntotspout I    number of total tracer species written to output file
c   ntotspdep I    number of total tracer deposition species written to output file
c   notimespc I    number of "regular" tracers (not including timing tracers)
c   ntrcls    I    number of tracer classes actually used
c   idxcls    I    index into the master list of classes
c   idxipt    I    index into the list of classes actually used
c   iptcls    I    array of indexes into tracer list of the beginning of each
c                  class of concentration tracers
c   nptcls    I    array of indexes into tracer list of the end of each
c                 class of concentration tracers
c   ipttrc    I    array of indexes into tracer list of the beginning of each
c                 class of concentration tracers (for O3 this is O3V+O3N)
c   npttrc    I    array of indexes into tracer list of the end of each
c                 class of concentration tracers (for O3 this is O3V+O3N)
c   iemcls    I    array of indexes into tracer list of the beginning of each
c                 class of emissions tracers
c   nemcls    I    array of indexes into tracer list of the end of each
c                 class of emissions tracers
c   ipttim    I    index into tracer list of the beginning of the
c                 timing tracers
c   iemtim    I    index into tracer list of the begginning of the 
c                 timing emissions tracers
c   nbdic     I    number of boundary/initial conditions tracer species
c                 (depends on wether the boundary is stratified)
c   lsamap    I    dummy mapping of arrays species (lsamap(i) = i)
c   lsagas    L    true if tracer species is a gaseous species
c   trspmap   I    this is a sparse 2-D array that contains the coefficients
c                 for linear combination of modeled species that contribute
c                 to tracer classes -- for emissions and concentrations
c   fluxmap   I    this is a sparse 2-D array that contains the coefficients
c                 for linear combination of modeled species that contribute
c                 to tracer classes -- for advection fluxes
c   yhratmap  I    high-NOx yield rate map for each modeled species in each class
c   ylratmap  I    low-NOx yield rate map for each modeled species in each class
c   lusespc   L    flag to indicate if species is used in probing tools
c   clsnam    C    names of the tracer classes
c   ptop_fac  R    factor to be applied to the tracers to create a new
c                  type of top boundary condition 
c
c   sa_mole_weight R   "molecular weight" of SA species - used for output 
c                     in ug/m^3
c
      character*10, allocatable, dimension(:)   :: ptname
      character*14, allocatable, dimension(:)   :: ptnameout
      character*10, allocatable, dimension(:)   :: ptdepname
      integer,      allocatable, dimension(:)   :: lsamap
      logical,      allocatable, dimension(:)   :: lsagas
      real,         allocatable, dimension(:,:) :: trspmap
      real,         allocatable, dimension(:,:) :: fluxmap
      real,         allocatable, dimension(:,:) :: yhratmap
      real,         allocatable, dimension(:,:) :: ylratmap
      real,         allocatable, dimension(:)   :: sa_mole_weight
      logical,      allocatable, dimension(:)   :: lusespc
      real,         allocatable, dimension(:)   :: ptop_fac
c
      character*3  clsnam(MXALCLS)
      integer      iptcls(MXALCLS)
      integer      nptcls(MXALCLS)
      integer      ipttrc(MXALCLS)
      integer      npttrc(MXALCLS)
      integer      idxcls(MXALCLS)
      integer      idxipt(MXALCLS)
      integer      iemcls(MXALCLS)
      integer      nemcls(MXALCLS)
      integer      nsaspc
      integer      ntrcls
      integer      npttim
      integer      ntotsp
      integer      ntotspout
      integer      ntotspdep
      integer      notimespc
      integer      ipttim 
      integer      iemtim 
      integer      nbdic
c
      common /ptnchr/ clsnam
      common /ptndat/ nsaspc, ntrcls, iptcls, nptcls, ipttrc, npttrc, 
     &                idxcls, idxipt, npttim, ntotsp, ntotspout, ntotspdep,
     &                iemcls, nemcls, ipttim, iemtim, nbdic, notimespc
c
c-----------------------------------------------------------------------
c    Parameters for user options and flags:
c-----------------------------------------------------------------------
c
c   SA     C   code for SA technology
c   DDM    C   code for DDM technology
c   HDDM   C   code for HDDM technology
c   RTRAC  C   code for RTRAC technology
c   RTCMC  C   code for RTCMC technology
c
      character*10 SA
      character*10 DDM
      character*10 HDDM
      character*10 RTRAC
      character*10 RTCMC
c
      parameter( SA    = 'SA        ')
      parameter( DDM   = 'DDM       ')
      parameter( HDDM  = 'HDDM      ')
      parameter( RTRAC = 'RTRAC     ')
      parameter( RTCMC = 'RTCMC     ')
c
c-----------------------------------------------------------------------
c    Parameters for average output type:
c-----------------------------------------------------------------------
c
c    ALLOUT    C  string for indicating all species are output
c    SUMMARY   C  string for indicating that just summary species
c                 are output to the average file
c
      character*10 ALLOUT
      character*10 SUMMARY
c
      parameter( ALLOUT  = 'ALL' )
      parameter( SUMMARY = 'SUMMARY' )
c
c-----------------------------------------------------------------------
c    Parameters for types of SA:
c-----------------------------------------------------------------------
c
c    SULFATE   C  string for indicating Sulfate species
c    NITRATE   C  string for indicating Nitrate species
c    SOA       C  string for indicating secondary aerosal species
c    PRIMARY   C  string for indicating primary aerosal species
c    MERCURY   C  string for indicating Mercury species
c    OZONE     C  string for indicating Ozone species
c
      character*10 SULFATE
      character*10 NITRATE
      character*10 SOA
      character*10 PRIMARY
      character*10 MERCURY
      character*10 OZONE
c
      parameter( SULFATE = 'SULFATE' )
      parameter( NITRATE = 'NITRATE' )
      parameter( SOA     = 'SOA' )
      parameter( PRIMARY = 'PRIMARY' )
      parameter( MERCURY = 'MERCURY' )
      parameter( OZONE   = 'OZONE' )
c
c-----------------------------------------------------------------------
c    Variables for turning parts of PSAT:
c-----------------------------------------------------------------------
c
c   lsulfate   L   .TRUE. if doing SULFATE species
c   lnitrate   L   .TRUE. if doing NITRATE species
c   lsoa       L   .TRUE. if doing SOA species
c   lprimary   L   .TRUE. if doing PRIMARY species
c   lmercury   L   .TRUE. if doing MERCURY species
c   lozone     L   .TRUE. if doing Ozone species
c   lapca      L   .TRUE. if using APCA
c   lapcapt    L   .TRUE. if using APCA with optional point source override
c   lmineral   L   .TRUE. if MINERAL species are to be added to PRIMARY group
c   ldmschm    L   .TRUE. if DMS is to be added to SULFATE group
c
      logical lsulfate
      logical lnitrate
      logical lsoa
      logical lprimary
      logical lmercury
      logical lozone
      logical lapca
      logical lapcapt
      logical lmineral
      logical ldmschm
c
      common /psdat/ lsulfate, lnitrate, lsoa, lprimary, lmercury,
     &               lozone, lapca, lapcapt, lmineral, ldmschm
c
c-----------------------------------------------------------------------
c    Variables for user options and flags:
c-----------------------------------------------------------------------
c
c   ltrace   L   flag for determining if the passive tracer algorithm 
c                should be used
c   lddm     L   flag for determining if DDM is being used 
c   lhddm    L   flag for determining if HDDM is being used
c   lddmcalc L   flag for determining if (H)DDM sens are being calculated
c                in a grid
c   lrestrt  L   flag for determining if the simulation is a first day
c   leftovr_area  L   flag for determining if the left-over group should be 
c                used
c   lbndry   L   flag to determine if the boundary conditions should
c                be stratified by edge
c   ngroup   I   number of source groupings for emissions
c   tectyp   C   flag for determining which type of technology will
c                be performed.
c   lallout  L   flag for determining if all species are output to the 
c                average file, or just the summary
c   loutsa   L   flag for determining if the species should
c                be output to average file
c   lsa_3davrg L .TRUE. means write all layers to average file
c
      logical, allocatable, dimension(:) :: loutsa
      logical, allocatable, dimension(:) :: lddmcalc
c
      integer      ngroup
      integer      nchar
      logical      ltrace
      logical      lddm
      logical      lhddm
      logical      lrestrt
      logical      leftovr_area
      logical      lbndry
      character*10 tectyp
      logical      lallout
      logical      lsa_3davrg
      logical      luse_points
c
      common /usrchr/ tectyp
      common /usrdat/ ngroup, nchar, ltrace, lrestrt, leftovr_area, lbndry,
     &                lddm, lhddm, lallout, lsa_3davrg, luse_points
c
c-----------------------------------------------------------------------
c    Variables for region mapping:
c-----------------------------------------------------------------------
c
c   nregin   I   number of source regions 
c   nxcell   I   number of cells in X-direction
c   nycell   I   number of cells in X-direction
c   igrmap   I   grid that maps cell to source region
c   frcmap   R   fraction of cell in each region (for partial source
c                area map)
c
      integer, allocatable, dimension(:)         :: nxcell
      integer, allocatable, dimension(:)         :: nycell
      integer, allocatable, dimension(:,:,:,:,:) :: igrmap
      real,    allocatable, dimension(:,:,:,:,:) :: frcmap
c
      integer   nregin
c
      common /mapdat/ nregin
c
c-----------------------------------------------------------------------
c    Parameter for Bi-Di NH3 source group:
c-----------------------------------------------------------------------
c
c   IGRPBIDI I   source group to allocate Bi-Di NH3 re-emissions
c
      integer, parameter :: IGRPBIDI = 1
c
c-----------------------------------------------------------------------
c    Parmeters for receptor variables:
c-----------------------------------------------------------------------
c
c   MXRECP   I   maximum number of receptor locations
c   MXCELR   I   maximum number of cells in CELL AVERAGE type of receptor
c   CDPNT    C   string for indicating POINT type of recptor
c   IDPNT    I   id code for POINT type of receptor
c   CDCEL    C   string for indicating SINGLE CELL type of recptor
c   IDCEL    I   id code for SINGLE CELL type of receptor
c   CDAVG    C   string for indicating CELL AVERAGE type of recptor
c   IDAVG    I   id code for CELL AVERAGE type of receptor
c   CDWAL    I   string indicating WALL OF CELLS type of receptors
c   IDWAL    I   id code for WALL OF CELLS type of receptors
c
      character*20 CDPNT
      character*20 CDCEL
      character*20 CDAVG
      character*20 CDWAL
      integer      MXRECP
      integer      MXCELR
      integer      IDPNT
      integer      IDCEL
      integer      IDAVG
      integer      IDWAL
c
      parameter( CDPNT  = 'POINT          ' )
      parameter( CDCEL  = 'SINGLE CELL    ' )
      parameter( CDAVG  = 'CELL AVERAGE   ' )
      parameter( CDWAL  = 'WALL OF CELLS  ' )
      parameter( MXRECP = 100 )
      parameter( MXCELR = 70  )
      parameter( IDPNT  = 1   )
      parameter( IDCEL  = 2   ) 
      parameter( IDAVG  = 3   ) 
      parameter( IDWAL  = 4   ) 
c
c-----------------------------------------------------------------------
c    Variables for receptor data:
c-----------------------------------------------------------------------
c
c   rcpnam    C   names of the receptors
c   nrecep    I   number of receptors
c   idrcp     I   array of id codes for each receptor
c   igrdrcp   I   grid index for each receptor
c   irecep    I   array of I-cell locations for CELL type of receptors
c   jrecep    I   array of I-cell locations for CELL type of receptors
c   nclrcp    I   number of cells to average for CELL type of receptors
c   recepx    R   array of X-ccordinates for POINT type of receptors
c   recepy    R   array of Y-ccordinates for POINT type of receptors
c   conrcp    R   array of tracer surface concentrations at each receptor
c   volrcp    R   array of total volume for WALL OF CELLS receptor
c   ipekcl    I   index into gridded array of the hourly peak cell
c   iwalbg    I   beginning column for WALL OF CELLS receptor
c   iwalnd    I   ending column for WALL OF CELLS receptor
c   jwalbg    I   beginning row for WALL OF CELLS receptor
c   jwalnd    I   ending row for WALL OF CELLS receptor
c   kwalbg    I   beginning layer for WALL OF CELLS receptor
c   kwalnd    I   ending layer for WALL OF CELLS receptor
c   lwalls    L   flag to determine if any WALL OF CELLS receptors 
c                 were specified
c
      real, allocatable, dimension(:,:) :: conrcp
      real, allocatable, dimension(:)   :: volrcp
c
      character*10 rcpnam(MXRECP)
      integer      idrcp(MXRECP)
      integer      igrdrcp(MXRECP)
      integer      irecep(MXRECP,MXCELR)
      integer      jrecep(MXRECP,MXCELR)
      integer      nclrcp(MXRECP)
      integer      iwalbg(MXRECP)
      integer      iwalnd(MXRECP)
      integer      jwalbg(MXRECP)
      integer      jwalnd(MXRECP)
      integer      kwalbg(MXRECP)
      integer      kwalnd(MXRECP)
      real         recepx(MXRECP)
      real         recepy(MXRECP)
      integer      nrecep
      logical      lwalls
c
      common /rcpchr/ rcpnam
      common /rcpdat/ nrecep, idrcp, igrdrcp, irecep, jrecep, nclrcp, 
     &                recepx, recepy, iwalbg, iwalnd, 
     &                jwalbg, jwalnd, kwalbg, kwalnd, lwalls
c
c-----------------------------------------------------------------------
c    Parameters for boundary conditions:
c-----------------------------------------------------------------------
c
c   IDXBWS   I   index of the WEST boundary in arrays
c   IDXBES   I   index of the EAST boundary in arrays
c   IDXBST   I   index of the SOUTH boundary in arrays
c   IDXBNT   I   index of the NORTH boundary in arrays
c   IDXBTP   I   index of the TOP boundary in arrays
c
      integer   IDXBWS
      integer   IDXBES
      integer   IDXBST
      integer   IDXBNT
      integer   IDXBTP
c
      parameter( IDXBWS = 1 )
      parameter( IDXBES = 2 )
      parameter( IDXBST = 3 )
      parameter( IDXBNT = 4 )
      parameter( IDXBTP = 5 )
c
c-----------------------------------------------------------------------
c    Variables for gridded deposition fields:
c-----------------------------------------------------------------------
c
c   ptdryfld  R   gridded array of dry deposition mass
c   ptwetfld  R   gridded array of wet deposition mass
c   lptdepout L   flag to determine if probing tools depostion is on
c
      real, allocatable, dimension(:) :: ptdryfld
      real, allocatable, dimension(:) :: ptwetfld
c
      logical lptdepout
c
      common /depoutdat/ lptdepout
c
c-----------------------------------------------------------------------
c    Variables for species order and species flags:
c-----------------------------------------------------------------------
c
c  NOTE:  In all arrays, position zero is for the regular emissions
c         files.
c
c   lvocsp   L   flag to determine if species is VOC species
c   lnoxsp   L   flag to determine if species is NOx species
c   lo3sp    L   flag to determine if species is O3 species
c   lvocsoa  L   flag to determine if species is SOA-producing VOC
c   lhrvoc   L   flag to determine if species is HRVOC species
c   crbnum   R   carbon number of each species 
c   mwspec   R   molecular weight for each model species
c   idxems   I   index of each species in each emissions file into
c                species list arrays
c   idxpts   I   index of each species in each point source file
c                into species list arrays
c   idx_point_in_list I index of point source in master list from base 
c                     inventory
c   nspcem   I   number of species in each emissions file
c   nspcpt   I   number of species in each point source file
c   rkohrt   R   the reactivity fraction of each species
c   rmirrt   R   the MIR reactivity for each VOC species
c   wtkoh    R   the kOH reactivity weighting factor for each tracer
c                species (really only needed for VOC species, but easier
c                to code this way)
c   wtmir    R   the MIR reactivity weighting factor for each tracer
c                species (really only needed for VOC species, but easier
c                to code this way)
c   yhrates  R   high-NOx yield rates for each of the tracer species
c   ylrates  R   low-NOx yield rates for each of the tracer species
c
c   nptsrc_safile   number of point sources in each file
c   idx_start_sapts index of starting point into global point source array
c 
      integer, allocatable, dimension(:,:,:,:) :: idxems
      integer, allocatable, dimension(:,:,:)   :: nspcem
      integer, allocatable, dimension(:,:,:)   :: idxpts
      integer, allocatable, dimension(:,:,:)   :: idx_point_in_list
      integer, allocatable, dimension(:,:)     :: nspcpt
      integer, allocatable, dimension(:,:)     :: nptsrc_safile
      integer, allocatable, dimension(:,:)     :: idx_start_sapts
      logical, allocatable, dimension(:)       :: lvocsp
      logical, allocatable, dimension(:)       :: lnoxsp
      logical, allocatable, dimension(:)       :: lo3sp
      logical, allocatable, dimension(:)       :: lvocsoa
      logical, allocatable, dimension(:)       :: lhrvoc
      real,    allocatable, dimension(:)       :: crbnum
      real,    allocatable, dimension(:)       :: mwspec
      real,    allocatable, dimension(:)       :: rkohrt
      real,    allocatable, dimension(:)       :: rmirrt
      real,    allocatable, dimension(:)       :: wtkoh
      real,    allocatable, dimension(:)       :: wtmir
      real,    allocatable, dimension(:)       :: yhrates
      real,    allocatable, dimension(:)       :: ylrates
c
c-----------------------------------------------------------------------
c  Parameters for DDM flags:
c-----------------------------------------------------------------------
c
c    NAMVOC  C   character string for name of VOC species
c    NAMNOX  C   character string for name of NOX species
c    NAMALL  C   character string for name of ALL species
c    NAMHRV  C   character string for name of HRVOC species
c
c    IDVOC   I   species ID for VOC species
c    IDNOX   I   species ID for NOX species
c    IDALL   I   species ID for ALL species
c    IDHRV   I   species ID for HRVOC species
c
      character*10 NAMVOC
      character*10 NAMNOX
      character*10 NAMALL
      character*10 NAMHRV
      integer      IDVOC
      integer      IDNOX
      integer      IDALL
      integer      IDHRV
c
      parameter( NAMVOC = 'VOC       ')
      parameter( NAMNOX = 'NOX       ')
      parameter( NAMALL = 'ALL       ')
      parameter( NAMHRV = 'HRVOC     ')
      parameter( IDVOC  = -1 )
      parameter( IDNOX  = -2 )
      parameter( IDALL  = -3 )
      parameter( IDHRV  = -4 )
c
c-----------------------------------------------------------------------
c  Tracer species map for DDM
c-----------------------------------------------------------------------
c
c   icddmsp  C  species names for initial conditions treated by DDM
c   bcddmsp  C  species names for boundary condition treated by DDM
c   emddmsp  C  species names for emissions treated by DDM
c   rateddm  C  names for rate constant groups treated by DDM
c   termddm  C  names for rate term groups treated by DDM
c   hddmsp   C  1st-order sens parameter names to which HDDM sens is
c               computed
c   nicddm   I  number of initial conditions groups in DDM
c   nbcddm   I  number of boundary conditions groups in DDM
c   nemddm   I  number of emissions groups in DDM
c   nrateddm I  number of rate constant groups in DDM
c   ntermddm I  number of rate term groups in DDM
c   nhddm    I  number of HDDM sensitivity groups
c   iprate   I  index of rxns in rate constant groups;
c               iprate(0,i) stores the number of rxns in the i-th group
c   ipterm   I  index of rxn terms in rate term groups;
c               ipterm(0,j,i) = 1 if i-th group includes j-th rxn; otherwise, 0
c               ipterm(1,j,i) = # of reactants of j-th rxn included in the i-th group
c               ipterm(1+k,j,i) = pointer to k-th reactant of j-th rxn included ...
c               ipterm(2+ipt(1,j,i),j,i) = # of products of j-th rxn included ...
c               ipterm(2+ipt(1,j,i)+k,j,i) = pointer to k-th product of j-th rxn included ...
c   wfterm   R  weighting factors of reactant term sens in net product formation sens formula
c               same array structure as ipterm except for the 1st dimension (no 0-th element)
c   iphddm   I  index of the 1st-order sens parameters to which HDDM
c               sens is computed
c   ddmdesc  C  array of descriptions of DDM sensitivities for NCF file
c   cpadesc  C  array of descriptions for the CPA species
c   cpaunit  C  array of descriptions for the CPA species
c 
      character*10, allocatable, dimension(:)   :: icddmsp
      character*10, allocatable, dimension(:)   :: bcddmsp
      character*10, allocatable, dimension(:)   :: emddmsp
      character*10, allocatable, dimension(:)   :: rateddm
      character*10, allocatable, dimension(:)   :: termddm
      character*10, allocatable, dimension(:,:) :: hddmsp
      integer,      allocatable, dimension(:,:) :: iprate
      integer,      allocatable, dimension(:,:,:) :: ipterm
      real,         allocatable, dimension(:,:,:) :: wfterm
      integer,      allocatable, dimension(:,:) :: iphddm
      character*60, allocatable, dimension(:)   :: ddmdesc
      character*60, allocatable, dimension(:)   :: cpadesc
      character*20, allocatable, dimension(:)   :: cpaunit
c
      integer      nicddm
      integer      nbcddm
      integer      nemddm
      integer      nrateddm
      integer      ntermddm
      integer      nhddm
c
      common /ddmdat/ nicddm, nbcddm, nemddm, nrateddm, ntermddm, nhddm
c
c-----------------------------------------------------------------------
c  Temporary space for DDM sensitivities:
c-----------------------------------------------------------------------
c
c   sns   R  used for storing calculations to update DDM sensitivities
c      
      real, allocatable, dimension(:,:,:,:) :: sns
c
c-----------------------------------------------------------------------
c  Species list and pointers into arrays for DDM species:
c-----------------------------------------------------------------------
c
c   ptlong   C  long names of the DDM species
c   iptddm   I  index into the gridded arrays of the DDM families
c   nddmsp   I  total number of DDM parameters
c
      character*14, allocatable, dimension(:) :: ptlong
      integer,      allocatable, dimension(:) :: iptddm
c
      integer*4    nddmsp
c
      common /dspdat/ nddmsp
c
c-----------------------------------------------------------------------
c  Parameters for DDM-PM:
c-----------------------------------------------------------------------
c
c   indices to map species between ptconc and sddm arrays for DDM-PM
c
c   NSDDMPM  I  number of species used in DDM-PM
c   NGDDMPM  I  number of gas species used in DDM-PM
c   jdpmap   I  pointer into ptconc array of the DDM-PM species
c
      integer jdpSO2 , jdpH2O2, jdpO3  , jdpFOA , jdpMHP ,
     &        jdpOHP , jdpPAA , jdpOPA , jdpN2O5, jdpHNO3,
     &        jdpNH3 , jdpSULF, jdpHCL , jdpCG1 , jdpCG2 ,
     &        jdpCG3 , jdpCG4 , jdpGLY , jdpMGLY, jdpGLYD,
     &        jdpHO  , jdpNA  , jdpPSO4, jdpPNH4, jdpPNO3,
     &        jdpPCL , jdpSOA1, jdpSOA2, jdpSOA3, jdpSOA4,
     &        jdpSOPA, jdpSOPB, jdpPOA , jdpPH2O, jdpFPRM,
     &        jdpFCRS
      integer nsddmpm, ngddmpm
      parameter( jdpSO2  =  1 )
      parameter( jdpH2O2 =  2 )
      parameter( jdpO3   =  3 )
      parameter( jdpFOA  =  4 )
      parameter( jdpMHP  =  5 )
      parameter( jdpOHP  =  6 )
      parameter( jdpPAA  =  7 )
      parameter( jdpOPA  =  8 )
      parameter( jdpN2O5 =  9 )
      parameter( jdpHNO3 = 10 )
      parameter( jdpNH3  = 11 )
      parameter( jdpSULF = 12 )
      parameter( jdpHCL  = 13 )
      parameter( jdpCG1  = 14 )
      parameter( jdpCG2  = 15 )
      parameter( jdpCG3  = 16 )
      parameter( jdpCG4  = 17 )
      parameter( jdpGLY  = 18 )
      parameter( jdpMGLY = 19 )
      parameter( jdpGLYD = 20 )
      parameter( jdpHO   = 21 )
      parameter( jdpNA   = 22 )
      parameter( jdpPSO4 = 23 )
      parameter( jdpPNH4 = 24 )
      parameter( jdpPNO3 = 25 )
      parameter( jdpPCL  = 26 )
      parameter( jdpSOA1 = 27 )
      parameter( jdpSOA2 = 28 )
      parameter( jdpSOA3 = 29 )
      parameter( jdpSOA4 = 30 )
      parameter( jdpSOPA = 31 )
      parameter( jdpSOPB = 32 )
      parameter( jdpPOA  = 33 )
      parameter( jdpPH2O = 34 )
      parameter( jdpFPRM = 35 )
      parameter( jdpFCRS = 36 )
      parameter( NSDDMPM = 36 )
      parameter( NGDDMPM = 21 )
      integer jdpmap(NSDDMPM)
      common /ddmpm/ jdpmap
c
c-----------------------------------------------------------------------
c  Variables for Rtrac/PiG sampling grids:
c-----------------------------------------------------------------------
c
c   lsmptrc  L  include Rtrac concs in sampling grid output
c   iprtsmp  I  pointer array for 2-D RTRAC sampling grid variables
c   nrtsmpcels I  number of total sampling grids (size of rtsmpgrd array)
c   rtsmpcnc R  time-averaged concentration array on sampling grid
c               (gasses=ppm, PM=ug/m3)
c   puffrt   R  PiG mass (gasses=umol,PM=ug)
c
      integer,      allocatable, dimension(:)     ::  iprtsmp
      real,         allocatable, dimension(:,:,:) ::  puffrt
      real,         allocatable, dimension(:)     ::  rtsmpcnc
c
      integer nrtsmpcels
      logical lsmptrc
      common /rtracsmp/ nrtsmpcels, lsmptrc
c
c-----------------------------------------------------------------------
c     Variables for boundary condition concentrastions
c
c      bndsa  -- SA boundary conditions on edge cells
c      bndddm -- sensitivities on the edge cells
c      bndrt  -- RTRAC boundary conditions on edge cells
c
c-----------------------------------------------------------------------
c
      real, allocatable, dimension(:,:,:,:) :: bndsa
      real, allocatable, dimension(:,:,:,:) :: bndddm
      real, allocatable, dimension(:,:,:,:) :: bndrt
c
c-----------------------------------------------------------------------
c Variables for MPI
c-----------------------------------------------------------------------
c
      integer  mvecscr_pt
      real, allocatable, dimension(:) :: scr1_pt
      integer  :: pt_identifier = 411   ! assume there are at most ngrid=10
      integer  :: dry_identifier = 611
      integer  :: wet_identifier = 811
      integer  :: rtsol_identifier = 1011
      integer  :: rtveg_identifier = 1211
c
c-----------------------------------------------------------------------
c Parameters for Partial source maps
c-----------------------------------------------------------------------
c
c   MXPARTIAL I  maximum number of partial source maps allowed
c
      integer MXPARTIAL
c
      parameter( MXPARTIAL = 10 )
c
c-----------------------------------------------------------------------
c Variables for Partial source maps
c-----------------------------------------------------------------------
c
c   lpartial    L   flag to determine if using partial source area maps
c   npartial    I   number of source area maps (will be 1 if not using
c                   partial source area map)
c
      logical lpartial

      integer, allocatable, dimension(:,:) :: npartial

      common /partdat/ lpartial
