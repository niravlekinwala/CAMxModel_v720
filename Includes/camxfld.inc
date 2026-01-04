c----CAMx v7.20 220430
c
c     CAMxFLD.COM contains all multidimensional fields for CAMx
c 
c      Copyright 1996 - 2022
c     Ramboll
c
c     Modifications:
c        1/10/03    Added deposition field array
c        1/23/03    New cloud parameters from cloud/rain file
c       10/12/04    Added time-dependent fields for water vapor and
c                   vertical diffusivity
c        6/21/05    Added fields for rain, snow, and graupel cloud water
c       10/24/05    Added cloud pH field
c        4/27/06    Added topo field
c       04/30/13    Added surface model
c       08/08/14    Added new snow cover input fields
c       09/02/14    Added subgrid convective model
c       10/13/17    Added aerosol pH field (APH)
c       10/16/17    Added time-weighted NO2 photolysis rate field (AJNO2)
c       07/23/18    Added emiss flux of Bi-Di NH3 drydep (eflxnh3)
c
c-----------------------------------------------------------------------
c     Variables for 2-dimensional fields:
c
c     cellon --  array of longitudes for the cell centroids (deg)
c     cellat --  array of lattitudes for the cell centroids (deg)
c     mapscl --  array of map-scale factors at cell centroids
c     topo   --  array of topographic elevation (m MSL)
c     lai    --  array of leaf area index
c     lrdlai --  flag indicating that LAI was read
c     tsurf  --  array for surface temperature field (K)
c     pspt   --  array for time-rate change of surface temperature (K/s) 
c     sfcz0  --  array of surface roughness (m)
c     snow   --  array of snow cover water equivalent (m)
c     snowage--  array of snow cover age (hr)
c     snowrat--  array of snow water accumulation rate (m/hr)
c     snowfrc--  array of net snow cover fraction
c     snowalb--  array of snow surface albedo
c     albedo --  array of net surface albedo
c-----------------------------------------------------------------------
c
      real, allocatable, dimension(:) :: cellon
      real, allocatable, dimension(:) :: cellat
      real, allocatable, dimension(:) :: mapscl
      real, allocatable, dimension(:) :: topo
      real, allocatable, dimension(:) :: lai
      real, allocatable, dimension(:) :: tsurf
      real, allocatable, dimension(:) :: pspt
      real, allocatable, dimension(:) :: sfcz0
      real, allocatable, dimension(:) :: snow
      real, allocatable, dimension(:) :: snowage
      real, allocatable, dimension(:) :: snowrat
      real, allocatable, dimension(:) :: snowfrc
      real, allocatable, dimension(:) :: snowalb
      real, allocatable, dimension(:) :: albedo
      logical, allocatable, dimension(:) :: lrdlai
c
c-----------------------------------------------------------------------
c     Variables for 3-dimensional fields:
c
c      windu  --  U-component of the wind field (m/s)
c      windv  --  V-component of the wind field (m/s)
c      pupt   --  time-rate change of U-component wind (m/s2)
c      pvpt   --  time-rate change of V-component wind (m/s2)
c      tempk  --  temperature field (K) 
c      ptpt   --  time-rate change of temperature (K/s)
c      press  --  pressure field (mb)
c      pppt   --  time-rate change of pressure (mb/s)
c      height --  layer interface height field (m)
c      phpt   --  time-rate change of layer interface height (m/s)
c      rkv    --  vertical diffusion coefficient field (m2/s)
c      pkpt   --  time-rate change of Kv (m2/s2)
c      water  --  water vapor concentration field (ppm)
c      pwpt   --  time-rate change of water vapor (ppm/s)
c      fcloud --  vertically-accumulating cloud coverage field (fraction)
c      depth  --  layer depth field (m)
c      rkx    --  horizontal diffusion coefficient in X direction (m2/s)
c      rky    --  horizontal diffusion coefficient in Y direction (m2/s)
c      aremis --  surface emissions (moles/hour or g/hour)
c      ctop   --  top concentrations
c      cwc    --  cloud water content (g/m3)
c      pwr    --  rain water content (g/m3)
c      pws    --  snow water content (g/m3)
c      pwg    --  graupel water content (g/m3)
c      cod    --  cloud optical depth
c      cldtrns -  cloud energy transmission coefficient (fraction)
c      cph    --  cloud water pH
c      aph    --  aerosol pH
c      ajno2  --  time-weighted accumulated NO2 photolysis rate ([1/hr]*[hr])
c      cigfrc --  subgrid cloud fraction (unitless)
c      cigtim --  subgrid cloud timescale (s)
c      cigwtr --  subgrid cloud water (g/m3)
c      cigph  --  subgrid cloud pH
c      cigpcp --  subgrid cloud precip (g/m3)
c      cigent --  subgrid cloud entrainment rate (kg/m2/s)
c      cigdet --  subgrid cloud detrainment rate (kg/m2/s)
c      cigmas --  subgrid cloud transport matrix (unitless)
c      fluxtmp -- temporary gridded variable used to store fluxes for diffusion
c
c-----------------------------------------------------------------------
c
      real, allocatable, dimension(:) :: windu
      real, allocatable, dimension(:) :: windv
      real, allocatable, dimension(:) :: pupt
      real, allocatable, dimension(:) :: pvpt
      real, allocatable, dimension(:) :: tempk
      real, allocatable, dimension(:) :: ptpt
      real, allocatable, dimension(:) :: press
      real, allocatable, dimension(:) :: pppt
      real, allocatable, dimension(:) :: height
      real, allocatable, dimension(:) :: phpt
      real, allocatable, dimension(:) :: rkv
      real, allocatable, dimension(:) :: pkpt
      real, allocatable, dimension(:) :: water
      real, allocatable, dimension(:) :: pwpt
cae      real, allocatable, dimension(:) :: fcloud
      real, allocatable, dimension(:) :: depth
      real, allocatable, dimension(:) :: rkx
      real, allocatable, dimension(:) :: rky
      real, allocatable, dimension(:) :: aremis
      real, allocatable, dimension(:) :: ctop
      real, allocatable, dimension(:) :: cwc
      real, allocatable, dimension(:) :: pwr
      real, allocatable, dimension(:) :: pws
      real, allocatable, dimension(:) :: pwg
      real, allocatable, dimension(:) :: cod
      real, allocatable, dimension(:) :: cldtrns
      real, allocatable, dimension(:) :: cph
      real, allocatable, dimension(:) :: aph
      real, allocatable, dimension(:) :: ajno2
      real, allocatable, dimension(:) :: cigfrc
      real, allocatable, dimension(:) :: cigtim
      real, allocatable, dimension(:) :: cigwtr
      real, allocatable, dimension(:) :: cigph
      real, allocatable, dimension(:) :: cigpcp
      real, allocatable, dimension(:) :: cigent
      real, allocatable, dimension(:) :: cigdet
      real, allocatable, dimension(:) :: cigmas
      real, allocatable, dimension(:,:,:,:) :: fluxtmp
c
c-----------------------------------------------------------------------
c     Variables for 4-dimensional fields:
c
c      conc   --  species concentrations field (umol/m3)
c      avcnc  --  average species concentration (gas=ppm,other=ug/m3)
c
c-----------------------------------------------------------------------
c
      real, allocatable, dimension(:) :: conc
      real, allocatable, dimension(:) :: avcnc
c
c-----------------------------------------------------------------------
c     Variables for Surface Model fields:
c
c      solmas --  surface mass in soil (mol/ha)
c      vegmas --  surface mass in vegetation (mol/ha)
c      reemis --  re-emitted surface mass (mol/s)
c
c-----------------------------------------------------------------------
c
      real, allocatable, dimension(:) :: solmas
      real, allocatable, dimension(:) :: vegmas
      real, allocatable, dimension(:) :: reemis
c
c-----------------------------------------------------------------------
c     Variables for boundary condition concentrations
c
c      bndconc -- concentrations on the edge cells
c
c-----------------------------------------------------------------------
c
      real, allocatable, dimension(:,:,:,:) :: bndconc
c
c-----------------------------------------------------------------------
c     Variables for calculating vertical transport:
c     NOTE:  These fields are over-written when processing each grid.
c
c     entrn  --  entrainment rate (m/s)
c     dilut  --  dilution rate (m/s)
c-----------------------------------------------------------------------
c   
      real, allocatable, dimension(:) :: entrn
      real, allocatable, dimension(:) :: dilut
c
c-----------------------------------------------------------------------
c     Variables for calculating depostion rates:
c
c     vdep   --  species-dependent deposition velocity field (m/s)
c     fsurf  --  fractional landuse cover field (fraction)
c     depfld --  2-D array containing dry, wet dep mass (mol/ha, g/ha) and
c                precip liquid concentration (mol/l, g/l)
c     eflxnh3 -- 2-D array for emiss flux of Bi-Di NH3 drydep (umol/m2-s)
c-----------------------------------------------------------------------
c
      real, allocatable, dimension(:) :: vdep
      real, allocatable, dimension(:) :: fsurf
      real, allocatable, dimension(:) :: depfld
      real, allocatable, dimension(:) :: eflxnh3
c
      integer, parameter :: EX_SCRATCH=100
      integer  mvecscr
      integer  mvecscr_dp
      real, allocatable, dimension(:) :: scr1
      real, allocatable, dimension(:) :: scr1_dp
      integer  :: dp_identifier = 111 
      integer  :: rad_identifier = 711
      integer  :: smsol_identifier = 811 
      integer  :: smveg_identifier = 911
c
c-----------------------------------------------------------------------
c     Variables for mass flux calculations:
c
c     xmass   -- current total grid mass by species (moles)
c     xmass0  -- initial total grid mass by species (moles)
c     armass  -- total area emissions mass by species (moles)
c     ptmass  -- total point emissions mass by species (moles)
c     fluxes  -- array of mass tranport by species (moles)
c     xmschem -- mass change due to chemistry by species (moles)
c     xmsold  -- total grid mass by species before last process (moles)
c     resid   -- residual (error) mass by species (moles)
c     xmsfin  -- current total fine grid mass by species (moles)
c     xmstmp  -- temporary total fine grid mass by species (moles)
c     pigdump -- total mass transferred from PiG to grid by species (umol)
c     pgmserr -- error in mass due to aborted pig dumping (umol)
c
c-----------------------------------------------------------------------
c
      real*8, allocatable, dimension(:,:) :: xmass
      real*8, allocatable, dimension(:,:) :: xmass0
      real*8, allocatable, dimension(:,:) :: armass
      real*8, allocatable, dimension(:,:) :: ptmass
      real*8, allocatable, dimension(:,:) :: fluxes
      real*8, allocatable, dimension(:,:) :: xmschem
      real*8, allocatable, dimension(:,:) :: xmsold
      real*8, allocatable, dimension(:,:) :: resid
      real*8, allocatable, dimension(:,:) :: xmsfin
      real*8, allocatable, dimension(:,:) :: xmstmp
      real*8, allocatable, dimension(:,:) :: pigdump
      real*8, allocatable, dimension(:,:) :: pgmserr
c
c-----------------------------------------------------------------------
c     Variables for mass balance calculations for the extent of the
c     simulation:
c
c     tarmass  -- total area emissions mass by species (moles)
c     tptmass  -- total point emissions mass by species (moles)
c     tfluxes  -- array of mass tranport by species (moles)
c     txmschem -- mass change due to chemistry by species (moles)
c     tresid   -- residual (error) mass by species (moles)
c     txmsfin  -- current total fine grid mass by species (moles)
c-----------------------------------------------------------------------
c
      real*8, allocatable, dimension(:,:) :: tarmass
      real*8, allocatable, dimension(:,:) :: tptmass
      real*8, allocatable, dimension(:,:) :: tfluxes
      real*8, allocatable, dimension(:,:) :: tresid
      real*8, allocatable, dimension(:,:) :: txmschem
      real*8, allocatable, dimension(:,:) :: txmsfin
c
c-----------------------------------------------------------------------
c     Variables for storing the next hours meteorological fields:
c
c     hnxt   --  layer interface height at next hour (m)
c     pnxt   --  pressure field at next hour (mb)
c     unxt   --  U-component of the wind field at next hour (m/s)
c     vnxt   --  V-component of the wind field at next hour (m/s)
c     tnxt   --  temperature field at next hour (K) 
c     tsnxt  --  surface temperature field at next hour (K)
c     knxt   --  vertical diffusivity field at next hour (m2/s)
c     wnxt   --  water vapor field at next hour (ppm)
c-----------------------------------------------------------------------
c
      real, allocatable, dimension(:) :: hnxt
      real, allocatable, dimension(:) :: pnxt
      real, allocatable, dimension(:) :: unxt
      real, allocatable, dimension(:) :: vnxt
      real, allocatable, dimension(:) :: tnxt
      real, allocatable, dimension(:) :: tsnxt
      real, allocatable, dimension(:) :: knxt
      real, allocatable, dimension(:) :: wnxt
