      subroutine initnml()
c
c-----------------------------------------------------------------------
c    Modules used:
c-----------------------------------------------------------------------
c
      use o3colmap
      use camxfld
      use camxcom
      use pigsty
c
      implicit none
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c         This routine initializes the namelist variables
c         to default values.
c-----------------------------------------------------------------------
c
c    Argument descriptions:
c     Input:  
c     Output:  
c
c    Called by:
c       STARTUP
c
c      Copyright 1996 - 2022
c     Ramboll
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c       03/15/09     Added code for deposition output for tracers 
c       10/29/09     Added code for RTRAC surface model
c       01/04/11     Revised for new met input format
c       05/07/12     Added flexi-nesting flag
c       04/30/13     Added surface model
c       09/02/14     Added subgrid convective model
c       04/10/15     Added WRF polar and mercator projections
c       09/11/15     Revised for SA v3
c       03/01/16     Added partial source area map
c       05/13/16     Added in-line Ix emissions flag
c       11/09/16     Added Baker APCA point source override option
c       07/23/18     Added Bi-Di NH3 drydep flag
c       08/09/18     Added variables for rate term sensitivity
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'flags.inc'
      include 'namelist.inc'
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      integer :: i
      integer :: n
      integer :: m
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
      Flexi_Nest              = .false.
      Restart                 = .false.
      Chemistry               = .false.
      Wet_Deposition          = .false.
      ACM2_Diffusion          = .false.
      Gridded_Emissions       = .false.
      Point_Emissions         = .false.
      Ignore_Emission_Dates   = .false.
      Diagnostic_Error_Check  = .false.
      PiG_Sampling_Grid       = .false.
      Sample_Background       = .false.
      Average_Output_3D       = .false.
      Output_Gas_Concs_PPM    = .true.
      Super_Stepping          = .true.
      Surface_Model           = .false.
      Subgrid_Convection      = .false.
      Bidi_NH3_Drydep         = .false.
      Strat_Ozone_Profile     = .false.
      NetCDF_Format_Output    = .false.
      do i = 1,MXNAM
        Output_3D_Grid(i)    = .false.
      enddo

      Root_Output_Name     = ' '
      Chemistry_Parameters = ' '
      Photolyis_Rates      = ' '
      Photolysis_Rates     = ' '
      Ozone_Column         = ' '
      Initial_Conditions   = ' '
      Boundary_Conditions  = ' '
      Top_Concentrations   = ' '
      Master_Grid_Restart  = ' '
      Nested_Grid_Restart  = ' '
      PiG_Restart          = ' '
      do i = 1,MXNAM
        Point_Sources(i)   = ' '
        Surface_Grid(i)    = ' '
        Met2D_Grid(i)      = ' '
        Met3D_Grid(i)      = ' '
        Cloud_Grid(i)      = ' '
        Vdiff_Grid(i)      = ' '
        Srfmod_Grid(i)     = ' '
      enddo
      do i = 1,MXGRID
       do n = 1,MXNAM
            Emiss_Grid(i,n) = ' '
        enddo
      enddo
      Run_Message          = ' '
      Map_Projection       = ' '
      Advection_Solver     = ' '
      Vadvection_Solver    = ' '
      Chemistry_Solver     = ' '
      PiG_Submodel         = ' '
      Probing_Tool         = ' '
      Inline_Ix_Emissions  = ' '
      do i = 1,MXNAM
        Output_Species_Names(i) = ' '
      enddo

      Time_Zone = 0
      do i = 1,4
        Start_Date_Hour(i)     = 0
        End_Date_Hour(i)       = 0
      enddo
      UTM_Zone                 = 0
      Number_of_Grids          = 0
      Master_Grid_Columns      = 0
      Master_Grid_Rows         = 0
      Number_of_Layers         = 0
      do i = 1,MXNAM
        Nest_Meshing_Factor(i) = 0
        Nest_Beg_I_Index(i)    = 0
        Nest_End_I_Index(i)    = 0
        Nest_Beg_J_Index(i)    = 0
        Nest_End_J_Index(i)    = 0
      enddo

      Maximum_Timestep         = 15.
      Met_Input_Frequency      = 60.
      Ems_Input_Frequency      = 60.
      Output_Frequency         = 60.
      Longitude_Pole           = 0.
      Latitude_Pole            = 0.
      True_Latitude1           = 0.
      True_Latitude2           = 0.
      Master_SW_XCoord         = 0.
      Master_SW_YCoord         = 0.
      Master_Cell_XSize        = 0.
      Master_Cell_YSize        = 0.

      Number_of_Sampling_Grids     = 0
      do n = 1,MXNAM
        SG_Beg_I_Index(n)          = 0
        SG_End_I_Index(n)          = 0
        SG_Beg_J_Index(n)          = 0
        SG_End_J_Index(n)          = 0
        SG_Mesh_Factor(n)          = 0.
      enddo
c
c  --- Defaults for NetCDF chunking parameters ---
c
      NetCDF_Use_Compression = .FALSE.
c
c======================== Probing Tool Begin ===========================
c
      SA_Master_Sfc_Output  = .false.
      SA_Nested_Sfc_Output  = .false.
      SA_Deposition_Output  = .false.
      SA_Stratify_Boundary  = .false.
      Use_Leftover_Group    = .false.
      Use_Gridded_Leftover_Group    = .false.
      SA_Summary_Output     = .false.
      SA_PT_Override        = .false.
      SA_3D_Average         = .false.

      DDM_Master_Sfc_Output = .false.
      DDM_Nested_Sfc_Output = .false.
      DDM_Stratify_Boundary = .false.
      DDM_PT_Override       = .false.

      SA_File_Root                     = ' '
      SA_Receptor_Definitions          = ' '
      SA_Initial_Conditions            = ' '
      SA_Boundary_Conditions           = ' '
      SA_Top_Concentrations            = ' '
      do n = 1,MXNAM
        SA_Source_Area_Map(n)          = ' '
      enddo
      SA_Master_Restart                = ' '
      SA_Nested_Restart                = ' '
      do n = 1,MXNAM
        do i = 1,MXFILES
           SA_Points_Group(n,i)        = ' '
        enddo
      enddo
      do n = 1,MXGRID
        do m = 1,MXNAM
          do i = 1,MXFILES
            SA_Emiss_Group_Grid(m,n,i) = ' '
          enddo
        enddo
      enddo
      SA_Treat_SULFATE_Class = .false.
      SA_Treat_NITRATE_Class = .false.
      SA_Treat_SOA_Class     = .false.
      SA_Treat_PRIMARY_Class = .false.
      SA_Treat_MERCURY_Class = .false.
      SA_Treat_OZONE_Class   = .false.
      SA_Use_APCA            = .false.
      SA_Use_APCA_Ptoverride = .false.
      SA_Use_Partial_SourceMap = .false.
      Partial_Source_Area_Map  = ' '

      DDM_File_Root                     = ' '
      DDM_Receptor_Definitions          = ' '
      do n = 1,MXNAM
        DDM_Source_Area_Map(n)          = ' '
      enddo
      DDM_Initial_Conditions            = ' '
      DDM_Boundary_Conditions           = ' '
      DDM_Top_Concentrations            = ' '
      DDM_Master_Restart                = ' '
      DDM_Nested_Restart                = ' '
      do n = 1,MXNAM
        do i = 1,MXFILES
           DDM_Points_Group(n,i)        = ' '
        enddo
        do m = 1,MXNAM
          HDDM_parameters(m,n)          = ' '
          Rate_Term_Groups(m,n)     = ' '
        enddo
        IC_Species_Groups(n)            = ' '
        BC_Species_Groups(n)            = ' '
        Emis_Species_Groups(n)          = ' '
        Rate_Const_Groups(n)            = ' '
        DDM_Calc_Grid(n)                = .TRUE.
      enddo
      do n = 1,MXGRID
        do m = 1,MXNAM
          do i = 1,MXFILES
            DDM_Emiss_Group_Grid(m,n,i) = ' '
          enddo
        enddo
      enddo

      RT_File_Root                  = ' '
      RT_Initial_Conditions         = ' '
      RT_Boundary_Conditions        = ' '
      RT_Top_Concentrations         = ' '
      RT_Master_Restart             = ' '
      RT_Nested_Restart             = ' '
      RT_Chemistry_Parameters       = ' '
      RT_Receptor_Definitions       = ' '
      RT_Point_Sources              = ' '
      do m = 1,MXNAM
        RT_Emiss_Grid(m)            = ' '
        RT_Srfmod_Grid(m)           = ' '
      enddo
      RT_PiG_Sample                 = .false.
      RT_Surface_Model              = .false.
      RT_Partitioning               = .false.

      PA_File_Root                  = ' '

      SA_Number_of_Source_Regions  = 0
      SA_Number_of_Source_Groups   = 0
      Number_of_Timing_Releases    = 0

      DDM_Number_of_Source_Regions = 0
      DDM_Number_of_Source_Groups  = 0
      Number_of_IC_Species_Groups  = 0
      Number_of_BC_Species_Groups  = 0
      Number_of_EM_Species_Groups  = 0
      Number_of_Rate_Const_Groups  = 0
      Number_of_HDDM_Sens_Groups   = 0
      Number_of_Rate_Term_Groups   = 0

      Number_of_PA_Domains         = 0
      do n = 1,MXNAM
        Within_CAMx_Grid(n)        = 0
        PA_Beg_I_Index(n)          = 0
        PA_End_I_Index(n)          = 0
        PA_Beg_J_Index(n)          = 0
        PA_End_J_Index(n)          = 0
        PA_Beg_K_Index(n)          = 0
        PA_End_K_Index(n)          = 0
      enddo
c
c======================== Probing Tool End =============================
c
      end
