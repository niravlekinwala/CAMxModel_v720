c
c----CAMx v7.20 220430
c
c     NAMELIST.COM contains all input namelist variables, for both
c     the core model and Probing Tools
c                           
c      Copyright 1996 - 2022
c     Ramboll
c           
c     Modifications: 
c        7/20/05       Moved PiG sampling grid inputs to main namelist, added
c                      new options for sampling standard species
c       03/15/09       Added code for deposition output for tracers
c        07/16/07 -bkoo-     Added variables for HDDM
c        06/11/08 -bkoo-     Added variables for rate constant sensitivity
c        07/16/08 -bkoo-     Added DDM turn-off flag
c        10/10/08 -cemery-   Added ACM2 option
c        10/29/09 -cemery-   Added RTRAC surface model option
c        10/10/08       Added ACM2 option
c        7/14/10        Added in-line TUV option
c        1/04/11        Revised for new met input format
c        03/29/11       Support in-line TUV with aerosol optical depth
c        08/31/11  -gwilson- Added easter egg parameter for APCA treatment
c                            with PSAT
c        04/02/12      Removed RADM cloud adjustment option, cloud/aerosol
c                      adjustments now always done with in-line TUV; AHO
c                      file is now just ozone column; number of output
c                      species now determined internally from user list
c        05/07/12      Added flexi-nesting flag
c        04/30/13      Added surface model
c        09/02/14      Added subgrid convective model
c        09/11/15      Revised for SA v3
c        03/01/16      Added partial source area map
c        05/13/16      Added in-line Ix emissions flag
c        11/09/16      Added Baker APCA point source override option
c        02/28/17      Added NetCDF output flag
c        07/23/18      Added Bi-Di NH3 drydep flag
c        08/09/18      Added variables for rate term sensitivity
c        07/10/19      Added stratospheric ozone treatment flag
c        04/15/21      Added flag for output gas species in PPM
c        01/05/22      Added flag for PPM vertical advection
c
c-----------------------------------------------------------------------
c
      integer MXNAM
      parameter( MXNAM = 999 )
      integer MXFILES
      parameter( MXFILES = 100 )
c
      logical       Restart
      logical       Flexi_Nest
      logical       Chemistry
      logical       Wet_Deposition
      logical       Gridded_Emissions
      logical       Point_Emissions
      logical       Ignore_Emission_Dates
      logical       Diagnostic_Error_Check
      logical       Average_Output_3D
      logical       Output_Gas_Concs_PPM
      logical       Output_3D_Grid(MXNAM)
      logical       NetCDF_Format_Output
      logical       PiG_Sampling_Grid
      logical       Sample_Background
      logical       Super_Stepping
      logical       ACM2_Diffusion
      logical       Surface_Model
      logical       Subgrid_Convection
      logical       Bidi_NH3_Drydep
      logical       Strat_Ozone_Profile
c
      character*200 Root_Output_Name
      character*200 Chemistry_Parameters
      character*200 Photolyis_Rates
      character*200 Photolysis_Rates
      character*200 Initial_Conditions
      character*200 Boundary_Conditions
      character*200 Top_Concentrations
      character*200 Ozone_Column
      character*200 Point_Sources(MXNAM)
      character*200 Master_Grid_Restart
      character*200 Nested_Grid_Restart
      character*200 PiG_Restart
      character*200 Surface_Grid(MXNAM)
      character*200 Met2D_Grid(MXNAM)
      character*200 Met3D_Grid(MXNAM)
      character*200 Cloud_Grid(MXNAM)
      character*200 Vdiff_Grid(MXNAM)
      character*200 Emiss_Grid(MXGRID,MXNAM)
      character*200 Srfmod_Grid(MXNAM)
c
      character*60  Run_Message
c
      character*10  Map_Projection
      character*10  Advection_Solver
      character*10  Vadvection_Solver
      character*10  Chemistry_Solver
      character*10  Drydep_Model
      character*10  PiG_Submodel
      character*10  Probing_Tool
      character*10  Inline_Ix_Emissions
      character*10  Output_Species_Names(MXNAM)
c
      integer       Time_Zone
      integer       Start_Date_Hour(4)
      integer       End_Date_Hour(4)
      integer       UTM_Zone
      integer       Number_of_Grids
      integer       Master_Grid_Columns
      integer       Master_Grid_Rows
      integer       Number_of_Layers
      integer       Nest_Meshing_Factor(MXNAM)
      integer       Nest_Beg_I_Index(MXNAM)
      integer       Nest_End_I_Index(MXNAM)
      integer       Nest_Beg_J_Index(MXNAM)
      integer       Nest_End_J_Index(MXNAM)
      integer       Number_of_Sampling_Grids
      integer       SG_Beg_I_Index(MXNAM)
      integer       SG_End_I_Index(MXNAM)
      integer       SG_Beg_J_Index(MXNAM)
      integer       SG_End_J_Index(MXNAM)
c
      real          Maximum_Timestep
      real          Met_Input_Frequency
      real          Ems_Input_Frequency
      real          Output_Frequency
      real          Longitude_Pole
      real          Latitude_Pole
      real          True_Latitude1
      real          True_Latitude2
      real          Master_SW_XCoord
      real          Master_SW_YCoord
      real          Master_Cell_XSize
      real          Master_Cell_YSize
      real          SG_Mesh_Factor(MXNAM)
c
      common /chr_camx_cntrl/ 
     &             Root_Output_Name, Chemistry_Parameters,
     &             Photolyis_Rates,Photolysis_Rates, Initial_Conditions,
     &             Boundary_Conditions, Top_Concentrations,
     &             Ozone_Column, Point_Sources,
     &             Master_Grid_Restart, Nested_Grid_Restart,
     &             PiG_Restart, Surface_Grid, Met2D_Grid, Met3D_Grid,
     &             Cloud_Grid, Vdiff_Grid, Emiss_Grid, Srfmod_Grid,
     &             Run_Message, Map_Projection, Advection_Solver,
     &             Vadvection_Solver,Chemistry_Solver, Drydep_Model,
     &             PiG_Submodel,Probing_Tool,Inline_Ix_Emissions,
     &             Output_Species_Names
c
      common /dat_camx_cntrl/ 
     &      Restart, Flexi_Nest, Chemistry, Wet_Deposition,
     &      Gridded_Emissions, Point_Emissions, Ignore_Emission_Dates,
     &      Diagnostic_Error_Check, Average_Output_3D, Output_3D_Grid,
     &      Output_Gas_Concs_PPM, NetCDF_Format_Output, Super_Stepping, 
     &      ACM2_Diffusion, Surface_Model, Subgrid_Convection, 
     &      Bidi_NH3_Drydep, Strat_Ozone_Profile,
     &      Time_Zone, Start_Date_Hour, End_Date_Hour,
     &      UTM_Zone, Number_of_Grids, Master_Grid_Columns, 
     &      Master_Grid_Rows,Number_of_Layers, Nest_Meshing_Factor,
     &      Nest_Beg_I_Index, Nest_End_I_Index,
     &      Nest_Beg_J_Index, Nest_End_J_Index,
     &      Maximum_Timestep, Met_Input_Frequency, Ems_Input_Frequency,
     &      Output_Frequency, Longitude_Pole, Latitude_Pole,
     &      True_Latitude1, True_Latitude2, Master_SW_XCoord,
     &      Master_SW_YCoord, Master_Cell_XSize, Master_Cell_YSize,
     &      PiG_Sampling_Grid, Sample_Background, Sampling_Grid, 
     &      Number_of_Sampling_Grids, 
     &      SG_Beg_I_Index, SG_End_I_Index, SG_Beg_J_Index, 
     &      SG_End_J_Index, SG_Mesh_Factor
c
c  --- Variables for NETCDF ---
c
      logical NetCDF_Use_Compression

      common /dat_ncf_chunk/ NetCDF_Use_Compression
c
c======================== Probing Tool Begin ===========================
c
      logical       SA_Master_Sfc_Output
      logical       SA_Nested_Sfc_Output
      logical       SA_Deposition_Output
      logical       SA_Stratify_Boundary
      logical       SA_PT_Override
      logical       SA_3D_Average
      logical       Use_Leftover_Group
      logical       Use_Gridded_Leftover_Group
      logical       SA_Summary_Output
      logical       DDM_Master_Sfc_Output
      logical       DDM_Nested_Sfc_Output
      logical       DDM_Stratify_Boundary
      logical       DDM_Calc_Grid(MXNAM)
      logical       DDM_PT_Override
      logical       Sampling_Grid
      logical       SA_Treat_SULFATE_Class
      logical       SA_Treat_NITRATE_Class
      logical       SA_Treat_SOA_Class
      logical       SA_Treat_PRIMARY_Class
      logical       SA_Treat_MERCURY_Class
      logical       SA_Treat_OZONE_Class
      logical       SA_Use_APCA
      logical       SA_Use_APCA_Ptoverride
      logical       SA_Use_Partial_SourceMap
      logical       RT_Surface_Model
      logical       RT_Partitioning
c
      character*200 SA_File_Root
      character*200 SA_Receptor_Definitions
      character*200 SA_Source_Area_Map(MXNAM)
      character*200 Partial_Source_Area_Map(MXNAM,MXGRID)
      character*200 SA_Master_Restart
      character*200 SA_Nested_Restart
      character*200 SA_Points_Group(0:MXNAM,MXFILES)
      character*200 SA_Emiss_Group_Grid(0:MXNAM,MXGRID,MXFILES)
      character*200 SA_Initial_Conditions
      character*200 SA_Boundary_Conditions
      character*200 SA_Top_Concentrations
      character*200 DDM_File_Root
      character*200 DDM_Receptor_Definitions
      character*200 DDM_Source_Area_Map(MXNAM)
      character*200 DDM_Initial_Conditions
      character*200 DDM_Boundary_Conditions
      character*200 DDM_Top_Concentrations
      character*200 DDM_Master_Restart
      character*200 DDM_Nested_Restart
      character*200 DDM_Points_Group(MXNAM,MXFILES)
      character*200 DDM_Emiss_Group_Grid(MXNAM,MXGRID,MXFILES)
      character*200 RT_File_Root
      character*200 RT_Initial_Conditions
      character*200 RT_Boundary_Conditions
      character*200 RT_Top_Concentrations
      character*200 RT_Master_Restart
      character*200 RT_Nested_Restart
      character*200 RT_Chemistry_Parameters
      character*200 RT_Receptor_Definitions
      character*200 RT_Point_Sources
      character*200 RT_Emiss_Grid(MXNAM)
      character*200 RT_Srfmod_Grid(MXNAM)
      character*200 PA_File_Root
      character*200 Rate_Const_Groups(MXNAM)
      character*200 Rate_Term_Groups(MXNAM,MXNAM)
c
      character*10  IC_Species_Groups(MXNAM)
      character*10  BC_Species_Groups(MXNAM)
      character*10  Emis_Species_Groups(MXNAM)
      character*10  HDDM_parameters(MXNAM,MXNAM)
c
      integer       SA_Number_of_Source_Regions
      integer       SA_Number_of_Source_Groups
      integer       Number_of_Timing_Releases
      integer       DDM_Number_of_Source_Regions
      integer       DDM_Number_of_Source_Groups
      integer       Number_of_IC_Species_Groups
      integer       Number_of_BC_Species_Groups
      integer       Number_of_EM_Species_Groups
      integer       Number_of_HDDM_Sens_Groups
      integer       Number_of_Rate_Const_Groups
      integer       Number_of_Rate_Term_Groups
      integer       Number_of_PA_Domains
      integer       Within_CAMx_Grid(MXNAM)
      integer       PA_Beg_I_Index(MXNAM)
      integer       PA_End_I_Index(MXNAM)
      integer       PA_Beg_J_Index(MXNAM)
      integer       PA_End_J_Index(MXNAM)
      integer       PA_Beg_K_Index(MXNAM)
      integer       PA_End_K_Index(MXNAM)
      logical       RT_PiG_Sample
c
      common /chr_sa_cntrl/
     &        SA_File_Root, SA_Receptor_Definitions, SA_Source_Area_Map,
     &        SA_Master_Restart, SA_Nested_Restart, SA_Points_Group,
     &        SA_Emiss_Group_Grid, Partial_Source_Area_Map,
     &        SA_Initial_Conditions, SA_Boundary_Conditions,
     &        SA_Top_Concentrations
c
      common /dat_sa_cntrl/
     &        SA_Master_Sfc_Output, SA_Nested_Sfc_Output,
     &        SA_Deposition_Output,
     &        SA_Stratify_Boundary, Use_Leftover_Group,
     &        Use_Gridded_Leftover_Group,
     &        SA_Summary_Output, SA_Treat_SULFATE_Class,
     &        SA_Treat_NITRATE_Class, SA_Treat_SOA_Class,
     &        SA_Treat_PRIMARY_Class, SA_Treat_MERCURY_Class,
     &        SA_Treat_OZONE_Class, SA_Number_of_Source_Regions,
     &        SA_Number_of_Source_Groups, Number_of_Timing_Releases,
     &        SA_Use_APCA, SA_Use_APCA_Ptoverride, 
     &        SA_Use_Partial_SourceMap, SA_PT_Override, SA_3D_Average
c
      common /chr_ddm_cntrl/
     &        DDM_File_Root, DDM_Receptor_Definitions,
     &        DDM_Source_Area_Map, DDM_Initial_Conditions,
     &        DDM_Boundary_Conditions, DDM_Top_Concentrations,
     &        DDM_Master_Restart, DDM_Nested_Restart,
     &        DDM_Points_Group, DDM_Emiss_Group_Grid,
     &        IC_Species_Groups, BC_Species_Groups, Emis_Species_Groups,
     &        HDDM_parameters,Rate_Const_Groups, Rate_Term_Groups
c
      common /dat_ddm_cntrl/
     &        DDM_Master_Sfc_Output, DDM_Nested_Sfc_Output,
     &        DDM_Stratify_Boundary, DDM_Number_of_Source_Regions,
     &        DDM_Number_of_Source_Groups, Number_of_IC_Species_Groups,
     &        Number_of_BC_Species_Groups, Number_of_EM_Species_Groups,
     &        Number_of_HDDM_Sens_Groups,Number_of_Rate_Const_Groups,
     &        Number_of_Rate_Term_Groups, DDM_Calc_Grid, DDM_PT_Override
c
      common /chr_rt_cntrl/
     &        RT_File_Root, RT_Initial_Conditions,
     &        RT_Boundary_Conditions, RT_Top_Concentrations,
     &        RT_Master_Restart, RT_Nested_Restart,
     &        RT_Chemistry_Parameters, RT_Receptor_Definitions,
     &        RT_Point_Sources, RT_Emiss_Grid, RT_PiG_Sample,
     &        RT_Surface_Model,RT_Partitioning,
     &        RT_Srfmod_Grid
c
      common /chr_pa_cntrl/
     &        PA_File_Root
c
      common /dat_pa_cntrl/
     &        Number_of_PA_Domains, Within_CAMx_Grid,
     &        PA_Beg_I_Index, PA_End_I_Index,
     &        PA_Beg_J_Index, PA_End_J_Index,
     &        PA_Beg_K_Index, PA_End_K_Index
c
c======================== Probing Tool End ===========================
c
