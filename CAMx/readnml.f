      subroutine readnml(version,enddate,endtim,numprocs)
      use filunit
      use grid
      use chmstry
      use o3colmap
      use camxfld
      use camxcom
      use pigsty
      use procan
      use tracer
      use node_mod
      use rtracchm
c
      implicit none
c
c----CAMx v7.20 220430
c
c     READNML opens and reads the CAMx input namelist file called "CAMx.in"
c     that defines user inputs.  All namelist variables are mapped to 
c     internal variables and checked for appropriate/consistent values.
c                          
c      Copyright 1996 - 2022
c     Ramboll 
c          
c     Modifications:
c        7/20/05       Moved PiG sampling grid inputs to main namelist, added
c                      new options for sampling standard species
c        7/11/07       Added new RTCMC Probing Tool option
c        07/16/07 -bkoo-     Revised for HDDM
c        04/24/08 -gyarwood- Added EBI chemistry solver option
c        06/11/08 -bkoo-     Added rate constant sensitivity
c        10/09/08      Added ACM2 option
c        07/16/08 -bkoo-     Added DDM turn-off flag
c        01/30/09 -bkoo-     Removed the CMC fast solver
c        10/29/09      Added code for RTRAC surface model
c        07/14/10      Added in-line TUV option
c        03/29/11      Added option for inert PM with gas-phase chemistry
c                      and support in-line TUV with aerosol optical depth
c        04/02/12      Removed RADM cloud adjustment option, cloud/aerosol
c                      adjustments now always done with in-line TUV; AHO
c                      file is now just ozone column; number of output
c                      species now determined internally from user list
c        05/07/12      Added flexi-nesting flag
c        04/30/13      Added Surface Model
c        12/05/13      Added "ALL" option for output average species
c        09/02/14      Added subgrid convective model
c        4/10/15       Added WRF polar and mercator projections
c        06/29/15      Added "ALLR" option to include radicals
c        09/11/15      Revised for SA v3
c        03/01/16      Added partial source area map
c        05/13/16      Added in-line Ix emissions flag
c        05/31/16      Extended ACM2 to H/DDM
c        11/09/16      Added Baker APCA point source override option
c        07/23/18      Added Bi-Di NH3 drydep flag
c                      Added ALLOC_LDDMCALC call
c        08/09/18      Added variables for rate term sensitivity
c        02/21/20      Added Mechnism 1 (CB6r5) -rlb-
c        03/18/21      Initialized lddmcalc for all probing tools
c
c     Input arguments:
c        version             model version character string
c
c     Output arguments:
c        none
c
c     Routines Called:
c        CVTDATE
c        JSTLFT
c        ISTRLN
c        TOUPPER
c        INIPTR
c        OPENFILS
c
c     Called by:
c        STARTUP
c
      include 'camx.prm'
      include 'flags.inc'
      include 'namelist.inc'
      include 'deposit.inc'
      include 'mpif.h'
      include 'rtracsrf.inc'
      include 'vbs.inc'
c
      logical   Dry_Deposition
      namelist /CAMx_Control/
     & Run_Message,Time_Zone,Restart,Start_Date_Hour,End_Date_Hour,
     & Maximum_Timestep,Met_Input_Frequency,Ems_Input_Frequency,
     & Output_Frequency,Map_Projection,UTM_Zone,Longitude_Pole,
     & Latitude_Pole,True_Latitude1,True_Latitude2,Number_of_Grids,
     & Master_SW_XCoord,Master_SW_YCoord,Master_Cell_XSize,
     & Master_Cell_YSize,Master_Grid_Columns,Master_Grid_Rows,
     & Number_of_Layers,Nest_Meshing_Factor,Nest_Beg_I_Index,
     & Nest_End_I_Index,Nest_Beg_J_Index,Nest_End_J_Index,
     & Diagnostic_Error_Check,Advection_Solver,Vadvection_Solver,
     & Chemistry_Solver,
     & Drydep_Model,PiG_Submodel,Probing_Tool,Chemistry,Super_Stepping,
     & Wet_Deposition,ACM2_Diffusion,Gridded_Emissions,Point_Emissions,
     & Ignore_Emission_Dates,Surface_Model,Subgrid_Convection,
     & Inline_Ix_Emissions,Bidi_NH3_Drydep,Strat_Ozone_Profile,
     & Root_Output_Name,Average_Output_3D,Output_Gas_Concs_PPM,
     & Output_3D_Grid,NetCDF_Format_Output,NetCDF_Use_Compression,
     & Output_Species_Names,PiG_Sampling_Grid,
     & Sample_Background,Number_of_Sampling_Grids,SG_Beg_I_Index,
     & SG_End_I_Index,SG_Beg_J_Index,SG_End_J_Index,SG_Mesh_Factor,
     & Flexi_Nest,Chemistry_Parameters,Photolyis_Rates,Photolysis_Rates,
     & Initial_Conditions,Boundary_Conditions,Top_Concentrations,
     & Ozone_Column,Point_Sources,Master_Grid_Restart,
     & Nested_Grid_Restart,PiG_Restart,Surface_Grid,Met2D_Grid,
     & Met3D_Grid,Cloud_Grid,Vdiff_Grid,Emiss_Grid,Srfmod_Grid,
     & Dry_Deposition

c
c======================== Probing Tool Begin ===========================
c
      namelist /SA_Control/
     & SA_File_Root,SA_Master_Sfc_Output,SA_Nested_Sfc_Output,
     & SA_Deposition_Output,
     & SA_Stratify_Boundary,SA_Number_of_Source_Regions,
     & SA_Number_of_Source_Groups,Use_Leftover_Group,
     & Use_Gridded_Leftover_Group,
     & Number_of_Timing_Releases,SA_Receptor_Definitions,
     & SA_Source_Area_Map,SA_Master_Restart,SA_Nested_Restart,
     & SA_Points_Group,SA_Emiss_Group_Grid,SA_Summary_Output,
     & SA_Treat_SULFATE_Class,SA_Treat_NITRATE_Class,
     & SA_Treat_SOA_Class,SA_Treat_PRIMARY_Class,
     & SA_Treat_MERCURY_Class,SA_Treat_OZONE_Class,SA_Use_APCA,
     & SA_Use_APCA_Ptoverride,SA_Use_Partial_SourceMap, 
     & Partial_Source_Area_Map,SA_PT_Override,SA_3D_Average,
     & SA_Initial_Conditions, SA_Boundary_Conditions, SA_Top_Concentrations
c
      namelist /DDM_Control/
     & DDM_File_Root,DDM_Master_Sfc_Output,DDM_Nested_Sfc_Output,
     & DDM_Stratify_Boundary,DDM_Number_of_Source_Regions,
     & DDM_Number_of_Source_Groups,Number_of_IC_Species_Groups,
     & IC_Species_Groups,Number_of_BC_Species_Groups,BC_Species_Groups,
     & Number_of_EM_Species_Groups,Emis_Species_Groups,
     & Number_of_Rate_Const_Groups,Rate_Const_Groups,
     & Number_of_Rate_Term_Groups,Rate_Term_Groups,
     & Number_of_HDDM_Sens_Groups,HDDM_parameters,
     & DDM_Receptor_Definitions,DDM_Source_Area_Map,
     & DDM_Initial_Conditions,DDM_Boundary_Conditions,
     & DDM_Top_Concentrations,DDM_Master_Restart,DDM_Nested_Restart,
     & DDM_Points_Group,DDM_Emiss_Group_Grid,DDM_Calc_Grid,
     & DDM_PT_Override
c
      namelist /RT_Control/
     & RT_File_Root,RT_Initial_Conditions,RT_Boundary_Conditions,
     & RT_Top_Concentrations,RT_Master_Restart,RT_Nested_Restart,
     & RT_Chemistry_Parameters,RT_Receptor_Definitions,RT_Point_Sources,
     & RT_Emiss_Grid,RT_PiG_Sample,RT_Surface_Model,RT_Partitioning,
     & RT_Srfmod_Grid
c
      namelist /PA_Control/
     & PA_File_Root,Number_of_PA_Domains,Within_CAMx_Grid,
     & PA_Beg_I_Index,PA_End_I_Index,PA_Beg_J_Index,PA_End_J_Index,
     & PA_Beg_K_Index,PA_End_K_Index
c
c======================== Probing Tool End =============================
c
      character*20  version
      integer       enddate
      real          endtim
      integer       numprocs
c      
      integer istrln,i,jj,l,lav,n,ibyr,ibmo,ibdy,ieyr,iemo,iedy
      integer ng
      character*200 ctlfil,filtmp
      character*100 action
      character*30  keyword
      character*20  namegrp
      integer   nemiss,cbdate,cedate,inp,ii,iifroot,nopen,igrd
      logical   lexist, lacross, lmech1_ok, lmech3_ok, lmech4_ok
      logical   lmech5_ok, lmech6_ok, lmech7_ok, lmech10_ok
      logical   lixemspc, is_any_grid_emiss
c
      data inp /3/
      data ctlfil /'CAMx.in'/
c
c-----Entry point
c
c-----Open user control file and read core model parameters
c
      Emiss_Grid = ' '
      Dry_Deposition = .FALSE.
      inquire(file=ctlfil,exist=lexist)
      if( .NOT. lexist ) goto 7007
      open(unit=inp,file=ctlfil,STATUS='UNKNOWN',ERR=7005)
      namegrp = 'CAMx_Control'
      action = ' '
      read(inp,CAMx_Control,END=7100,ERR=7102)
c
      runmsg = Run_Message
      filroot = Root_Output_Name
      call jstlft( filroot )
      ii = istrln( filroot )
      if (filroot .EQ. ' ') goto 7006
c
c-----Open ASCII output, diagnostic, and mass reporting files
c
      filtmp = filroot
      filtmp(ii+1:) = '.out'
      call getunit(iout)
      action = 'Opening OUT message file.'
      open(unit=iout,file=filtmp(1:ii+4),status='UNKNOWN',ERR=7000)
c
      filtmp(ii+1:) = '.diag'
      call getunit(idiag)
      action = 'Opening DIAG diagnostic file.'
      open(unit=idiag,file=filtmp(1:ii+5),status='UNKNOWN',ERR=7000)
c
      filtmp(ii+1:) = '.mass'
      call getunit(imass)
      action = 'Opening MASS summary file.'
      open(unit=imass,file=filtmp(1:ii+5),status='UNKNOWN',ERR=7000)
      iifroot = ii
      nopen = 3
c
c-----Write model version to output and diagnostic files
c
      write(iout,8000) version(:istrln(version))
      write(idiag,8000) version(:istrln(version))
 8000 format(//,30x,20('*'),/,30x,a,/,30x,20('*'),//) 
      write(iout,8001) runmsg(:istrln(runmsg))
      write(idiag,8001) runmsg(:istrln(runmsg))
 8001 format(/,a,/)
c
c-----Set internal variables from namelist parameters and check inputs
c
c-----Clock management
c
      itzon = Time_Zone
      ibyr = Start_Date_Hour(1)
      ibmo = Start_Date_Hour(2)
      ibdy = Start_Date_Hour(3)
      begtim = float(Start_Date_hour(4))
      if (ibyr.eq.0 .or. ibmo.eq.0 .or. ibdy.eq.0) then
        write(iout,'(//,a)') 'ERROR in READNML:'
        write(iout,'(a)') 'Simulation start date contains zeros '
        write(iout,'(a,i5)') 'Year:  ',ibyr
        write(iout,'(a,i5)') 'Month: ',ibmo
        write(iout,'(a,i5)') 'Day:   ',ibdy
        call camxerr()
      endif
      call cvtdate(ibyr,ibmo,ibdy,begtim,cbdate,begdate)
c
      ieyr = End_Date_Hour(1)
      iemo = End_Date_Hour(2)
      iedy = End_Date_Hour(3)
      endtim = float(End_Date_Hour(4))
      if (ieyr.eq.0 .or. iemo.eq.0 .or. iedy.eq.0)  then
        write(iout,'(//,a)') 'ERROR in READNML:'
        write(iout,'(a)') 'Simulation end date contains zeros '
        write(iout,'(a,i5)') 'Year:  ',ieyr
        write(iout,'(a,i5)') 'Month: ',iemo
        write(iout,'(a,i5)') 'Day:   ',iedy
        call camxerr()
      endif
      call cvtdate(ieyr,iemo,iedy,endtim,cedate,enddate)
c
      if( enddate .LT. begdate ) then
        write(iout,'(//,a)') 'ERROR in READNML:'
        write(iout,'(a,a)') 'Simulation end date is less than ',
     &                                 'simulation start date.'
        write(iout,'(a,i10.5)') 'Simulation start: ',begdate
        write(iout,'(a,i10.5)') 'Simulation end:   ',enddate
        call camxerr()
      elseif( enddate .EQ. begdate ) then
        if( endtim .LE. begtim ) then
          write(iout,'(//,a)') 'ERROR in READNML:'
          write(iout,'(a,a)') 'Simulation end time is less than ',
     &                                   'simulation start time.'
          write(iout,'(a,f10.0)') 'Simulation start: ',begtim
          write(iout,'(a,f10.0)') 'Simulation end:   ',endtim
          call camxerr()
        endif
      endif
c
c-----Max timestep and I/O frequencies
c
      dtmax = Maximum_Timestep
      dtinp = Met_Input_Frequency
      dtems = Ems_Input_Frequency
      dtout = Output_Frequency
c
      if( dtmax .LE. 1 .OR. dtinp .LE. 1. .or. dtems .LE. 1.
     &                                    .or. dtout .LE. 1.) then
        write(iout,'(//,a)') 'ERROR in READNML:'
        write(iout,'(/,a,/)')
     &            'I/O frequencies must be greater than 1 minute'
        write(iout,'(a,f6.3)')'Maximum time step   (DTMAX): ',dtmax
        write(iout,'(a,f6.3)')'Input interval      (DTINP): ',dtinp
        write(iout,'(a,f6.3)')'Emissions interval  (DTEMS): ',dtems
        write(iout,'(a,f6.3)')'Output interval     (DTOUT): ',dtout
        call camxerr()
      endif
      if( dtmax .GT. MAXDT ) then
        write(iout,'(//,A)') 'WARNING in READNML:'
        write(iout,'(A)') 
     &               'Invalid value specified for maximum time step.'
        write(iout,'(A)') 'A default value will be used instead.'
        write(iout,'(a,f6.0)')'Maximum time step (DTMAX)     : ',dtmax
        write(iout,'(a,f6.0)')'Default value to be used      : ',MAXDT
        dtmax = MAXDT
      endif
      if( (dtmax .GT. 60. .AND. amod(dtmax,60.) .GT. 0. ) .OR.
     &    (dtinp .GT. 60. .AND. amod(dtinp,60.) .GT. 0. ) .OR.
     &    (dtems .GT. 60. .AND. amod(dtems,60.) .GT. 0. ) .OR.
     &    (dtout .GT. 60. .AND. amod(dtout,60.) .GT. 0. )) then
        write(iout,'(//,a)') 'ERROR in READNML:'
        write(iout,'(a)')'An Input/Output interval is > 60 minutes.'
        write(iout,'(a)')'It must be an integer multiple of 60 minutes.'
        write(iout,'(a,f6.0)')'Maximum time step   (DTMAX): ',dtmax
        write(iout,'(a,f6.0)')'Input interval      (DTINP): ',dtinp
        write(iout,'(a,f6.0)')'Emissions interval  (DTEMS): ',dtems
        write(iout,'(a,f6.0)')'Output interval     (DTOUT): ',dtout
        call camxerr()
      endif
      if( (dtmax .LT. 60. .AND. amod(60.,dtmax) .GT. 0. ) .OR.
     &    (dtinp .LT. 60. .AND. amod(60.,dtinp) .GT. 0. ) .OR.
     &    (dtems .LT. 60. .AND. amod(60.,dtems) .GT. 0. ) .OR.
     &    (dtout .LT. 60. .AND. amod(60.,dtout) .GT. 0. )) then
        write(iout,'(//,a)') 'ERROR in READNML:'
        write(iout,'(a)')'An Input/Output interval is < 60 minutes.'
        write(iout,'(a)')'It must divide 60 minutes evenly.'
        write(iout,'(a,f6.0)')'Maximum time step   (DTMAX): ',dtmax
        write(iout,'(a,f6.0)')'Input interval      (DTINP): ',dtinp
        write(iout,'(a,f6.0)')'Emissions interval  (DTEMS): ',dtems
        write(iout,'(a,f6.0)')'Output interval     (DTOUT): ',dtout
        call camxerr()
      endif
      if( amod(amax1(dtinp,dtems),amin1(dtinp,dtems)) .GT. 0. .OR.
     &    amod(amax1(dtinp,dtout),amin1(dtinp,dtout)) .GT. 0. .OR.
     &    amod(amax1(dtems,dtout),amin1(dtems,dtout)) .GT. 0. ) then
        write(iout,'(//,a)') 'ERROR in READNML:'
        write(iout,'(a)')'Input/Output intervals must be even multiples'
        write(iout,'(a)')'of each other.'
        write(iout,'(a,f6.0)')'Input interval      (DTINP): ',dtinp
        write(iout,'(a,f6.0)')'Emissions interval  (DTEMS): ',dtems
        write(iout,'(a,f6.0)')'Output interval     (DTOUT): ',dtout
        call camxerr()
      endif
c
c-----Projection parameters
c
      llatlon = .FALSE.
      lutm    = .FALSE.
      lrpolar = .FALSE.
      lambrt  = .FALSE.
      lpolar  = .FALSE.
      lmerc   = .FALSE.
      call jstlft( Map_Projection )
      call toupper( Map_Projection )
      if( Map_Projection .EQ. 'LATLON    ' ) then
         llatlon = .TRUE.
      elseif( Map_Projection .EQ. 'UTM       ' ) then
         lutm = .TRUE.
      elseif( Map_Projection .EQ. 'RPOLAR    ' ) then
        lrpolar = .TRUE.
      elseif( Map_Projection .EQ. 'LAMBERT   ' ) then
        lambrt = .TRUE.
      elseif( Map_Projection .EQ. 'POLAR     ' ) then
        lpolar = .TRUE.
      elseif( Map_Projection .EQ. 'MERCATOR  ' ) then
        lmerc  = .TRUE.
      else
        write(iout,'(//,a)') 'ERROR in READNML:'
        write(iout,'(3A)') 'Incorrect coordinate ID specified in ',
     &                     'control file: ',Map_Projection
        write(iout,'(1X,A)') 'Acceptable options are:'
        write(iout,'(10X,A)') 'LATLON'
        write(iout,'(10X,A)') 'UTM'
        write(iout,'(10X,A)') 'RPOLAR'
        write(iout,'(10X,A)') 'LAMBERT'
        write(iout,'(10X,A)') 'POLAR'
        write(iout,'(10X,A)') 'MERCATOR'
        call camxerr()
      endif
c
      if( lutm ) then
         iuzon = UTM_Zone
         if( iuzon .eq. 0 ) then
           write(iout,'(//,A)')'ERROR in READNML:'
           write(iout,'(A)')   '  The UTM zone can not be set to zero'
           write(iout,'(2A)')  '  Use +60 for the northern or -60 for',
     &                         ' the southern hemisphere'
           call camxerr()
         endif
      elseif( lrpolar .or. lambrt .or. lpolar .or. lmerc ) then
         polelon = Longitude_Pole
         polelat = Latitude_Pole
         if( lambrt .or. lpolar .or. lmerc ) then
           tlat1 = True_Latitude1
           if ( lambrt ) tlat2 = True_Latitude2
         endif
      endif
c
c-----Number of Grids
c
      ngrid = Number_of_Grids
c
c----Calculate number of emissions files for each grid ---
c
      is_any_grid_emiss = .FALSE.
      allocate(nemiss_files(ngrid))
      do igrd=1,ngrid
         nemiss_files(igrd) = 0
      enddo
      do igrd = 1,ngrid
        do i=1,MXNAM
          if( istrln(Emiss_Grid(igrd,i)) .NE. 0 ) then
            nemiss_files(igrd) = i
            is_any_grid_emiss = .TRUE.
          endif
        enddo
        do i=1,MXNAM
          if( istrln(Point_Sources(i)) .NE. 0 ) npoint_files = i
        enddo
      enddo
c
c----Call routine to allocate arrays for file units ---
c
      call alloc_filunit(ngrid,nemiss_files,MAX(npoint_files,1))
c
      nnest = ngrid - 1
      if( ngrid .LT. 1 ) then
        write(iout,'(//,a)') 'ERROR in READNML:'
        write(iout,'(a)') 'Number of grids must be 1 or more'
        call camxerr()
      endif
c
c----Call routine to allocate arrays for grid definitions ---
c
      call alloc_grid()
c
c-----Master grid parameters
c
      xorg     = Master_SW_XCoord
      yorg     = Master_SW_YCoord
      delx     = Master_Cell_XSize
      dely     = Master_Cell_YSize
      ncol(1)  = Master_Grid_Columns
      nrow(1)  = Master_Grid_Rows
      nlay(1)  = Number_of_Layers
      inst1(1) = 2
      inst2(1) = ncol(1)-1
      jnst1(1) = 2
      jnst2(1) = nrow(1)-1
c
c-----Nested grid parameters
c
      if (ngrid.gt.1) then
        do n = 2,ngrid
          meshold(n) = Nest_Meshing_Factor(n)
          inst1(n)   = Nest_Beg_I_Index(n)
          inst2(n)   = Nest_End_I_Index(n)
          jnst1(n)   = Nest_Beg_J_Index(n)
          jnst2(n)   = Nest_End_J_Index(n)
          nlay(n)    = Number_of_Layers
        enddo
c
        do n = 2,ngrid
          if( inst1(n) .LE. 1 .OR. inst1(n) .GE. ncol(1) .OR.
     &        inst2(n) .LE. 1 .OR. inst2(n) .GE. ncol(1) .OR.
     &        jnst1(n) .LE. 1 .OR. jnst1(n) .GE. nrow(1) .OR.
     &        jnst2(n) .LE. 1 .OR. jnst2(n) .GE. nrow(1) ) then
              write(iout,'(//,a)') 'ERROR in READNML:'
              write(iout,'(a,i5)') 'For grid # ',n
              write(iout,'(2a)') 'Starting/ending indices exceed ',
     &                           'extent of coarse grid.'
              write(iout,'(a,2i5)') 'Nest column range   :',
     &                              inst1(n),inst2(n)
              write(iout,'(a,2i5)') 'Coarse column range :',
     &                              1,ncol(1)
              write(iout,'(a,2i5)') 'Nest row range      :',
     &                              jnst1(n),jnst2(n)
              write(iout,'(a,2i5)') 'Coarse row range    :',
     &                              1,nrow(1)
              call camxerr()
          endif
c
c-----Calculate dimensions for the grid and check for array overflow
c
          ncol(n) = (inst2(n) - inst1(n) + 1)*meshold(n) + 2
          nrow(n) = (jnst2(n) - jnst1(n) + 1)*meshold(n) + 2
        enddo
      endif
c
c-----Model options/flags
c
      call jstlft( Advection_Solver )
      call toupper( Advection_Solver )
      if( Advection_Solver .EQ. 'BOTT      ' ) then
        iadvct = 2
      elseif( Advection_Solver .EQ. 'PPM       ' ) then
        iadvct = 3
      else
        write(iout,'(//,a)') 'ERROR in READNML:' 
        write(iout,'(3A)') 'Incorrect horizontal advection solver ',
     &                'specified in control file: ',Advection_Solver
        write(iout,'(1X,A)')'Acceptable options are:'
        write(iout,'(10X,A)') 'BOTT'
        write(iout,'(10X,A)') 'PPM'
        call camxerr()
      endif
c
      call jstlft( Vadvection_Solver )
      call toupper( Vadvection_Solver )
      if( Vadvection_Solver .EQ. 'IMPLICIT  ' ) then
        lzppm = .false.
      elseif( Vadvection_Solver .EQ. 'PPM       ' ) then
        lzppm = .true.
      else
        write(iout,'(//,a)') 'ERROR in READNML:'
        write(iout,'(3A)') 'Incorrect vertical advection solver ',
     &                'specified in control file: ',Vadvection_Solver
        write(iout,'(1X,A)')'Acceptable options are:'
        write(iout,'(10X,A)') 'IMPLICIT'
        write(iout,'(10X,A)') 'PPM'
        call camxerr()
      endif
c
      call jstlft( Chemistry_Solver )
      call toupper( Chemistry_Solver )
      if( Chemistry_Solver .EQ. CDCMC ) then
          write(iout,'(//,a)') 'ERROR in READNML:'
          write(iout,'(/,1X,2A)')
     &                      'The CMC solver is no longer supported. ',
     &                      'Use the EBI solver instead.'
          call camxerr()
      else if( Chemistry_Solver .EQ. CDEBI ) then
          idsolv = IDEBI
      else if( Chemistry_Solver .EQ. CDLSOD ) then
          idsolv = IDLSOD
      else
          write(iout,'(//,a)') 'ERROR in READNML:'
          write(iout,'(/,1X,6A)') 'Invalid chemistry solver specified ',
     &                              'in control file: ',Chemistry_Solver
          write(iout,'(1X,A)') 'Acceptable options are: '
          write(iout,'(10X,A)') CDEBI
          write(iout,'(10X,A)') CDLSOD
          call camxerr()
      endif
      if( Dry_Deposition ) then
          write(iout,'(//,a)') 'ERROR in READNML:'
          write(iout,'(2A)') 'NOTE: This version uses a different ',
     &               'designation for the dry deposition option.'
          write(iout,'(6X,2A)') 'It is now a character ',
     &                         'variable called: Drydep_Model'
          write(iout,'(10X,A)') 'Acceptable options are: '
          write(iout,'(15X,A)') 'NONE'
          write(iout,'(15X,A)') 'WESELY89 - original scheme'
          write(iout,'(15X,A)') 'ZHANG03  - new scheme'
          write(iout,'(//)')
          call camxerr()
      endif
      call jstlft( Drydep_Model )
      call toupper( Drydep_Model )
      nlu = NLUZ03
      if( Drydep_Model .EQ. 'ZHANG03   ') then
          idrydep = 2
          ldry = .true.
      else if( Drydep_Model .EQ. 'WESELY89  ') then
          idrydep = 1
          ldry = .true.
          nlu = NLUW89
      else if( Drydep_Model .EQ. 'NONE      ' ) then
          idrydep = 0
          ldry = .false.
      else
          write(iout,'(//,a)') 'ERROR in READNML:'
          write(iout,'(2A)') 'NOTE: This version uses a different ',
     &               'designation for the dry deposition option.'
          write(iout,'(6X,2A)') 'It is now a character ',
     &                         'variable called: Drydep_Model'
          write(iout,'(/,1X,6A)') 'Invalid dry deposition option ',
     &                            'specified in control file: ',
     &                            Drydep_Model
          write(iout,'(1X,A)') 'Acceptable options are: '
          write(iout,'(10X,A)') 'NONE'
          write(iout,'(10X,A)') 'ZHANG03'
          write(iout,'(10X,A)') 'WESELY89'
          call camxerr()
      endif
c
      call jstlft( PiG_Submodel )
      call toupper( PiG_Submodel )
      if (PiG_Submodel .EQ. 'GREASD    ') then
        ipigflg = GRESPIG
      elseif (PiG_Submodel .EQ. 'IRON      ') then
        ipigflg = IRONPIG
      elseif (PiG_Submodel .EQ. 'NONE      ') then
        ipigflg = 0
      else
        write(iout,'(//,a)') 'ERROR in READNML:'
        write(iout,'(3A)') 'Incorrect PiG option specified in ',
     &                             'control file: ',PiG_Submodel
        write(iout,'(1X,A)')'Acceptable options are: '
        write(iout,'(10X,A)') 'NONE'
        write(iout,'(10X,A)') 'GREASD'
        write(iout,'(10X,A)') 'IRON'
        call camxerr()
      endif
      if ( Surface_Model .AND. Pig_Submodel .NE. 'NONE      ') then
          write(iout,'(//,a)') 'ERROR in READNML:'
          write(iout,'(2A)') 'The PiG sub-model cannot be used with',
     &                       ' the Surface Model'
          call camxerr()
      endif
c
c======================== Probing Tool Begin ===========================
c
      call jstlft( Probing_Tool )
      call toupper( Probing_Tool )
      if( Probing_Tool .NE. 'NONE      ' .AND.
     &    Probing_Tool .NE. SA     .AND.
     &    Probing_Tool .NE. RTRAC  .AND.
     &    Probing_Tool .NE. RTCMC  .AND.
     &    Probing_Tool .NE. DDM    .AND.
     &    Probing_Tool .NE. HDDM   .AND.
     &    Probing_Tool .NE. STRPA  .AND.
     &    Probing_Tool .NE. STRIPR .AND.
     &    Probing_Tool .NE. STRIRR) then
          write(iout,'(//,A)') 'ERROR in READNML:'
          write(iout,'(/,1X,6A)') 'Invalid technology type ',
     &             'specified in control file: ',Probing_Tool
          write(iout,'(1X,A)') 'Acceptable options are: '
          write(iout,'(10X,A)')'NONE'
          write(iout,'(10X,A)') SA
          write(iout,'(10X,A)') RTRAC
          write(iout,'(10X,A)') RTCMC
          write(iout,'(10X,A)') DDM
          write(iout,'(10X,A)') HDDM
          write(iout,'(10X,A)') STRPA
          write(iout,'(10X,A)') STRIPR
          write(iout,'(10X,A)') STRIRR
          call camxerr()
       endif
c
c  --- check that switches are compatible ---
c
       if( ipigflg .EQ. GRESPIG .AND. 
     &   Probing_Tool .NE. 'NONE' .AND.
     &   Probing_Tool .NE. SA     .AND. Probing_Tool .NE. STRPA .AND.
     &   Probing_Tool .NE. STRIPR .AND. Probing_Tool .NE. STRIRR .AND.
     &   Probing_Tool .NE. RTRAC ) then
          write(iout,'(//,a)') 'ERROR in READNML:'
          write(iout,'(2a)') 'GREASD PiG only works with the following',
     &                       ' Probing Tools:'
          write(iout,'(a)') 'SA,PA,IPR,IRR'
          write(iout,'(a)') 'Turn off PiG switch and try again.'
          call camxerr()
       endif
       if( ipigflg .EQ. IRONPIG .AND.
     &     Probing_Tool .NE. 'NONE' .AND. Probing_Tool .NE. RTRAC .AND.
     &     Probing_Tool .NE. RTCMC) then
          write(iout,'(//,a)') 'ERROR in READNML:'
          write(iout,'(2a)') 'IRON PiG only works with the following',
     &                       ' Probing Tools:'
          write(iout,'(a)') 'RTRAC'
          write(iout,'(a)') 'RTCMC'
          write(iout,'(a)') 'Turn off PiG switch and try again.'
          call camxerr()
       endif
       if( ipigflg .NE. 0 .AND. Probing_Tool .NE. 'NONE' 
     &                                          .AND. LVISPIG ) then
          write(iout,'(//,a)') 'ERROR in READNML:'
          write(iout,'(a)') 'You are running PiG with Probing Tools.'
          write(iout,'(2a)') 'PiG visualization in the average file ',
     &                          'is not supported with Probing Tools.' 
          write(iout,'(2a,/,a)') 'Please set the LVISPIG parameter in ',
     &                        'the camx.prm file to FALSE to run PiG ',
     &                                             'with Probing Tools.'
          write(iout,'(2a)') 'See the CAMx Users Guide for a ',
     &                  'description of the PiG visualization feature.'
          call camxerr()
        endif
cbk       if( iadvct .ne. 2 .AND. Probing_Tool .EQ. DDM ) then
cbk          write(iout,'(//,a)') 'ERROR in READNML:'
cbk          write(iout,'(2a)') 'Advection solver used with DDM ',
cbk     &                                          'must be BOTT.'
cbk          write(iout,*)'Change the control file and try again.'
cbk          call camxerr()
cbk       endif
       if( idsolv .NE. IDEBI .AND. Probing_Tool .EQ. DDM ) then
          write(iout,'(//,a)') 'ERROR in READNML:'
          write(iout,'(2a)') 'Chemistry solver used with DDM ',
     &                                           'must be EBI.'
          write(iout,*)'Change the control file and try again.'
          call camxerr()
       endif
       if( ACM2_Diffusion .AND. ( Probing_Tool .EQ. STRIPR .OR.
     &                            Probing_Tool .EQ. STRPA )) then
          write(iout,'(//,a)') 'ERROR in READNML:'
          write(iout,'(a)') 'PA/IPR does not work with ACM2 diffusion.'
          write(iout,*)'Change the control file and try again.'
          call camxerr()
       endif
       if( Subgrid_Convection .AND. Probing_Tool .NE. 'NONE') then
          write(iout,'(2a)') 'Probing Tools do not work with the',
     &                       ' subgrid convection option.'
          write(iout,*)'Change the control file and try again.'
          call camxerr()
       endif
       if( Subgrid_Convection .AND. Flexi_Nest) then
          write(iout,'(2a)') 'Flexi-nesting is not allowed with the',
     &                       ' subgrid convection option.'
          write(iout,*)'Change the control file and try again.'
          call camxerr()
       endif
       if( Strat_Ozone_Profile .AND. Probing_Tool .NE. 'NONE' .AND.
     &                               Probing_Tool .NE. SA) then
          write(iout,'(2a)') 'The stratospheric ozone profile scheme',
     &                       ' only works with the SA Probing Tool.'
          write(iout,*)'Change the control file and try again.'
          call camxerr()
       endif
c
c======================== Probing Tool End =============================
c 
      lflexi = Flexi_Nest
      ldiag  = Diagnostic_Error_Check
      lrstrt = Restart
      lchem  = Chemistry
      lwet   = Wet_Deposition
      lacm2  = ACM2_Diffusion
      le1day = Ignore_Emission_Dates
      larsrc = Gridded_Emissions
      lptsrc = Point_Emissions
      lsuper = Super_Stepping
      lsrfmod = Surface_Model
      lbidinh3 = Bidi_NH3_Drydep
      lstrato3 = Strat_Ozone_Profile
      do i = 1,ngrid
        lcig(i) = Subgrid_Convection
      enddo
      if( ipigflg .NE. 0 .AND. .NOT.lptsrc ) then
        write(iout,'(//,a)') 'ERROR in READNML:'
        write(iout,'(a)')
     &         'Point source emissions are required for PiG module.'
        write(iout,'(2a)')'Either turn off PiG module or turn on ',
     &                                    'point sources treatment.'
        call camxerr()
      endif
      lixemis   = .false.
      lixbypass = .false.
      call jstlft( Inline_Ix_Emissions )
      call toupper( Inline_Ix_Emissions )
      if (Inline_Ix_Emissions .EQ. 'TRUE      ') then
        lixemis = .true.
      elseif (Inline_Ix_Emissions .EQ. 'FALSE     ') then
        lixemis = .false.
      elseif (Inline_Ix_Emissions .EQ. 'BYPASS    ') then
        lixbypass = .true.
      else
        write(iout,'(//,a)') 'ERROR in READNML:'
        write(iout,'(3A)') 'Incorrect Ix emission option specified in ',
     &                             'control file: ',Inline_Ix_Emissions
        write(iout,'(1X,A)')'Acceptable options are: '
        write(iout,'(10X,A)') 'TRUE'
        write(iout,'(10X,A)') 'FALSE'
        write(iout,'(10X,A)') 'BYPASS'
        call camxerr()
      endif
c
c-----Output options and species
c
      do i=1,ngrid
         l3davg(i) = Average_Output_3D
         if( Output_3D_Grid(i) ) l3davg(i) = .TRUE.
      enddo
      lcdfout = NetCDF_Format_Output
      ncf_compress = NetCDF_Use_Compression
      call ncf_set_compress_flag()
      lgas_out_ppm = Output_Gas_Concs_PPM
c
c-----Call routine to  set the NetCDF cache parameters ---
c
      if( lcdfout .AND. ncf_compress ) call ncf_set_cache()
c
      do n = MXNAM,1,-1
        if (Output_Species_Names(n).ne.' ') then
          navspc = n
          call alloc_chmstry_avg(navspc)
          do l = 1,navspc
            if (Output_Species_Names(n).eq.' ') then
              write(iout,'(//,a)') 'ERROR in READNML:'
              write(iout,'(a)')'An output species name is blank'
              write(iout,'(a,i5)')'at species number: ',l
              call camxerr()
            endif 
            spavg(l) = Output_Species_Names(l)
            if ((spavg(l).eq.'ALL' .or. spavg(l).eq.'ALLR') .and.
     &           navspc.ne.1) then
              write(iout,'(//,a)') 'ERROR in READNML:'
              write(iout,'(3a)')'If using "ALL" or "ALLR" as an',
     &              ' Output_Species_Name it must be the the',
     &              ' only name listed.'
              call camxerr()
            endif
          enddo
          goto 50
        endif
      enddo
      write(iout,'(//,a)') 'ERROR in READNML:'
      write(iout,'(a)')'Number of average species is zero.'
      write(iout,'(a)')'At least 1 species needs to be output.'
      call camxerr()
 50   continue
c
c-----PiG sampling grids
c
      nsample = 0
      lsample = PiG_Sampling_Grid
      if( .NOT. lsample ) goto 100
      if( ipigflg .EQ. 0 ) then
        write(iout,'(//,a)') 'ERROR in READNML:'
        write(iout,'(/,1X,A,A)')'PiG flag must be turned on',
     &                           ' to use sampling grids.'
        call camxerr()
      endif
      lbckgrd = Sample_Background
c
      if( lmpi .AND. lbckgrd ) then
        write(iout,'(//,a)') 'ERROR in READNML:'
        write(iout,'(/,1X,2A)') 'Adding background to the sampling ',
     &                                     'grid is disabled with MPI.'
        write(iout,'(/,1X,2A)') 'You must either set Sample_Background ',
     &                                    'to false or run without MPI.'
        call camxerr()
      endif
c
      nsample = Number_of_Sampling_Grids
      if (nsample .EQ. 0) then
        write(iout,'(//,a)') 'ERROR in READNML:'
        write(iout,'(/,1X,A,A)')'Sampling grids are turned on ',
     &                'but number of sampling grids is set to zero'
        call camxerr()
      endif
c
c----Call routine to allocate arrays for sampling file units ---
c
      call alloc_filunit_sample(nsample)
c
c   --- allocate the arrays that depend on number of samples ---
c
      call alloc_pigsty_sample()
c
      do n = 1,nsample
        ismp1(n)   = SG_Beg_I_Index(n)
        ismp2(n)   = SG_End_I_Index(n)
        jsmp1(n)   = SG_Beg_J_Index(n)
        jsmp2(n)   = SG_End_J_Index(n)
        meshsmp(n) = SG_Mesh_Factor(n)
      enddo
      do n = 1,nsample
        if (ismp1(n) .LE. 1 .OR. ismp1(n) .GE. ncol(1) .OR.
     &      ismp2(n) .LE. 1 .OR. ismp2(n) .GE. ncol(1) .OR.
     &      jsmp1(n) .LE. 1 .OR. jsmp1(n) .GE. nrow(1) .OR.
     &      jsmp2(n) .LE. 1 .OR. jsmp2(n) .GE. nrow(1)) then
          write(iout,'(//,a)') 'ERROR in READNML:'
          write(iout,'(a,i3)') 'For sampling grid # ',n
          write(iout,'(2a)') 'Starting/ending indices exceed extent',
     &                       ' of coarse grid.'
          write(iout,'(a,2i5)')
     &                   'Sample column range   :',ismp1(n),ismp2(n)
          write(iout,'(a,2i5)')
     &                   'Sample row range      :',jsmp1(n),jsmp2(n)
          write(iout,'(a,2i5)') 'Coarse column range   :',1,ncol(1)
          write(iout,'(a,2i5)') 'Coarse row range      :',1,nrow(1)
          call camxerr()
        endif
        ncolsmp(n) = (ismp2(n) - ismp1(n) + 1)*meshsmp(n)
        nrowsmp(n) = (jsmp2(n) - jsmp1(n) + 1)*meshsmp(n)
        xorgsmp(n) = xorg + (ismp1(n) - 1)*delx
        yorgsmp(n) = yorg + (jsmp1(n) - 1)*dely
        if( ncolsmp(n) .LE. 0 .OR. nrowsmp(n) .LE. 0 ) then
          write(iout,'(//,a)') 'ERROR in READNML:'
          write(iout,'(a,i3)') 'For sampling grid # ',n
          write(iout,'(2a)') 'Starting/ending indices are not valid.'
          write(iout,'(a,2i5)')
     &                   'Sample column range   :',ismp1(n),ismp2(n)
          write(iout,'(a,2i5)')
     &                   'Sample row range      :',jsmp1(n),jsmp2(n)
          call camxerr()
        endif
      enddo
      do n = 1,nsample
        ismpgrd(n) = 1
        do 101 ng = 2,ngrid
c
c  --- check for across the seam of this grid ---
c
          lacross = .FALSE.
          if( ismp1(n) .LT. inst1(ng) .AND. 
     &                    ismp2(n) .GT. inst1(ng) ) lacross = .TRUE. 
          if( ismp1(n) .LT. inst2(ng) .AND. 
     &                    ismp2(n) .GT. inst2(ng) ) lacross = .TRUE. 
          if( jsmp1(n) .LT. jnst1(ng) .AND. 
     &                    jsmp2(n) .GT. jnst1(ng) ) lacross = .TRUE. 
          if( jsmp1(n) .LT. jnst2(ng) .AND. 
     &                    jsmp2(n) .GT. jnst2(ng) ) lacross = .TRUE. 
          if( lacross ) then
             write(iout,'(//,a)') 'ERROR in READNML:'
             write(iout,'(a,i3)') 'For sampling grid # ',n
             write(iout,'(2a,I3)') 'Sampling grid crosses the seam of ',
     &                       'grid # ',ng
             write(iout,'(a,2i5)')
     &                   'Sample column range   :',ismp1(n),ismp2(n)
             write(iout,'(a,2i5)')
     &                   'Sample row range      :',jsmp1(n),jsmp2(n)
             write(iout,'(a,2i5)') 'Grid column range   :',
     &                                             inst1(ng),inst2(ng)
             write(iout,'(a,2i5)') 'Grid row range      :',
     &                                             jnst1(ng),jnst2(ng)
             call camxerr()
          endif
          if (inst1(ng).le.ismp1(n) .and. inst2(ng).ge.ismp2(n) .and.
     &        jnst1(ng).le.jsmp1(n) .and. jnst2(ng).ge.jsmp2(n)) then
            ismpgrd(n) = ng
c
c  --- make sure the sampling grid is finer than current grid ---
c
            if( meshold(ng) .GT. meshsmp(n) ) then
             write(iout,'(//,a,/)') 'ERROR in READNML:'
             write(iout,'(a,i3)') 'For sampling grid # ',n
             write(iout,'(2a)') 'The modeling grid is finer than the ',
     &                                                   'sampling grid.'
             write(iout,'(2a,/)') 'Aggregation for sampling grids is ',
     &                                                  'not supported.'
             write(iout,'(a,i3,a,i3)') 'Sampling grid: ',n,
     &                             ' is contained in modeling grid: ',ng
             write(iout,'(a,i3,a,i3)') 'Sampling grid: ',n,
     &                               ' has meshing factor: ',meshsmp(n)
             write(iout,'(a,i3,a,i3)') 'Modeling grid: ',ng,
     &                               ' has meshing factor: ',meshold(ng)
             write(iout,'(2a,/,a)') 'Set the meshing factor for your ',
     &                              'sampling grid to be at least as ',
     &                       'large as the value for the modeling grid.'
             call camxerr()
            endif
          endif
 101    continue
      enddo
 100  continue
c
c-----Echo run control parameters
c
      write(idiag,'(A,F10.0,I10.6,I10.5)')
     &         'Simulation start time/date : ',begtim,cbdate,begdate
      write(idiag,'(A,F10.0,I10.6,I10.5)')
     &         'Simulation end time/date   : ',endtim,cedate,enddate
      write(idiag,'(A,I10)')
     &         'Time zone                  : ',itzon
      write(idiag,'(A,F10.0)')
     &         'Max timestep (min)         : ',dtmax
      write(idiag,'(A,F10.0)')
     &         'Met Input interval (min)   : ',dtinp
      write(idiag,'(A,F10.0)')
     &         'Emiss Input interval (min) : ',dtems
      write(idiag,'(A,F10.0)')
     &         'Output interval (min)      : ',dtout
      write(idiag,'(A,A)')
     &         'Grid Projection Type       : ',Map_Projection
      if (lutm) then
        write(idiag,'(A,I10)')
     &         'UTM: zone number           : ',iuzon
      elseif (lrpolar) then
        write(idiag,'(A,2F10.3)')
     &         'RPOLAR: pole lon/lat        : ',polelon,polelat
      elseif (lambrt) then
        write(idiag,'(A,4F10.3)')
     &         'LAMBERT: pole lon/lat      : ',polelon,polelat
        write(idiag,'(A,4F10.3)')
     &         'LAMBERT: true latitudes    : ',tlat1,tlat2
      elseif (lpolar) then
        write(idiag,'(A,4F10.3)')
     &         'POLAR: pole lon/lat        : ',polelon,polelat
        write(idiag,'(A,4F10.3)')
     &         'POLAR: true latitude       : ',tlat1
      elseif (lmerc) then
        write(idiag,'(A,4F10.3)')
     &         'MERCATOR: pole lon/lat     : ',polelon,polelat
        write(idiag,'(A,4F10.3)')
     &         'MERCATOR: true latitude    : ',tlat1
      endif
      write(idiag,'(A,4F10.3)')
     &         'Master Grid SW Corner X/Y  : ',xorg,yorg
      write(idiag,'(A,4F10.3)')
     &         'Master Grid cell size      : ',delx,dely
      write(idiag,'(A,3I10)')
     &         'Master Grid NCOL NROW NLAY : ',ncol(1),nrow(1),nlay(1)
      write(idiag,*)
      write(idiag,'(A)')'CAMx control flags'
      write(idiag,'(A,L10)') 
     &         'Stop after diagnostic check: ',ldiag
      write(idiag,'(A,A)')
     &         'Horizontal Advection Solver: ',Advection_Solver
      write(idiag,'(A,A)')
     &         'Vertical Advection Solver  : ',Vadvection_Solver
      write(idiag,'(A,A)')
     &         'Chemistry Solver           : ',Chemistry_Solver
      write(idiag,'(A,A)')
     &         'Dry Deposition             : ',Drydep_Model
      write(idiag,'(A,A)')
     &         'PiG submodel               : ',PiG_Submodel
      write(idiag,'(A,L10)')
     &         'PiG Sampling Grid          : ',PiG_Sampling_Grid
      write(idiag,'(A,A)')
     &         'Probing Tools              : ',Probing_Tool
      write(idiag,'(A,L10)') 
     &         'Flexi-Nest                 : ',lflexi 
      write(idiag,'(A,L10)') 
     &         'Restart                    : ',lrstrt 
      write(idiag,'(A,L10)') 
     &         'Chemistry                  : ',lchem
      write(idiag,'(A,L10)') 
     &         'Dry deposition             : ',ldry 
      write(idiag,'(A,L10)') 
     &         'Wet deposition             : ',lwet 
      write(idiag,'(A,L10)')
     &         'ACM2 Vertical Diffusion    : ',lacm2
      write(idiag,'(A,L10)') 
     &         'Super stepping advection   : ',lsuper
      write(idiag,'(A,L10)') 
     &         'Surface Model              : ',lsrfmod
      write(idiag,'(A,L10)') 
     &         'Subgrid Convection         : ',Subgrid_Convection
      write(idiag,'(A,L10)') 
     &         'Inline Ix Emissions        : ',lixemis
      write(idiag,'(A,L10)') 
     &         'Bypass Ix Emissions        : ',lixbypass
      write(idiag,'(A,L10)') 
     &         'Bidi NH3 Drydep            : ',lbidinh3
      write(idiag,'(A,L10)') 
     &         'Strat Ozone Profile        : ',lstrato3
      write(idiag,'(A,L10)') 
     &         'Area sources               : ',larsrc
      write(idiag,'(A,L10)') 
     &         'Point sources              : ',lptsrc
      write(idiag,'(A,L10)') 
     &         'Date-insensitive emissions : ',le1day
      do i=1,ngrid
         write(idiag,'(A,I3,L10)') 
     &         '3-D average file - Grid #  : ',i,l3davg(i)
      enddo
      write(idiag,'(A,L10)') 
     &         'NetCDF output format          : ',lcdfout
      write(idiag,*)
c
      if( nnest .GT. 0 ) then
        write(idiag,'(A,I3)')'Number of nested fine grids: ',nnest
        write(idiag,*)' Nest        x-range       ncol        y-range',
     &                '       nrow   mesh factor  nlay'      
        do n = 2,ngrid
          write(idiag,'(I5,2(5X,I5,1X,I5,5X,I5),5X,I5,5X,I5)')
     &          n,inst1(n),inst2(n),ncol(n),jnst1(n),jnst2(n),
     &          nrow(n),meshold(n),nlay(n)
        enddo
        write(idiag,*)
        write(idiag,*) '|',('-',i=1,74),'|'
        write(idiag,*) '|',(' ',i=1,74),'|'
        write(idiag,'(2A)') ' | NOTE:  The nest order listed above ',
     &       'is the original order specified in    |'
        write(idiag,'(2A)') ' |        the CAMx.in file.           ',
     &       '                                      |'
        write(idiag,'(2A)') ' |        CAMx may re-order the nests.',
     &      ' See the internal nest order provided |'
        write(idiag,'(2A)') ' |        in the table below.         ',
     &       '                                      |'
        write(idiag,*) '|',(' ',i=1,74),'|'
        write(idiag,*) '|',('-',i=1,74),'|'
      else
        write(idiag,*)'Fine grid nests are not specified.'
      endif
      write(idiag,*)
c
      if( lsample ) then
        write(idiag,*)
        write(idiag,'(A,L10)') 
     &         'Sample grid includes bckgrd: ',lbckgrd
        write(idiag,'(A,I3)')'Number PiG Sampling grids  : ',nsample
        write(idiag,*) ' ID     x-range     ncol     y-range',
     &                '     nrow   mesh factor  delx    dely'
        do n = 1,nsample
          write(idiag,'(1X,I2,2(3X,I4,3X,I4,3X,I4),6X,I4,4X,
     &                  2(2X,F6.3))')
     &          n,ismp1(n),ismp2(n),ncolsmp(n),jsmp1(n),jsmp2(n),
     &          nrowsmp(n),meshsmp(n),delx/float(meshsmp(n)),
     &          dely/float(meshsmp(n))
        enddo
        write(idiag,*)
      endif
c
c-----Open the remaining I/O files
c
      call openfils(iifroot,nopen)
c
c-----Read chemistry parameters
c
      time = begtim
      date = begdate
      call readchm()
c
c-----Check number of average output species
c
      if( navspc .GT. nspec ) then
        write(iout,'(//,A)') 'ERROR in READNML:'
        write(iout,*) 'Number of average species to be output ',
     &               'is greater than number of species to model'
        write(iout,*) 'Average species (NAVSPC): ',navspc
        write(iout,*) 'Model species (NSPEC)   : ',nspec
        call camxerr()
      endif
c
c-----Map average species to model species
c
      write(idiag,*)
      ndepspc = 0
      if (spavg(1).eq.'ALL') then
        navspc = nspec - nrad
        lav = 0
        do l = nrad+1,nspec
          lav = lav + 1
          lavmap(lav) = l
          ndepspc = ndepspc + 1
          ldepmap(ndepspc) = l
          write(idiag,'(2(A,I5,2X,A))')
     &                   'Average species ',lav,spname(l),
     &                   ' mapped to modeled species ',l,spname(l)
        enddo
      elseif (spavg(1).eq.'ALLR') then
        navspc = nspec
        lav = 0
        do l = 1,nspec
          lav = lav + 1
          lavmap(lav) = l
          ndepspc = ndepspc + 1
          ldepmap(ndepspc) = l
          write(idiag,'(2(A,I5,2X,A))')
     &                   'Average species ',lav,spname(l),
     &                   ' mapped to modeled species ',l,spname(l)
        enddo
      else
        do lav = 1,navspc
          lavmap(lav) = 0
          do l = 1,nspec
            if( spavg(lav) .EQ. spname(l) ) then
              lavmap(lav) = l
              ndepspc = ndepspc + 1
              ldepmap(ndepspc) = l
              write(idiag,'(2(A,I5,2X,A))')
     &                   'Average species ',lav,spavg(lav),
     &                   ' mapped to modeled species ',l,spname(l)
            endif
          enddo
          if( lavmap(lav) .EQ. 0 ) then
            write(iout,'(//,a)') 'ERROR in STARTUP:'
            write(iout,*) 'Did not find average species: ',
     &            spavg(lav)(:istrln(spavg(lav))),' in chemparam list.'
            write(iout,*)'Either remove this name from the ',
     &                                    'average species list in your'
            write(iout,*)'CAMx control file, or use the appropriate ',
     &                                               'chemparam file.'
            if( aeropt .EQ. 'CMU' ) then
              write(iout,*) 'If you are using the PM sectional model ',
     &                                     'you have to include the '
              write(iout,*) 'section number in the species name.'
            endif
            call camxerr()
          endif
        enddo
      endif
      write(idiag,*)
c
c  --- display warning if running CMU, in case this has OMP ---
c
      if( aeropt.EQ.'CMU' ) then
          write(*,'(//,36(2A))') ('/\\',i=1,36)
          write(*,'(72A)') '/',(' ',i=1,70),'\\'
          write(*,'(A,25X,A,26X,A)') '/','WARNING in READNML:','\\'
          write(*,'(72A)') '/',(' ',i=1,70),'\\'
          write(*,'(A,2X,2A)') '/',
     &    'You have selected CMU chemistry. CAMx does not support OMP with     ','\\'
          write(*,'(A,2X,2A)') '/',
     &    'the CMU module. Please make sure you are not running OMP with CMU.  ','\\'
          write(*,'(72A)') '/',(' ',i=1,70),'\\'
          write(*,'(36(2A),/)') ('/\\',i=1,36)
          write(idiag,'(//,36(2A))') ('/\\',i=1,36)
          write(idiag,'(72A)') '/',(' ',i=1,70),'\\'
          write(idiag,'(A,25X,A,26X,A)') '/','WARNING in READNML:','\\'
          write(idiag,'(72A)') '/',(' ',i=1,70),'\\'
          write(idiag,'(A,2X,2A)') '/',
     &    'You have selected CMU chemistry. CAMx does not support OMP with     ','\\'
          write(idiag,'(A,2X,2A)') '/',
     &    'the CMU module. Please make sure you are not running OMP with CMU.  ','\\'
          write(idiag,'(72A)') '/',(' ',i=1,70),'\\'
          write(idiag,'(36(2A),/)') ('/\\',i=1,36)
      endif
c
c  --- read Probing Tool namelists ---
c
c
c   --- allocate the arrays for PiG --
c
      if( ipigflg .NE. 0 .AND. .NOT. lrstrt )
     &    call alloc_pigsty(nspec,nreactr,ngrid,
     &                                  MAX(1,numprocs-1),ipigflg)
      if( ipigflg .NE. 0 ) 
     &         call alloc_pigsty_vpconc(ncol,nrow,nlay,nspec,ngrid)
c
c  --- setup some pointers that depend on number of species ---
c
      if( lsample ) then
        ipsmp(1) = 1
        do i = 2,nsample
          ipsmp(i) = ipsmp(i-1) + ncolsmp(i-1)*nrowsmp(i-1)*navspc
        enddo
c
c   --- allocate the arrays for sampling that depend on species ---
c
         call alloc_pigsty_smpgrd(navspc,nsample,
     &                                    nrowsmp,ncolsmp,nsmpcels)
      endif
c
c======================== Probing Tool Begin ===========================
c
c  --- set Probing Tool switches according to technology type ---
c
      call alloc_lddmcalc(ngrid)
      ltrace = .FALSE.
      lddm   = .FALSE.
      lhddm  = .FALSE.
      lproca = .FALSE.
      lipr   = .FALSE.
      lirr   = .FALSE.
      lpartial = .FALSE.
      lptdepout = .FALSE.
      lsfcfl = .FALSE.
      tectyp = Probing_Tool
      lmech1_ok = .FALSE.
      lmech3_ok = .FALSE.
      lmech4_ok = .FALSE.
      lmech5_ok = .FALSE.
      lmech6_ok = .FALSE.
      lmech7_ok = .FALSE.
      lmech10_ok = .FALSE.
      if( tectyp .EQ. 'NONE      '  ) then
c
c---- call routines to allocate the arrays that need a pointer ---
c
         call alloc_tracer_null(nspec,lirr,ngrid,ncol,nrow)
         call alloc_procan_null()
         goto 8888
      elseif( tectyp .EQ. SA    ) then
         ltrace = .true.
         lsfcfl = .TRUE.
         verson = VERSA
         call alloc_procan_null()
         lmech1_ok = .TRUE.
         lmech3_ok = .TRUE.
         lmech4_ok = .TRUE.
         lmech5_ok = .TRUE.
         lmech6_ok = .TRUE.
         lmech7_ok = .TRUE.
         if (lvbs) goto 7009
         lmineral = .FALSE.
         ldmschm = .FALSE.
         if (kpfe.lt.nspec+1 .and. kpmn.lt.nspec+1 .and.
     &       kpmg.lt.nspec+1 .and. kpk.lt.nspec+1 .and.
     &       kpca.lt.nspec+1 .and. kpal.lt.nspec+1 .and.
     &       kpsi.lt.nspec+1 .and. kpti.lt.nspec+1) lmineral = .TRUE.
         if (kdms.lt.nspec+1) ldmschm = .TRUE.
         lddmcalc = .TRUE.
      elseif( tectyp .EQ. RTRAC ) then
         ltrace = .true.
         lsfcfl = .TRUE.
         verson = VERRTRAC
         call alloc_procan_null()
         lmech1_ok = .TRUE.
         lmech3_ok = .TRUE.
         lmech4_ok = .TRUE.
         lmech5_ok = .TRUE.
         lmech6_ok = .TRUE.
         lmech7_ok = .TRUE.
         lddmcalc = .TRUE.
      elseif( tectyp .EQ. RTCMC ) then
         ltrace = .true.
         lsfcfl = .TRUE.
         verson = VERRTCMC
         call alloc_procan_null()
         lmech1_ok = .TRUE.
         lmech3_ok = .TRUE.
         lmech4_ok = .TRUE.
         lmech5_ok = .TRUE.
         lmech6_ok = .TRUE.
         lmech7_ok = .TRUE.
         lddmcalc = .TRUE.
      elseif( tectyp .EQ. DDM  ) then
         lddm = .TRUE.
         lsfcfl = .TRUE.
         verson = VERDDM
         call alloc_procan_null()
         lmech1_ok = .TRUE.
         lmech3_ok = .TRUE.
         lmech4_ok = .TRUE.
         lmech5_ok = .TRUE.
         lmech6_ok = .TRUE.
         lmech7_ok = .TRUE.
         if (lvbs) goto 7009
      elseif( tectyp .EQ. HDDM  ) then
         lhddm = .TRUE.
         lsfcfl = .TRUE.
         verson = VERHDDM
         call alloc_procan_null()
         lmech1_ok = .TRUE.
         lmech3_ok = .TRUE.
         lmech4_ok = .TRUE.
         lmech5_ok = .TRUE.
         lmech6_ok = .TRUE.
         lmech7_ok = .TRUE.
         if (lvbs) goto 7009
      elseif( tectyp .EQ. STRPA ) then
         lproca = .TRUE.
         lipr = .TRUE.
         lirr = .TRUE.
         lsfcfl = .TRUE.
         verson = VERPA
         call alloc_procan_init(ngrid,MXCPA)
         call alloc_ptwet_null()
         lmech1_ok = .TRUE.
         lmech3_ok = .TRUE.
         lmech4_ok = .TRUE.
         lmech5_ok = .TRUE.
         lmech6_ok = .TRUE.
         lmech7_ok = .TRUE.
         if (lvbs) goto 7009
         lddmcalc = .TRUE.
      elseif( tectyp .EQ. STRIPR ) then
         lproca = .TRUE.
         lipr = .TRUE.
         verson = VERPA
         call alloc_tracer_null(nspec,lirr,ngrid,ncol,nrow)
         call alloc_procan_init(ngrid,MXCPA)
         call alloc_ptwet_null()
         lmech1_ok = .TRUE.
         lmech3_ok = .TRUE.
         lmech4_ok = .TRUE.
         lmech5_ok = .TRUE.
         lmech6_ok = .TRUE.
         lmech7_ok = .TRUE.
         lmech10_ok = .TRUE.
         lddmcalc = .TRUE.
      elseif( tectyp .EQ. STRIRR ) then
         lproca = .TRUE.
         lirr = .TRUE.
         lsfcfl = .TRUE.
         verson = VERPA
         call alloc_procan_init(ngrid,MXCPA)
         call alloc_ptwet_null()
         lmech1_ok = .TRUE.
         lmech3_ok = .TRUE.
         lmech4_ok = .TRUE.
         lmech5_ok = .TRUE.
         lmech6_ok = .TRUE.
         lmech7_ok = .TRUE.
         if (lvbs) goto 7009
         lddmcalc = .TRUE.
      endif
      if( idmech.EQ.1 .AND. .NOT. lmech1_ok ) then
         write(iout,'(//,a)') 'ERROR in READNML:'
         write(iout,'(/,2A)') ' You have selected a probing tool that is ',
     &               'currently not supported with Mechanism 1.'
         write(iout,'(/,2A)') ' Turn off probing tools or select ',
     &                                         'another mechanism.'
         call camxerr()
      endif
      if( idmech.EQ.3 .AND. .NOT. lmech3_ok ) then
         write(iout,'(//,a)') 'ERROR in READNML:'
         write(iout,'(/,2A)') ' You have selected a probing tool that is ',
     &               'currently not supported with Mechanism 3.'
         write(iout,'(/,2A)') ' Turn off probing tools or select ',
     &                                         'another mechanism.'
         call camxerr()
      endif
      if( idmech.EQ.4 .AND. .NOT. lmech4_ok ) then
         write(iout,'(//,a)') 'ERROR in READNML:'
         write(iout,'(/,2A)') ' You have selected a probing tool that is ',
     &               'currently not supported with Mechanism 4.'
         write(iout,'(/,2A)') ' Turn off probing tools or select ',
     &                                         'another mechanism.'
         call camxerr()
      endif
      if( idmech.EQ.5 .AND. .NOT. lmech5_ok ) then
         write(iout,'(//,a)') 'ERROR in READNML:'
         write(iout,'(/,2A)') ' You have selected a probing tool that is ',
     &               'currently not supported with Mechanism 5.'
         write(iout,'(/,2A)') ' Turn off probing tools or select ',
     &                                         'another mechanism.'
         call camxerr()
      endif
      if( idmech.EQ.6 .AND. .NOT. lmech6_ok ) then
         write(iout,'(//,a)') 'ERROR in READNML:'
         write(iout,'(/,2A)') ' You have selected a probing tool that is ',
     &               'currently not supported with Mechanism 6.'
         write(iout,'(/,2A)') ' Turn off probing tools or select ',
     &                                         'another mechanism.'
         call camxerr()
      endif
      if( idmech.EQ.7 .AND. .NOT. lmech7_ok ) then
         write(iout,'(//,a)') 'ERROR in READNML:'
         write(iout,'(/,2A)') ' You have selected a probing tool that is ',
     &               'currently not supported with Mechanism 7.'
         write(iout,'(/,2A)') ' Turn off probing tools or select ',
     &                                         'another mechanism.'
         call camxerr()
      endif
      if( idmech.eq.10 .AND. .NOT. lmech10_ok ) then
         write(iout,'(//,a)') 'ERROR in READNML:'
         write(iout,'(/,2A)') ' You have selected a probing tool that is ',
     &               'currently not supported with Mechanism 10.'
         write(iout,'(/,2A)') ' Turn off probing tools or select ',
     &                                         'another mechanism.'
        call camxerr()
      endif
c
c  --- CMU chemistry and probing tools is not allowed ---
c
      if( (aeropt.EQ.'INERT' .OR. aeropt.EQ.'CMU' .OR.
     &     lvbs) .AND. tectyp.NE.'NONE      ') then
          write(iout,'(//,a)') 'ERROR in READNML:'
          write(iout,'(/,2A)') ' Probing Tools are ',
     &     'not supported with the INERT, VBS or CMU Aerosol Schemes.'
          write(iout,'(2A)')' Either set Probing Tools = NONE or ',
     &                'choose the CF Aerosol Scheme with SOAP chemistry.'
          call camxerr()
      endif
c
c  --- DDM/HDDM is not allowed with EQSAM ---
c
      if( leqsam .AND. (lddm .OR. lhddm) ) then
          write(iout,'(//,a)') 'ERROR in READNML:'
          write(iout,'(/,2A)') ' DDM/HDDM is ',
     &     'not supported with the EQSAM aerosol option.'
          call camxerr()
      endif
c
      if( ltrace ) then
         if( tectyp .EQ. RTRAC .OR. tectyp .EQ. RTCMC ) then
            namegrp = 'RT_Control'
            action =  'You must have the '//namegrp(:istrln(namegrp))//
     &          ' namelist group when probing tools is set to '//
     &                                    tectyp(:istrln(tectyp))//'.'
            read(inp,RT_Control,END=7100,ERR=7101)
            flrtsa = RT_File_Root
         else
            namegrp = 'SA_Control'
            action =  'You must have the '//namegrp(:istrln(namegrp))//
     &          ' namelist group when probing tools is set to '//
     &                                    tectyp(:istrln(tectyp))//'.'
            read(inp,SA_Control,END=7100,ERR=7103)
            flrtsa = SA_File_Root
            lptdepout = SA_Deposition_Output
            if( .NOT. ldry .AND. .NOT. lwet ) lptdepout = .FALSE.
         endif 
      elseif( lddm .OR. lhddm ) then
         namegrp = 'DDM_Control'
         action =  'You must have the '//namegrp(:istrln(namegrp))//
     &          ' namelist group when probing tools is set to '//
     &                                    tectyp(:istrln(tectyp))//'.'
         read(inp,DDM_Control,END=7100,ERR=7101)
         flrtsa = DDM_File_Root
      elseif( lproca ) then
         namegrp = 'PA_Control'
         action =  'You must have the '//namegrp(:istrln(namegrp))//
     &          ' namelist group when probing tools is set to '//
     &                                    tectyp(:istrln(tectyp))//'.'
         read(inp,PA_Control,END=7100,ERR=7101)
         flrtsa = PA_File_Root
      endif
c
c======================== Probing Tool End =============================
c 
c
c======================== Source Apportion Begin =======================
c
      if( ltrace ) then
c
c  --- call routines to get the rest of the options ---
c
c
         if( tectyp .EQ. RTRAC .OR. tectyp .EQ. RTCMC) then
            lsrfmodrt = RT_Surface_Model
            lparttn = RT_Partitioning
            if (lsrfmodrt .AND. idrydep .NE. 1) then
              write(iout,*)'The RTRAC Surface Model must be used',
     &                     ' with the WESELY89 deposition'
              write(iout,*)'Stopping'
              call camxerr()
            endif
            if (lsrfmodrt .AND. ipigflg.ne.0) then
              write(iout,*)'The RTRAC Surface Model cannot be used',
     &                     ' with the PiG model'
              write(iout,*)'Stopping'
              call camxerr()
            endif
            if (lparttn .AND. (aeropt.ne.'CF' .or. lvbs)) then
              write(iout,*)'RTRAC gas-aerosol partitioning must',
     &                     ' use the CF/SOAP aerosol scheme'
              write(iout,*)'Stopping'
              call camxerr()
            endif
            write(idiag,*)
            write(idiag,'(A,L10)')
     &                  'RTRAC Surface Model        : ',lsrfmodrt
            write(idiag,'(A,L10)')
     &                  'RTRAC Gas-PM Partitioning  : ',lparttn
c
            if( lsample ) then
               lsmptrc = RT_PiG_Sample
               write(idiag,*)
               write(idiag,'(A,L10)')
     &                  'Sample grid includes RTRAC : ',lsmptrc
               if( lsmptrc ) then
c
c  --- call routine to allocate the variables for RTRAC sampling ---
c
                  call alloc_tracer_sample_io(nsample,0)
                  iprtsmp(1) = 1
                  do i=2,nsample
                    iprtsmp(i) = iprtsmp(i-1) +
     &                           ncolsmp(i-1)*nrowsmp(i-1)*ntotsp
                  enddo
               endif
            endif
            write(idiag,*)
            call startrt(iout,nopen)
         else
            call rdoptsa()
            lptoverride = SA_PT_Override
            nemiss = ngroup
            if( leftovr_area ) nemiss = ngroup + 1
            write(idiag,'(A,L10)') 
     &           'Output Deposition Fields            : ',lptdepout
            write(idiag,'(A,L10)') 
     &           'Stratify Boundary                   : ',lbndry
            write(idiag,'(A,I10)') 
     &           'Number of source areas              : ',nregin
            write(idiag,'(A,I10)') 
     &           'Number of source groups             : ',nemiss
            write(idiag,'(A,L10)') 
     &           'Leftover group for griddded sources : ',leftovr_area
            write(idiag,'(A,I10)') 
     &           'Number of time releases             : ',ntrtim
            write(idiag,'(A,L10)')
     &           'Using Partial Source Area map       : ',lpartial
            write(idiag,*)
            call startsa(iout,nopen)
         endif
      endif
c
c======================== Source Apportion End =========================
c
c
c======================== DDM Begin ====================================
c
      if( lddm .OR. lhddm ) then
         call alloc_ddm_species(nspec)
         call rdoptddm()
         lptoverride = DDM_PT_Override
         write(idiag,'(A,I10)') 
     &            'Number of IC groups                 : ',nicddm
         if( nicddm .GT. 0 ) then
            do i = 1,nicddm
              write(idiag,'(24X,i3,'': '',A)') i,icddmsp(i)
            enddo
         else
            write(idiag,'(29X,A)') 'No IC species modeled.'
         endif
         write(idiag,'(A,I10)') 
     &            'Number of BC groups                 : ',nbcddm
         if( nbcddm .GT. 0 ) then
            do i = 1,nbcddm
              write(idiag,'(24X,i3,'': '',A)') i,bcddmsp(i)
            enddo
         else
            write(idiag,'(29X,A)') 'No BC species modeled.'
         endif
         write(idiag,'(A,I10)') 
     &             'Number of EM groups                : ',nemddm
         if( nemddm .GT. 0 ) then
            lixemspc = .FALSE.
            do i = 1,nemddm
              write(idiag,'(24X,i3,'': '',A)') i,emddmsp(i)
              if ( emddmsp(i) .EQ. 'I2        ' .OR.
     &             emddmsp(i) .EQ. 'HOI       ' ) lixemspc = .TRUE.
            enddo
            if ( lixemis ) then
              write(idiag,'(2A)') 'WARNING: In-line Ix emissions are',
     &                            ' invoked, but DDM currently does'
              write(idiag,'(2A)') '      not account for sensitivities',
     &                            ' to the in-line emissions.'
              if ( lixemspc ) then
                write(iout,'(//,A)') 'ERROR in READNML:'
                write(iout,'(2A)') 'In-line Ix emissions are invoked,',
     &                            ' but DDM currently does not support'
                write(iout,'(2A)') 'in-line emissions. You need to',
     &                          ' either turn off in-line Ix emissions'
                write(iout,'(2A)') 'or remove I2 and HOI from the DDM',
     &                             ' emission species group list.'
                call camxerr()
              endif
            endif
         else
            write(idiag,'(29X,A)') 'No emission species modeled.'
         endif
         write(idiag,'(A,I10)') 
     &             'Number of rate const groups        : ',nrateddm
         if( nrateddm .GT. 0 ) then
            if( lddm ) then
              write(iout,'(//,a)') 'ERROR in READNML:'
              write(iout,'(2a)') 'Rate constant sensitivity cannot ',
     &                                           'be used with DDM.'
              call camxerr()
            endif
            do i = 1,nrateddm
              write(idiag,'(24X,i3,'': '',A)') i,rateddm(i)
              write(idiag,'(29X,99I4)') (iprate(n,i),n=1,iprate(0,i))
            enddo
         else
            write(idiag,'(29X,A)') 'No rate constant sens modeled.'
         endif
         write(idiag,'(A,I10)') 'Number of rate term groups : ',ntermddm
         if( ntermddm .GT. 0 ) then
            if( lddm ) then
              write(iout,'(//,a)') 'ERROR in READNML:'
              write(iout,'(2a)') 'Rate term sensitivity cannot ',
     &                                           'be used with DDM.'
              call camxerr()
            endif
            do i = 1,ntermddm
              write(idiag,'(24X,i3,'': '',A)') i,termddm(i)
              do n = 1, nreact
                if (ipterm(0,n,i).eq.1) then
                  write(idiag,'(29X,A,I4)') 'Rxn:',n
                  do l = 1, ipterm(1,n,i)
                    write(idiag,'(29X,2A,$)') '   [R]: ',
     &                    spname(ipterm(1+l,n,i))
                    if (wfterm(1+l,n,i).eq.1.0) then
                      write(idiag,*)
                    else
                      write(idiag,*) wfterm(1+l,n,i)
                    endif
                  enddo
                  do l = 1, ipterm(2+ipterm(1,n,i),n,i)
                  write(idiag,'(29X,2A)') '   [P]: ',
     &                    spname(ipterm(2+ipterm(1,n,i)+l,n,i))
                  enddo
                endif
              enddo
            enddo
         else
            write(idiag,'(29X,A)') 'No rate term sens modeled.'
         endif
         write(idiag,'(A,I10)') 
     &             'Number of HDDM sens groups         : ',nhddm
         if( nhddm .GT. 0 ) then
            do i = 1,nhddm
              write(idiag,'(24X,i3,'': '',A,'', '',A)') i,hddmsp(1,i),
     &                                                    hddmsp(2,i)
            enddo
         else
            write(idiag,'(29X,A)') 'No HDDM sensitivity modeled.'
         endif
         write(idiag,'(2A)') 'DDM turn-off flags are shown below ',
     &                       'after the nested-grid maps.'
         write(idiag,*)
         call startddm(iout,nopen)
      endif
c
c========================= DDM End =====================================
c
c
c======================== Process Analysis Begin =======================
c
c   --- echo the irmb specific flags ---
c
      if( lproca ) then
         call rdoptpa
         write(idiag,'(A,L10)') 'Integrated Process Rates   : ',lipr
         write(idiag,'(A,L10)') 'Integrated Reaction Rates  : ',lirr
         write(idiag,*)
      endif
c
      if( lipr ) then
         jj = istrln(flrtsa)
         call getunit(ipr_unit)
         filtmp = flrtsa
         filtmp(jj+1:) = '.ipr'
         nopen = nopen + 1
         open(unit=ipr_unit,file=filtmp(1:jj+4),form='UNFORMATTED',
     &                                  status= 'UNKNOWN',ERR=7000)
         write(iout,9000)
     &   'Cell Specific Process output         (unit):',ipr_unit
         write(iout,9002) filtmp(1:jj+4)
      endif
c
      if( lirr ) then
         jj = istrln(flrtsa)
         call getunit(irr_unit)
         filtmp = flrtsa
         filtmp(jj+1:) = '.irr'
         nopen = nopen + 1
         open(unit=irr_unit,file=filtmp(1:jj+4),form='UNFORMATTED',
     &                                  status= 'UNKNOWN',ERR=7000)
         write(iout,9000)
     &   'Cell Specific Rates output           (unit):',irr_unit
         write(iout,9002) filtmp(1:jj+4)
c
         if( lsfcfl ) then
            do n = 1,ngrid
               filtmp = flrtsa
               if( .NOT. lcdfout ) then
                  write(filtmp(jj+1:),'(a,i2.2)') '.cpa.grd',n
                  sfcfil(n) = filtmp
                  nopen = nopen + 1
                  call getunit(iowsfc(n))
                  open(unit=iowsfc(n),file=sfcfil(n),
     &                  form='UNFORMATTED',status= 'UNKNOWN',ERR=7000)
                  write(iout,9000)
     &                  'Gridded CPA file                     (unit):',
     &                                                     iowsfc(n)
                  write(iout,9002) filtmp(1:jj+10)
               else
                  write(filtmp(jj+1:),'(a,i2.2,a)') '.cpa.grd',n,'.nc'
                  sfcfil(n) = filtmp
                  nopen = nopen + 1
                  action = 'Opening output CPA file'
                  call ncf_createfile(filtmp,action,iowsfc(n))
                  write(iout,9000) 'Gridded CPA file:'
                  write(iout,9002) filtmp(1:jj+13)
               endif
            enddo
         endif
      endif
c
c========================= Process Analysis End ========================
c
c-----Everything worked correctly, return to calling routine
c
 8888 continue
      if( .NOT. ltrace .AND. lptsrc .AND.  npoint_files .LE. 0 ) goto 7104
      if( larsrc .AND. .NOT. is_any_grid_emiss ) goto 7105
      write(iout,*)
      write(iout,*)'Finished reading control file.'
      write(iout,'(A,I3)')'Number of input files opened: ',nopen
      write(iout,*)
      write(iout,*)
      if( ntrtim .GT. 0 ) then
        write(iout,*) '   ---------------------------------------------'
        write(iout,*) '   |                                           |'
        write(iout,*) '   |    NOTE:                                  |'
        write(iout,*) '   |    You have chosen to include timing      |'
        write(iout,*) '   |    tracers. Be aware that this will       |'
        write(iout,*) '   |    add a significant amount of memory     |'
        write(iout,*) '   |    to your application. The memory        |'
        write(iout,*) '   |    requirements may exceed available      |'
        write(iout,*) '   |    resources.                             |'
        write(iout,*) '   |                                           |'
        write(iout,*) '   ---------------------------------------------'
        write(iout,*)
      endif
      call flush(iout)
      call flush(idiag)
      return
c
 9000 format(/,A,I6)
 9002 format(2A)
c
c-----Error messages
c
 7000 continue
      write(*,'(//,A)') 'ERROR in READNML:'
      write(*,'(A)') action(:istrln(action))
      write(*,'(2A,//)') ' Could not open file: ',
     &                                 filtmp(:istrln(filtmp))
      call exit(1)
      stop
c
 7005 continue
      write(*,'(//,A)') 'ERROR in READNML:'
      write(*,'(2A,//)') ' Could not open control file: ',
     &                                   ctlfil(:istrln(ctlfil))
      call exit(1)
      stop
c
 7006 continue
      write(*,'(//,A)') 'ERROR in READNML:'
      write(*,'(A)') ' Root output file name is blank.'
      write(*,'(A)') ' A valid path/filename must be provided in '
      write(*,'(A,//)') ' the namelist variable: Root_Output_Name'
      call exit(1)
      stop
c
 7007 continue
      write(*,'(//,a)') 'ERROR in READNML:'
      write(*,'(2A,//)') ' Cannot find control file: ',
     &                     ctlfil(:istrln(ctlfil))
      call exit(1)
      stop
c
 7009 continue
      write(*,'(//,a)') 'ERROR in READNML:'
      write(*,'(2A)') ' VBS organic aerosol chemistry treatment cannot',
     &                ' be used with the chosen Probing Tool'
      write(*,'(A)') 'Set Probing_Tool = .FALSE.'
      call exit(1)
      stop
c
 7100 continue
      write(*,'(//,a)') 'ERROR in READNML:'
      write(*,'(2A)') ' Cannot find namelist group: ',
     &                                  namegrp(:istrln(namegrp))
      write(*,'(1X,A)') action(:istrln(action))
      write(*,'(3A)') ' Please make sure the CAMx control file has ',
     &        'a section beginning with &',namegrp(:istrln(namegrp))
      write(*,'(A,//)') ' and terminated with a & character.'
      call exit(1)
      stop
c
 7101 continue
      write(*,'(//,a)') 'ERROR in READNML:'
      write(*,'(2A)') ' Error reading the namelist group: ',
     &                                  namegrp(:istrln(namegrp))
      if( lddm .OR. lhddm ) then
         write(*,'(2A,/,2A,/)') ' This version allows for multiple gridded ',
     &                'emissions files.',' The namelist variable should ',
     &                'look like:  DDM_Emiss_Group_Grid(1,3,1)'
      endif
      write(*,'(2A)') ' It is probably a problem with ',
     &                         'syntax.  Please see the user''s'
      write(*,'(2A)') ' guide for examples ',
     &                'of possible problems with namelist syntax.'
      write(*,'(A,//)') ' Then review your namelist carefully.'
      call exit(1)
      stop
c
 7102 continue
      write(*,'(//,a)') 'ERROR in READNML:'
      write(*,'(2A)') ' Error reading the namelist group: ',
     &                                  namegrp(:istrln(namegrp))
      write(*,'(2A)') ' It is probably a problem with ',
     &                         'syntax.  Please see the users'
      write(*,'(2A)') ' guide for examples ',
     &                'of possible problems with namelist syntax.'
      write(*,'(A,/)') ' Then review your namelist carefully.'
      write(*,'(2A)') ' Be aware that in this version some of the',
     &                ' namelist variables have changed.'
      write(*,'(/,2A,/,2A)') ' This version allows for multiple ',
     &                'emissions files.',' The namelist variables should ',
     &                'look like: '
      write(*,'(10x,A)') 'Emiss_Grid(1,1)'
      write(*,'(10x,A)') 'Point_Sources(1)'
      write(*,'(/,2A,/,2A)') ' The Inline_Ix_Emissions name list variable has changed.',
     &                ' It is now a text variable',' instead of a logical variable.',
     &                ' These settings will NOT work.'
      write(*,'(10X,A)') '.true.'
      write(*,'(10X,A)') '.false.'
      write(*,'(1X,A)')'Acceptable options are now.'
      write(*,'(10X,A)') "'TRUE'"
      write(*,'(10X,A)') "'FALSE'"
      write(*,'(10X,A)') "'BYPASS'"
      write(*,'(/,A,/,A,//)') 
     &  ' Please refer to the Release Notes and the template included ',
     &                               ' in the source code distribution.'
      call exit(1)
      stop
c
 7103 continue
      write(*,'(//,a)') 'ERROR in READNML:'
      write(*,'(2A)') ' Error reading the namelist group: ',
     &                                  namegrp(:istrln(namegrp))
      write(*,'(2A)') ' It is probably a problem with ',
     &                         'syntax.  Please see the user''s'
      write(*,'(2A)') ' guide for examples ',
     &                'of possible problems with namelist syntax.'
      write(*,'(A,//)') ' Then review your namelist carefully.'
      write(*,'(1X,2A)') 'With this version the SA_Control namelist ',
     &                   'variables have changed'
      write(*,'(/,2A,/,2A)') ' This version allows for multiple ',
     &                'emissions files.',' The namelist variables should ',
     &                'look like: '
      write(*,'(10x,A)') 'SA_Emiss_Group_Grid(1,1,1)'
      write(*,'(10x,A)') 'SA_Points_Group(4,1)'
      write(*,'(1X,2A)') "See the User's Guide or the namelist template ",
     &           'for details.'
      call exit(1)
      stop
c
 7104 continue
      write(*,'(//,a)') 'ERROR in READNML:'
      write(*,'(2A)') ' No point source files provided.'
      write(*,'(2A,//)') ' You requested point source emissions treatment ',
     &            'but did not supply any point source files.'
      call exit(1)
c
 7105 continue
      write(*,'(//,a)') 'ERROR in READNML:'
      write(*,'(2A)') ' No gridded emissions files provided.'
      write(*,'(2A,//)') ' You requested gridded emissions treatment ',
     &            'but did not supply any gridded emissions files.'
      call exit(1)
c
c-----Return point
c
      end
