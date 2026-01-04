c**** NCF_SET_GLOBAL
c
      subroutine ncf_set_global(inname,igrd,begin_date,begin_time,
     &                       ending_date,ending_time,nlays,nspcs)
      use camxcom
      use grid
      use pigsty
      use filunit
      implicit none
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routine sets the global file attributes for the NetCDF file
c
c      Copyright 1996 - 2022
c     Ramboll
c      Argument description:
c       Inputs:
c         inname      C filetype of this kind of file
c         igrd        I grid index
c         begin_date  I model begin date (YYJJJ)
c         begin_time  R model begin time
c         ending_date I model end date (YYJJJ)
c         ending_time R model end time
c         nlays       I number of layers in this file
c         nspcs       I number of species in this file
c       Outputs:
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     02/20/17   --gwilson--    Original development
c     07/23/18   --bkoo--       Added ncf_bidi_nh3_drydep
c     07/23/19   Staggered wind flag is now grid-specific
c     04/28/21   --cemery--     Added ncf_strat_o3_profile
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'flags.inc'
      include 'ncf_iodat.inc'
      include 'namelist.inc'
      include 'chmdat.inc'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      character*10 inname
      integer      igrd
      integer      begin_date
      real         begin_time
      integer      ending_date
      real         ending_time
      integer      nlays
      integer      nspcs
c
c-----------------------------------------------------------------------
c    External functions:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      real xorg_grid, yorg_grid, delx_grid, dely_grid
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c
c --- get current date/time ---
c
      call getime(ncf_cdate,ncf_ctime)
      ncf_wdate = ncf_cdate
      ncf_wtime = ncf_ctime
c
c --- domain definition attributes ---
c
      xorg_grid = xorg
      yorg_grid = yorg
      delx_grid = delx
      dely_grid = dely
      if( igrd .GT. 1 ) then
        delx_grid   = delx/FLOAT( meshold(igrd) )
        dely_grid   = delx/FLOAT( meshold(igrd) )
        xorg_grid = xorg + delx*(inst1(igrd)-1)
        yorg_grid = yorg + dely*(jnst1(igrd)-1)
      endif
      ncf_iutm = iuzon
      ncf_istag = 0
      if (lstagw(igrd)) ncf_istag = 1
      ncf_xorig = xorg_grid*1000.
      ncf_yorig = yorg_grid*1000.
      ncf_xcell = delx_grid*1000.
      ncf_ycell = dely_grid*1000.
      if( llatlon ) then
         ncf_cproj = 0
         ncf_gdtyp = 1
         ncf_xorig = xorg_grid
         ncf_yorig = yorg_grid
         ncf_xcell = delx_grid
         ncf_ycell = dely_grid
      else if( lutm ) then
         ncf_cproj = 1
         ncf_gdtyp = 5
      else if( lambrt ) then
         ncf_cproj = 2
         ncf_gdtyp = 2
      else if( lrpolar ) then
         ncf_cproj = 3
         ncf_gdtyp = 4
      else if( lpolar ) then
         ncf_cproj = 4
         ncf_gdtyp = 6
      else if( lmerc ) then
         ncf_cproj = 5
         ncf_gdtyp = 7
      endif
      ncf_xcent = polelon
      ncf_ycent = polelat
      if (tlat1.gt.0.) then
        ncf_p_alp = min(tlat1,tlat2)
        ncf_p_bet = max(tlat1,tlat2)
      else
        ncf_p_alp = max(tlat1,tlat2)
        ncf_p_bet = min(tlat1,tlat2)
      endif
      ncf_p_gam = polelon
      ncf_nlays = nlays
      ncf_nrows = nrow(igrd)
      ncf_ncols = ncol(igrd)
      if( igrd .GT. 1 )  then
         ncf_nrows = nrow(igrd) - 2
         ncf_ncols = ncol(igrd) - 2
      endif
      ncf_nthik = 1
      ncf_nvars = nspcs + 1
c
c --- file description attributes ---
c
      ncf_name = inname
      ncf_note = runmsg
      ncf_itzon = itzon
      ncf_ftype = 1
      ncf_vgtyp = 6
      ncf_vgtop = 10000.
      ncf_vglvls = 0.
      ncf_gdnam = "CAMx v7.20"
      ncf_gdnam = "CAMx v7.20"
      ncf_ioapi_ver = "IOAPI-CAMx"
      ncf_execid = "CAMx"
      ncf_filedesc = inname
      ncf_sdate = Start_Date_Hour(1)*1000 + MOD(begin_date,1000)
      ncf_stime = begin_time*100
      if( dtout .LT. 60. ) then
        ncf_tstep = INT(dtout)*100
      else
        ncf_tstep = INT(dtout/60.)*10000
      endif
      ncf_grid_id = igrd
      ncf_i_grid_start = 1
      ncf_i_grid_end = ncol(1)
      ncf_j_grid_start = 1
      ncf_j_grid_end = nrow(1)
      if( igrd .GT. 1 ) then
         ncf_i_grid_start = inst1(igrd)
         ncf_i_grid_end = inst2(igrd)
         ncf_j_grid_start = jnst1(igrd)
         ncf_j_grid_end = jnst2(igrd)
      endif
      ncf_grid_mesh_factor = meshold(igrd)
c
c --- simulation description attributes ---
c
      ncf_flexi_nest = 0
      if( lflexi ) ncf_flexi_nest = 1
      ncf_advection = Advection_Solver
      ncf_chem_solver = Chemistry_Solver
      ncf_pig = PiG_Submodel
      ncf_probing_tool = Probing_Tool
      ncf_total_species = nspec
      ncf_radical_species = nrad
      ncf_gas_species = ngas
      ncf_pm_species = naero
      ncf_reactions = nreact
      ncf_drydep = Drydep_Model
      ncf_wetdep = 0
      if( lwet ) ncf_wetdep = 1
      ncf_acm2 = 0
      if( lacm2 ) ncf_acm2 = 1
      ncf_cig_model = 0
      if( Subgrid_Convection ) ncf_cig_model = 1
      ncf_surface_model = 0
      if( lsrfmod ) ncf_surface_model = 1
      ncf_inline_ix_emiss = 0
      if( lixemis ) ncf_inline_ix_emiss = 1
      ncf_bidi_nh3_drydep = 0
      if( lbidinh3 ) ncf_bidi_nh3_drydep = 1
      ncf_strat_o3_profile = 0
      if( lstrato3 ) ncf_strat_o3_profile = 1
      ncf_super_stepping = 0
      if( lsuper ) ncf_super_stepping = 1
      ncf_gridded_emiss = 0
      if( larsrc ) ncf_gridded_emiss = 1
      ncf_point_emiss = 0
      if( lptsrc ) ncf_point_emiss = 1
      ncf_ignore_emiss_dates = 0
      if( le1day ) ncf_ignore_emiss_dates = 1
      ncf_output_3d = 0
      if( l3davg(igrd) ) ncf_output_3d = 1
      ncf_pig_sample_grid = 0
      if( lsample ) ncf_pig_sample_grid = 1
      ncf_pig_sample_bckgnd = 0
      if( lbckgrd ) ncf_pig_sample_bckgnd = 1
      ncf_conventions = "CF-1.6"
      ncf_history = 'Generated by '//PRMVERSION
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
 
