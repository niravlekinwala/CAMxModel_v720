c**** NCF_WRT_GLOBAL
c
      subroutine ncf_wrt_global(action,iounit,nspcs,spcname,is_sample)
      use ncf_iomod
      use grid
      use filunit
      implicit none
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routine writes the Global attributes to the NetCDF file
c
c      Copyright 1996 - 2022
c     Ramboll
c      Argument description:
c       Inputs:
c           action    C  name of file to open
c           iounit    I  NetCDF file ID of file
c           nspcs     I  numver of species in file
c           spcname   C name of each species
c           is_sample L .TRUE. of this is a sampling grif file
c       Outputs:
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     02/20/17   --gwilson--    Original development
c     07/23/18   --bkoo--       Added ncf_bidi_nh3_drydep
c     04/28/21   --cemery--     Added ncf_strat_o3_profile
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'netcdf.inc'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      character*(*) action
      integer       iounit
      integer       nspcs
      character*(*) spcname(nspcs)
      logical       is_sample
c
c-----------------------------------------------------------------------
c    External functions:
c-----------------------------------------------------------------------
c
      integer istrln
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      character*80                     ioapi, exec_id
      character*((MXSPEC*4+MXTRSP)*16) string
      character*16                     varname(MXSPEC*4+2*MXTRSP)
      integer                          ierr, i, string_length
c
      data varname(1:NCF_BASE_VARS) 
     &                  /'X               ','Y               ',
     &                   'layer           ','TFLAG           ',
     &                   'ETFLAG          ','longitude       ',
     &                   'latitude        ','topo            ',
     &                   'z               '/
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'FTYPE', NF_INT,
     &                                                    1, ncf_ftype)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'CDATE', NF_INT,
     &                                                    1, ncf_cdate)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'CTIME', NF_INT,
     &                                                    1, ncf_ctime)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'WDATE', NF_INT,
     &                                                    1, ncf_wdate)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'WTIME', NF_INT,
     &                                                    1, ncf_wtime)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'SDATE', NF_INT,
     &                                                    1, ncf_sdate)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'STIME', NF_INT,
     &                                                    1, ncf_stime)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'TSTEP', NF_INT,
     &                                                    1, ncf_tstep)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'NTHIK', NF_INT,
     &                                                 1, ncf_nthik)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'NCOLS', NF_INT,
     &                                                 1, ncf_ncols)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'NROWS', NF_INT,
     &                                                 1, ncf_nrows)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'NLAYS', NF_INT,
     &                                                 1, ncf_nlays)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'NVARS', NF_INT,
     &                                                 1, ncf_nvars)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'GDTYP', NF_INT,
     &                                                   1, ncf_gdtyp)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_double(iounit, NF_GLOBAL, 'P_ALP', NF_DOUBLE,
     &                                           1, DBLE(ncf_p_alp) )
      if( ierr .NE. NF_NOERR ) goto 7000
c
       ierr = nf_put_att_double(iounit, NF_GLOBAL, 'P_BET', NF_DOUBLE,
     &                                           1, DBLE(ncf_p_bet) )
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_double(iounit, NF_GLOBAL, 'P_GAM', NF_DOUBLE,
     &                                           1, DBLE(ncf_p_gam) )
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_double(iounit, NF_GLOBAL, 'XCENT', NF_DOUBLE,
     &                                           1, DBLE(ncf_xcent) )
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_double(iounit, NF_GLOBAL, 'YCENT', NF_DOUBLE,
     &                                           1, DBLE(ncf_ycent) )
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_double(iounit, NF_GLOBAL, 'XORIG', NF_DOUBLE, 
     &                                           1, DBLE(ncf_xorig) )
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_double(iounit, NF_GLOBAL, 'YORIG', NF_DOUBLE,
     &                                           1, DBLE(ncf_yorig) )
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_double(iounit, NF_GLOBAL, 'XCELL', NF_DOUBLE,
     &                                           1, DBLE(ncf_xcell) )
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_double(iounit, NF_GLOBAL, 'YCELL', NF_DOUBLE,
     &                                           1, DBLE(ncf_ycell) )
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'VGTYP', NF_INT,
     &                                                    1, ncf_vgtyp)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_real(iounit, NF_GLOBAL, 'VGTOP', NF_FLOAT, 
     &                                                    1, ncf_vgtop)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'VGLVLS', NF_INT,
     &                                            nlay(1)+1, ncf_vglvls)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_text(iounit, NF_GLOBAL, 'GDNAM',
     &                 istrln(ncf_gdnam), ncf_gdnam(:istrln(ncf_gdnam)) )
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_text(iounit, NF_GLOBAL, 'UPNAM', 
     &                 istrln(ncf_upnam), ncf_upnam(:istrln(ncf_upnam)) )
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_text(iounit, NF_GLOBAL, 'IOAPI_VERSION', 
     &                 istrln(ncf_ioapi_ver), ncf_ioapi_ver(:istrln(ncf_ioapi_ver)) )
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_text(iounit, NF_GLOBAL, 'EXEC_ID', 
     &                 istrln(ncf_execid), ncf_execid(:istrln(ncf_execid)) )
      if( ierr .NE. NF_NOERR ) goto 7000
c
      varname(1) = 'z               '
      do i=1,nspcs
        call jstlft( spcname(i) )
        varname(i+1) = spcname(i)
      enddo
      string = varname(1)
      string_length = 16
      do i=2,nspcs+1
        string_length = string_length + 16
        string = string(1:(i-1)*16) // TRIM(varname(i))
      enddo
      ierr = nf_put_att_text(iounit, NF_GLOBAL, 'VAR-LIST',
     &                              string_length, string(:string_length) )
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_text(iounit, NF_GLOBAL, 'FILEDESC', 
     &        istrln(ncf_filedesc), ncf_filedesc(:istrln(ncf_filedesc)) )
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_text(iounit, NF_GLOBAL, 'HISTORY', 
     &             istrln(ncf_history), ncf_history(:istrln(ncf_history)) )
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'IUTM', NF_INT,
     &                                                   1, ncf_iutm)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'ISTAG', NF_INT,
     &                                                   1, ncf_istag)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'CPROJ', NF_INT,
     &                                                   1, ncf_cproj)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'NSTEPS', NF_INT,
     &                                                 1, ncf_nsteps)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_text(iounit, NF_GLOBAL, 'CAMx_NAME',
     &                   istrln(ncf_name), ncf_name(:istrln(ncf_name)) )
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_text(iounit, NF_GLOBAL, 'NOTE',
     &                   istrln(ncf_note), ncf_note(:istrln(ncf_note)) )
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'ITZON', NF_INT,
     &                                                    1, ncf_itzon)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_text(iounit, NF_GLOBAL, 'UPDSC', 
     &                 istrln(ncf_upnam), ncf_upnam(:istrln(ncf_upnam)) )
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'GRID_ID', NF_INT,
     &                                                    1, ncf_grid_id)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'I_GRID_START', NF_INT,
     &                                            1, ncf_i_grid_start)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'I_GRID_END', NF_INT,
     &                                            1, ncf_i_grid_end)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'J_GRID_START', NF_INT,
     &                                            1, ncf_j_grid_start)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'J_GRID_END', NF_INT,
     &                                            1, ncf_j_grid_end)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'GRID_MESH_FACTOR',
     &                                        NF_INT, 1, ncf_grid_mesh_factor)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'FLEXI_NEST', NF_INT,
     &                                              1, ncf_flexi_nest)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_text(iounit, NF_GLOBAL, 'ADVECTION', 
     &      istrln(ncf_advection), ncf_advection(:istrln(ncf_advection)) )
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_text(iounit, NF_GLOBAL, 'CHEM_SOLVER', 
     &  istrln(ncf_chem_solver), ncf_chem_solver(:istrln(ncf_chem_solver)) )
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_text(iounit, NF_GLOBAL, 'PIG', 
     &                      istrln(ncf_pig), ncf_pig(:istrln(ncf_pig)) )
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_text(iounit, NF_GLOBAL, 'PROBING_TOOL', 
     &  istrln(ncf_probing_tool), ncf_probing_tool(:istrln(ncf_probing_tool)) )
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_text(iounit, NF_GLOBAL, 'CHEMISTRY', 
     &         istrln(ncf_chemistry), ncf_chemistry(:istrln(ncf_chemistry)) )
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'TOTAL_SPECIES', NF_INT,
     &                                              1, ncf_total_species)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'RADICAL_SPECIES', NF_INT,
     &                                            1, ncf_radical_species)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'GAS_SPECIES', NF_INT,
     &                                            1, ncf_gas_species)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'PM_SPECIES', NF_INT,
     &                                            1, ncf_pm_species)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'REACTIONS', NF_INT,
     &                                            1, ncf_reactions)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_text(iounit, NF_GLOBAL, 'DRYDEP', 
     &            istrln(ncf_drydep), ncf_drydep(:istrln(ncf_drydep)) )
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'WETDEP', NF_INT,
     &                                              1, ncf_wetdep)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'ACM2', NF_INT, 1, ncf_acm2)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'CIG_MODEL', NF_INT,
     &                                              1, ncf_cig_model)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'SURFACE_MODEL', NF_INT,
     &                                          1, ncf_surface_model)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'INLINE_IX_EMISS', NF_INT,
     &                                        1, ncf_inline_ix_emiss)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'BIDI_NH3_DRYDEP', NF_INT,
     &                                        1, ncf_bidi_nh3_drydep)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'STRAT_O3_PROFILE', NF_INT,
     &                                        1, ncf_strat_o3_profile)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'SUPER_STEPPING', NF_INT,
     &                                        1, ncf_super_stepping)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'GRIDDED_EMISS', NF_INT,
     &                                        1, ncf_gridded_emiss)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'POINT_EMISS', NF_INT,
     &                                        1, ncf_point_emiss)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'IGNORE_EMISS_DATES', NF_INT,
     &                                    1, ncf_ignore_emiss_dates)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'OUTPUT_3D', NF_INT,
     &                                            1, ncf_output_3d)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'PIG_SAMPLE_GRID', NF_INT,
     &                                       1, ncf_pig_sample_grid)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      if( is_sample ) then
          ierr = nf_put_att_int(iounit, NF_GLOBAL, 'PIG_SAMPLE_GRID_ID', 
     &                               NF_INT, 1, ncf_pig_sample_grid_id)
          if( ierr .NE. NF_NOERR ) goto 7000
          ierr = nf_put_att_int(iounit, NF_GLOBAL, 'I_SAMPLE_START', 
     &                                  NF_INT, 1, ncf_i_sample_start)
          if( ierr .NE. NF_NOERR ) goto 7000
          ierr = nf_put_att_int(iounit, NF_GLOBAL, 'I_SAMPLE_END', 
     &                                  NF_INT, 1, ncf_i_sample_end)
          if( ierr .NE. NF_NOERR ) goto 7000
          ierr = nf_put_att_int(iounit, NF_GLOBAL, 'J_SAMPLE_START', 
     &                                  NF_INT, 1, ncf_j_sample_start)
          if( ierr .NE. NF_NOERR ) goto 7000
          ierr = nf_put_att_int(iounit, NF_GLOBAL, 'J_SAMPLE_END', 
     &                                  NF_INT, 1, ncf_j_sample_end)
          if( ierr .NE. NF_NOERR ) goto 7000
      endif
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'PIG_SAMPLE_BCKGND', NF_INT,
     &                                     1, ncf_pig_sample_bckgnd)
      if( ierr .NE. NF_NOERR ) goto 7000
c
cgwilson      ierr = nf_put_att_text(iounit, NF_GLOBAL, 'Conventions', 
cgwilson     &   istrln(ncf_conventions), ncf_conventions(:istrln(ncf_conventions)) )
cgwilson      if( ierr .NE. NF_NOERR ) goto 7000
c
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,a)') 'ERROR in NCF_WRT_GLOBAL:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(A)') 'Cannot write global atttributes to file.'
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
 
