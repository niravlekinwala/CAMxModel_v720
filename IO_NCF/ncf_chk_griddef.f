c**** NCF_CHK_GRIDDEF
c
      subroutine ncf_chk_griddef(iounit,action,igrd,lchk_layers,
     &        lchk_cells,l_allow_buffer_cells,is_met_file,buffer_offset)
      use filunit
      use grid
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
c         iounit               I NCF file ID
c         action               C string that describes file being read
c         igrd                 I grid number
c         lchk_layers          L .TRUE. if number of layers is to be checked
c         lchk_cells           L .TRUE. if number of cells is to be checked
c                                       (this is for the point source exception)
c         l_allow_buffer_cells L .TRUE. if it OK to  have buffer cells
c       Outputs:
c         buffer_offset        I number of buffer to add( 0 or 1)
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     02/20/17   --gwilson--    Original development
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'flags.inc'
      include 'ncf_iodat.inc'
      include 'netcdf.inc'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      integer       iounit
      character*(*) action
      integer       igrd
      logical       lchk_layers
      logical       lchk_cells
      logical       l_allow_buffer_cells
      logical       is_met_file
      integer       buffer_offset
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
      character*10 this_var
c
      integer ierr, itzon_in, iutm_in, ncols_in, nrows_in, nlays_in
      integer cproj_in, cproj_mod, gdtyp_in, gdtyp_mod
      integer this_dimid, ncols_match, nrows_match
      real*8  dble_xcent, dble_ycent, dble_xorig, dble_yorig
      real*8  dble_xcell, dble_ycell, dble_p_alp, dble_p_bet
      real    xcent_in, ycent_in, xorig_in, yorig_in, xcell_in, ycell_in
      real    p_alp_in, p_bet_in, xorig_grid, yorig_grid
      real    xorig_mod, yorig_mod, xcell_mod, ycell_mod
      real    delx_grid, dely_grid, delx_mod, dely_mod
      real    xorig_match, yorig_match
      logical lcheck_ok
c
c-----------------------------------------------------------------------
c    Parameters:
c-----------------------------------------------------------------------
c
      real FUZZ
      parameter( FUZZ = 0.01 )
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c  --- get global variables for domain definition ----
c
      this_var = 'ITZON'
      ierr = nf_get_att_int(iounit, NF_GLOBAL, 'ITZON', itzon_in)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      this_var = 'XCENT'
      ierr = nf_get_att_double(iounit, NF_GLOBAL, 'XCENT', dble_xcent)
      if( ierr .NE. NF_NOERR ) goto 7000
      xcent_in = REAL(dble_xcent)
c
      this_var = 'YCENT'
      ierr = nf_get_att_double(iounit, NF_GLOBAL, 'YCENT', dble_ycent)
      if( ierr .NE. NF_NOERR ) goto 7000
      ycent_in = REAL(dble_ycent)
c
      this_var = 'IUTM'
      ierr = nf_get_att_int(iounit, NF_GLOBAL, 'IUTM', iutm_in)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      this_var = 'XORIG_BUF'
      ierr = nf_get_att_double(iounit, NF_GLOBAL, 'XORIG_BUF', dble_xorig)
      if( ierr .NE. NF_NOERR ) then
         this_var = 'XORIG'
         ierr = nf_get_att_double(iounit, NF_GLOBAL, 'XORIG', dble_xorig)
         if( ierr .NE. NF_NOERR ) goto 7000
      endif
      xorig_in = REAL(dble_xorig)
c
      this_var = 'YORIG_BUF'
      ierr = nf_get_att_double(iounit, NF_GLOBAL, 'YORIG_BUF', dble_yorig)
      if( ierr .NE. NF_NOERR ) then
         this_var = 'YORIG'
         ierr = nf_get_att_double(iounit, NF_GLOBAL, 'YORIG', dble_yorig)
         if( ierr .NE. NF_NOERR ) goto 7000
      endif
      yorig_in = REAL(dble_yorig)
c
      this_var = 'XCELL'
      ierr = nf_get_att_double(iounit, NF_GLOBAL, 'XCELL', dble_xcell)
      if( ierr .NE. NF_NOERR ) goto 7000
      xcell_in = REAL(dble_xcell)
c
      this_var = 'YCELL'
      ierr = nf_get_att_double(iounit, NF_GLOBAL, 'YCELL', dble_ycell)
      if( ierr .NE. NF_NOERR ) goto 7000
      ycell_in = REAL(dble_ycell)
c
      if( .NOT. is_met_file ) then
          this_var = 'COL'
          ierr = nf_inq_dimid(iounit, "COL", this_dimid  )
          if( ierr .NE. NF_NOERR ) goto 7006
          ierr = nf_inq_dimlen(iounit,this_dimid,ncols_in)
          if( ierr .NE. NF_NOERR ) goto 7006
      else
          this_var = 'NCOLS_BUF'
          ierr = nf_get_att_int(iounit, NF_GLOBAL, 'NCOLS_BUF', ncols_in)
          if( ierr .NE. NF_NOERR ) then
            this_var = 'NCOLS'
            ierr = nf_get_att_int(iounit, NF_GLOBAL, 'NCOLS', ncols_in)
            if( ierr .NE. NF_NOERR ) goto 7000
          endif
      endif
c
      if( .NOT. is_met_file ) then
         this_var = 'ROW'
         ierr = nf_inq_dimid(iounit, "ROW", this_dimid  )
         if( ierr .NE. NF_NOERR ) goto 7006
         ierr = nf_inq_dimlen(iounit,this_dimid,nrows_in)
         if( ierr .NE. NF_NOERR ) goto 7006
      else
         this_var = 'NROWS_BUF'
         ierr = nf_get_att_int(iounit, NF_GLOBAL, 'NROWS_BUF', nrows_in)
         if( ierr .NE. NF_NOERR ) then
           this_var = 'NROWS'
           ierr = nf_get_att_int(iounit, NF_GLOBAL, 'NROWS', nrows_in)
           if( ierr .NE. NF_NOERR ) goto 7000
         endif
      endif
c
      this_var = 'LAY'
      ierr = nf_inq_dimid(iounit, "LAY", this_dimid  )
      if( ierr .NE. NF_NOERR ) goto 7006
      ierr = nf_inq_dimlen(iounit,this_dimid,nlays_in)
      if( ierr .NE. NF_NOERR ) goto 7006
c
      this_var = 'CPROJ'
      ierr = nf_get_att_int(iounit, NF_GLOBAL, 'CPROJ', cproj_in)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      this_var = 'GDTYP'
      ierr = nf_get_att_int(iounit, NF_GLOBAL, 'GDTYP', gdtyp_in)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      this_var = 'P_ALP'
      ierr = nf_get_att_double(iounit, NF_GLOBAL, 'P_ALP', dble_p_alp)
      if( ierr .NE. NF_NOERR ) goto 7000
      p_alp_in = REAL(dble_p_alp)
c
      this_var = 'P_BET'
      ierr = nf_get_att_double(iounit, NF_GLOBAL, 'P_BET', dble_p_bet)
      if( ierr .NE. NF_NOERR ) goto 7000
      p_bet_in = REAL(dble_p_bet)
c
c  --- check that they all match user definition for this grid ---
c
      xorig_grid = xorg
      yorig_grid = yorg
      delx_grid = delx
      dely_grid = dely
      if( igrd .GT. 1 ) then
        delx_grid = delx/FLOAT( meshold(igrd) )
        dely_grid = delx/FLOAT( meshold(igrd) )
        xorig_grid = xorg + delx*(inst1(igrd)-1) - delx_grid
        yorig_grid = yorg + dely*(jnst1(igrd)-1) - dely_grid
      endif
      xorig_mod = xorig_grid*1000.
      yorig_mod = yorig_grid*1000.
      delx_mod = delx_grid*1000.
      dely_mod = dely_grid*1000.
      if( llatlon ) then
         cproj_mod = 0
         gdtyp_mod = 1
         xorig_mod = xorig_grid
         yorig_mod = yorig_grid
         xcell_mod = delx_grid
         ycell_mod = dely_grid
      else if( lutm ) then
         cproj_mod = 1
      else if( lambrt ) then
         cproj_mod = 2
         gdtyp_mod = 2
      else if( lrpolar ) then
         cproj_mod = 3
         gdtyp_mod = 4
      else if( lpolar ) then
         cproj_mod = 4
         gdtyp_mod = 6
      else if( lmerc ) then
         cproj_mod = 5
         gdtyp_mod = 7
      endif
c
c   --- stop if domain definition does not match ---
c
      ncols_match = ncol(igrd)
      nrows_match = nrow(igrd)
      xorig_match = xorig_mod
      yorig_match = yorig_mod
      lcheck_ok = .TRUE.
      if( ABS(xorig_match-xorig_in) .GT. FUZZ ) lcheck_ok = .FALSE.
      if( ABS(yorig_match-yorig_in) .GT. FUZZ ) lcheck_ok = .FALSE.
      if( l_allow_buffer_cells .AND. .NOT. lcheck_ok ) then
         lcheck_ok = .TRUE.
         xorig_match = xorig_mod+delx_mod
         yorig_match = yorig_mod+dely_mod
         if( ABS(xorig_match-xorig_in) .GT. FUZZ ) lcheck_ok = .FALSE.
         if( ABS(yorig_match-yorig_in) .GT. FUZZ ) lcheck_ok = .FALSE.
      endif
      if( ABS(delx_mod-xcell_in) .GT. FUZZ ) lcheck_ok = .FALSE.
      if( ABS(dely_mod-ycell_in) .GT. FUZZ ) lcheck_ok = .FALSE.
      if( lchk_cells .AND. ncols_in .NE. ncols_match ) lcheck_ok = .FALSE.
      if( lchk_cells .AND. nrows_in .NE. nrows_match ) lcheck_ok = .FALSE.
      if( l_allow_buffer_cells .AND. .NOT. lcheck_ok ) then
         ncols_match = ncols_match - 2
         nrows_match = nrows_match - 2
         lcheck_ok = .TRUE.
         if( lchk_cells .AND. ncols_in .NE. ncols_match ) lcheck_ok = .FALSE.
         if( lchk_cells .AND. nrows_in .NE. nrows_match ) lcheck_ok = .FALSE.
      endif
      if( lchk_layers .AND. nlays_in .NE. nlay(igrd) ) lcheck_ok = .FALSE.
      if( .NOT. lcheck_ok ) goto 7001
      buffer_offset = 0
      if( ncols_in+2 .EQ. ncol(igrd) 
     &          .AND. nrows_in+2 .EQ. nrow(igrd)) buffer_offset = 1
c
c   --- stop if projection type does not match ---
c
      lcheck_ok = .TRUE.
      if( cproj_in .NE. cproj_mod ) goto 7002
c
c   --- stop if projection parameters do not match ---
c
      if( itzon_in .NE. itzon ) goto 7003
      if(  lutm ) then
         if( iutm_in .NE. iuzon ) goto 7004
      else
         lcheck_ok = .TRUE.
         if( ABS(xcent_in-polelon) .GT. FUZZ ) lcheck_ok = .FALSE.
         if( ABS(ycent_in-polelat) .GT. FUZZ ) lcheck_ok = .FALSE.
         if( .NOT. lcheck_ok ) goto 7005
         lcheck_ok = .TRUE.
         if( ABS(p_alp_in-tlat1) .GT. FUZZ ) lcheck_ok = .FALSE.
         if( ABS(p_bet_in-tlat2) .GT. FUZZ 
     &                       .AND. .NOT. lpolar ) lcheck_ok = .FALSE.
c
c  --- order of the poles doesn't matter ---
c
         if( .NOT. lcheck_ok .AND. .NOT. lpolar ) then
            if( ABS(p_alp_in-tlat2) .GT. FUZZ ) goto 7005
            if( ABS(p_bet_in-tlat1) .GT. FUZZ ) goto 7005
         endif
      endif
c 
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,a)') 'ERROR in NCF_CHK_GRIDDEF:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot find global attribute: ',
     &                                      this_var(:istrln(this_var))
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7001 continue
      write(iout,'(//,a)') 'ERROR in NCF_CHK_GRIDDEF:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(A)') 'Domain definition does not match: '
      write(iout,'(/,10X,A)') 'User supplied:'
      write(iout,'(A,F15.2,A,F15.2,A)') 'Grid origin: (',
     &                                     xorig_mod,',',yorig_mod,')'
      write(iout,'(A,F15.2,A,F15.2,A)') 'Cell Width: (',
     &                                     delx_mod,',',dely_mod,')'
      if( lchk_cells) write(iout,'(A,I4,A,I4,A)') 'Number of cells: (',
     &                                 ncols_match,',',nrows_match,')'
      if( lchk_layers ) write(iout,'(A,I4)') 'Number of layers: ',nlay(igrd)
      write(iout,'(/,10X,A)') 'Values in file:'
      write(iout,'(A,F15.2,A,F15.2,A)') 'Grid origin: (',
     &                                     xorig_in,',',yorig_in,')'
      write(iout,'(A,F15.2,A,F15.2,A)') 'Cell Width: (',
     &                                     ycell_in,',',ycell_in,')'
      if( lchk_cells) write(iout,'(A,I4,A,I4,A)') 'Number of cells: (',
     &                                 ncols_in,',',nrows_in,')'
      if( lchk_layers ) write(iout,'(A,I4)') 'Number of layers: ',nlays_in
      call camxerr()
c
 7002 continue
      write(iout,'(//,a)') 'ERROR in NCF_CHK_GRIDDEF:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(A)') 'Projection type does not match.'
      write(iout,'(A,I5)') 'User supplied: ',cproj_mod
      write(iout,'(A,I5)') 'Value in file: ',cproj_in
      call camxerr()
c
 7003 continue
      write(iout,'(//,a)') 'ERROR in NCF_CHK_GRIDDEF:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(A)') 'Time zone does not match.'
      write(iout,'(A,I5)') 'User supplied: ',itzon
      write(iout,'(A,I5)') 'Value in file: ',itzon_in
      call camxerr()
c
 7004 continue
      write(iout,'(//,a)') 'ERROR in NCF_CHK_GRIDDEF:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(A)') 'UTM zone does not match.'
      write(iout,'(A,I5)') 'User supplied: ',iuzon
      write(iout,'(A,I5)') 'Value in file: ',iutm_in
      call camxerr()
c
 7005 continue
      write(iout,'(//,a)') 'ERROR in NCF_CHK_GRIDDEF:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(A)') 'Projection parameters do not match.'
      write(iout,'(/,10X,A)') 'User supplied:'
      write(iout,'(A,F10.4)') 'Longitude of projection pole:',polelon
      write(iout,'(A,F10.4)') 'Latitude of projection pole:',polelat
      write(iout,'(A,F10.4)') '1st true latitude:',tlat1
      if( .NOT. lpolar ) write(iout,'(A,F10.4)') '2nd true latitude:',tlat2
      write(iout,'(/,10X,A)') 'Values in file:'
      write(iout,'(A,F10.4)') 'Longitude of projection pole:',xcent_in
      write(iout,'(A,F10.4)') 'Latitude of projection pole:',ycent_in
      write(iout,'(A,F10.4)') '1st true latitude:',p_alp_in
      if( .NOT. lpolar ) write(iout,'(A,F10.4)') '2nd true latitude:',p_bet_in
      call camxerr()
c
 7006 continue
      write(iout,'(//,a)') 'ERROR in NCF_CHK_GRIDDEF:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(A)') 'Cannot find the dimension id for dimension:',
     &                            this_var(:istrln(this_var))
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
