      subroutine startup(version,inptim,inpdate,emstim,emsdate,
     &                   ozntim,ozndate,bndtim,bnddate,toptim,topdate,
     &                   wrttim,wrtdate,endtim,enddate,numprocs)
      use filunit
      use grid
      use chmstry
      use o3colmap
      use bndary
      use camxfld
      use camxcom
      use ptemiss
      use pigsty
      use procan
      use tracer
      use rtracchm
     
      implicit none
c
c----CAMx v7.20 220430
c
c     STARTUP is the main initialization and setup routine for CAMx.
c     It performs the following tasks:
c        - initializes certain vector/array and scalar variables 
c        - reads and checks the user input file
c        - sets model simulation clock
c        - opens all I/O files
c        - reads all time-invariant files and look-up tables
c        - reads/writes headers from/to UAM-formatted I/O files
c        - calculates grid parameters
c        - initializes the PiG submodel
c                          
c      Copyright 1996 - 2022
c     Ramboll 
c          
c     Modifications:
c        1/29/99   Added diagnostic printout of layer nesting, and error
c                  messages if layer nesting calculation is invalid
c        2/12/99   Removed assignment of negative longitude for xorg
c        4/26/99   Added Piecewise Parabolic Method for horizontal advection
c        5/24/99   Fixed bug in filling idfin array
c        10/20/00  Added CAMx version as first record on control file
c        11/06/01  Added Y2K flag; conversion of simulation date to Julian 
c                   now done immediately after reading from control file
c        1/21/02   Added RTRAC tech type to OSAT
c        1/25/02   Revised input I/O frequencies and max time step to minutes,
c                  and improved checks on values
c        7/5/02    Added code to handle new IRON-PiG option
c        1/10/03   Added prep of deposition output files 
c        03/21/03  Removed the OSAT technology type OPPAT
c        05/01/03  Added snow cover and other options to AHO file
c        11/10/03  Added sampling grid setup for RTRAC+IRONPIG
c        10/06/04  Restructured for namelist input
c        10/12/04  Added read of water vapor and vertical diffusivity as
c                  time-interpolated variables
c        10/14/04  Special checks added for mechanism 10
c        7/29/05   Added sampling grid setup for IRONPIG for average concs
c        12/15/08  Added code to handle averaging of radicals
c        03/15/09  Added code for deposition output for tracers
c        07/16/07 -bkoo-     Revised for HDDM
c        07/16/08 -bkoo-     Added code to print/check DDM turn-off flags
c        11/04/09  Removed setting nest grid boundary winds from parent
c        11/08/09  Removed ptconc from argument to rdinstsa
c        11/08/09  Fixed bug in checking parameter against nrad
c        01/04/11  Revised for new met input format
c        04/02/12  Removed drought stress and snow flag; AHO
c                  file is now just ozone column; replaced haze
c                  dimension with terrain height in photo file
c        11/05/12  Removed vertical nesting
c        04/30/13  Added surface model
c        08/23/13 -bkoo-     Enabled DDM with PM
c        4/07/14   Added top con file
c        02/21/20 -rlb- Added mechanism 1 (CB6r5)
c        04/01/21  Added routine to write TOPO to SA files
c        08/26/21  RTRAC surface file doubles as (default) dep output
c                  or (flagged) surface model output
c
c     Input arguments:
c        version             model version character string
c
c     Output arguments:
c        inptim              next time to read environmental fields (HHMM)
c        inpdate             next date to read environmental fields (YYJJJ)
c        emstim              next time to read emissions (HHMM)
c        emsdate             next date to read emissions (YYJJJ)
c        ozntim              next update time for ozone map (HHMM)  
c        ozndate             next update date for ozone map (YYJJJ) 
c        bndtim              next update time for boundary conditions (HHMM)  
c        bnddate             next update date for boundary conditions (YYJJJ) 
c        toptim              next update time for top concentrations (HHMM)  
c        topdate             next update date for top concentrations (YYJJJ) 
c        wrttim              next time to output concentrations (HHMM)
c        wrtdate             next date to output concentrations (YYJJJ)
c        endtim              model end time (HHMM)
c        enddate             model end date (YYJJJ)
c        numprocs            number of processsors for MPI
c
c     Routines Called:
c        READNML,  READCHM,  READPHT,  O3COLPREP, BNDPREP,  GRDPREP,
c        NSTPREP,  IASSGN2D, INTERP2D, METINIT,  VNMSHCAL, SRFPREP, 
c        LUASSGN,  PNTPREP,  AREAPREP, CNCPREP,  PIGPREP,
c        FINWIND,  RASSGN3D, DEPPREP
c
c     Called by:
c        CAMx
c
      include 'camx.prm'
      include 'flags.inc'
      include 'deposit.inc'
      include 'rtracsrf.inc'
      include 'vbs.inc'
c
      logical lddmtmp(MXGRID)
c
      integer istrln
c
      character*200 action
      character*20  version
      character*10  cparm
      character*20  spec_units(4*MXSPEC+2+2*MXTRSP)
      character*20  spec_long_name(4*MXSPEC+2+2*MXTRSP)
      character*60  spec_desc(4*MXSPEC+2+2*MXTRSP)
      character*60  spec_coords(4*MXSPEC+2+2*MXTRSP)
      integer       inpdate, emsdate, wrtdate, enddate, ozndate
      integer       bnddate, topdate, numprocs
      integer       numspcs, nlayer, lav, i, j, l, n, ip, ic 
      integer       igrd, ig, kg, kg1, kp
      integer       icheck, iparm, ism, ifile, numpts_lst
      real          inptim, emstim, wrttim, endtim, ozntim
      real          bndtim, toptim, dxmod, dymod, whr, wmn
      real          orgx,orgy
      logical       lfirst_call, li2
      logical       l_exist_ivoc, l_exist_poa_xx, l_exist_pap

c
c-----Entry point
c
      lfirst_call = .TRUE.
      icur_unit = 6
      call initnml( )
      call readnml(version,enddate,endtim,numprocs)
c
c-----Initialize simulation clock
c
      time = begtim
      date = begdate
c
c-----Call routine to set up the pointers for grid vectors
c
      call iniptr(ncol,nrow)
c
c-----Allocate the arrays for Ozone Column
c
      call alloc_o3col(ngrid,ncol,nrow)
c
c----make sure arrays are large enough
c
      icheck = nspec+1
      iparm = MXSPEC
      cparm = 'MXSPEC'
      if( icheck .GT. iparm ) goto 7000
c  
c-----Call the pointer routine again to update for deposition arrays
c
      call iniptr(ncol,nrow)
c
c-----Allocate the arrays for time steps 
c
      call alloc_camxcom(ngrid)
c
c-----Allocate the arrays for gridded fields
c
      call alloc_camxfld(ngrid,ncol,nrow,nlay,nspec,navspc,
     &                                    ndepspc,nsmspc,l3davg)
c
c-----Read photolysis rates lookup table
c
      if (lchem .and. idmech.ne.10) call readpht()
c
c-----Read ozone column file header
c
      if (io3col .NE. 0 .AND. idmech .NE. 10) call o3colprep()
c
c-----Call routine to allocate arrays for boundary data ---
c
      call alloc_bndary(ngrid,ncol,nrow,nlay,nspec)
c
c-----Call routine to write the headers of all sampling grids
c
      if (lsample) then
        call smpprep(.false.,endtim,enddate)
      endif
c
c-----Call routine to allocate the arrays by row
c
      call alloc_grid_row(nrow)
c
c-----Calculate grid parameters for coarse grid
c
      call grdprep(ncol(1),nrow(1),cellon(1),cellat(1),mapscl(1),
     &                        xorg,yorg,delx,dely,deltax(1,1),deltay(1) )
c
c-----Call routine to allocate arrays by layer --
c
      call alloc_grid_lay()
c
c-----Call routine to allocate gridded 2D arrays
c
      call alloc_grid_2d(ncol, nrow)
c
c-----Calculate nested grid mapping parameters
c
      if (ngrid.gt.1) then
        call nstprep()
      else
        do j=1,nrow(1)
          do i=1,ncol(1)
            n = (j-1)*ncol(1) + i
            idfin(n) = 0
          enddo
        enddo
        mapgrd(1) = 1
        nchdrn(1) = 0
        meshold(1) = 1
        nmesh(1) = 1
      endif
c
c======================== Source Apportion Begin =======================
c
c   --- make sure that the parameters are large enough ---
c
      icheck = 0
      do i=1,ngrid
        icheck = MAX(icheck,nrow(i),ncol(i))
      enddo
      iparm = MXCELLS
      cparm = 'MXCELLS'
      if( icheck .GT. iparm ) goto 7000
c
      if( ltrace ) then
c
c   --- call routine to initialize the source apportionment
c       data structures ---
c
          call initsa(version,ncol,nrow,
     &                             begdate,begtim,enddate,endtim)
c
c   --- make sure that the parameters are large enough ---
c
          icheck = ntotsp
          iparm = MXTRSP
          cparm = 'MXTRSP'
          if( icheck .GT. iparm ) goto 7000
          icheck = nlay(1)
          iparm = MXLAYER
          cparm = 'MXLAYER'
          if( icheck .GT. iparm ) goto 7000
c
c   --- call the pointer routine again to update for tracer arrays ---
c
          call iniptr(ncol,nrow)
          if( lrstrt ) then
             call rdinstsa(begdate,time,ncol(1),nrow(1),nlay(1),ntotsp)
          endif
c
c  ---- call routine to calculate the average reactivity
c       of boundary conditions ---
c
          if( tectyp .NE. RTRAC .AND. tectyp .NE. RTCMC ) then
              call clcbwt(begdate,begtim,enddate,endtim,
     &                                        ncol(1),nrow(1),nlay(1))
          endif
c
c   --- call routine to write the header of the coarse grid
c       average surface tracer concentrations file ---
c
          if( lsfcfl ) then
             if( tectyp .EQ. RTRAC .OR. tectyp .EQ. RTCMC ) then
                call ncf_set_specatt_rtrac(spec_units,spec_long_name,spec_desc,
     &                                                           spec_coords)
             else
                call ncf_set_specatt_sa(spec_units,spec_long_name,spec_desc,
     &                                                         spec_coords)
             endif
             call hdrwsa(iowsfc,sfcfil,'AVERAGE   ',ntotsp,1,begdate,begtim,
     &                              enddate,endtim,spec_units,spec_long_name,
     &                                                  spec_desc,spec_coords)
          endif
c
c   --- call routine to write the header of the RTRAC
c       surface model mass file ---
c
          if( tectyp .EQ. RTRAC .OR. tectyp .EQ. RTCMC ) then
             call ncf_set_specatt_rtrac_dep(spec_units,spec_long_name,spec_desc,
     &                                                           spec_coords)
             call hdrwsrf(lsrfmodrt,iowrtsrf,rtsrfout,ntotsp,1,
     &               begdate,begtim,enddate,endtim,spec_units,spec_long_name,
     &                                                  spec_desc,spec_coords)
          endif
c
c   --- call routine to write the header of the coarse grid
c       average surface tracer depositions file ---
c
          call ncf_set_specatt_sadep(spec_units,spec_long_name,spec_desc,
     &                                                         spec_coords)
          if( lptdepout ) call hdrdepsa(iowptdep,ptdepfil,notimespc,
     &                          begdate,begtim,enddate,endtim,spec_units,
     &                                 spec_long_name,spec_desc,spec_coords)
c
c   --- if doing ozone or nitrate species in SA, initialize the
c       number of IRR reactions ---
c
          if( lozone .OR. lnitrate ) nirrrxn =  nreact
c
c   --- call routine to write the headers of all sampling grids
c
          if ((tectyp .EQ. RTRAC .OR. tectyp .EQ. RTCMC) .AND. lsample
     &                                             .AND. lsmptrc) then
             call smpprep(.TRUE.,endtim,enddate)
          endif
      endif
c
c========================= Source Apportion End ========================
c
c
c======================== DDM Begin ====================================
c
      if( lddm .OR. lhddm ) then
c
c   --- call routine to initialize the source apportionment
c       data structures ---
c
         call initsa(version,ncol,nrow,begdate,begtim,enddate,endtim)
c
c   --- if this is a restart, call routine to read the instantaneous files ---
c
         icheck = nddmsp*nspec
         iparm = MXTRSP
         cparm = 'MXTRSP'
         if( icheck .GT. iparm ) goto 7000
         icheck = ntotsp
         iparm = MXTRSP
         cparm = 'MXTRSP'
         if( icheck .GT. iparm ) goto 7000
c
c   --- call the pointer routine again to update for tracer arrays ---
c
         call iniptr(ncol,nrow)
         if( lrstrt ) then
            if( lddmcalc(1) )
     &       call rdinstsa(begdate,time,ncol(1),nrow(1),nlay(1),ntotsp)
         endif
c
c   --- call routine to write the header of the coarse grid
c       average surface tracer concentrations file ---
c
         if( lsfcfl ) then
             if( lcdfout ) then
                call ncf_set_specatt_ddm(spec_units,spec_long_name,spec_desc,
     &                                                         spec_coords)
             endif
             call hdrwsa(iowsfc,sfcfil,'AVERAGE   ',ntotsp,1,begdate,begtim,
     &                              enddate,endtim,spec_units,spec_long_name,
     &                                                  spec_desc,spec_coords)
         endif
      endif
c
c======================== DDM End ====================================
c
c-----Read BC/TC file header and set irregular boundary cells
c
      if( is_netcdf_ibc ) then
          call ncf_bndprep(begtim,begdate,endtim,enddate)
      else
          call bndprep(begtim,begdate,endtim,enddate)
      endif
      if (ltopcon) then
         if( is_netcdf_itc ) then
             call ncf_topprep(begtim,begdate,endtim,enddate)
         else
             call topprep(begtim,begdate,endtim,enddate)
         endif
      endif
c
c------Calculate grid parameters for fine grids
c
      do ip = 1,ngrid
        do ic = 1,nchdrn(ip)
          igrd = idchdrn(ic,ip)
          call interp2d(ncol(ip),nrow(ip),1,i1(igrd),j1(igrd),
     &                  nmesh(igrd),ncol(igrd),nrow(igrd),
     &                  cellat(iptr2d(ip)),cellat(iptr2d(igrd)) )
          call interp2d(ncol(ip),nrow(ip),1,i1(igrd),j1(igrd),
     &                  nmesh(igrd),ncol(igrd),nrow(igrd),
     &                  cellon(iptr2d(ip)),cellon(iptr2d(igrd)) )
          call interp2d(ncol(ip),nrow(ip),1,i1(igrd),j1(igrd),
     &                  nmesh(igrd),ncol(igrd),nrow(igrd),
     &                  mapscl(iptr2d(ip)),mapscl(iptr2d(igrd)) )
        enddo
      enddo
c
c-----Read met fields that are to be time-interpolated, for all grids,
c     to current time
c
      do igrd = 1,ngrid
        icheck = nlay(igrd)
        iparm = MXLAYER
        cparm = 'MXLAYER'
        if( icheck .GT. iparm ) goto 7000

        orgx = xorg
        orgy = yorg
        dxmod = delx
        dymod = dely
        if (igrd.gt.1) then
          dxmod = delx/meshold(igrd)
          dymod = dely/meshold(igrd)
          orgx = xorg + (inst1(igrd) - 1)*delx - dxmod
          orgy = yorg + (jnst1(igrd) - 1)*dely - dymod
        endif
        if( is_netcdf_i3dmet(igrd) ) then
          call ncf_metprep(igrd,begtim,begdate,endtim,enddate)
          call ncf_metinit_3d(igrd,ncol(igrd),nrow(igrd),nlay(igrd),
     &               height(iptr3d(igrd)),press(iptr3d(igrd)),
     &               depth(iptr3d(igrd)),windu(iptr3d(igrd)),
     &               windv(iptr3d(igrd)),tempk(iptr3d(igrd)),
     &               water(iptr3d(igrd)))
          call ncf_metinit_2d(igrd,ncol(igrd),nrow(igrd),
     &               tsurf(iptr2d(igrd)),snow(iptr2d(igrd)),
     &               snowage(iptr2d(igrd)),snowalb(iptr2d(igrd)))
          call ncf_metinit_kv(igrd,ncol(igrd),nrow(igrd),nlay(igrd),
     &               rkv(iptr3d(igrd)))
        else
           call metinit(igrd,ncol(igrd),nrow(igrd),nlay(igrd),
     &               endtim,enddate,orgx,orgy,dxmod,dymod,
     &               height(iptr3d(igrd)),press(iptr3d(igrd)),
     &               depth(iptr3d(igrd)),windu(iptr3d(igrd)),
     &               windv(iptr3d(igrd)),tempk(iptr3d(igrd)),
     &               tsurf(iptr2d(igrd)),water(iptr3d(igrd)),
     &               rkv(iptr3d(igrd)),snow(iptr2d(igrd)),
     &               snowage(iptr2d(igrd)),snowalb(iptr2d(igrd)))
        endif
      enddo
c
c-----Loop over all grids to identify missing input fields
c
      do ip = 1,ngrid
        do ic = 1,nchdrn(ip)
          ig = idchdrn(ic,ip)

          icheck = nlay(ip)
          iparm = MXLAYER
          cparm = 'MXLAYER'
          if( icheck .GT. iparm ) goto 7000

          icheck = nlay(ig)
          iparm = MXLAYER
          cparm = 'MXLAYER'
          if( icheck .GT. iparm ) goto 7000
c
c-----Compare vertical layer structures to ensure consistency
c
          if( i3dmet(ig) .GT. 0 ) then
            call vnmshcal(ig,ncol(ip),nrow(ip),nlay(ip),i1(ig),j1(ig),
     &                    ncol(ig),nrow(ig),nlay(ig),height(iptr3d(ip)),
     &                    height(iptr3d(ig)))
          else
c
c-----Assign vertical grid structures that were not read to each grid
c
            write(iout,'(a40,f7.0,i8.5,a,i3)')
     &                 'Assigning heights from parent grid',
     &                             time, date,' grid',ig
            call rassgn3d(ncol(ip),nrow(ip),nlay(ip),
     &           i1(ig),j1(ig),nmesh(ig),ncol(ig),nrow(ig),
     &                        height(iptr3d(ip)),height(iptr3d(ig)) )
            call rassgn3d(ncol(ip),nrow(ip),nlay(ip),
     &           i1(ig),j1(ig),nmesh(ig),ncol(ig),nrow(ig),
     &                          depth(iptr3d(ip)),depth(iptr3d(ig)) )
c
c-----Interpolate wind fields that were not read to each grid
c
            write(iout,'(a40,f7.0,i8.5,a,i3)')
     &                 'Interpolating winds from parent grid',
     &                             time, date,' grid',ig
            call finwind(ncol(ip),nrow(ip),nlay(ip),i1(ig),
     &                 j1(ig),nmesh(ig),ncol(ig),nrow(ig),nlay(ig),
     &                 windu(iptr3d(ip)),windv(iptr3d(ip)),
     &                 windu(iptr3d(ig)),windv(iptr3d(ig)))
c 
c-----Interpolate pressure fields that were not read to each grid
c 
            write(iout,'(a40,f7.0,i8.5,a,i3)')
     &                 'Interpolating pressure from parent grid',
     &                             time, date,' grid',ig
            call interp2d(ncol(ip),nrow(ip),nlay(ip),
     &                    i1(ig),j1(ig),nmesh(ig),ncol(ig),nrow(ig),
     &                            press(iptr3d(ip)),press(iptr3d(ig)) )
c 
c-----Interpolate temperature fields that were not read to each grid
c 
            write(iout,'(a40,f7.0,i8.5,a,i3)')
     &                 'Interpolating temps from parent grid',
     &                             time, date,' grid',ig
            call interp2d(ncol(ip),nrow(ip),nlay(ip),i1(ig),j1(ig),
     &                    nmesh(ig),ncol(ig),nrow(ig),
     &                    tempk(iptr3d(ip)),
     &                    tempk(iptr3d(ig)) )
c
c-----Interpolate water vapor fields that were not read to each grid
c
            write(iout,'(a40,f7.0,i8.5,a,i3)')
     &                 'Interpolating humidity from parent grid',
     &                             time, date,' grid',ig
            call interp2d(ncol(ip),nrow(ip),nlay(ip),i1(ig),j1(ig),
     &                    nmesh(ig),ncol(ig),nrow(ig),
     &                    water(iptr3d(ip)),
     &                    water(iptr3d(ig)) )
          endif
c
c-----Interpolate diffusivity fields that were not read to each grid
c
          if (ikv(ig).eq. 0) then
            write(iout,'(a40,f7.0,i8.5,a,i3)')
     &                 'Interpolating VDiff from parent grid',
     &                             time, date,' grid',ig
            call interp2d(ncol(ip),nrow(ip),nlay(ip),i1(ig),j1(ig),
     &                    nmesh(ig),ncol(ig),nrow(ig),
     &                    rkv(iptr3d(ip)),rkv(iptr3d(ig)) )
          endif
c
c-----Interpolate 2D surface temperature fields that were not read to each grid
c
          if (i2dmet(ig) .eq. 0) then
            write(iout,'(a40,f7.0,i8.5,a,i3)')
     &                 'Interpolating sfc temp from parent grid',
     &                             time, date,' grid',ig
            call interp2d(ncol(ip),nrow(ip),1,i1(ig),j1(ig), 
     &                    nmesh(ig),ncol(ig),nrow(ig),
     &                    tsurf(iptr2d(ip)),tsurf(iptr2d(ig)) ) 
          endif
        enddo
      enddo
c
c-----Read surface files and initialize arrays
c
      do igrd = 1,ngrid
        orgx = xorg
        orgy = yorg
        dxmod = delx
        dymod = dely
        if (igrd.gt.1) then
          dxmod = delx/meshold(igrd)
          dymod = dely/meshold(igrd)
          orgx = xorg + (inst1(igrd) - 1)*delx - dxmod
          orgy = yorg + (jnst1(igrd) - 1)*dely - dymod
        endif
        if (isurf(igrd).ne.0) lrdocn(igrd) = .true.
        if( is_netcdf_isurf(igrd) ) then
          call ncf_luseprep(igrd,begtim,begdate,endtim,enddate)
          call ncf_readinp_lu(igrd,ncol(igrd),nrow(igrd),nlu,
     &               fsurf(iptrlu(igrd)),topo(iptr2d(igrd)),
     &               lai(iptr2d(igrd)),lrdlai(igrd),
     &               icdocn(iptr2d(igrd)))
        else
           call srfprep(igrd,ncol(igrd),nrow(igrd),orgx,orgy,dxmod,dymod,
     &               fsurf(iptrlu(igrd)),topo(iptr2d(igrd)),
     &               lai(iptr2d(igrd)),lrdlai(igrd),
     &               icdocn(iptr2d(igrd)))
        endif
      enddo
c
c-----Assign fine grid landuse fractions and ocean mask
c
      do ip = 1,ngrid
        do ic = 1,nchdrn(ip)
          igrd = idchdrn(ic,ip)
          if (isurf(igrd).eq.0) then
            call luassgn(ncol(ip),nrow(ip),nlu,i1(igrd),j1(igrd),
     &                   nmesh(igrd),ncol(igrd),nrow(igrd),
     &                   fsurf(iptrlu(ip)),fsurf(iptrlu(igrd)) )
            call interp2d(ncol(ip),nrow(ip),1,i1(igrd),j1(igrd), 
     &                    nmesh(igrd),ncol(igrd),nrow(igrd),
     &                    topo(iptr2d(ip)),topo(iptr2d(igrd)) ) 
            call iassgn2d(ncol(ip),nrow(ip),i1(igrd),j1(igrd),
     &                    nmesh(igrd),ncol(igrd),nrow(igrd),
     &                    icdocn(iptr2d(ip)),icdocn(iptr2d(igrd)))
            lrdocn(igrd) = .true.
            if (lrdlai(ip)) then
              call luassgn(ncol(ip),nrow(ip),1,i1(igrd),j1(igrd),
     &                     nmesh(igrd),ncol(igrd),nrow(igrd),
     &                     lai(iptr2d(ip)),lai(iptr2d(igrd)) )
              lrdlai(igrd) = .true.
            endif
          elseif (.NOT.lrdlai(igrd) .AND. lrdlai(ip)) then
            call luassgn(ncol(ip),nrow(ip),1,i1(igrd),j1(igrd),
     &                   nmesh(igrd),ncol(igrd),nrow(igrd),
     &                   lai(iptr2d(ip)),lai(iptr2d(igrd)) )
            lrdlai(igrd) = .true.
          endif
        enddo
      enddo
c
c-----Read emission file headers
c
      if (lptsrc) then
         call pntprep(begtim,begdate,endtim,enddate)
         call alloc_tracer_ptsrce(nptsrc)
         call alloc_tracer_sapnts(nptsrc,ntotsp)
      else
         call alloc_grid_ptsrc(MAX(1,nptsrc))
         call alloc_ptemiss_null(nspec)
         call alloc_tracer_ptsrce(1)
         call alloc_tracer_sapnts(1,1)
      endif
      if (larsrc) then
        do igrd = 1,ngrid
          dxmod = delx
          dymod = dely
          if (igrd.gt.1) then
            dxmod = delx/meshold(igrd)
            dymod = dely/meshold(igrd)
          endif
          do ifile=1,nemiss_files(igrd)
             if( is_netcdf_iarem(igrd,ifile) ) then
                write(action,'(2A,I3,A,I3)') 'Reading the NetCDF ',
     &                              'gridded emissions file. Grid: ',
     &                                           igrd,' File: ',ifile
                call ncf_areaprep(action,iarem(igrd,ifile),igrd,
     &                              begtim,begdate,endtim,enddate,
     &                                  buffer_offset_iarem(igrd,ifile))
             else
                call areaprep(igrd,ifile,begtim,begdate,endtim,enddate,
     &                  iarem(igrd,ifile),iout,idiag,dxmod,dymod)
             endif
             li2 = .FALSE.
             do i=1,nemspc
               if( emspcname(i)(1:3) .EQ. 'I2 ' .OR. 
     &             emspcname(i)(1:4) .EQ. 'HOI ' ) li2 = .TRUE.
             enddo
          enddo
          if (.not.li2 .AND. .not.lixemis .AND. .not.lixbypass .AND. 
     &        (idmech.eq.1 .OR. idmech.eq.3 .OR. idmech.eq.4 .OR. idmech.eq.7)) then
            write(iout,'(//,A)') 'ERROR in STARTUP:'
            write(iout,'(A)') 'You are running with iodine chemistry,'
            write(iout,'(A)') 'but in-line Ix emissions are not invoked,'
            write(iout,'(2A)')'and I2 and/or HOI are not found',
     &                    ' in an input gridded emissions file.'
            write(iout,'(A)') 'Either turn on in-line Ix emissions,'
            write(iout,'(2A)')'or add I2 and HOI emissions to',
     &                        ' the gridded emissions file,'
            write(iout,'(A)') 'or set the Ix emission flag to BYPASS'
            call camxerr()
          endif
        enddo
      endif
c
c   --- if doing VBS check that emissions are compatible ---
c
      l_exist_ivoc = .FALSE.
      l_exist_poa_xx = .FALSE.
      if( lvbs .AND. LVBSPREPROC ) then
         numspcs = nemspc
         do i=1,numspcs
           if( emspcname(i) .eq. 'IVOG      ' ) l_exist_ivoc = .TRUE.
           if( emspcname(i) .eq. 'IVOD      ' ) l_exist_ivoc = .TRUE.
           if( emspcname(i) .eq. 'IVOA      ' ) l_exist_ivoc = .TRUE.
           if( emspcname(i) .eq. 'IVOB      ' ) l_exist_ivoc = .TRUE.
           if( emspcname(i) .eq. 'POA_OP    ' ) then
              l_exist_poa_xx = .TRUE.
              l_exist_pap = .TRUE.
           endif
           if( emspcname(i) .eq. 'POA_GV    ' ) then
              l_exist_poa_xx = .TRUE.
              l_exist_pap = .TRUE.
           endif
           if( emspcname(i) .eq. 'POA_DV    ' ) then
              l_exist_poa_xx = .TRUE.
              l_exist_pap = .TRUE.
           endif
           if( emspcname(i) .eq. 'POA_MC    ' ) then
              l_exist_poa_xx = .TRUE.
              do j=0,NVOLBIN
                nemspc = nemspc + 1
                emspcname(nemspc) = spname(kpcp_c(j))
              enddo
           endif
           if( emspcname(i) .eq. 'POA_BB    ' ) then
              l_exist_poa_xx = .TRUE.
              do j=0,NVOLBIN
                nemspc = nemspc + 1
                emspcname(nemspc) = spname(kpfp_c(j))
              enddo
           endif
         enddo
         if( .NOT. l_exist_ivoc .OR. .NOT. l_exist_poa_xx ) then
            write(iout,'(//,a)') 'ERROR in STARTUP:'
            write(iout,'(2A)')'Emissions source species are not compatible',
     &                                                       ' with VBS.'
            write(iout,'(2A)')'VBS requires IVOC and sector-specific',
     &                                                   ' POA emissions.'
            call camxerr()
         endif
         if( l_exist_pap ) then
            do j=0,NVOLBIN
              nemspc = nemspc + 1
              emspcname(nemspc) = spname(kpap_c(j))
            enddo
         endif
      endif
c
c   --- call the pointer routine again to update for emissions arrays ---
c
      call iniptr(ncol,nrow)
c
c  --- call routine to allocate the emissions arrray ---
c
      call alloc_emiss_array(ngrid,ncol,nrow,nlayers_ems,MAX(1,nemspc))
c
c======================== Source Apportion Begin =======================
c
c
c   --- if this is a restart, call routine to read the
c       instantaneous files for fine grids ---
c
      if( (ltrace.OR.lddm.OR.lhddm) .AND. lrstrt .AND. ngrid .GT. 1 )
     &                               call rdfgsa(begdate,time)
c
c======================== Source Apportion End =======================
c
c
c======================== DDM Begin ====================================
c
      if( lddm .OR. lhddm ) then
c
c   --- get each emissions files to the correct place,
c       NOTE:  Done here for DDM because we need the point locations
c       from regular model ----
c
          do i=1,ngrid
            call emprepsa(begdate,begtim,i,numpts_lst)
            call ncf_emprepddm(i,begdate,begtim,enddate,endtim,numpts_lst)
          enddo
          do i=1,nptsrc
            lpigsa(i) = .FALSE.
          enddo
          if( ipigflg .EQ. GRESPIG ) then
             do i=1,nptsrc
               lpigsa(i) = lpiglet(i)
             enddo
          endif
      endif
c
c======================== DDM End ====================================
c
c-----Read IC or restart files headers and write output concentration
c     file headers 
c 
      if( is_netcdf_iic ) then
         call ncf_cncprep(begtim,begdate,endtim,enddate,lfirst_call)
      else
         call cncprep(endtim,enddate,lfirst_call) 
      endif
c
c-----Write deposition output file headers
c
      if (ldry .or. lwet) call depprep(endtim,enddate)
c
c======================== Source Apportion Begin =======================
c
      if( ltrace ) then
          if( lsfcfl .AND. lcdfout ) then
             do igrd=1,ngrid
                if( lddmcalc(igrd) ) then
                    write(action,'(2A,I3)') 'Writing TOPO to ',
     &                                   'SA average file. Grid: ',igrd
                    call ncf_wrt_topo(action,igrd,iowsfc(igrd),
     &                      ncol(igrd),nrow(igrd),topo(iptr2d(igrd)))
                    write(action,'(2A,I3)') 'Writing TOPO to ',
     &                                'SA deposition file. Grid: ',igrd
                    if( lptdepout ) call ncf_wrt_topo(action,igrd,
     &                           iowptdep(igrd),ncol(igrd),nrow(igrd),
     &                                              topo(iptr2d(igrd)))
                 endif
             enddo
          endif
          if( tectyp .EQ. RTRAC .OR. tectyp .EQ. RTCMC) then
c
c   --- get each emissions file to the correct place ----
c
              do i=1,nptsrc
                 lpigsa(i) = .FALSE.
              enddo
              do i=1,ngrid
                call empreprt(begdate,begtim,i)
                call ncf_empreprt(i,begdate,begtim,enddate,endtim,
     &                                        buffer_offset_iortem(i,1,1))
              enddo
c
c   ---- call routine to read the receptor definition file ---
c
              call rdrcprt()
c
c  --- call routine to write the header of receptor average file ---
c
              call hdrcprt(begdate,begtim,enddate,endtim)
c
c  ---- call routine to calculate the average reactivity
c       of initial conditions ---
c
          else
              if( .NOT. is_netcdf_iic ) then
                 call clciwt(begdate,begtim,enddate,endtim,
     &                    ncol(1),nrow(1),nlay(1),height(iptr3d(1)))
                 call cncprep(endtim,enddate,lfirst_call) 
              else
                 call ncf_clciwt(begdate,begtim,enddate,endtim,
     &                    ncol(1),nrow(1),nlay(1),height(iptr3d(1)))
                  call ncf_cncprep(begtim,begdate,endtim,enddate,lfirst_call)
              endif
c
c   --- get each emissions files to the correct place ----
c
              icheck = ntotsp
              iparm = MXTRSP
              cparm = 'MXTRSP'
              if( icheck .GT. iparm ) goto 7000
              do i=1,ngrid
                  call emprepsa(begdate,begtim,i,numpts_lst)
                  call ncf_emprepsa(i,begdate,begtim,enddate,endtim)
              enddo
              do i=1,nptsrc
                 lpigsa(i) = .FALSE.
              enddo
              if( ipigflg .EQ. GRESPIG ) then
                do i=1,nptsrc
                   lpigsa(i) = lpiglet(i)
                enddo
             endif
c
c  ---- call routine to calculate the average reactivity
c       of emissions ---
c
              call clcewt(enddate,endtim)
c
c   --- get each emissions files BACK to the correct place ----
c
              do i=1,ngrid
                  call emprepsa(begdate,begtim,i,numpts_lst)
              enddo
c
c   ---- call routine to read the receptor definition file ---
c
              call rercp()
c
c  --- call routine to write the header of receptor average file ---
c
              call hdrrcp(begdate,begtim,enddate,endtim)
c
c  --- call routine to get metfiles synced up again ---
c
              if( .NOT. is_netcdf_i3dmet(1) ) then
                 call metinit(1,ncol(1),nrow(1),nlay(1),
     &               endtim,enddate,xorg,yorg,delx,dely,
     &               height(iptr3d(1)),press(iptr3d(1)),
     &               depth(iptr3d(1)),windu(iptr3d(1)),
     &               windv(iptr3d(1)),tempk(iptr3d(1)),
     &               tsurf(iptr2d(1)),water(iptr3d(1)),
     &               rkv(iptr3d(1)),snow(iptr2d(1)),
     &               snowage(iptr2d(1)),snowalb(iptr2d(1)))
              endif
          endif
      endif
c
c-----Initialize point source locations and PiG submodel
c
      call pnt_update()
      if( ipigflg .NE. 0 ) call pigprep(begdate,begtim,numprocs)
c
c========================= Source Apportion End ========================
c
c======================== DDM Begin ====================================
c
          if( lddm .OR. lhddm ) then
c
c  --- print DDM turn-off flags (done here because we need mapgrd)
c
            do i=1,ngrid
              lddmtmp(i) = lddmcalc(i)
            enddo
            do i=1,ngrid
              lddmcalc(i) = lddmtmp( mapgrd(i) )
            enddo
            write(idiag,*) '    DDM turn-off flag table'
            write(idiag,*) '----------------------------------------'
            write(idiag,*) ' Internal  Original         Calculate'
            write(idiag,*) '   grid      grid              DDM?'
            write(idiag,*) '----------------------------------------'
            do i=1,ngrid
              write(idiag,'(1x,i6,4x,i6,9x,l10)') i,mapgrd(i),
     &                                                 lddmcalc(i)
            enddo
            write(idiag,*) '--------------------------------'

            do ip=1,ngrid
              if ( lddmcalc(ip) ) then
                do ic = 1,nchdrn(ip)
                  igrd = idchdrn(ic,ip)
                  if ( lddmcalc(igrd) ) CYCLE
                  write(iout,'(//,a)') 'ERROR in STARTUP:'
                  write(iout,*)
     &                     'DDM cannot be turned off in a nested grid',
     &                     'if its parent grid does DDM calculation.'
                  call camxerr()
                enddo
              endif
            enddo
c
c   ---- call routine to read the receptor definition file ---
c
          call rercp()
c
c  --- call routine to write the header of receptor average file ---
c
          call hdrrcpddm(begdate,begtim,enddate,endtim)
      endif
c
c======================== DDM End ====================================
c
c
c=================== Process Analysis Begin ==========================
c
c-----Call routine to initialize Process Analysis subdomains
c
      if( lproca ) then
         call pagrids()
c
c-----Call routine to write the header to the otuput files ---
c
         if( lipr ) call wrtiprhdr(begdate,begtim,enddate,endtim)
         if( lirr ) then
           call pasetup()
           icheck = ntotsp
           iparm = MXTRSP
           cparm = 'MXTRSP'
           if( icheck .GT. iparm ) goto 7000
           call wrtirrhdr(begdate,begtim,enddate,endtim)
           if( lsfcfl ) then
              if( l3davg(1) ) then
                  nlayer = nlay(1) 
              else
                  nlayer = 1
              endif
              call ncf_set_specatt_cpa(spec_units,spec_long_name,spec_desc,
     &                                                            spec_coords)
              call hdrwsa(iowsfc,sfcfil,'AVERAGE   ',ntotsp,nlayer,begdate,
     &                    begtim,enddate,endtim,spec_units,spec_long_name,
     &                                                  spec_desc,spec_coords)
           endif
           call alloc_procan_irr(ngrid,ncol,nrow,nlay,ntotsp,iout)
         endif
c
c----Need to allocate arrays that are used as arguments
c
      else
        call alloc_procan_ipa(ngrid,ncol,nrow,nlay,tectyp,0)
      endif
c
c===================== Process Analysis End ============================
c
c
c-----Determine times and dates for next inputs/emissions/write
c
      inptim = time
      inpdate = date
      emstim = time 
      emsdate = date
      bndtim = time 
      bnddate = date
      toptim = time 
      topdate = date
      ozntim = time 
      ozndate = date

      whr = aint(time/100.)
      wmn = amod(time,100.)
      wrttim = 100.*(whr + aint((wmn + dtout)/60.)) +
     &             amod((wmn + dtout),60.)
      wrtdate = date
      if (wrttim.ge.2400.) then 
        wrttim = wrttim - 2400. 
        wrtdate = wrtdate + 1 
        if( MOD(wrtdate,1000) .GT. 365 ) then
           if( MOD(INT(wrtdate/1000),4) .EQ. 0 ) then
              if( MOD(wrtdate,1000) .EQ. 367 )
     &                   wrtdate = (INT(wrtdate/1000)+1)*1000 + 1
           else
              wrtdate = (INT(wrtdate/1000)+1)*1000 + 1
           endif
        endif
      endif
c
c----Make sure all of the parameters used for local arrays are large enough ---
c
      do igrd=1,ngrid
         icheck = ncol(igrd)
         iparm = MXCELLS
         cparm = 'MXCELLS'
         if( icheck .GT. iparm ) goto 7000
         icheck = nrow(igrd)
         iparm = MXCELLS
         cparm = 'MXCELLS'
         if( icheck .GT. iparm ) goto 7000
         icheck = nlay(igrd)
         iparm = MXLAYER
         cparm = 'MXLAYER'
         if( icheck .GT. iparm ) goto 7000
         do ifile=1,nemiss_files(igrd)
            icheck = narspc(igrd,ifile)
         enddo
         iparm = MXSPEC
         cparm = 'MXSPEC'
         if( icheck .GT. iparm ) goto 7000
      enddo
      icheck = ntotsp
      iparm = MXTRSP
      cparm = 'MXTRSP'
      if( icheck .GT. iparm ) goto 7000
      icheck = ntrcls
      iparm = MXALCLS
      cparm = 'MXALCLS'
      if( icheck .GT. iparm ) goto 7000
      icheck = ntrcls
      iparm = MXTRCLS
      cparm = 'MXTRCLS'
      if( icheck .GT. iparm ) goto 7000
      do ism = 1,nsample
         icheck = ncolsmp(ism)
         iparm = MXCOLSMP
         cparm = 'MXCOLSMP'
         if( icheck .GT. iparm ) goto 7000
         icheck = nrowsmp(ism)
         iparm = MXROWSMP
         cparm = 'MXROWSMP'
         if( icheck .GT. iparm ) goto 7000
      enddo
      icheck = ngrid
      iparm = MXGRID
      cparm = 'MXGRID'
      if( icheck .GT. iparm ) goto 7000
      icheck = nptsrc
      iparm = MXPTSRC
      cparm = 'MXPTSRC'
      if( icheck .GT. iparm ) goto 7000
      icheck = nreact
      iparm = MXREACT
      cparm = 'MXREACT'
      if( icheck .GT. iparm ) goto 7000
      icheck = nreactr
      iparm = MXREACT
      cparm = 'MXREACT'
      if( icheck .GT. iparm ) goto 7000
      icheck = nspec+1
      iparm = MXSPEC
      cparm = 'MXSPEC'
      if( icheck .GT. iparm ) goto 7000
      icheck = navspc
      iparm = MXSPEC
      cparm = 'MXSPEC'
      if( icheck .GT. iparm ) goto 7000
      icheck = nbcspc
      iparm = MXSPEC
      cparm = 'MXSPEC'
      if( icheck .GT. iparm ) goto 7000
      icheck = nicspc
      iparm = MXSPEC
      cparm = 'MXSPEC'
      if( icheck .GT. iparm ) goto 7000
      do i=1,npoint_files
         icheck = nptspc(i)
         iparm = MXSPEC
         cparm = 'MXSPEC'
         if( icheck .GT. iparm ) goto 7000
      enddo
      icheck = nrtrac
      iparm = MXSPEC
      cparm = 'MXSPEC'
      if( icheck .GT. iparm ) goto 7000
c
c-----Everything worked correctly, return to calling routine
c
      call flush(iout)
      call flush(idiag)
      return
c
c-----Error message for array bounds
c
 7000 continue
      write(iout,'(//,A)') 'ERROR in STARTUP:'
      write(iout,*) 'A parameter in the camx.prm is not ',
     &                                        'sufficiently large.'
      write(iout,*) 'Please change the value for parameter: ',cparm
      write(iout,*) 'It should be set to a value of at least: ',icheck
      call flush(iout)
      call camxerr()
c
c-----Return point
c

      end
