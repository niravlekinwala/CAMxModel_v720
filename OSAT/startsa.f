c**** STARTSA
c
      subroutine startsa(iounit,nopen)
      use grid
      use chmstry
      use filunit
      use ptemiss
      use procan
      use tracer
      implicit none
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c     This routine loads all of the files needed for the source 
c     apportionment algorithm.  The output files are opened, the input
c     files are opened as needed.
c
c      Copyright 1996 - 2022
c     Ramboll
c
c     Argument description:
c        Inputs:
c          iounit     I   unit number for output
c        Outputs:   
c          nopen    I   number of files opened
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     12/16/99   --gwilson--   Fixed a bug which caused the model to try
c                              and read a restart file for each find grid
c     07/19/02   --gwilson--   Added code for source area map for nests
c     10/06/04   --cemery --   Restructured for namelist input
c      8/23/06   --cemery--    Instantaneous restart files reduced to 1
c                              per grid type
c      8/25/06   --cemery--    Surface output files now all UAM format,
c                              one file per grid
c     03/15/09   --gwilson--   Added code for deposition output for tracers
c     05/07/12   --cemery--    Added flexi-nesting flag
c     03/01/16   --gwilson--   Added partial source area map
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'flags.inc'
      include 'namelist.inc'
      include 'netcdf.inc'
c
c-----------------------------------------------------------------------
c    Argument declaration:
c-----------------------------------------------------------------------
c
      integer   iounit
      integer   nopen
c
c-----------------------------------------------------------------------
c    External functions:
c-----------------------------------------------------------------------
c
      integer istrln
      integer ncf_chkfile
      integer ncf_get_nlayers
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      character*200 fname
      character*80 action
      character*10 cname
      integer      igrd, igrp, i, j, n, idxfile, ierr, nlays_in
      logical      lexist
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c-----if doing SA then point source files are read from SA_control section ----
c
      if( tectyp .EQ. SA .AND. luse_points ) then
         action = 'Reading point source filenames.'
         goto 7009
      endif
c
c   --- if option is turned off, just return ---
c
      if( .NOT. ltrace ) goto 9999
      nchar = istrln( flrtsa )
c
c  -- call routine to allocate arrays ----
c
      call alloc_tracer(ngroup,ngrid,ncol,nrow,nlay,nspec)
      call alloc_tracer_full(ngroup,ngrid,ncol,nrow)
      mapfil = ' '
c
c   --- intialize first time through to true ---
c
      lfirst = .TRUE.
      write(iounit,'(/,A,/)') 
     &            '           **** Source Apportionment Files ****'
c
c   --- source area mapping file ---
c
      if( nregin .GT. 0 ) then
         do 10 igrd = 1,ngrid
            write(action,'(A,I2)')
     &                        'Opening Source Area Map for grid: ',igrd
            mapfil(0,igrd) = SA_Source_Area_Map(igrd)
            fname = mapfil(0,igrd)
            if( fname .EQ. ' ' ) then
               if( igrd .EQ. 1 ) goto 7004
               lmapfl(0,igrd) = .FALSE.
               write(iounit,9000)
     &               'No Source Area Map file for Grid #         :',igrd
               goto 10
            endif
            inquire(file=mapfil(0,igrd),exist=lexist)
            if( .NOT. lexist ) goto 7000
            nopen =  nopen + 1
            call getunit(iormap(0,igrd))
            write(iounit,9001)
     &               'Source area map file for grid # ',igrd,'   (unit):',
     &                                                      iormap(0,igrd)
            write(iounit,9002) mapfil(0,igrd)(:istrln(mapfil(0,igrd)))
            lmapfl(0,igrd) = .TRUE.
   10    continue
         if( lpartial ) then
            do igrp=1,ngroup
               do 50 igrd = 1,ngrid
                 write(action,'(2A,I2,A,I2)') 'Opening Partial Source Area ',
     &                             'Map for grid: ',igrd,' group: ',igrp
                 mapfil(igrp,igrd) = Partial_Source_Area_Map(igrp,igrd)
                 fname = mapfil(igrp,igrd)
                 if( fname .EQ. ' ' ) then
                    lmapfl(igrp,igrd) = .FALSE.
                    write(iounit,9000)
     &                 'No Partial Source Area Map file for Grid # :',
     &                                             igrd,' group: ',igrp
                    goto 50
                 endif
                 inquire(file=mapfil(igrp,igrd),exist=lexist)
                 if( .NOT. lexist ) goto 7000
                 nopen =  nopen + 1
                 call getunit(iormap(igrp,igrd))
                 write(iounit,9003)
     &               'Partial Source area map file for grid # ',igrd,
     &                  ' group: ',igrp,'   (unit):',iormap(igrp,igrd)
                 write(iounit,9002) mapfil(igrp,igrd)
     &                              (:istrln(mapfil(igrp,igrd)))
                 lmapfl(igrp,igrd) = .TRUE.
   50          continue
            enddo
         endif
      endif
c
c   --- receptor definition file ----
c
      rcpfil = SA_Receptor_Definitions
      if( rcpfil .NE. ' ' ) then
         action = 'Opening SA Receptor Definition file.'
         fname = rcpfil
         inquire(file=rcpfil,exist=lexist)
         if( .NOT. lexist ) goto 7000
         nopen =  nopen + 1
         call getunit(iorrcp)
         write(iounit,9000)'SA Receptor definition file          (unit):',
     &                                                          iorrcp
         write(iounit,9002) rcpfil(:istrln(rcpfil))
         lrcpfil = .TRUE.
      else
         write(iounit,9000)'No SA Receptor definition file provided.'
         lrcpfil = .FALSE.
      endif
c
c   --- instantaneous file used for initialization ---
c
      if( lrstrt ) then
         action = 'Opening SA master grid Restart file.'
         inifil(IDXCRS) = SA_Master_Restart
         fname = inifil(IDXCRS)
         inquire(file=inifil(IDXCRS),exist=lexist)
         if( .NOT. lexist ) goto 7000
         nopen =  nopen + 1
         call getunit(iorini(IDXCRS))
         write(iounit,9000)
     &                'SA master grid Restart file          (unit):',
     &                                                  iorini(IDXCRS)
         write(iounit,9002) inifil(IDXCRS)(:istrln(inifil(IDXCRS)))
c
         if( ngrid .GT. 1 ) then
            action = 'Opening SA nested grid Restart file.'
            inifil(IDXFIN) = SA_Nested_Restart
            if( inifil(IDXFIN) .NE. ' ' ) then
                fname = inifil(IDXFIN)
                inquire(file=inifil(IDXFIN),exist=lexist)
                if( .NOT. lexist ) goto 7000
                nopen =  nopen + 1
                call getunit(iorini(IDXFIN))
                write(iounit,9000)
     &                  'SA nested grid Restart file          (unit):',
     &                                                    iorini(IDXFIN)
                write(iounit,9002) inifil(IDXFIN)(:istrln(inifil(IDXFIN)))
            elseif (.not.lflexi) then
                goto 7003
            endif
         endif
      endif
c
c   --- figure out how many files in each grid/group ---
c
      num_iortem = 0
      do j=1,ngrid
         num_iortem(j,0) = nemiss_files(j)
         if( nemiss_files(j) .GT. MXFILES ) goto 7006
         do idxfile=1,nemiss_files(j)
            SA_Emiss_Group_Grid(0,j,idxfile) = Emiss_Grid(j,idxfile)
         enddo
      enddo
      if( npoint_files .GT. MXFILES ) goto 7006
      do idxfile=1,npoint_files
         SA_Points_Group(0,idxfile) = Point_Sources(idxfile)
      enddo
      do i=0,ngroup
         do j=1,ngrid
            do idxfile=1,MXFILES
              if( istrln(SA_Emiss_Group_Grid(i,j,idxfile)) .NE. 0 )  
     &                                     num_iortem(j,i) = idxfile
            enddo
         enddo
      enddo
      num_iortpt = 0
      do i=0,ngroup
        do idxfile=1,MXFILES
          if( istrln(SA_Points_Group(i,idxfile)) .NE. 0 )
     &                                      num_iortpt(i) = idxfile
        enddo
      enddo
      call alloc_tracer_emfiles(ngroup,ngrid,nemiss_files,npoint_files,
     &                                      num_iortem,num_iortpt,nspec)
c
c   --- emissions files for the source groupings ---
c
      do 20 i = 1,ngroup
c
c   --- surface emissions file --- 
c
         do 30 j = 1,ngrid
            action = 'Opening SA Gridded Emissions file.'
            if( num_iortem(j,i) .EQ. 0 ) then
                write(iounit,9003)
     &                 'Gridded Emissions file for grid#/group#    :',
     &                                               j,i,' Not supplied'
            else
                do idxfile=1,num_iortem(j,i)
                   if( .NOT. larsrc ) cycle
                   temfil(j,i,idxfile) = ' '
                   ltemfl(j,i,idxfile) = .FALSE.
                   is_netcdf_iortem(j,i,idxfile) = .FALSE.
                   temfil(j,i,idxfile) = SA_Emiss_Group_Grid(i,j,idxfile)
                   ltemfl(j,i,idxfile) = .TRUE.
                   fname = temfil(j,i,idxfile)
                   if( fname .EQ. ' ' ) goto 7007
                   inquire(file=temfil(j,i,idxfile),exist=lexist)
                   if( .NOT. lexist ) goto 7000
                   nopen = nopen + 1
                   call getunit(iortem(j,i,idxfile))
                   cname = 'EMISSIONS '
                   ierr = ncf_chkfile(iortem(j,i,idxfile),fname,action,cname)
                   if( ierr .EQ. ISUCES ) then
                     is_netcdf_iortem(j,i,idxfile) = .FALSE.
                     open(unit=iortem(j,i,idxfile),file=fname,
     &                    ERR=7001,form='UNFORMATTED',status='UNKNOWN')
                   else
                     is_netcdf_iortem(j,i,idxfile) = .TRUE.
                     ierr = nf_open(fname,NF_NOWRITE,iortem(j,i,idxfile))
                     if( ierr .NE. NF_NOERR ) goto 7001
                     nlayers_ems = MAX(nlayers_ems,ncf_get_nlayers(iortem(j,i,idxfile),fname,action))
                   endif
                   write(iounit,9005)
     &                'Gridded Emissions for grid#/group#/file#   (unit):',
     &                                            j,i,idxfile,iortem(j,i,idxfile)
                   write(iounit,9002) temfil(j,i,idxfile)(:istrln(temfil(j,i,idxfile)))
                enddo
            endif
c
   30    continue
c
c   --- elevated point source emissions ----
c
         if( num_iortpt(i) .EQ. 0 ) then
            write(iounit,9004)
     &               'Point Source Emissions file for group#     :',
     &                                                 i,' Not supplied'
         else
            do idxfile=1,num_iortpt(i)
               if( .NOT. lptsrc ) cycle
               tptfil(i,idxfile) = ' '
               ltptfl(i,idxfile) = .FALSE.
               is_netcdf_iortpt(i,idxfile) = .FALSE.
               tptfil(i,idxfile) = SA_Points_Group(i,idxfile)
               fname = tptfil(i,idxfile)
               if( fname .EQ. ' ' ) goto 7008
               action = 'Opening SA Point Source file.'
               if( fname .EQ. ' ' ) then
                  tptfil(i,idxfile) = ' '
                  ltptfl(i,idxfile) = .FALSE.
                  is_netcdf_iortpt(i,idxfile) = .FALSE.
               else
                  ltptfl(i,idxfile) = .TRUE.
                  inquire(file=tptfil(i,idxfile),exist=lexist)
                  if( .NOT. lexist ) goto 7000
                  nopen = nopen + 1
                  call getunit(iortpt(i,idxfile))
                  cname = 'PTSOURCE  '
                  ierr = ncf_chkfile(iortpt(i,idxfile),tptfil(i,idxfile),action,cname)
                  if( ierr .EQ. ISUCES ) then
                    is_netcdf_iortpt(i,idxfile) = .FALSE.
                    open(unit=iortpt(i,idxfile),file=tptfil(i,idxfile),ERR=7001,
     &                        form='UNFORMATTED',status='UNKNOWN')
                  else
                    is_netcdf_iortpt(i,idxfile) = .TRUE.
                    ierr = nf_open(tptfil(i,idxfile),NF_NOWRITE,iortpt(i,idxfile))
                    if( ierr .NE. NF_NOERR ) goto 7001
                  endif
                  write(iounit,9007)
     &                 'Point Source Emissions for group#    (unit):',
     &                                                       i,iortpt(i,idxfile)
                  write(iounit,9002) tptfil(i,idxfile)(:istrln(tptfil(i,idxfile)))
               endif
            enddo
         endif
   20 continue
c
c   --- initial conditions file ---
c
      icfil = SA_Initial_Conditions
      if( icfil .NE. ' ' ) then
         action = 'Opening SA Initial Condition file.'
         fname = icfil
         inquire(file=icfil,exist=lexist)
         if( .NOT. lexist ) goto 7000
         nopen =  nopen + 1
         call getunit(ioric)
         cname = 'AIRQUALITY'
         ierr = ncf_chkfile(ioric,icfil,action,cname)
         if( ierr .EQ. ISUCES ) then
            is_netcdf_ioric = .FALSE.
            open(unit=ioric,file=icfil,ERR=7002,
     &                           form='UNFORMATTED',status='UNKNOWN')
          else
            is_netcdf_ioric = .TRUE.
            ierr = nf_open(icfil,NF_NOWRITE,ioric)
            if( ierr .NE. NF_NOERR ) goto 7001
          endif
          lsa_ioric = .TRUE.
          write(iounit,9000)'SA Initial Condition file            (unit):',
     &                                                          ioric
          write(iounit,9002) icfil(:istrln(icfil))
      else
         write(iounit,9000)'No SA Initial Condition file provided.'
         lsa_ioric = .FALSE.
      endif
c
c   --- boundary conditions file ---
c
      bcfil = SA_Boundary_Conditions
      if( bcfil .NE. ' ' ) then
         action = 'Opening SA Boundary Condition file.'
         fname = bcfil
         inquire(file=bcfil,exist=lexist)
         if( .NOT. lexist ) goto 7000
         call getunit(iorbc)
         nopen =  nopen + 1
         cname = 'BOUNDARY  '
         ierr = ncf_chkfile(iorbc,bcfil,action,cname)
         if( ierr .EQ. ISUCES ) then
            is_netcdf_iorbc = .FALSE.
            open(unit=iorbc,file=bcfil,ERR=7002,
     &                            form='UNFORMATTED',status='UNKNOWN')
         else
            is_netcdf_iorbc = .TRUE.
            ierr = nf_open(bcfil,NF_NOWRITE,iorbc)
            if( ierr .NE. NF_NOERR ) goto 7001
         endif
         write(iounit,9000)'SA Boundary Condition file           (unit):',
     &                                                             iorbc
         write(iounit,9002) bcfil(:istrln(bcfil))
         lsa_iorbc = .TRUE.
      else
         write(iounit,9000)'No SA Boundary Condition file provided.'
         lsa_iorbc = .FALSE.
      endif
c
c
c   --- boundary conditions file ---
c
      tcfil = SA_Top_Concentrations
      if( tcfil .NE. ' ' ) then
         action = 'Opening SA Top Concentrations file.'
         if( .NOT. lsa_iorbc ) goto 7005
         lsa_iortc = .TRUE.
         fname = tcfil
         inquire(file=tcfil,exist=lexist)
         if( .NOT. lexist ) goto 7000
         nopen =  nopen + 1
         if( is_netcdf_iorbc ) then
            ierr = nf_open(tcfil,NF_NOWRITE,iortc)
            if( ierr .NE. NF_NOERR ) goto 7001
         else
            call getunit(iortc)
            write(iounit,9000)'SA Top Concentrations file           (unit):',
     &                                                          iortc
            write(iounit,9002) tcfil(:istrln(tcfil))
            open(unit=iortc,file=tcfil,ERR=7002,
     &                           form='UNFORMATTED',status='UNKNOWN')
         endif
      else
         write(iounit,9000)'No SA Top Concentrations file provided.'
         lsa_iortc = .FALSE.
      endif
c
c   --- output filenames ---
c
      do 40 i=IDXCRS,IDXFIN
         if( i .EQ. IDXFIN .AND. ngrid .EQ. 1 ) goto 40
c
c   --- output instantaneous file ---
c
         confil(i) = flrtsa
         if( i .EQ. IDXCRS ) then
            confil(i)(nchar+1:) = '.sa.inst'
         else
            confil(i)(nchar+1:) = '.sa.finst'
         endif
         fname = confil(i)
         nopen = nopen + 1
         call getunit(iowcon(i))
         open(unit=iowcon(i),file=confil(i),ERR=7002,
     &                           form='UNFORMATTED',status='UNKNOWN')
         if( i .EQ. IDXCRS ) then
            write(iounit,9000)
     &                 'SA INST file for master grid         (unit):',
     &                                                         iowcon(i)
            write(iounit,9002) confil(i)(:istrln(confil(i)))
         else
            write(iounit,9000)
     &                 'SA INST file for nested grids        (unit):',
     &                                                         iowcon(i)
            write(iounit,9002) confil(i)(:istrln(confil(i)))
         endif
   40 continue
c
c    ---- surface concentrations file ----
c
      if( .NOT. lsfcfl ) then
         write(iounit,9006)
     &                 'SA Surface file                      (unit):',
     &                                                  ' Not supplied'
      else
         do n = 1,ngrid
            fname = flrtsa
            if( .NOT. lcdfout ) then
               write(fname(nchar+1:),'(a,i2.2)') '.sa.grd',n
               sfcfil(n) = fname
               nopen = nopen + 1
               call getunit(iowsfc(n))
               open(unit=iowsfc(n),file=sfcfil(n),ERR=7002,
     &                            form='UNFORMATTED',status='UNKNOWN')
               write(iounit,9000)
     &                     'SA Surface file                      (unit):',
     &                                                            iowsfc(n)
               write(iounit,9002) sfcfil(n)(:istrln(sfcfil(n)))
            else
               write(fname(nchar+1:),'(a,i2.2,a)') '.sa.grd',n,'.nc'
               sfcfil(n) = fname
               nopen = nopen + 1
               action = 'Opening SA surface output file.'
               call ncf_createfile(fname,action,iowsfc(n))
               write(iounit,9000) 'SA Surface file:'
               write(iounit,9002) sfcfil(n)(:istrln(sfcfil(n)))
            endif
         enddo
c
c  --- if doing deposition output, open the files ---
c
         if( lptdepout ) then
            do n = 1,ngrid
               fname = flrtsa
               if( .NOT. lcdfout ) then
                  write(fname(nchar+1:),'(a,i2.2)') '.sa.depn.grd',n
                  ptdepfil(n) = fname
                  nopen = nopen + 1
                  call getunit(iowptdep(n))
                  open(unit=iowptdep(n),file=ptdepfil(n),ERR=7002,
     &                            form='UNFORMATTED',status='UNKNOWN')
                  write(iounit,9000)
     &                     'SA Deposition file                   (unit):',
     &                                                      iowptdep(n)
                  write(iounit,9002) ptdepfil(n)(:istrln(ptdepfil(n)))
               else
                  write(fname(nchar+1:),'(a,i2.2,a)') '.sa.depn.grd',n,'.nc'
                  ptdepfil(n) = fname
                  nopen = nopen + 1
                  call ncf_createfile(fname,action,iowptdep(n))
                  action = 'Opening SA deposition output file.'
                  write(iounit,9000) 'SA Deposition file:'
                  write(iounit,9002) ptdepfil(n)(:istrln(ptdepfil(n)))
               endif
            enddo
         endif
      endif
c
c    ---- tracer receptor file ----
c
      if( lrcpfil ) then
         avgfil(1:) = flrtsa
         avgfil(nchar+1:) = '.sa.receptor'
         fname = avgfil
         nopen = nopen + 1
         call getunit(iowrcp)
         open(unit=iowrcp,file=avgfil,ERR=7002,form='FORMATTED',
     &                                               status='UNKNOWN')
         write(iounit,9000)
     &            'SA Receptor concentration file       (unit):',iowrcp
         write(iounit,9002) avgfil(:istrln(avgfil))
       endif
c
c  --- return to calling routine ---
c
      goto 9999
c
c-----------------------------------------------------------------------
c    Format statements:
c-----------------------------------------------------------------------
c
 9000 format(/,A,I8)
 9001 format(/,A,I2,A,I8)
 9002 format(2A)
 9003 format(/,A,I2,1X,I2,A,I8)
 9004 format(/,A,I2,A)
 9005 format(/,A,I2,1X,I2,1X,I2,1X,I8)
 9006 format(/,2A)
 9007 format(/,A,I2,1X,I8)
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iounit,'(//,A)') 'ERROR in STARTSA:'
      write(iounit,'(A)') action(:istrln(action))
      write(iounit,'(/,1X,2A)') 'Input file does not exist: ',
     &                                             fname(:istrln(fname))
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7001 continue
      write(iounit,'(//,A)') 'ERROR in STARTSA:'
      write(iounit,'(A)') action(:istrln(action))
      write(iounit,'(/,1X,2A)') 'Cannot open input file: ',
     &                                             fname(:istrln(fname))
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7002 continue
      write(iounit,'(//,A)') 'ERROR in STARTSA:'
      write(iounit,'(/,1X,2A)') 'Cannot open output file: ',
     &                                             fname(:istrln(fname))
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7003 continue
      write(iounit,'(//,A)') 'ERROR in STARTSA:'
      write(iounit,'(/,1X,3A,/,A)') 'Nested grid restart file for ',
     &                     tectyp(:istrln(tectyp)),' not supplied.',
     &                  ' Set Flexi_Nest = .true. to ignore this file.'
      call camxerr()
c
 7004 continue
      write(iounit,'(//,A)') 'ERROR in STARTSA:'
      write(iounit,'(/,2A)') 'A source area mapping file must be ',
     &                'supplied for the master grid, Grid #1'
      call camxerr()
c
 7005 continue
      write(iounit,'(//,A)') 'ERROR in STARTSA:'
      write(iounit,'(/,2A,/,A)') 'You supplied a file for ',
     &                            'SA_Top_Concentrations but not',
     &                                      'for SA_Boundary_Conditions.'
      write(iounit,'(/,2A,/,A)') 'Top concentrations for source ',
     &                       'apportionment are not',
     &                    'allowed without a boundary conditions file.'
      call camxerr()
c
 7006 continue
      write(iounit,'(//,A)') 'ERROR in STARTSA:'
      write(iounit,'(1X,A)') 'Number of SA emissions files exceeds max.'
      write(iounit,'(1X,2A)') 'Increase parameter MXFILES and recompile.'
      write(iounit,'(1X,2A)') 'It must be at least: ',
     &                              MAX(npoint_files,nemiss_files(j))
      call camxerr()
c
 7007 continue
      write(iounit,'(//,A)') 'ERROR in STARTSA:'
      write(iounit,'(A)') action(:istrln(action))
      write(iounit,'(1X,2A)') 'You supplied a blank filename in the middle ',
     &                  'of the list of gridded emissions '
      write(iounit,'(1X,A,I3)') 'filenames for grid: ',i
      write(iounit,'(1X,2A)') 'Check the numbering of the list of files.'
      call camxerr()
c
 7008 continue
      write(iounit,'(//,A)') 'ERROR in STARTSA:'
      write(iounit,'(A)') action(:istrln(action))
      write(iounit,'(1X,2A)') 'You supplied a blank filename in the middle ',
     &                  'of the list of point source '
      write(iounit,'(1X,A,I3)') 'filenames for grid: ',i
      write(iounit,'(1X,2A)') 'Check the numbering of the list of files.'
      call camxerr()
c
 7009 continue
      write(iounit,'(//,a)') 'ERROR in STARTSA:'
      write(iounit,'(A)') action(:istrln(action))
      write(iounit,'(1X,2A)') 'In this version the point source inventory ',
     &           'for a simulation using the SA '
      write(iounit,'(1X,2A)') 'probing tool is provided in the SA_Control ',
     &            'section of the CAMx control file.'
      write(iounit,'(1X,2A)') 'You have supplied filenames in the ',
     &                                         'CAMx_control section.'
      write(iounit,'(1X,2A)') 'Replace these filenames with blank strings ',
     &                              'and include the point source files'
      write(iounit,'(1X,2A)') 'in the SA_Control section.'
      call camxerr()
c
 7010 continue
      write(iounit,'(//,a)') 'ERROR in STARTSA:'
      write(iounit,'(A)') action(:istrln(action))
      write(iounit,'(1X,2A,/,1X,A)') 'In this version the leftover group is ',
     &           'not allowed when using point sources','in any source group.'
      write(iounit,'(1X,2A)') 'Your source groups contain point sources ',
     &            'and the leftover flag is set to .TRUE.'
      write(iounit,'(1X,2A,/,1X,A)') 'When using point sources in any SA emissions ',
     &                 'group you must supply the','entire inventory as emissions groups.'
      write(iounit,'(1X,2A,/,1X,A)') 'Make sure all point sources are included ',
     &         'in emissions groups and turn off','the leftover group flag.'
      call camxerr()
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
