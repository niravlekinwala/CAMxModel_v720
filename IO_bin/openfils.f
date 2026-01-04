      subroutine openfils(ii,nopen)
      use camxcom
      use camxfld
      use filunit
      use chmstry
      use o3colmap
      use grid
      use pigsty
      use tracer
      use procan
      implicit none
c
c----CAMx v7.20 220430
c
c     OPENFILS opens Fortran I/O files
c                          
c      Copyright 1996 - 2022
c     Ramboll
c          
c     Modifications:
c        7/5/02    Changed to account for new type of PiG flag
c        8/30/02   Modified to read combined cloud/rain file, and now
c                  water vapor and cloud/rain files can be provided for each
c                  nest
c        01/30/02  Added code for RTRAC probing tool
c        1/10/03   Added open of deposition output file
c        4/2/03    Removed option for UAM-V type cloud file
c        3/3/04    Added checks for reading water files
c                  if either chemistry, dry dep or wet dep is on
c        3/19/04   Added checks for reading AHO file if
c                  if either chemistry or dry dep is on
c        10/4/04   Restructured for namelist input
c        6/21/05   Cloud/rain header modified for new structure
c        8/02/05   Added open of general PiG sampling grid files
c        8/25/05   PiG restart file is now optional
c        8/23/06   Instantaneous restart files reduced to 1 per grid type
c        8/25/06   Average and deposition output files now all UAM format,
c                  one file per grid
c        11/4/09   Removed input top concentrations
c        1/04/11   Revised for new met input format
c        04/02/12  Removed drought stress and snow flag; AHO
c                  file is now just ozone column
c        05/07/12  Added flexi-nesting flag
c        04/30/13  Added surface model
c        4/07/14   Added top con file
c        09/02/14  Added subgrid convective model
c
c     Input arguments:
c        ii                  Length of fileroot string 
c        nopen               number of files opened
c
c     Output arguments:
c        none
c
c     Routines Called:
c        none
c
c     Called by:
c        READNML
c
      include 'camx.prm'
      include 'flags.inc'
      include 'namelist.inc'
      include 'netcdf.inc'
c
c-----External functions
c
      integer istrln
      integer ncf_chkfile
      integer ncf_get_nlayers
c
      character*200 filtmp
      character*80  action
      character*10  cname
      character*4   name(10)
      logical       lexist
      integer       nopen,nfils,ii,n,ifile,ierr
c
c-----Entry point
c
c-----Open text output files
c
      nlayers_ems = 1
      filtmp = ' '
      filtmp = filroot
      write(iout,9000) 'Output OUT message file              (unit):',
     &                                                             iout
      write(iout,9003) filtmp(1:ii),'.out'
      write(iout,9000) 'Output DIAG diagnostic file          (unit):',
     &                                                             idiag
      write(iout,9003) filtmp(1:ii),'.diag'
      write(iout,9000) 'Output MASS summary file             (unit):',
     &                                                             imass
      write(iout,9003) filtmp(1:ii),'.mass'
c
c-----Open master grid instantaneous concentration output file
c
      filtmp(ii+1:) = '.inst'
      nopen = nopen + 1
      call getunit(iconc)
      action = 'Opening output INST file for master grid.'
      open(unit=iconc,file=filtmp(1:ii+5),form='UNFORMATTED',
     &                                    status= 'UNKNOWN',ERR=7000)
      write(iout,9000) 'Output INST file for master grid     (unit):',
     &                                                           iconc
      write(iout,9002) filtmp(1:ii+5)
c
c-----Open fine grid instantaneous concentration output file
c
      nfils = 4
      if (ngrid.gt.1) then
        filtmp(ii+1:) = '.finst'
        nopen = nopen + 1
        call getunit(ifconc)
        action = 'Opening output INST file for nest grids.'
        open(unit=ifconc,file=filtmp(1:ii+6),form='UNFORMATTED',
     &                                    status= 'UNKNOWN',ERR=7000)
        write(iout,9000) 'Output INST file for nest grids      (unit):',
     &                                                            ifconc
        write(iout,9002) filtmp(1:ii+6)
        nfils = nfils + 1
      endif
c
c-----Open average concentration output file(s)
c
      if (navspc.gt.0) then
        filtmp(ii+1:) = '.avrg'
        do n = 1,ngrid
          if( .NOT. lcdfout ) then
             write(filtmp(ii+6:),'(a,i2.2)') '.grd',n
             nopen = nopen + 1
             call getunit(iavg(n))
             action = 'Opening output AVERAGE file'
             open(unit=iavg(n),file=filtmp(1:ii+11),form='UNFORMATTED',
     &                                    status='UNKNOWN',ERR=7000)
             write(iout,9000) 
     &           'Output AVERAGE file                  (unit):',iavg(n)
             write(iout,9002) filtmp(1:ii+11)
             nfils = nfils + 1
          else
             write(filtmp(ii+6:),'(a,i2.2,a)') '.grd',n,'.nc'
             action = 'Opening output AVERAGE file for NCF'
             call ncf_createfile(filtmp,action,iavg(n))
             nopen = nopen + 1
             write(iout,9000) 'Output AVERAGE file:'
             write(iout,9002) filtmp(1:ii+14)
             nfils = nfils + 1
          endif
        enddo
      endif
c
c-----Open deposition output file(s)
c
      if (navspc.gt.0 .and. (ldry .or. lwet)) then
        filtmp(ii+1:) = '.depn'
        do n = 1,ngrid
          if( .NOT. lcdfout ) then
             write(filtmp(ii+6:),'(a,i2.2)') '.grd',n
             nopen = nopen + 1
             call getunit(idep(n))
             action = 'Opening output DEPOSITION file'
             open(unit=idep(n),file=filtmp(1:ii+11),form='UNFORMATTED',
     &                                   status= 'UNKNOWN',ERR=7000)
             write(iout,9000) 
     &           'Output DEPOSITION file               (unit):',idep(n)
             write(iout,9002) filtmp(1:ii+11)
             nfils = nfils + 1
          else
             write(filtmp(ii+6:),'(a,i2.2,a)') '.grd',n,'.nc'
             nopen = nopen + 1
             action = 'Opening output DEPOSITION file'
             call ncf_createfile(filtmp,action,idep(n))
             write(iout,9000) 'Output DEPOSITION file:'
             write(iout,9002) filtmp(1:ii+14)
             nfils = nfils + 1
          endif
        enddo
      endif
c
c-----Open surface model output file(s)
c
      if (lsrfmod) then
        filtmp(ii+1:) = '.srf'
        do n = 1,ngrid
          if( .NOT. lcdfout ) then
             write(filtmp(ii+5:),'(a,i2.2)') '.grd',n
             nopen = nopen + 1
             call getunit(ismout(n))
             action = 'Opening output SURFACE MODEL file'
             open(unit=ismout(n),file=filtmp(1:ii+11),form='UNFORMATTED',
     &                                   status= 'UNKNOWN',ERR=7000)
             write(iout,9000)
     &           'Output SURFACE MODEL file:           (unit):',ismout(n)
             write(iout,9002) filtmp(1:ii+11)
             nfils = nfils + 1
           else
             write(filtmp(ii+5:),'(a,i2.2,a)') '.grd',n,'.nc'
             nopen = nopen + 1
             action = 'Opening output SURFACE MODEL file'
             call ncf_createfile(filtmp,action,ismout(n))
             write(iout,9000) 'Output SURFACE MODEL file:'
             write(iout,9002) filtmp(1:ii+14)
             nfils = nfils + 1
           endif
        enddo
      endif
c
c-----Open PiG output file
c
      if( ipigflg .NE. 0 ) then
        filtmp(ii+1:) = '.pig'
        nopen = nopen + 1
        call getunit(ipig)
        action = 'Opening output PiG file.'
        open(unit=ipig,file=filtmp(1:ii+4),form='UNFORMATTED',
     &                             status= 'UNKNOWN',ERR=7000)
        write(iout,9000) 'Output PiG file                      (unit):',
     &                                                              ipig
        write(iout,9002) filtmp(1:ii+4)
        nfils = nfils + 1
c
c-----Open PiG sampling grid output file(s)
c
        if (lsample) then
          filtmp(ii+1:) = '.smp'
          do n = 1,nsample
            if( .NOT. lcdfout ) then
               write(filtmp(ii+5:),'(i2.2,a)') n
               nopen = nopen + 1
               call getunit(isample(n))
               action = 'Opening output PiG sampling grid file.'
               open(unit=isample(n),file=filtmp(1:ii+6),form='UNFORMATTED',
     &              status='UNKNOWN',ERR=7000)
               write(iout,9000)
     &                   'PiG Sampling Grid output file        (unit):',
     &                                                        isample(n)
               write(iout,9002) filtmp(1:ii+6)
             else
               write(filtmp(ii+5:),'(i2.2,a)') n,'.nc'
               nopen = nopen + 1
               call getunit(isample(n))
               action = 'Opening output PiG sampling grid file.'
               call ncf_createfile(filtmp,action,isample(n))
               write(iout,9000) 'PiG Sampling Grid output file:'
               write(iout,9002) filtmp(1:ii+9)
             endif
          enddo
        endif
      endif
c
c-----Open chemistry parameters input file
c
      filtmp = Chemistry_Parameters
      action = 'Opening Chemistry Parameters file.'
      inquire(file=filtmp,exist=lexist)
      if( .NOT. lexist ) goto 7002
      nopen = nopen + 1
      call getunit(ichem)
      open(unit=ichem,file=filtmp,status='OLD',ERR=7000)
      write(iout,9000) 'Chemistry Parameters file            (unit):',
     &                                                           ichem
      write(iout,9002) filtmp(:istrln(filtmp))
      nfils = nfils + 1
c
c-----Open photoloysis rates input file
c
      if( lchem ) then
        filtmp = Photolyis_Rates
        action = 'Opening Photolysis Rates file.'
        inquire(file=filtmp,exist=lexist)
        if( .NOT. lexist .OR. filtmp .EQ. ' ' ) then
          filtmp = Photolysis_Rates
          inquire(file=filtmp,exist=lexist)
          if( .NOT. lexist ) goto 7002
        endif
        nopen = nopen + 1
        call getunit(iphot)
        open(unit=iphot,file=filtmp,status='OLD',ERR=7000) 
        write(iout,9000) 'Photolysis Rates file                (unit):',
     &                                                             iphot
        write(iout,9002) filtmp(:istrln(filtmp))
        nfils = nfils + 1
      else
        iphot = 0
        write(iout,9000) 'Photolysis Rates file                      :'
        write(iout,9002) '   Ignored.' 
      endif
c
c-----Open initial conditions input file
c
      if( .NOT. lrstrt ) then
        filtmp = Initial_Conditions
        action = 'Opening Initial Conditions file.'
        is_netcdf_iic = .FALSE.
        inquire(file=filtmp,exist=lexist)
        if( .NOT. lexist ) goto 7002
        nopen = nopen + 1
        call getunit(iic)
        cname = 'AIRQUALITY'
        ierr = ncf_chkfile(iic,filtmp,action,cname)
        if( ierr .EQ. ISUCES ) then
           open(unit=iic,file=filtmp,form='UNFORMATTED',
     &                                          status='OLD',ERR=7000)
        else
           is_netcdf_iic = .TRUE.
           ierr = nf_open(filtmp,NF_NOWRITE,iic)
           if( ierr .NE. NF_NOERR ) goto 7000
        endif
        write(iout,9000) 'Initial Conditions file              (unit):',
     &                                                               iic
        write(iout,9002) filtmp(:istrln(filtmp))
        nfils = nfils + 1
      else
        iic = 0
        write(iout,9000) 'Initial Conditions file                    :'
        write(iout,9002) '   Ignored.'
      endif
c
c-----Open boundary conditions input file
c
      filtmp = Boundary_Conditions
      action = 'Opening Boundary Conditions file.'
      is_netcdf_ibc = .FALSE.
      inquire(file=filtmp,exist=lexist)
      if( .NOT. lexist ) goto 7002
      nopen = nopen + 1
      call getunit(ibc)
      cname = 'BOUNDARY  '
      ierr = ncf_chkfile(ibc,filtmp,action,cname)
      if( ierr .EQ. ISUCES ) then
           open(unit=ibc,file=filtmp,form='UNFORMATTED',
     &                                          status='OLD',ERR=7000)
      else
         is_netcdf_ibc = .TRUE.
         ierr = nf_open(filtmp,NF_NOWRITE,ibc)
         if( ierr .NE. NF_NOERR ) goto 7000
      endif
      write(iout,9000) 'Boundary Conditions file             (unit):',
     &                                                             ibc
      write(iout,9002) filtmp(:istrln(filtmp))
      nfils = nfils + 1
c
c-----Open top concentration input file
c
      filtmp = Top_Concentrations
      if( filtmp .NE. ' ' ) then
        action = 'Opening Top Concentrations file.'
        is_netcdf_itc = .FALSE.
        inquire(file=filtmp,exist=lexist)
        if( .NOT. lexist ) goto 7002
        nopen = nopen + 1
        call getunit(itc)
        cname = 'AIRQUALITY'
        ierr = ncf_chkfile(itc,filtmp,action,cname)
        if( ierr .EQ. ISUCES ) then
           open(unit=itc,file=filtmp,form='UNFORMATTED',
     &                                          status='OLD',ERR=7000)
        else
           is_netcdf_itc = .TRUE.
           ierr = nf_open(filtmp,NF_NOWRITE,itc)
           if( ierr .NE. NF_NOERR ) goto 7000
        endif
        write(iout,9000) 'Top Concentrations file              (unit):',
     &                                                             itc
        write(iout,9002) filtmp(:istrln(filtmp))
        ltopcon = .true.
        nfils = nfils + 1
      else
        write(iout,9000) 'Top Concentrations file                    :'
        write(iout,9002) '   Ignored.'
        ltopcon = .false.
        if( lstrato3 ) goto 7011
      endif
c
c-----Open ozone column input file
c
      filtmp = Ozone_Column
      if( filtmp .NE. ' ' ) then
        action = 'Opening ozone column file.'
        inquire(file=filtmp,exist=lexist)
        if( .NOT. lexist ) goto 7002
        nopen = nopen + 1
        call getunit(io3col)
        open(unit=io3col,file=filtmp,status='OLD',ERR=7000)
        write(iout,9000) 'Ozone column file                    (unit):',
     &                                                            io3col
        write(iout,9002) filtmp(:istrln(filtmp))
        nfils = nfils + 1
      else
        io3col = 0
        if( lchem ) goto 7010
        write(iout,9000) 'Ozone column file                          :'
        write(iout,9002) '   Ignored.'
      endif
c
c-----Open point source emissions input file
c
      luse_points = .FALSE.
      if( lptsrc ) then
        do ifile=1,npoint_files
           filtmp = Point_Sources(ifile)
           action = 'Opening Point Source Emissions file.'
           if( filtmp .EQ. ' ' ) goto 7008
           inquire(file=filtmp,exist=lexist)
           if( .NOT. lexist ) goto 7002
           nopen = nopen + 1
           call getunit(iptem(ifile))
           cname = 'PTSOURCE  '
           ierr = ncf_chkfile(iptem(ifile),filtmp,action,cname)
           if( ierr .EQ. ISUCES ) then
              is_netcdf_iptem(ifile) = .FALSE.
              open(unit=iptem(ifile),file=filtmp,
     &                        form='UNFORMATTED',status='OLD',ERR=7000)
           else
              is_netcdf_iptem(ifile) = .TRUE.
              ierr = nf_open(filtmp,NF_NOWRITE,iptem(ifile))
              if( ierr .NE. NF_NOERR ) goto 7000
           endif
           write(iout,9000) 'Point Source Emissions file          (unit):',
     &                                                           iptem(ifile)
           write(iout,9002) filtmp(:istrln(filtmp))
           nfils = nfils + 1
           luse_points = .TRUE.
        enddo
      else
        npoint_files = 0
        iptem(1) = 0
        write(iout,9000) 'Point Source Emissions file                :'
        write(iout,9002) '   Ignored.'
      endif
c
c-----if doing SA then point source files are read from SA_control section ----
c
      if( tectyp .EQ. SA .AND. luse_points ) then
         action = 'Reading point source filenames.'
         goto 7006
      endif
c
c-----Open master grid restart input file
c
      if( lrstrt ) then
        filtmp = Master_Grid_Restart
        action = 'Opening Restart file for master grid.'
        inquire(file=filtmp,exist=lexist)
        if( .NOT. lexist ) goto 7002
        nopen = nopen + 1
        call getunit(irstc)
        open(unit=irstc,file=filtmp,form='UNFORMATTED',
     &                                       status='OLD',ERR=7005)
        write(iout,9000) 'Master grid Restart file             (unit):',
     &                                                             irstc
        write(iout,9002) filtmp(:istrln(filtmp))
        nfils = nfils + 1
c
c-----Open fine grid restart input file
c
        if( ngrid .GT. 1 ) then
          filtmp = Nested_Grid_Restart
          action = 'Opening Restart file for nest grids.'
          if( filtmp .EQ. ' ' ) then
             if (lflexi) then
               irstf = 0
               write(iout,9001)
     &                  'Nest grid Restart file                     :'
               write(iout,9002) '   Ignored.'
             else
               goto 7003
             endif
          else
             inquire(file=filtmp,exist=lexist)
             if( .NOT. lexist ) goto 7003
             nopen = nopen + 1
             call getunit(irstf)
             open(unit=irstf,file=filtmp,form='UNFORMATTED',
     &                                         status='OLD',ERR=7005)
             write(iout,9000)
     &              'Nest grid Restart file               (unit):',irstf
             write(iout,9002) filtmp(:istrln(filtmp))
             nfils = nfils + 1
          endif
        endif
c
c-----Open PiG restart input file
c
        if( ipigflg .NE. 0 ) then
          filtmp = PiG_Restart
          if( filtmp .EQ. ' ' ) then
             irstp = 0
             write(iout,9001)
     &                   'PiG Restart file                           :'
             write(iout,9002) '   Ignored.'
          else
             action = 'Opening PiG Restart file.'
             inquire(file=filtmp,exist=lexist)
             if( .NOT. lexist ) goto 7002
             nopen = nopen + 1
             call getunit(irstp)
             open(unit=irstp,file=filtmp,form='UNFORMATTED',
     &                                status='OLD',ERR=7005)
             write(iout,9000)
     &         'PiG Restart file                     (unit):',irstp
             write(iout,9002) filtmp(:istrln(filtmp))
             nfils = nfils + 1
          endif
        endif
      endif
c
c-----Open master grid 2D Surface input file
c
      filtmp = Surface_Grid(1)
      write(action,'(A,I4)') 'Opening 2D Surface file for grid:',1
      is_netcdf_isurf(1) = .FALSE.
      inquire(file=filtmp,exist=lexist)
      if( .NOT. lexist ) goto 7002
      nopen = nopen + 1
      call getunit(isurf(1))
      cname = 'AVERAGE   '
      ierr = ncf_chkfile(isurf(1),filtmp,action,cname)
      if( ierr .EQ. ISUCES ) then
         open(unit=isurf(1),file=filtmp,form='UNFORMATTED',
     &                                          status='OLD',ERR=7000)
      else
         is_netcdf_isurf(1) = .TRUE.
         ierr = nf_open(filtmp,NF_NOWRITE,isurf(1))
         if( ierr .NE. NF_NOERR ) goto 7000
      endif
      write(iout,9001)
     &        '2D Surface file for grid # ',1,'        (unit):',isurf(1)
      write(iout,9002) filtmp(:istrln(filtmp))
      nfils = nfils + 1
c
c-----Open nested grid 2D Surface input file(s)
c
      do n = 2,ngrid
        filtmp = Surface_Grid(n)
        write(action,'(A,I4)') 'Opening 2D Surface file for grid:',n
        is_netcdf_isurf(n) = .FALSE.
        if( filtmp .EQ. ' ')  then
           if (lflexi) then
             isurf(n) = 0
             write(iout,9001)
     &                 '2D Surface file for grid # ',n,'              :'
             write(iout,9002) '   Ignored.'
           else
             goto 7003
           endif
        else
           inquire(file=filtmp,exist=lexist)
           if( .NOT. lexist ) goto 7003
           nopen = nopen + 1
           call getunit(isurf(n))
           cname = 'AVERAGE   '
           ierr = ncf_chkfile(isurf(n),filtmp,action,cname)
           if( ierr .EQ. ISUCES ) then
              open(unit=isurf(n),file=filtmp,form='UNFORMATTED',
     &                                          status='OLD',ERR=7000)
           else
              is_netcdf_isurf(n) = .TRUE.
              ierr = nf_open(filtmp,NF_NOWRITE,isurf(n))
              if( ierr .NE. NF_NOERR ) goto 7000
           endif
           write(iout,9001)
     &             '2D Surface file for grid # ',n,'        (unit):',
     &                                                         isurf(n)
           write(iout,9002) filtmp(:istrln(filtmp))
           nfils = nfils + 1
        endif
      enddo
c
c-----Open master grid 3D Met input file
c
      filtmp = Met3d_Grid(1)
      is_netcdf_i3dmet(1) = .FALSE.
      write(action,'(A,I4)') 'Opening 3D Met file for grid:',1
      inquire(file=filtmp,exist=lexist)
      if( .NOT. lexist ) goto 7002
      nopen = nopen + 1
      call getunit(i3dmet(1))
      cname = 'AVERAGE   '
      ierr = ncf_chkfile(i3dmet(1),filtmp,action,cname)
      if( ierr .EQ. ISUCES ) then
         open(unit=i3dmet(1),file=filtmp,form='UNFORMATTED',
     &                                          status='OLD',ERR=7000)
      else
         is_netcdf_i3dmet(1) = .TRUE.
         ierr = nf_open(filtmp,NF_NOWRITE,i3dmet(1))
         if( ierr .NE. NF_NOERR ) goto 7000
      endif
      write(iout,9001)'3D Met file for grid # ',1,'            (unit):',
     &                                                         i3dmet(1)
      write(iout,9002) filtmp(:istrln(filtmp))
      nfils = nfils + 1
c
c-----Open nested grid 3D Met input file(s)
c
      do n = 2,ngrid
        filtmp = Met3D_Grid(n)
        write(action,'(A,I4)') 'Opening 3D Met file for grid:',n
        is_netcdf_i3dmet(n) = .FALSE.
        if( filtmp .EQ. ' ' ) then
          if (lflexi) then
            i3dmet(n) = 0
            write(iout,9001)
     &                 '3D Met file for grid # ',n,'                  :'
            write(iout,9002) '   Ignored.'
          else
            goto 7003
          endif
        else
          inquire(file=filtmp,exist=lexist)
          if( .NOT. lexist ) goto 7003
          nopen = nopen + 1
          call getunit(i3dmet(n))
          cname = 'AVERAGE   '
          ierr = ncf_chkfile(i3dmet(n),filtmp,action,cname)
          if( ierr .EQ. ISUCES ) then
             open(unit=i3dmet(n),file=filtmp,form='UNFORMATTED',
     &                                          status='OLD',ERR=7000)
          else
             is_netcdf_i3dmet(n) = .TRUE.
             ierr = nf_open(filtmp,NF_NOWRITE,i3dmet(n))
             if( ierr .NE. NF_NOERR ) goto 7000
          endif
          write(iout,9001)
     &                '3D Met file for grid # ',n,'            (unit):',
     &                                                         i3dmet(n)
          write(iout,9002) filtmp(:istrln(filtmp))
          nfils = nfils + 1
        endif
      enddo
c
c-----Open master grid 2D Met input file
c
      filtmp = Met2D_Grid(1)
      write(action,'(A,I4)') 'Opening 2D Met file for grid:',1
      is_netcdf_i2dmet(1) = .FALSE.
      inquire(file=filtmp,exist=lexist)
      if( .NOT. lexist ) goto 7002
      nopen = nopen + 1
      call getunit(i2dmet(1))
      cname = 'AVERAGE   '
      ierr = ncf_chkfile(i2dmet(1),filtmp,action,cname)
      if( ierr .EQ. ISUCES ) then
           open(unit=i2dmet(1),file=filtmp,form='UNFORMATTED',
     &                                          status='OLD',ERR=7000)
      else
         is_netcdf_i2dmet(1) = .TRUE.
         ierr = nf_open(filtmp,NF_NOWRITE,i2dmet(1))
         if( ierr .NE. NF_NOERR ) goto 7000
      endif
      write(iout,9001)
     &              '2D Met file for grid # ',1,'            (unit):',
     &                                                       i2dmet(1)
      write(iout,9002) filtmp(:istrln(filtmp))
      nfils = nfils + 1
c
c-----Open nested grid 2D Met input file(s)
c
      do n = 2,ngrid
        filtmp = Met2D_Grid(n)
        write(action,'(A,I4)') 'Opening 2D Met file for grid:',n
        is_netcdf_i2dmet(n) = .FALSE.
        if( filtmp .EQ. ' ' ) then
           if (lflexi) then
             i2dmet(n) = 0
             write(iout,9001)
     &                 '2D Met file for grid # ',n,'                  :'
             write(iout,9002) '   Ignored.'
           else
             goto 7003
           endif
         else
           inquire(file=filtmp,exist=lexist)
           if( .NOT. lexist ) goto 7003
           nopen = nopen + 1
           call getunit(i2dmet(n))
           cname = 'AVERAGE   '
           ierr = ncf_chkfile(i2dmet(n),filtmp,action,cname)
           if( ierr .EQ. ISUCES ) then
              open(unit=i2dmet(n),file=filtmp,form='UNFORMATTED',
     &                                          status='OLD',ERR=7000)
           else
              is_netcdf_i2dmet(n) = .TRUE.
              ierr = nf_open(filtmp,NF_NOWRITE,i2dmet(n))
              if( ierr .NE. NF_NOERR ) goto 7000
           endif
           write(iout,9001)
     &               '2D Met file for grid # ',n,'            (unit):',
     &                                                        i2dmet(n)
           write(iout,9002) filtmp(:istrln(filtmp))
           nfils = nfils + 1
         endif
      enddo
c
c-----Open master grid 3D VDiff input file
c
      filtmp = Vdiff_Grid(1)
      write(action,'(A,I4)')
     &                'Opening 3D Vertical Diffusivity file for grid:',1
      is_netcdf_ikv(1) = .FALSE.
      inquire(file=filtmp,exist=lexist)
      if( .NOT. lexist ) goto 7002
      nopen = nopen + 1
      call getunit(ikv(1))
      cname = 'AVERAGE   '
      ierr = ncf_chkfile(ikv(1),filtmp,action,cname)
      if( ierr .EQ. ISUCES ) then
         open(unit=ikv(1),file=filtmp,form='UNFORMATTED',
     &                                          status='OLD',ERR=7000)
      else
         is_netcdf_ikv(1) = .TRUE.
         ierr = nf_open(filtmp,NF_NOWRITE,ikv(1))
         if( ierr .NE. NF_NOERR ) goto 7000
      endif
      write(iout,9001)
     &             '3D VDiff file for grid # ',1,'             (unit):',
     &                                                          ikv(1)
      write(iout,9002) filtmp(:istrln(filtmp))
      nfils = nfils + 1
c
c-----Open nested grid 3D VDiff input file(s)
c
      do n = 2,ngrid
        filtmp = Vdiff_Grid(n)
        write(action,'(A,I4)')
     &                'Opening 3D Vertical Diffusivity file for grid:',n
        is_netcdf_ikv(n) = .FALSE.
        if( filtmp .EQ. ' ' ) then
           if (lflexi) then
             ikv(n) = 0
             write(iout,9001)
     &              '3D VDiff file for grid # ',n,'                   :'
             write(iout,9002) '   Ignored.'
           else
             goto 7003
           endif
        else
           inquire(file=filtmp,exist=lexist)
           if( .NOT. lexist ) goto 7003
           nopen = nopen + 1
           call getunit(ikv(n))
           cname = 'AVERAGE   '
           ierr = ncf_chkfile(ikv(n),filtmp,action,cname)
           if( ierr .EQ. ISUCES ) then
              open(unit=ikv(n),file=filtmp,form='UNFORMATTED',
     &                                          status='OLD',ERR=7000)
           else
              is_netcdf_ikv(n) = .TRUE.
              ierr = nf_open(filtmp,NF_NOWRITE,ikv(n))
              if( ierr .NE. NF_NOERR ) goto 7000
           endif
           write(iout,9001)
     &             '3D VDiff file for grid # ',n,'             (unit):',
     &                                                           ikv(n)
           write(iout,9002) filtmp(:istrln(filtmp))
           nfils = nfils + 1
        endif
      enddo
c
c-----Open master grid 3D Cloud/rain input file
c
      icld = 0
      if(  .NOT. is_netcdf_i3dmet(1) ) then
         filtmp = Cloud_Grid(1)
         write(action,'(A,I4)') 'Opening 3D Cloud/Rain file for grid:',1
         if( filtmp .EQ. ' ' ) then
            if( lwet .or. lcig(1) ) goto 7002
            icld(1) = 0
            write(iout,9001)
     &                 '3D Cloud/Rain file for grid # ',1,'           :'
            write(iout,9002) '   Ignored.' 
         else
            inquire(file=filtmp,exist=lexist)
            if( .NOT. lexist ) goto 7002
            nopen = nopen + 1
            call getunit(icld(1))
            open(unit=icld(1),file=filtmp,form='UNFORMATTED',
     &                                            status='OLD',ERR=7000)
            write(iout,9001)
     &               '3D Cloud/Rain file for grid # ',1,'     (unit):',
     &                                                          icld(1)
            write(iout,9002) filtmp(:istrln(filtmp))
            nfils = nfils + 1
c
c-----Open nested grid 3D Cloud/rain input file(s)
c
            do n = 2,ngrid
              if( is_netcdf_i3dmet(n) ) cycle
              filtmp = Cloud_Grid(n)
              write(action,'(A,I4)')
     &                          'Opening 3D Cloud/Rain file for grid:',n
              if( filtmp .EQ. ' ' ) then
                 if (lcig(n)) goto 7002
                 if (lflexi .or. .not.lwet) then
                   icld(n) = 0
                   write(iout,9001)
     &                 '3D Cloud/Rain file for grid # ',n,'           :'
                   write(iout,9002) '   Ignored.'
                 else
                   goto 7003
                 endif
              else
                 inquire(file=filtmp,exist=lexist)
                 if( .NOT. lexist ) goto 7003
                 nopen = nopen + 1
                 call getunit(icld(n))
                 open(unit=icld(n),file=filtmp,form='UNFORMATTED',
     &                                           status='OLD',ERR=7000)
                 write(iout,9001)
     &               '3D Cloud/Rain file for grid # ',n,'     (unit):',
     &                                                          icld(n)
                 write(iout,9002) filtmp(:istrln(filtmp))
                 nfils = nfils + 1
c
              endif
            enddo
         endif
       endif
c
c-----Open master grid 2D Emissions input file
c
      if( larsrc ) then
         do ifile=1,nemiss_files(1)
            filtmp = Emiss_Grid(1,ifile)
            write(action,'(A,I4)') 
     &                     'Opening Gridded Emissions file for grid:',1
            if( filtmp .EQ. ' ' ) goto 7007
            inquire(file=filtmp,exist=lexist)
            if( .NOT. lexist ) goto 7002
            nopen = nopen + 1
            call getunit(iarem(1,ifile))
            cname = 'EMISSIONS '
            ierr = ncf_chkfile(iarem(1,ifile),filtmp,action,cname)
            if( ierr .EQ. ISUCES ) then
                 is_netcdf_iarem(1,ifile) = .FALSE.
                 open(unit=iarem(1,ifile),file=filtmp,form='UNFORMATTED',
     &                                             status='OLD',ERR=7000)
            else
               is_netcdf_iarem(1,ifile) = .TRUE.
               ierr = nf_open(filtmp,NF_NOWRITE,iarem(1,ifile))
               if( ierr .NE. NF_NOERR ) goto 7000
               nlayers_ems = MAX(nlayers_ems,ncf_get_nlayers(iarem(1,ifile),filtmp,action))
            endif
            write(iout,9001)
     &              'Gridded Emissions file for grid # ',1,' (unit):',
     &                                                         iarem(1,ifile)
            write(iout,9002) filtmp(:istrln(filtmp))
            nfils = nfils + 1
         enddo
c
c-----Open nested grid 2D Emissions input file(s)
c
         do n = 2,ngrid
           do ifile=1,nemiss_files(n)
             filtmp = Emiss_Grid(n,ifile)
             write(action,'(A,I4)') 
     &                     'Opening Gridded Emissions file for grid:',n
             if( filtmp .EQ. ' ' ) then
                if (lflexi) then
                  iarem(n,ifile) = 0
                  write(iout,9001)
     &                'Gridded Emissions file for grid # ',n,'       :'
                  write(iout,9002) '   Ignored.' 
                else
                  goto 7003
                endif
             else
                inquire(file=filtmp,exist=lexist)
                if( .NOT. lexist ) goto 7003
                nopen = nopen + 1
                call getunit(iarem(n,ifile))
                cname = 'EMISSIONS '
                ierr = ncf_chkfile(iarem(n,ifile),filtmp,action,cname)
                if( ierr .EQ. ISUCES ) then
                      is_netcdf_iarem(n,ifile) = .FALSE.
                      open(unit=iarem(n,ifile),file=filtmp,
     &                        form='UNFORMATTED',status='OLD',ERR=7000)
                else
                   is_netcdf_iarem(n,ifile) = .TRUE.
                   ierr = nf_open(filtmp,NF_NOWRITE,iarem(n,ifile))
                   if( ierr .NE. NF_NOERR ) goto 7000
                   nlayers_ems = MAX(nlayers_ems,ncf_get_nlayers(iarem(n,ifile),filtmp,action))
                endif
                write(iout,9001)
     &                'Gridded Emissions file for grid # ',n,' (unit):',
     &                                                          iarem(n,ifile)
                write(iout,9002) filtmp(:istrln(filtmp))
                nfils = nfils + 1
             endif
           enddo
         enddo
      else
        do n = 1,ngrid
          nemiss_files(n) = 0
          iarem(n,1) = 0
          write(iout,9001)
     &              'Gridded Emissions file for grid # ',n,'       :'
          write(iout,9002) '   Ignored.' 
        enddo
      endif
c
c-----Open input surface model mass files for every grid 
c
      do n = 1,ngrid
        ismin(n) = 0
      enddo
      if (lsrfmod .and. lrstrt) then
        do n = 1,ngrid
          filtmp = Srfmod_Grid(n)
          write(action,'(A,I4)') 
     &                  'Opening Restart Surface Model file for grid:',n
          if( filtmp .EQ. ' ' ) then
             write(iout,9001)
     &             'Restart Surface Model file for grid # ',n,'       :'
             write(iout,9002) ' Not supplied'
          else
             inquire(file=filtmp,exist=lexist)
             if( .NOT. lexist ) goto 7003
             nopen = nopen + 1
             call getunit(ismin(n))
             open(unit=ismin(n),file=filtmp,form='UNFORMATTED',
     &                                          status='OLD',ERR=7000)
             write(iout,9001)
     &            'Restart Surface Model file for grid # ',n,' (unit):',
     &                                                          ismin(n)
             write(iout,9002) filtmp(:istrln(filtmp))
             nfils = nfils + 1
          endif
        enddo
      endif
c
      goto 9999
c
c  --- Format statements ---
c
 9000 format(/,A,I8)
 9001 format(/,A,I2,A,I8)
 9002 format(2A)
 9003 format(3A)
c
c  --- Error messages
c
 7000 continue
      write(iout,'(//,a)') 'ERROR in OPENFILS:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Could not open file: ',
     &                                   filtmp(:istrln(filtmp))
      write(iout,'(3A)') 'If this is a NCF file, check that its',
     &                   ' format is consistent with the NCF library',
     &                   ' used to build CAMx.'
      write(iout,'(2A)') 'The CAMx makefile supports netCDF3-Classic,',
     &                   ' and netCDF4 compressed or uncompressed'
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7002 continue
      write(iout,'(//,a)') 'ERROR in OPENFILS:'
      write(iout,'(A)') action(:istrln(action))
      if( filtmp .EQ. ' ' ) then
         write(iout,'(A)') 'Blank filename provided in control file'
      else
         write(iout,'(2A)') 'Input file does not exist: ',
     &                                       filtmp(:istrln(filtmp))
      endif
      call camxerr()
c
 7003 continue
      write(iout,'(//,a)') 'ERROR in OPENFILS:'
      write(iout,'(A)') action(:istrln(action))
      if( filtmp .EQ. ' ' ) then
         write(iout,'(A)') 'Blank filename provided in control file'
         write(iout,'(2A)') 'If this file is for a nested grid, set',
     &                      ' the namelist variable '
         write(iout,'(5X,A,/,A)')  'Flexi_Nest = .true.',
     &                             'to ignore this file.'
      else
         write(iout,'(2A)') 'Input file does not exist: ',
     &                                       filtmp(:istrln(filtmp))
      endif
      call camxerr()
c
 7005 continue
      write(iout,'(//,a)') 'ERROR in OPENFILS:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Could not open file: ',
     &                                   filtmp(:istrln(filtmp))
      write(iout,'(10X,2A)') 'Make sure the names of restart files ',
     &                       'are for the previous simulation period.'
      call camxerr()
c
 7006 continue
      write(iout,'(//,a)') 'ERROR in OPENFILS:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'In this version the point source inventory ',
     &           'for a simulation using the SA '
      write(iout,'(2A)') 'probing tool is provided in the SA_Control ',
     &            'section of the CAMx control file.'
      write(iout,'(2A)') 'You have supplied filenames in the ',
     &                                         'CAMx_control section.'
      write(iout,'(2A)') 'Replace these filenames with blank strings ',
     &        'and include the point source files in the SA_Control section.'
      write(iout,'(2A)') 'Replace these filenames with blank strings ',
     &                                'and include the point source '
      write(iout,'(2A)') 'files in the SA_Control section.'
      call camxerr()
c
 7007 continue
      write(iout,'(//,A)') 'ERROR in OPENFILS:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(1X,2A)') 'You supplied a blank filename in the middle ',
     &                  'of the list of gridded emissions filenames.'
      write(iout,'(1X,2A)') 'Check the numbering of the list of files.'
      call camxerr()
c
 7008 continue
      write(iout,'(//,A)') 'ERROR in OPENFILS:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(1X,2A)') 'You supplied a blank filename in the middle ',
     &                  'of the list of point source filenames.'
      write(iout,'(1X,2A)') 'Check the numbering of the list of files.'
      call camxerr()
c
 7010 continue
      write(iout,'(//,a)') 'ERROR in OPENFILS:'
      write(iout,'(2A)') 'The chemistry flag is set but no ',
     &                         'Ozone Column file is supplied.'
      write(iout,'(2A)')'Either supply an Ozone Column file ',
     &                         'or set the Chemistry flag to false.'
      write(iout,'(2A)')'If using Mechanism 10 just supply ',
     &                  'a dummy file for the Ozone Column file.'
      call camxerr()
c
 7011 continue
      write(iout,'(//,a)') 'ERROR in OPENFILS:'
      write(iout,'(2A)') 'The Stratospheric Ozone treatment is selected',
     &                   ' but no Top Concentration file is supplied.'
      write(iout,'(2A)')'Either supply a Top Concentration file',
     &                  ' containing Ozone'
      write(iout,'(A)') 'or set the Stratospheric Ozone flag to false.'
      call camxerr()
c
 9999 continue
      return
      end
