c**** STARTDDM.F
c
      subroutine startddm(iounit,nopen)
      use grid
      use chmstry
      use filunit
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
c     07/19/02   --gwilson--   Added code for source area map for nests
c     10/06/04   --cemery --   Restructured for namelist input
c      8/23/06   --cemery--    Instantaneous restart files reduced to 1
c                              per grid type
c      8/25/06   --cemery--    Surface output files now all UAM format,
c                              one file per grid
c     07/16/07   --bkoo--      Added check for HDDM
c     11/4/09    --cemery--    Removed input top concentrations
c     05/07/12   --cemery--    Added flexi-nesting flag
c     03/01/16   --gwilson--   Added partial source area map
c     11/28/16   --bkoo--      Added code for new top con
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
      integer      i, j, igrd, n, idxfile, ierr
      logical      lexist
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c   --- if option is turned off, just return ---
c
      if( .NOT. lddm .AND. .NOT. lhddm ) goto 9999
      nchar = istrln( flrtsa )
c
c  -- call routine to allocate arrays ----
c
      call alloc_tracer(ngroup,ngrid,ncol,nrow,nlay,nspec)
      call alloc_tracer_full(ngroup,ngrid,ncol,nrow)
c
c   --- intialize first time through to true ---
c
      lfirst = .TRUE.
      write(iounit,'(/,A,/)') 
     &            '           **** DDM Files ****'
c
c   --- source area mapping file ---
c
      do 10 igrd = 1,ngrid
         write(action,'(A,I2)') 
     &   'Opening Source Area Map for grid: ',igrd
                                mapfil(0,igrd) = DDM_Source_Area_Map(igrd)
         fname = ' '
         fname = mapfil(0,igrd)
         if( nregin .LE. 0 ) then
            lmapfl(0,igrd) = .FALSE.
            goto 10
         endif 
         if( fname .EQ. ' ' ) then
            if( igrd .EQ. 1 ) goto 7004
            lmapfl(0,igrd) = .FALSE.
            write(iounit,9001)
     &              'No Source Area Map file for Grid #         :',igrd
            goto 10
         endif
         inquire(file=mapfil(0,igrd),exist=lexist)
         if( .NOT. lexist ) goto 7000
         nopen =  nopen + 1
         call getunit(iormap(0,igrd))
         write(iounit,9001)  'Source area map file for grid # ',igrd,
     &                                        '   (unit):',iormap(0,igrd)
         write(iounit,9002) mapfil(0,igrd)(:istrln(mapfil(0,igrd)))
         lmapfl(0,igrd) = .TRUE.
   10 continue
c
c   --- receptor definition file ----
c
      rcpfil = ' '
      rcpfil = DDM_Receptor_Definitions
      if( rcpfil .NE. ' ' ) then
         action = 'Opening DDM Receptor Definition file.'
         fname = rcpfil
         inquire(file=rcpfil,exist=lexist)
         if( .NOT. lexist ) goto 7000
         nopen =  nopen + 1
         call getunit(iorrcp)
         write(iounit,9000) 
     &       'DDM Receptor definition file         (unit):',iorrcp
         write(iounit,9002) rcpfil(:istrln(rcpfil))
         lrcpfil = .TRUE.
      else
         write(iounit,9000) 'No DDM Receptor definition file provided.'
         lrcpfil = .FALSE.
      endif
c
c   --- instantaneous file used for initialization ---
c
      if( lrstrt ) then
         action = 'Opening DDM master grid Restart file.'
         inifil(IDXCRS) = ' '
         inifil(IDXCRS) = DDM_Master_Restart
         fname = inifil(IDXCRS)
         inquire(file=inifil(IDXCRS),exist=lexist)
         if( .NOT. lexist ) goto 7000
         nopen =  nopen + 1
         call getunit(iorini(IDXCRS))
         write(iounit,9000)'DDM master grid Restart file         (unit):',
     &                                                    iorini(IDXCRS)
         write(iounit,9002) inifil(IDXCRS)(:istrln(inifil(IDXCRS)))
c
         if( ngrid .GT. 1 ) then
            action = 'Opening DDM nested grid Restart file.'
            inifil(IDXFIN) = ' '
            inifil(IDXFIN) = DDM_Nested_Restart
            if( inifil(IDXFIN) .NE. ' ' ) then
               fname = inifil(IDXFIN)
               inquire(file=inifil(IDXFIN),exist=lexist)
               if( .NOT. lexist ) goto 7000
               nopen =  nopen + 1
               call getunit(iorini(IDXFIN))
               write(iounit,9000)
     &                'DDM nested grid Restart file         (unit):',
     &                                                    iorini(IDXFIN)
               write(iounit,9002) inifil(IDXFIN)(:istrln(inifil(IDXFIN)))
            elseif (.not.lflexi) then
               goto 7003
            endif
         endif
      endif
c
c   --- IC file for DDM ---
c
      if( nicddm .GT. 0 .AND. .NOT. lrstrt ) then
         action = 'Opening DDM Initial Conditions file.'
         icfil = ' '
         icfil = DDM_Initial_Conditions
         fname = icfil
         inquire(file=icfil,exist=lexist)
         if( .NOT. lexist ) goto 7000
         nopen =  nopen + 1
         call getunit(ioric)
         write(iounit,9000)'DDM Initial Conditions file          (unit):',
     &                                                             ioric
         write(iounit,9002) icfil(:istrln(icfil))
      endif
c
c   --- BC file for DDM ---
c
      if( nbcddm .GT. 0 ) then
         action = 'Opening DDM Boundary Conditions file.'
         bcfil = ' '
         bcfil = DDM_Boundary_Conditions
         fname = bcfil
         inquire(file=bcfil,exist=lexist)
         if( .NOT. lexist ) goto 7000
         nopen =  nopen + 1
         call getunit(iorbc)
         write(iounit,9000)'DDM Boundary Conditions file         (unit):',
     &                                                             iorbc
         write(iounit,9002) bcfil(:istrln(bcfil))

         if ( ltopcon ) then
           action = 'Opening DDM Top Concentrations file.'
           tcfil = ' '
           tcfil = DDM_Top_Concentrations
           fname = tcfil
           if( fname .EQ. ' ' ) goto 7007
           inquire(file=tcfil,exist=lexist)
           if( .NOT. lexist ) goto 7000
           nopen =  nopen + 1
           call getunit(iortc)
           write(iounit,9000)'DDM Top Concentrations file          (unit):',
     &                                                             iortc
           write(iounit,9002) tcfil(:istrln(tcfil))
         endif
      endif
c
c   --- figure out how many files in each grid/group ---
c
      num_iortem = 0
      do i=1,ngroup
         do j=1,ngrid
            do idxfile=1,MXFILES
              if( istrln(DDM_Emiss_Group_Grid(i,j,idxfile)) .NE. 0 )
     &                                      num_iortem(j,i) = idxfile
            enddo
         enddo
      enddo
      num_iortpt = 0
      do i=1,ngroup
        do idxfile=1,MXFILES
          if( istrln(DDM_Points_Group(i,idxfile)) .NE. 0 )
     &                                      num_iortpt(i) = idxfile
        enddo
      enddo
      call alloc_tracer_emfiles(ngroup,ngrid,nemiss_files,npoint_files,
     &                                      num_iortem,num_iortpt,nspec)
c
c   --- emissions files for the source groupings,
c       skip if not doing any emissions DDM species ---
c
      if( nemddm .EQ. 0 ) goto 111
      nlayers_ems = MAX(nlayers_ems,1)
      do 20 i=1,ngroup
c
c   --- surface emissions file --- 
c
         action = 'Opening DDM Gridded Emissions file.'
         do 30 j=1,ngrid
          if( num_iortem(j,i) .EQ. 0 ) then
               write(iounit,9003)
     &         'Gridded emissions file for grid#/group#    :',
     &                                     j,i,' Not supplied'
          else
           do idxfile=1,num_iortem(j,i)
             temfil(j,i,idxfile) = ' '
             temfil(j,i,idxfile) = DDM_Emiss_Group_Grid(i,j,idxfile)
             fname = temfil(j,i,idxfile)
             ltemfl(j,i,idxfile) = .TRUE.
             if( fname .EQ. ' ' ) goto 7005
             inquire(file=temfil(j,i,idxfile),exist=lexist)
             if( .NOT. lexist ) goto 7000
             nopen = nopen + 1
             cname = 'EMISSIONS '
             call getunit(iortem(j,i,idxfile))
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
     &         'Gridded emissions for grid#/group#   (unit):',
     &                                                 j,i,iortem(j,i,idxfile)
             write(iounit,9002) temfil(j,i,idxfile)(:istrln(temfil(j,i,idxfile)))
           enddo
          endif
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
                tptfil(i,idxfile) = ' '
                tptfil(i,idxfile) = DDM_Points_Group(i,idxfile)
                fname = tptfil(i,idxfile)
                action = 'Opening DDM Point Source file.'
                if( fname .EQ. ' ' ) goto 7006
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
     &                                                      i,iortpt(i,idxfile)
                write(iounit,9002) tptfil(i,idxfile)(:istrln(tptfil(i,idxfile)))
             enddo
          endif
   20 continue
c
c   --- output filenames ---
c
  111 continue
      do 40 i=IDXCRS,IDXFIN
         if( i .EQ. IDXFIN .AND. ngrid .EQ. 1 ) goto 40
c
c   --- output instantaneous file ---
c
         confil(i)(1:) =  ' '
         confil(i)(1:) = flrtsa(1:nchar)
         if( i .EQ. IDXCRS ) then
            confil(i)(nchar+1:) = '.ddm.inst'
         else
            confil(i)(nchar+1:) = '.ddm.finst'
         endif
         fname = confil(i)
         nopen = nopen + 1
         call getunit(iowcon(i))
         open(unit=iowcon(i),file=confil(i),ERR=7002,
     &                              form='UNFORMATTED',status='UNKNOWN')
         if( i .EQ. IDXCRS ) then
            write(iounit,9000)
     &          'DDM INST file for master grid        (unit):',iowcon(i)
            write(iounit,9002) confil(i)(:istrln(confil(i)))
         else
            write(iounit,9000)
     &          'DDM INST file for nested grids       (unit):',iowcon(i)
            write(iounit,9002) confil(i)(:istrln(confil(i)))
         endif
   40 continue
c
c    ---- surface concentrations file ----
c
      if( .not.lsfcfl ) then
         write(iounit,9006)
     &              'DDM Surface file                     (unit):',
     &                                                ' Not supplied'
      else
         do n = 1,ngrid
            if (.NOT. lddmcalc(n) ) cycle
            if( .NOT. lcdfout ) then
               sfcfil(n) = ' '
               fname = flrtsa(1:nchar)
               write(fname(nchar+1:),'(a,i2.2)') '.ddm.grd',n
               sfcfil(n)(1:nchar+10) = fname
               nopen = nopen + 1
               call getunit(iowsfc(n))
               open(unit=iowsfc(n),file=sfcfil(n),ERR=7002,
     &                           form='UNFORMATTED',status='UNKNOWN')
               write(iounit,9000)
     &                   'DDM Surface file                     (unit):',
     &                                                      iowsfc(n)
               write(iounit,9002) sfcfil(n)(:istrln(sfcfil(n)))
            else
               sfcfil(n) = ' '
               fname = flrtsa(1:nchar)
               write(fname(nchar+1:),'(a,i2.2,a)') '.ddm.grd',n,'.nc'
               sfcfil(n)(1:nchar+13) = fname
               nopen = nopen + 1
               action = 'Opening DDM surface output file.'
               call ncf_createfile(sfcfil(n),action,iowsfc(n))
               write(iounit,9000) 'DDM Surface file:'
               write(iounit,9002) sfcfil(n)(:istrln(sfcfil(n)))
            endif
         enddo
      endif
c
c    ---- receptor file ----
c
      if( lrcpfil ) then
         avgfil(1:) = ' '
         avgfil(1:) = flrtsa(1:nchar)
         avgfil(nchar+1:) = '.ddm.receptor'
         fname = avgfil
         nopen = nopen + 1
         call getunit(iowrcp)
         open(unit=iowrcp,file=avgfil,ERR=7002,form='FORMATTED',
     &                                      status='UNKNOWN')
         write(iounit,9000) 
     &     'DDM Receptor concentration file       (unit):',iowrcp
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
 9001 format(/,A,I6,A,I8)
 9002 format(2A)
 9003 format(/,A,I2,1X,I2,A)
 9004 format(/,A,I2,A)
 9005 format(/,A,I2,1X,I2,1X,I8)
 9006 format(/,2A)
 9007 format(/,A,I6,2I8)
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iounit,'(//,A)') 'ERROR in STARTDDM:'
      write(iounit,'(A)') action(:istrln(action))
      write(iounit,'(/,1X,2A)') 'Input file does not exist: ',fname
      call camxerr()
c
 7001 continue
      write(iounit,'(//,A)') 'ERROR in STARTDDM:'
      write(iounit,'(A)') action(:istrln(action))
      write(iounit,'(/,1X,2A)') 'Cannot open input file: ',fname
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7002 continue
      write(iounit,'(//,A)') 'ERROR in STARTDDM:'
      write(iounit,'(/,1X,2A)') 'Cannot open output file: ',fname
      call camxerr()
c
 7003 continue
      write(iounit,'(//,A)') 'ERROR in STARTDDM:'
      write(iounit,'(/,1X,3A,/,A)') 'Nested grid restart file for ',
     &                     tectyp(:istrln(tectyp)),' not supplied.',
     &                  ' Set Flexi_Nest = .true. to ignore this file.'
      call camxerr()
c
 7004 continue
      write(iounit,'(//,A)') 'ERROR in STARTDDM:'
      write(iounit,'(/,A)') 'A source area mapping file must be ',
     &                'supplied for the master grid, Grid #1'
      call camxerr()
c
 7005 continue
      write(iounit,'(//,A)') 'ERROR in STARTDDM:'
      write(iounit,'(A)') action(:istrln(action))
      write(iounit,'(1X,2A)') 'You supplied a blank filename in the middle ',
     &                  'of the list of gridded emissions '
      write(iounit,'(1X,A,I3)') 'filenames for grid: ',i
      write(iounit,'(1X,2A)') 'Check the numbering of the list of files.'
      call camxerr()
c
 7006 continue
      write(iounit,'(//,A)') 'ERROR in STARTDDM:'
      write(iounit,'(A)') action(:istrln(action))
      write(iounit,'(1X,2A)') 'You supplied a blank filename in the middle ',
     &                  'of the list of point source '
      write(iounit,'(1X,A,I3)') 'filenames for grid: ',i
      write(iounit,'(1X,2A)') 'Check the numbering of the list of files.'
      call camxerr()
c
 7007 continue
      write(iounit,'(//,A)') 'ERROR in STARTDDM:'
      write(iounit,'(A)') action(:istrln(action))
      write(iounit,'(1X,2A)') 'You supplied a blank filename for the DDM ',
     &                  'Top Conenctrations file.'
      write(iounit,'(1X,2A)') 'But you did provide a Top Concntration file ',
     &                  'in the CAMx_Control section.'
      write(iounit,'(1X,2A)') 'This is not allowed. You must supply both ',
     &                     'or neither.'
      call camxerr()

c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
