      subroutine rdsrfmod(igrd,nox,noy,nsurf,endtim,enddate,
     &                    soilmass,vegiemass)
      use filunit
      use grid
      use chmstry
      use camxcom
      use camxfld
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routine reads a surface model restart file and fills data into
c   the surface mass arrays.  It then writes the header of the new output
c   surface model file.
c
c      Copyright 1996 - 2022
c     Ramboll
c
c      Argument description:
c       Inputs:
c           igrd      I  grid index
c           nox       I  number of X cells in the grid
c           noy       I  number of Y cells in the grid
c           nsurf     I  number of surface model species
c           endtim    R  Simulation end time (HHMM)
c           enddate   I  Simulation end date (YYJJJ)
c       Outputs:
c           soilmass  R  surface soil mass (mol/ha, g/ha)
c           vegiemass R  surface veg mass (mol/ha, g/ha)
c       
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     05/14/13   --cemery--    Original development
c     09/08/21   --cemery--    Modified output species names for
c                              consistency with SAT and RTRAC dep output
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      implicit none
      include 'camx.prm'
      include 'flags.inc'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      integer   igrd
      integer   nox
      integer   noy
      integer   nsurf
      integer   enddate
      real      endtim
      real      soilmass(nox,noy,nsurf)
      real      vegiemass(nox,noy,nsurf)
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
      character*200 action
      character*20  spec_units(2*MXSPEC)
      character*20  spec_long_name(2*MXSPEC)
      character*60  spec_desc(2*MXSPEC)
      character*60  spec_coords(2*MXSPEC)
      character*10  avfil, cfile, cspec, srfmodname(2*MXSPEC)
      character*4   ifile(10), note(60), ispec(10), mspec(10,2*MXSMSPC)
      integer       iunit,izero,ione, numxcells, numycells
      integer       isp, i, j, n
      integer       itzn, nspc, idat1, idat2
      integer       iutm, nx, ny, nz, idum, imod
      integer       iproj,istag
      real          tim1, tim2, orgx, orgy, plon, plat, dx, dy
      real          t1, t2
      real          dxf, dyf
      real          zero, rdum
c
      real srfin(MXCELLS,MXCELLS)
c
c-----------------------------------------------------------------------
c    Data statements:
c-----------------------------------------------------------------------
c
      data avfil /'AVERAGE'/
      data izero,ione,zero /0,1,0./
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c  --- initailize to lower bound ---
c
      do isp = 1,nsurf
        do j = 1,noy
          do i = 1,nox
            soilmass(i,j,isp) = 0.
            vegiemass(i,j,isp) = 0.
          enddo
        enddo
      enddo
c
c  --- read 1st input header record and check inputs ---
c
      iunit = ismin(igrd)
      if (iunit.ne.0) then
        read(iunit,ERR=7002) ifile,note,itzn,nspc,idat1,tim1,idat2,tim2
        if (INT(tim2) .EQ. 24 ) then
          idat2 = idat2 + 1
          tim2 = 0.
          if( MOD(idat2,1000) .GT. 365 ) then
             if( MOD(INT(idat2/1000),4) .EQ. 0 ) then
                if( MOD(idat2,1000) .EQ. 367 )
     &                     idat2 = (INT(idat2/1000)+1)*1000 + 1
             else
                idat2 = (INT(idat2/1000)+1)*1000 + 1
             endif
          endif
        endif
        write(cfile,'(10A1)') (ifile(i),i=1,10)
        if( cfile .NE. avfil ) goto 7003
        tim1 = 100.*tim1
        tim2 = 100.*tim2
        if( idat1 .GT. begdate ) goto 7004
        if( idat1 .EQ. begdate .AND. tim1 .GT. begtim ) goto 7004
c
c  --- read 2nd header record and check inputs ---
c
        read(iunit,ERR=7002) plon,plat,iutm,orgx,orgy,dx,dy,nx,ny,nz,
     &                       iproj,istag,t1,t2,rdum
        if( nx .NE. nox .OR. ny .NE. noy .OR. nz .NE. 1 ) goto 7005
c
c  --- skip the next 2 records ---
c
        read(iunit,ERR=7002)
        read(iunit,ERR=7002) 
c
c  --- read the date and time and make sure it is the correct hour ---
c
  111   continue
        read(iunit,ERR=7007,END=7006) idat1, tim1, idat2, tim2
        tim1 = 100*tim1
        tim2 = 100*tim2
        if( INT(tim2) .EQ. 24 ) then
          idat2 = idat2 + 1
          tim2 = 0.
          if( MOD(idat2,1000) .GT. 365 ) then
             if( MOD(INT(idat2/1000),4) .EQ. 0 ) then
                if( MOD(idat2,1000) .EQ. 367 )
     &                       idat2 = (INT(idat2/1000)+1)*1000 + 1
             else
                idat2 = (INT(idat2/1000)+1)*1000 + 1
             endif
          endif
        endif
        write(iout,'(a40,f7.0,i8.5,a,i3)')
     &    'Read surface model file at ',tim1, idat1,' grid',igrd

        if( (idat1 .LT. begdate .OR. 
     &      (idat1 .EQ. begdate .AND. tim1 .LE. begtim)) .AND.
     &      (idat2 .GT. begdate .OR. 
     &      (idat2 .EQ. begdate .AND. tim2 .GT. begtim))) then
c
c  --- read the concentrations for each species ---
c
          call flush(6)
          do 20 isp = 1,nspc
            read(iunit) idum,(ispec(n),n=1,10),
     &                       ((srfin(i,j),i=1,nx),j=1,ny)
            write(cspec,'(10A1)') (ispec(n),n=1,10)
c
c  --- find this species in the surface model species list ---
c
            do imod = 1,nsurf
              if( cspec(3:) .EQ. smspc(imod) ) then
                do j = 2,ny-1
                  do i = 2,nx-1
                    if (cspec(1:1) .EQ. 'S') then
                      soilmass(i,j,imod) = srfin(i,j)
                    elseif (cspec(1:1) .EQ. 'V') then
                      vegiemass(i,j,imod) = srfin(i,j)
                    endif
                  enddo
                enddo
                goto 20
              endif
            enddo
            write(iout,'(//,A)') 'ERROR in RDSRFMOD:'
            write(iout,'(A)') 'Species in surface mass file: '
            write(iout,'(A)') cspec(:istrln(cspec))
            write(iout,'(2A)') 'not found in surface model species',
     &                        ' list .. Stopping.'
            call camxerr()
 20       continue
c
c  --- if not the right hour, read through it ---
c
        else
          do isp = 1,nspc
            read(iunit,ERR=7007) 
          enddo
          goto 111
        endif
c
c  --- close the inpu file ---
c
        close(iunit)
      endif
c
c --- Write output surface model file headers
c
      idat1 = begdate
      idat2 = enddate
      tim1 = begtim/100.
      tim2 = endtim/100.
      iutm = iuzon
      plon = polelon
      plat = polelat
      t1   = tlat1
      t2   = tlat2
      if (llatlon) then
        iproj = 0
        orgx = xorg
        orgy = yorg
        dx = delx
        dy = dely
      else
        orgx = 1000.*xorg
        orgy = 1000.*yorg
        dx = 1000.*delx
        dy = 1000.*dely
        if (lutm) then
          iproj = 1
        elseif (lambrt) then
          iproj = 2
        elseif (lrpolar) then
          iproj = 3
        elseif (lpolar) then
          iproj = 4
        elseif (lmerc) then
          iproj = 5
        endif
      endif
      if (igrd.gt.1) then
        dxf  = dx/float(meshold(igrd))
        dyf  = dy/float(meshold(igrd))
        orgx = orgx + dx*(inst1(igrd)-1)
        orgy = orgy + dy*(jnst1(igrd)-1)
        dx = dxf
        dy = dyf
      endif
      numxcells = nox - 2 
      numycells = noy - 2 

      read(runmsg(1:60),'(60a1)') (note(n),n=1,60)
      read(avfil,'(10a1)') (ifile(n),n=1,10)
      do isp = 1,nsurf
         mspec(1,isp) = 'S'
         mspec(2,isp) = '_'
         read(smspc(isp)(1:8),'(8A1)') (mspec(i,isp),i=3,10)
         srfmodname(isp) = 'S_'//smspc(isp)(1:8)
      enddo
      do isp = 1,nsurf
         mspec(1,isp+nsurf) = 'V'
         mspec(2,isp+nsurf) = '_'
         read(smspc(isp)(1:8),'(8A1)') (mspec(i,isp+nsurf),i=3,10)
         srfmodname(isp+nsurf) = 'V_'//smspc(isp)(1:8)
      enddo
c
      if( lcdfout) then
         call ncf_set_vars_base()
         call ncf_set_tstep(begdate,begtim,enddate,endtim)
         call ncf_set_specatt_srfmod(spec_units,spec_long_name,spec_desc,
     &                                     spec_coords,2*nsurf,srfmodname)
      endif
      iunit = ismout(igrd)
      if( .NOT. lcdfout ) then
         write(iunit) ifile,note,itzon,2*nsurf,idat1,tim1,idat2,tim2
         write(iunit) plon,plat,iutm,orgx,orgy,dx,dy,numxcells,numycells,
     &             ione,iproj,izero,t1,t2,zero
         write(iunit) ione,ione,numxcells,numycells
         write(iunit) ((mspec(n,isp),n=1,10),isp=1,2*nsurf)
      else
         write(action,'(2A,I2)') 'Writing surface model ',
     &                                    'output file for grid ',igrd
         call ncf_set_global('SURFACE   ',igrd,idat1,tim1,idat2,tim2,
     &                                                       ione,2*nsurf)
         call ncf_wrt_dim(action,iunit,igrd,nox,noy,ione,2*nsurf)
         call ncf_wrt_global(action,iunit,2*nsurf,srfmodname,.FALSE.)
         call ncf_wrt_vars_base(action,iunit,igrd)
         call ncf_wrt_vars_species(action,iunit,nox,noy,ione,
     &                    2*nsurf,srfmodname,spec_units,spec_long_name,
     &                                           spec_desc,spec_coords,4)
         call ncf_enddef_file(action,iunit)
         call ncf_wrt_data_grid(action,iunit,igrd,nox,noy,
     &                           orgx,orgy,dx,dy,ione,cellat(iptr2d(igrd)),
     &                                   cellon(iptr2d(igrd)),topo(iptr2d(igrd)))
      endif
 
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7002 continue
      write(iout,'(//,A)') 'ERROR in RDSRFMOD:'
      write(iout,'(2A,i3)')'ERROR: Reading header of Surface mass',
     &                     ' file for grid: ',igrd
      call camxerr()
c
 7003 continue
      write(iout,'(//,A)') 'ERROR in RDSRFMOD:'
      write(iout,'(2A)') 'ERROR:  Surface mass input file is not ',
     &                   'labelled AVERAGE.'
      call camxerr()
c
 7004 continue
      write(iout,'(//,A)') 'ERROR in RDSRFMOD:'
      write(iout,'(2A)') 'ERROR:  Surface mass input file is not for ',
     &                                          'correct time period.'
      write(iout,*) '   ***  Episode ***'
      write(iout,'(a,i10.5)') 'Date   : ',begdate
      write(iout,'(a,f10.1)') 'Time   : ',begtim
      write(iout,*) '   ***  Surface Model File ***'
      write(iout,'(a,i10.5)') 'Date   : ',idat1
      write(iout,'(a,f10.1)') 'Time   : ',tim1
      call camxerr()
c
 7005 continue
      write(iout,'(//,A)') 'ERROR in RDSRFMOD:'
      write(iout,'(2A)') 'ERROR: Surface mass input file does not',
     &                   ' appear to be for the correct grid.' 
      write(iout,*) '   ***  Grid ***'
      write(iout,*) 'No. of Cells : (',nox,',',noy,')'
      write(iout,*) 'No. of layers:  ',1
      write(iout,*) '   ***  RTRAC File ***'
      write(iout,*) 'No. of Cells : (',nx,',',ny,')'
      write(iout,*) 'No. of layers:  ',nz
      call camxerr()
c
 7006 continue
      write(iout,'(//,A)') 'ERROR in RDSRFMOD:'
      write(iout,'(2A)') 'ERROR:  Premature end-of-file reached in ',
     &                                    'surface mass input file.'
      write(iout,'(2A)') 'Make sure the file contains the correct',
     &                   ' date/time.'
      call camxerr()
c
 7007 continue
      write(iout,'(//,A)') 'ERROR in RDSRFMOD:'
      write(iout,'(2A,I3)') 'ERROR:  Reading surface mass input ',
     &                      ' file for grid: ',igrd
      call camxerr()
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
c
      return
      end
