      subroutine smpprep(lrtrac,endtim,enddate)
      use camxcom
      use camxfld
      use grid
      use chmstry
      use pigsty
      use filunit
      use tracer
      use rtracchm
      use ncf_iomod
c 
c----CAMx v7.20 220430
c 
c     SMPPREP writes headers to new SAMPLING GRID output files.
c 
c      Copyright 1996 - 2022
c     Ramboll
c           
c     Modifications: 
c        8/2/05       Generalized for regular model species or
c                     RTRAC species.
c        1/04/11      Revised for new header format
c 
c     Input arguments: 
c        endtim              model end time (HHMM)
c        enddate             model end date (YYJJJ)
c             
c     Output arguments: 
c        none
c             
c     Routines Called: 
c        none
c             
c     Called by: 
c        STARTUP 
c
      implicit none
      include 'flags.inc'
c
      integer MXOUT
      parameter (MXOUT = MAX(MXSPEC,MXTRSP))
c
      logical lrtrac
      integer enddate
      real endtim
      integer idat1,idat2,n,iunit,iutm,l,nsg,izero,ione,nav
      integer iproj
      real tim1,tim2,orgx,orgy,dx,dy,zero
      real plon,plat,t1,t2
      character*4 ifile(10),note(60)
      character*10 avfil
c
      character*4   smpspec(10,MXOUT)
      character*10  smpname(MXOUT), filedesc
      character*20  spec_units(MXSPEC)
      character*20  spec_long_name(MXSPEC)
      character*60  spec_desc(MXSPEC)
      character*60  spec_coords(MXSPEC)
      character*200 action
c
      integer istrln
c
      data avfil /'AVERAGE   '/
      data izero,ione /0,1/
      data zero /0./
c
c-----Entry point
c
      filedesc = avfil
      idat1 = begdate
      idat2 = enddate
      tim1 = begtim/100.
      tim2 = endtim/100.
      read(runmsg(1:60),'(60a1)') (note(n),n=1,60)
      read(avfil,'(10a1)') (ifile(n),n=1,10)
c
c-----Loop over sampling grids; prep grid and species info
c
      do nsg = 1,nsample
        iunit = isample(nsg)
        if (lrtrac) iunit = iowsmp(nsg)
c
        iutm = iuzon
        plon = polelon
        plat = polelat
        t1   = tlat1
        t2   = tlat2
        if (llatlon) then
          iproj = 0
          orgx = xorgsmp(nsg)
          orgy = yorgsmp(nsg)
          dx = delx/meshsmp(nsg)
          dy = dely/meshsmp(nsg)
        else
          orgx = 1000.*xorgsmp(nsg)
          orgy = 1000.*yorgsmp(nsg)
          dx = 1000.*delx/meshsmp(nsg)
          dy = 1000.*dely/meshsmp(nsg)
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
c
        if (lrtrac) then
          filedesc = 'RTRAC     '
          nav = nrtrac
          do l = 1,nrtrac
            read(ptname(l),'(10a1)') (smpspec(n,l),n=1,10)
          enddo
        else
          nav = navspc
          do l = 1,navspc
            read(spname(lavmap(l)),'(10a1)') (smpspec(n,l),n=1,10)
            smpname(l) = spname(lavmap(l))
            if( lgas(lavmap(l)) .AND. lgas_out_ppm ) then
               spec_units(l) = "ppmv"
            else
               spec_units(l) = "micrograms m-3"
            endif
            spec_long_name(l) = smpname(l)
            spec_desc(l) = smpname(l)(:istrln(smpname(l)))//" air concentration"
            spec_coords(l) = "latitude longitude"
          enddo
        endif
c
c-----Write header
c
        if( lcdfout ) then
          call ncf_set_vars_base()
          call ncf_set_tstep(begdate,begtim,enddate,endtim)
        endif
        if( .NOT. lcdfout ) then
           rewind(iunit)
           write(iunit) ifile,note,itzon,nav,idat1,tim1,idat2,tim2
           write(iunit) plon,plat,iutm,orgx,orgy,dx,dy,ncolsmp(nsg),
     &                  nrowsmp(nsg),ione,iproj,izero,t1,t2,zero
           write(iunit) ione,ione,ncolsmp(nsg),nrowsmp(nsg)
           write(iunit) ((smpspec(n,l),n=1,10),l=1,nav)
        else
           write(action,'(2A,I2)') 'Writing sampling ',
     &                                    'grid file for sampling grid: ',nsg
           call ncf_set_global_smp(filedesc,nsg,orgx,orgy,dx,dy,
     &                         idat1,tim1*100.,idat2,tim2,ione,navspc)
           call ncf_wrt_dim(action,iunit,nsg,ncolsmp(nsg),nrowsmp(nsg),ione,nav)
           call ncf_wrt_global(action,iunit,nav,smpname,.TRUE.)
           call ncf_wrt_vars_base(action,iunit,nsg)
           call ncf_wrt_vars_species(action,iunit,ncolsmp(nsg),nrowsmp(nsg),ione,
     &                                   navspc,smpname,spec_units,spec_long_name,
     &                                                     spec_desc,spec_coords,4)
           call ncf_enddef_file(action,iunit)
           call ncf_wrt_data_gridsmp(action,iunit,nsg,ncolsmp(nsg),nrowsmp(nsg),
     &                                                         orgx,orgy,dx,dy)
        endif
      enddo
c
      return
      end
