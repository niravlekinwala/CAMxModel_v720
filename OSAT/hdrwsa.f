c*** HDRWSA
c
      subroutine hdrwsa(iounit,fname,ftype,nspcs,nzcell,
     &             idate,btim,jdate,etim,spec_units,spec_long_name,
     &                                          spec_desc,spec_coords)
      use filunit
      use grid
      use camxcom
      use camxfld
      use tracer
      use procan
      use ncf_iomod
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c   Description:
c     This routine writes the header information to the output tracer
c     instantaneous files and tracer surface concentration files.
c     The data is written to the unit number in the argument list.
c
c      Copyright 1996 - 2022
c     Ramboll
c
c    Argument description:
c      iounit         I  unit number of file to write
c      fname          C  name of file being written
c      ftype          C  file type of file being written 
c      nspcs          I  number of species to write to file
c      nzcell         I  number of layers (will be 1 for surface file)
c      idate          I   beginning date of simulation (YYJJJ)
c      btim           R   beginning time of simulation
c      jdate          I   ending date of simulation (YYJJJ)
c      etim           R   ending time of simulation
c      spec_units     C array of units for this each species
c      spec_long_name C array of "long names" for each each species
c      spec_desc      C array of desciption for this each species
c      spec_coords    C array of coordinates for this each species

c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c      1/20/99   Grid cell size on file should be meters for all
c                cartesian projections (UTM, LCP, PSP)
c      10/24/01  Removed BSWAP and converted integer strings to character*4
c      11/06/01  Input dates are now Julian
c      8/25/06   Surface output files now all UAM format, one file per grid
c      1/04/11   Revised for new header format
c
c-----------------------------------------------------------------------
c   Include files:
c-----------------------------------------------------------------------
c
      include 'flags.inc'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      integer      iounit(*)
      character*(*) fname(*)
      character*10 ftype
      integer      nspcs
      integer      nzcell       
      integer      idate
      real         btim
      integer      jdate
      real         etim
      character*20 spec_units(*)
      character*20 spec_long_name(*)
      character*60 spec_desc(*)
      character*60 spec_coords(*)
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      character*200 action
      character*10 ncf_filtype
      character*4  ifile(10), note(60), ispec(10,MXTRSP)
      integer      i, j, ndate, ndlast, izero, ione, num_dims
      real         ttime, ttlast, zero, factr
      real         xorgf, yorgf, dxf, dyf
      logical      lwrite_file
c
c-----------------------------------------------------------------------
c    Data statements:
c-----------------------------------------------------------------------
c
      data zero /0.0/
      data izero /0/
      data ione  /1/
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
      if( tectyp .EQ. SA ) then
         ncf_filtype = "SA        "
      else if( tectyp .EQ. DDM .OR. tectyp .EQ. HDDM ) then
         ncf_filtype = "DDM       "
      else if( tectyp .EQ. RTRAC .OR. tectyp .EQ. RTCMC ) then
         ncf_filtype = "RTRAC     "
      else
         ncf_filtype = "CPA       "
      endif
      ndate = idate
      ndlast = jdate
      ttime  = AINT(ANINT(btim)/100.) + amod(ANINT(btim),100.)/60.
      ttlast = AINT(ANINT(etim)/100.) + amod(ANINT(etim),100.)/60.
c
c   --- set scaling factor for coordinates ---
c
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
c
c   ---- put species names into integer array ----
c 
      nspcout = 0
      ptnameout = ' '
      do 10 j=1,nspcs
        if( loutsa(j) .OR. ftype .EQ. 'AIRQUALITY' ) then
            nspcout = nspcout + 1
            read(ptname(j),'(10A1)') (ispec(i,nspcout),i=1,10)
            if( lcdfout ) then
                if( lddm .OR. lhddm ) then
                   ptnameout(nspcout) = ptlong(j)
                else
                   ptnameout(nspcout)(1:10) = ptname(j)
                endif
            endif
        endif
   10 continue
c
c   --- master grid file header ---
c
      lwrite_file = .TRUE.
      read(ftype,'(10A1)') (ifile(i),i=1,10)
      read(runmsg(1:60),'(60A1)') (note(i),i=1,60)
      if( (lddm .OR. lhddm) .AND. .NOT. lddmcalc(1) ) lwrite_file = .FALSE.
      if( lwrite_file ) then
          n = 1
          nz = 1
          num_dims = 4
          if (nzcell.gt.1) nz = nlay(1)
          if( ltrace .AND. lsa_3davrg ) nz = nlay(1)
          if( lirr ) then
             nz = 1
             if( l3davg(n) ) nz = nlay(n)
          endif
           ntotspout = nspcout
          if( lcdfout ) then
              call ncf_set_vars_base()
              call ncf_set_tstep(idate,btim,jdate,etim)
          endif
          if( .NOT. lcdfout .OR. ftype .EQ. 'AIRQUALITY' ) then
             write(iounit(1),ERR=7000) ifile,note,itzon,nspcout,ndate,ttime,
     &                                                    ndlast,ttlast
             write(iounit(1),ERR=7000) plon,plat,iutm,orgx,orgy, 
     &                                   dx,dy,ncol(1),nrow(1), 
     &                                   nz,iproj,izero,t1,t2,zero
             write(iounit(1),ERR=7000) ione,ione,ncol(1),nrow(1)
             write(iounit(1),ERR=7000) ((ispec(i,j),i=1,10),j=1,nspcout)
             if( ftype .EQ. 'AIRQUALITY' ) goto 9999
          else if( lcdfout .AND. ftype .NE. 'AIRQUALITY' ) then
              action = 'Writing tracer output file: '//fname(1)(:istrln(fname(1)))
              call ncf_set_global(ncf_filtype,1,ndate,ttime*100.,ndlast,ttlast,
     &                                                       nz,nspcout)
              call ncf_wrt_dim(action,iounit(1),1,ncol(1),nrow(1),nz,nspcout)
              call ncf_wrt_global(action,iounit(1),nspcout,ptnameout,.FALSE.)
              if( tectyp .EQ. SA ) then
                 call ncf_wrt_global_sa(action,iounit(1))
              else if( tectyp .EQ. DDM .OR. tectyp .EQ. 'HDDM' ) then
                 call ncf_wrt_global_ddm(action,iounit(1))
              endif
              call ncf_wrt_vars_base(action,iounit(1))
              call ncf_wrt_vars_species(action,iounit(1),ncol(1),nrow(1),nz,
     &                        nspcout,ptnameout,spec_units,spec_long_name,
     &                                       spec_desc,spec_coords,num_dims)
              call ncf_enddef_file(action,iounit(1))
              call ncf_wrt_data_grid(action,iounit(1),1,ncol(1),nrow(1),
     &                           orgx,orgy,dx,dy,nz,cellat(iptr2d(1)),
     &                                   cellon(iptr2d(1)),topo(iptr2d(1)))
          endif
      endif
c
c   --- nested grid headers ---
c
      do n = 2,ngrid
         if( (lddm .OR. lhddm) .AND. .NOT. lddmcalc(n) ) cycle
         dxf   = dx/float(meshold(n))
         dyf   = dy/float(meshold(n))
         xorgf = orgx + dx*(inst1(n)-1)
         yorgf = orgy + dy*(jnst1(n)-1)
         nz = 1
         num_dims = 4
         if (nzcell.gt.1) nz = nlay(n)
         if( ltrace .AND. lsa_3davrg ) nz = nlay(1)
         if( lirr ) then
            nz = 1
            if( l3davg(n) ) nz = nlay(n)
         endif
         if( .NOT. lcdfout .OR. ftype .EQ. 'AIRQUALITY' ) then
            write(iounit(n),ERR=7000) ifile,note,itzon,nspcout,ndate,ttime,
     &                                                       ndlast,ttlast
            write(iounit(n),ERR=7000) plon,plat,iutm,xorgf,yorgf,dxf,dyf,
     &                                   ncol(n)-2,nrow(n)-2,nz,iproj,izero,
     &                                                        t1,t2,zero
            write(iounit(n),ERR=7000) ione,ione,ncol(n)-2,nrow(n)-2
            write(iounit(n),ERR=7000) ((ispec(i,j),i=1,10),j=1,nspcout)
            if( ftype .EQ. 'AIRQUALITY' ) goto 111
         else
          action = 'Writing tracer output file: '//fname(n)(:istrln(fname(n)))
          call ncf_set_global(ncf_filtype,n,ndate,ttime*100.,ndlast,ttlast,
     &                                                       nz,nspcout)
          call ncf_wrt_dim(action,iounit(n),n,ncol(n),nrow(n),nz,nspcout)
          call ncf_wrt_global(action,iounit(n),nspcout,ptnameout,.FALSE.)
          if( tectyp .EQ. SA ) then
             call ncf_wrt_global_sa(action,iounit(n))
          else if( tectyp .EQ. DDM .OR. tectyp .EQ. 'HDDM' ) then
             call ncf_wrt_global_ddm(action,iounit(n))
          endif
          call ncf_wrt_vars_base(action,iounit(n))
          call ncf_wrt_vars_species(action,iounit(n),ncol(n),nrow(n),nz,
     &                        nspcout,ptnameout,spec_units,spec_long_name,
     &                                         spec_desc,spec_coords,num_dims)
          call ncf_enddef_file(action,iounit(n))
          call ncf_wrt_data_grid(action,iounit(n),n,ncol(n),nrow(n),
     &                           orgx,orgy,dx,dy,nz,cellat(iptr2d(n)),
     &                                   cellon(iptr2d(n)),topo(iptr2d(n)))
         endif
      enddo
c
c  --- reset the output species name ---
c
  111 continue
      ptnameout = ' '
      do j=1,nspcs
         if( lddm .OR. lhddm ) then
             ptnameout(j) = ptlong(j)
         else
             ptnameout(j)(1:10) = ptname(j)
         endif
      enddo
c
c  --- return to calling routine ---
c
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,A)') 'ERROR in HDRWSA:'
      write(iout,9000,ERR=9999)'Writing output tracer file: ',
     &                                     fname(n)(:istrln(fname(n)))
      call camxerr()
c
c-----------------------------------------------------------------------
c    Format statements:
c-----------------------------------------------------------------------
c
 9000 format(/,1X,2A)
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
c
