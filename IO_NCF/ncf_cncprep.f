C*** NCF_CNCPREP
c
      subroutine ncf_cncprep(timbeg,datebeg,timend,dateend,lfirst)
      use grid
      use chmstry
      use camxfld
      use filunit
      implicit none
c 
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c     NCF_CNCPREP reads the initial conditions NETCDF file and checks
c     the global variables for consistency with user defined values.
c     It also prepares mapping variables that map the IC species list
c     to the internal CAMx species list
c 
c      Copyright 1996 - 2022
c     Ramboll
c           
c     Input arguments: 
c         timbeg  R model begin time
c         begdat  I model begin date (YYJJJ)
c         timend  R model end time
c         dateend I model end date (YYJJJ)
c        lfirst  - .TRUE. if this is the first time this routine is called
c     Output arguments: 
c        none
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
      include 'camx.inc'
      include 'flags.inc'
      include 'ncf_iodat.inc'
      include 'netcdf.inc'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      real    timbeg
      integer datebeg
      real    timend
      integer dateend
      logical lfirst
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
      character*20  spec_units(MXSPEC)
      character*20  spec_long_name(MXSPEC)
      character*60  spec_desc(MXSPEC)
      character*60  spec_coords(MXSPEC)
      character*10  icfil, name_in, this_var, avfil
      character*4   note(60), ifile(10)
      character*4   icspec(10,MXSPEC)
      real          orgx, orgy, dx, dy, dxf, dyf, orgxf, orgyf, zero
      integer       ierr, nlayer, num_dims, l, iproj, ione, i, n, izero
      integer       itmp
c
c-----------------------------------------------------------------------
c    Data Statements:
c-----------------------------------------------------------------------
c
      data icfil /'INITIAL   '/
      data avfil /'AVERAGE   '/
      data ione /1/
      data zero /0./
      data izero /0/
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c --- set the string for error messages ---
c
      action = 'Reading the NetCDF initial conditions file.'
c
c --- get the type of file to make sure it is a boundary file ---
c
      this_var = 'CAMx_NAME'
      ierr = nf_get_att_text(iic, NF_GLOBAL, this_var, name_in)
      if( ierr .NE. NF_NOERR ) then
         this_var = 'NAME'
         ierr = nf_get_att_text(iic, NF_GLOBAL, this_var, name_in)
         if( ierr .NE. NF_NOERR ) goto 7002
      endif
      if( name_in(:7) .NE. icfil(:7) ) goto 7001
c
c --- call routine to make sure grid defintion is consistent ---
c
      call ncf_chk_griddef(iic,action,1,.TRUE.,.TRUE.,.FALSE.,.FALSE.,itmp)
c
c --- call routine to make sure file spans the episode ---
c
      call ncf_chk_tstep(iic,action,datebeg,timbeg,datebeg,timbeg,.FALSE.)
c
c --- call routine to setup species mappping array ---
c
      call ncf_set_species_mapping(iic,action,spname,nspec,nicspc,licmap)
c
c --- assign values for the output file ---
c
      read(runmsg,'(60A1)') note
      read(avfil,'(10a1)') (ifile(n),n=1,10)
      do l = 1,navspc
        if(lavmap(l) .GT. 0 ) then
          spavnam(l) = spname(lavmap(l))
          read(spname(lavmap(l)),'(10a1)') (icspec(n,l),n=1,10)
        endif
      enddo
      if( l3davg(1) ) then
        nlayer = nlay(1)
        num_dims = 4
      else
        nlayer = 1
        num_dims = 4
      endif
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
c --- Write average output concentration file headers ---
c
      if( lfirst ) then
          if( lcdfout ) then
            call ncf_set_vars_base()
            call ncf_set_tstep(datebeg,timbeg,dateend,timend)
            action = 'Writing master grid output average file.'
            call ncf_set_specatt_avrg(spec_units,spec_long_name,spec_desc,
     &                                                        spec_coords)
            call ncf_set_global('AVERAGE   ',1,datebeg,timbeg,dateend,timend,
     &                                                       nlayer,navspc)
            call ncf_wrt_dim(action,iavg(1),1,ncol(1),nrow(1),nlayer,navspc)
            call ncf_wrt_global(action,iavg(1),navspc,spavnam,.FALSE.)
            call ncf_wrt_vars_base(action,iavg(1))
            call ncf_wrt_vars_species(action,iavg(1),ncol(1),nrow(1),nlayer,
     &                              navspc,spavnam,spec_units,spec_long_name,
     &                                        spec_desc,spec_coords,num_dims)
            call ncf_enddef_file(action,iavg(1))
            call ncf_wrt_data_grid(action,iavg(1),1,ncol(1),nrow(1),
     &                           orgx,orgy,dx,dy,nlayer,cellat(iptr2d(1)),
     &                                   cellon(iptr2d(1)),topo(iptr2d(1)))
          else
            rewind(iavg(1))
            write(iavg(1)) ifile,note,itzon,navspc,datebeg,timbeg/100.,
     &                                                     dateend,timend/100.
            write(iavg(1)) polelon,polelat,iuzon,orgx,orgy,dx,dy,ncol(1),nrow(1),nlayer,
     &                                       iproj,izero,tlat1,tlat2,zero
            write(iavg(1)) ione,ione,ncol(1),nrow(1)
            write(iavg(1)) ((icspec(n,l),n=1,10),l=1,navspc)
          endif
      endif
      if( ngrid .EQ. 1 ) goto 9999
c
c-----Fine grid average headers
c
      do i = 2,ngrid 
        dxf   = dx/float(meshold(i))
        dyf   = dy/float(meshold(i))
        orgxf = orgx + dx*(inst1(i)-1)
        orgyf = orgy + dy*(jnst1(i)-1)
        if( l3davg(i) ) then
          nlayer = nlay(i)
        else
          nlayer = 1
        endif
        if( l3davg(i) ) then
          nlayer = nlay(i)
        else
          nlayer = 1
        endif
        if( lfirst ) then
          if( lcdfout ) then
             write(action,'(A,I2)') 'Writing output average file for grid: ',i
             call ncf_set_specatt_avrg(spec_units,spec_long_name,spec_desc,
     &                                                        spec_coords)
             call ncf_set_global('AVERAGE   ',i,datebeg,timbeg,dateend,timend,
     &                                                      nlayer,navspc)
             call ncf_wrt_dim(action,iavg(i),i,ncol(i),nrow(i),nlayer,navspc)
             call ncf_wrt_global(action,iavg(i),navspc,spavnam,.FALSE.)
             call ncf_wrt_vars_base(action,iavg(i))
             call ncf_wrt_vars_species(action,iavg(i),ncol(i),nrow(i),nlayer,
     &                            navspc,spavnam,spec_units,spec_long_name,
     &                                       spec_desc,spec_coords,num_dims)
             call ncf_enddef_file(action,iavg(i))
             call ncf_wrt_data_grid(action,iavg(i),i,ncol(i),nrow(i),
     &                           orgx,orgy,dx,dy,nlayer,cellat(iptr2d(i)),
     &                                  cellon(iptr2d(i)),topo(iptr2d(i)))
          else
             rewind(iavg(i))
             write(iavg(i)) ifile,note,itzon,navspc,datebeg,timbeg/100.,
     &                                                    dateend,timend/100.
             write(iavg(i)) polelon,polelat,iuzon,orgxf,orgyf,dxf,dyf,ncol(i)-2,
     &                       nrow(i)-2,nlayer,iproj,izero,tlat1,tlat2,zero
             write(iavg(i)) ione,ione,ncol(i)-2,nrow(i)-2
             write(iavg(i)) ((icspec(n,l),n=1,10),l=1,navspc)
          endif
        endif
      enddo 
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,a)') 'ERROR in NCF_CNCPREP:'
      write(iout,'(2A)',ERR=9999)'Reading global attributes from ',
     &                                          'initial conditions file.'
      write(iout,'(2A)',ERR=9999)'Cannot read: ',this_var(:istrln(this_var))
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7001 continue
      write(iout,'(//,a)') 'ERROR in NCF_CNCPREP:'
      write(iout,'(2A)',ERR=9999)'Reading the header of initial ',
     &                                        'conditions file.'
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7002 continue
      write(iout,'(//,A)') 'ERROR in NCF_CNCPREP:'
      write(iout,'(2A)',ERR=9999)'Reading global attributes from ',
     &                                          'initial conditions file.'
      write(iout,'(2A)',ERR=9999) 'Cannot find global variable for ',
     &                                       'type of file: CAMx_NAME'
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 9999 continue
      lfirst = .FALSE.
      return
      end
