      subroutine depprep(endtim,enddate)
      use camxcom
      use camxfld
      use filunit
      use grid
      use chmstry
c 
c----CAMx v7.20 220430
c 
c     DEPPREP generates new deposition output species names and writes 
c     headers to new DEPOSITION output files.
c 
c      Copyright 1996 - 2022
c     Ramboll
c           
c     Modifications: 
c        8/25/06   Dep output files now all UAM format, one file per grid
c        1/04/11   Revised for new header format
c        05/17/16  Added in-line Ix emissions to dep output
c        09/08/21  Modified output species names for consistency with
c                  SAT and RTRAC dep output
c 
c     Input arguments: 
c        endtim              model end time (HHMM)
c        enddate             model end date (YYJJJ)
c             
c     Output arguments: 
c        none
c             
c     Routines Called: 
c        ISTRLN
c             
c     Called by: 
c        STARTUP 
c
      include 'camx.prm'
      include 'flags.inc'
c
      character*200 action
      character*20 spec_units(4*MXSPEC+2)
      character*20 spec_long_name(4*MXSPEC+2)
      character*60 spec_desc(4*MXSPEC+2)
      character*60 spec_coords(4*MXSPEC+2)
      character*4 ifile(10),note(60)
      integer enddate
      character*10 avfil,tmpnam
c
      integer ndps, ninline
      character*4 dpspec(10,4*MXSPEC+2)
c
      integer istrln   !--- external function
c
      data avfil /'AVERAGE   '/
      data izero,ione /0,1/
      data zero /0./
c
c-----Entry point
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
      read(runmsg(1:60),'(60a1)') (note(n),n=1,60)
      read(avfil,'(10a1)') (ifile(n),n=1,10)
      do l = 1,ndepspc
        if( ldepmap(l) .GT. 0 ) then
           tmpnam(3:10) = spname(ldepmap(l))(1:8)
c
           tmpnam(1:2) = 'R_'
           depsp(l) = tmpnam
           read(tmpnam,'(10a1)') (dpspec(n,l),n=1,10)
c
           tmpnam(1:2) = 'D_'
           depsp(ndepspc+l) = tmpnam
           read(tmpnam,'(10a1)') (dpspec(n,ndepspc+l),n=1,10)
c
           tmpnam(1:2) = 'W_'
           depsp(2*ndepspc+l) = tmpnam
           read(tmpnam,'(10a1)') (dpspec(n,2*ndepspc+l),n=1,10)
c
           tmpnam(1:2) = 'L_'
           depsp(3*ndepspc+l) = tmpnam
           read(tmpnam,'(10a1)') (dpspec(n,3*ndepspc+l),n=1,10)
        endif
      enddo
c
c-----If in-line Ix emissions are active, add species I2 and HOI
c
      ndps = 4*ndepspc
      ninline = 0
      if (lixemis) then
        ninline = 2
        ndps = ndps + 2
        tmpnam = 'I2_EMIS   '
        depsp(4*ndepspc+1) = tmpnam
        read(tmpnam,'(10a1)') (dpspec(n,4*ndepspc+1),n=1,10)
        tmpnam = 'HOI_EMIS  '
        depsp(4*ndepspc+2) = tmpnam
        read(tmpnam,'(10a1)') (dpspec(n,4*ndepspc+2),n=1,10)
      endif
c
c-----Master grid header
c
      if( lcdfout) then
          call ncf_set_vars_base()
          call ncf_set_specatt_dep(spec_units,spec_long_name,spec_desc,
     &                                                        spec_coords)
      endif
      if( .NOT. lcdfout ) then
         write(idep(1)) ifile,note,itzon,ndps,idat1,tim1,idat2,tim2
         write(idep(1)) plon,plat,iutm,orgx,orgy,dx,dy,ncol(1),nrow(1),
     &                  ione,iproj,izero,t1,t2,zero
         write(idep(1)) ione,ione,ncol(1),nrow(1)
         write(idep(1)) ((dpspec(n,l),n=1,10),l=1,ndps)
      else
          action = 'Writing master grid output deposition file.'
          call ncf_set_global('DEPOSITION',1,begdate,begtim,enddate,endtim,
     &                                                       ione,ndps)
          call ncf_wrt_dim(action,idep(1),1,ncol(1),nrow(1),ione,ndps)
          call ncf_wrt_global(action,idep(1),ndps,depsp,.FALSE.)
          call ncf_wrt_vars_base(action,idep(1))
          call ncf_wrt_vars_species(action,idep(1),ncol(1),nrow(1),ione,
     &                 ndps,depsp,spec_units,spec_long_name,
     &                                             spec_desc,spec_coords,4)
          call ncf_enddef_file(action,idep(1))
          call ncf_wrt_data_grid(action,idep(1),1,ncol(1),nrow(1),
     &                           orgx,orgy,dx,dy,ione,cellat(iptr2d(1)),
     &                                   cellon(iptr2d(1)),topo(iptr2d(1)))
      endif
      if( ngrid .EQ. 1 ) goto 9999
c
c-----Fine grid header
c
      do i = 2,ngrid
        dxf   = dx/float(meshold(i))
        dyf   = dy/float(meshold(i))
        orgxf = orgx + dx*(inst1(i)-1)
        orgyf = orgy + dy*(jnst1(i)-1)
        if( .NOT. lcdfout ) then
            write(idep(i)) ifile,note,itzon,ndps,idat1,tim1,idat2,tim2
            write(idep(i)) plon,plat,iutm,orgxf,orgyf,dxf,dyf,ncol(i)-2,
     &                 nrow(i)-2,ione,iproj,izero,t1,t2,zero
            write(idep(i)) ione,ione,ncol(i)-2,nrow(i)-2
            write(idep(i)) ((dpspec(n,l),n=1,10),l=1,ndps)
        else
          write(action,'(A,I2)') 'Writing output deposition file for grid: ',i
          call ncf_set_global('DEPOSITION',i,begdate,begtim,enddate,endtim,
     &                                                       ione,ndps)
          call ncf_wrt_dim(action,idep(i),i,ncol(i),nrow(i),ione,ndps)
          call ncf_wrt_global(action,idep(i),ndps,depsp,.FALSE.)
          call ncf_wrt_vars_base(action,idep(i))
          call ncf_wrt_vars_species(action,idep(i),ncol(i),nrow(i),ione,
     &                        ndps,depsp,spec_units,spec_long_name,
     &                                              spec_desc,spec_coords,4)
          call ncf_enddef_file(action,idep(i))
          call ncf_wrt_data_grid(action,idep(i),i,ncol(i),nrow(i),
     &                           orgxf,orgyf,dxf,dyf,ione,cellat(iptr2d(i)),
     &                                   cellon(iptr2d(i)),topo(iptr2d(i)))
        endif
      enddo
c
 9999 continue
      return
      end
