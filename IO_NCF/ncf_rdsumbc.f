C**** NCF_RDSUMBC
c
      subroutine ncf_rdsumbc(idtnow,timnow,jdlast,ttlast,
     &                                   ncolx,nrowy,nlays,consum)
      use filunit
      use grid
      use chmstry
      use bndary
      use camxcom
      use tracer
      implicit none
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routines reads the boundary condition data for the entire
c   simulation period and sums up the BC tracers for purposes of
c   calculating koh reactivity for these tracers.
c
c      Copyright 1996 - 2022
c     Ramboll
c
c      Argument description:
c       Inputs:
c         idtnow  I date of beginning of simulation
c         timnow  R time of beginning of simulation
c         jdlast  I date of end of sumulation
c         ttlast  R time of end of sumulation
c         ncolx   I number of columns  in master grid 
c         nrowy   I number of rows in master grid 
c         nlays   I number of layers in master grid 
c       Inputs:
c         consum  R sum of tracer concentrations
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
      include 'netcdf.inc'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      integer idtnow
      real    timnow
      integer jdlast
      real    ttlast
      integer ncolx
      integer nrowy
      integer nlays
      real    consum(MXSPEC,0:IDXBTP)
c
c-----------------------------------------------------------------------
c    External functions:
c-----------------------------------------------------------------------
c
      integer istrln
      integer ncf_get_tstep
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      character*200 action
      character*4   iname(10)
      character*4   fname(10),note(60)
      character*10  filnote
      logical       lfound,lfrst
      integer       ncola, nrowa, nlaya, nvar, this_varid, ierr, ispc, j
      integer       this_time_tflag, this_time_etflag, this_tstep
      integer       data_start(4), data_count(4), ibgdhp, idum, izcl
      real          btimhp

      real, allocatable, dimension(:)     :: height
      real, allocatable, dimension(:)     :: depth
      real, allocatable, dimension(:,:,:) :: bctmp
c
      real conwst(MXLAYER,MXCELLS)
      real conest(MXLAYER,MXCELLS)
      real consth(MXLAYER,MXCELLS)
      real connth(MXLAYER,MXCELLS)
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c  --- allocate the local arrays ---
c
      ncola = maxval( ncol(1:ngrid) )
      nrowa = maxval( nrow(1:ngrid) )
      nlaya = maxval( nlay(1:ngrid) )
      allocate( height(ncola*nrowa*nspec) )
      allocate( depth(ncola*nrowa*nspec) )
      allocate( bctmp(ncol(1),nrow(1),nlay(1)) )
c
c   --- set the indexes for what to read ---
c
      data_start(1) = 1
      data_count(1) = ncol(1)
      data_start(2) = 1
      data_count(2) = nrow(1)
      data_start(3) = 1
      data_count(3) = nlay(1)
c
      ibgdhp = 0
      btimhp = 0 
c
c   --- initialize the array to zero ---
c
      consum = 0.
c
c   --- read the 3D met file to get the layer heights ----
c
      if( .NOT. is_netcdf_i3dmet(1) ) then
         rewind(i3dmet(1))
         read(i3dmet(1),ERR=7001) fname,note,idum,nvar
         write(filnote,'(10a1)') (note(j),j=1,10)
         if (filnote.ne.'3DMET     ') then
           write(iout,'(//,a)') 'ERROR in RDSUMBC:'
           write(iout,*) 'Mismatch in file type.'
           write(iout,*) 'File is:'
           write(iout,*) filnote
           write(iout,*) 'Expecting:'
           write(iout,*) '3DMET'
           write(iout,*)
           call camxerr()
         endif
         do j = 1,3
           read(i3dmet(1),ERR=7001)
         enddo
      endif
c
  333 continue
      if( is_netcdf_i3dmet(1) ) then
         call ncf_getdepth(nvar,ncolx,nrowy,nlays,ibgdhp,idtnow,btimhp,timnow,
     &              i3dmet(1),height,depth)
      else
         call getdepth(nvar,ncolx,nrowy,nlays,ibgdhp,idtnow,btimhp,timnow,
     &              i3dmet(1),height,depth)
      endif
c
c   --- get the time step for this hour ---
c
      action = 'Reading lateral boundary conditions file.'
c
c   --- get the timestep that encompasses this date/time ---
c
      this_tstep = ncf_get_tstep(ibc,action,idtnow,timnow*100.,
     &          bnddate_time_tflag,bnddate_time_etflag,.FALSE.,.TRUE.)
      data_start(4) = this_tstep
      data_count(4) = 1
c
c   --- read boundary conditions for this hour ----
c
      do ispc=1,nspec
c
c  --- skip if species not needed ---
c
         if( lbcmap(ispc) .LE. 0 .OR. .NOT. 
     &              (lvocsp(ispc) .OR. lvocsoa(ispc)) ) cycle
c
c  ---- load this species into local array ---
c
         this_varid = lbcmap(ispc)
         ierr = nf_get_vara_real(ibc,this_varid,data_start,data_count,bctmp)
         if( ierr .NE. NF_NOERR) then
           write(iout,'(//,a)') 'ERROR in NCF_RDSUMBC:'
           write(iout,*)'Cannot read boundary conditions data for speices: ',
     &                                   spname(ispc)(:istrln(spname(ispc)))
           write(iout,'(A,I5)') 'NetCDF error code: ',ierr
           call camxerr()
         endif
         do izcl=1,nlays
            do  j=1,nrowy
               conwst(izcl,j) = bctmp(1,j,izcl)
               conest(izcl,j) = bctmp(ncolx,j,izcl)
            enddo
            do  j=1,ncolx
               consth(izcl,j) = bctmp(j,1,izcl)
               connth(izcl,j) = bctmp(j,nrowy,izcl)
            enddo
         enddo
c
c   --- if the species is a not modelled or not a VOC species skip it ---
c
        call sumwt4(ispc,ncolx,nrowy,nlays,depth,conwst,conest,
     &                                          consth,connth,consum)
c
c   --- next species ---
c
      enddo
c
c   --- chack date and time, if it is still in the episode
c       go back to read the next hour ----
c
      timnow = timnow + 1.0
      if( timnow .GT. 24.0 ) then
          idtnow = idtnow + 1
          timnow = 0.0
          if( MOD(idtnow,1000) .GT. 365 ) then
             if( MOD(INT(idtnow/1000),4) .EQ. 0 ) then
                if( MOD(idtnow,1000) .EQ. 367 )
     &                     idtnow = (INT(idtnow/1000)+1)*1000 + 1
             else
                idtnow = (INT(idtnow/1000)+1)*1000 + 1
             endif
          endif
      endif
      if( idtnow .GT. jdlast .OR. 
     &          (idtnow .EQ. jdlast .AND. timnow .GE. ttlast) ) goto 444
c
c  --- process next hour ----
c
      goto 333
c
c  --- all concentrations are summed, return ---
c
  444 continue
c
c  --- deallocate the local arrays ---
c
      deallocate( height )
      deallocate( depth  )
      deallocate( bctmp  )
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7001 continue
      write(iout,'(//,a)') 'ERROR in NCF_RDSUMBC:'
      write(iout,'(/,1X,A)') 'Reading 3D Met file'
      call camxerr()
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
