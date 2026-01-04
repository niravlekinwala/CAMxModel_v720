C**** NCF_READAR
c
      subroutine ncf_readar(igrd,ifile,iunit,ncol,nrow,nlay,nlay_ems,
     &                                                aremis,numspcs)
      use chmstry
      use camxcom
      use filunit
      implicit none
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c     NCF_READAR reads the NetCDF gridded emissions file for 
c     the given grid and loads the data into global arrays for emissions
c     injection.
c
c      Copyright 1996 - 2022
c     Ramboll
c
c     Input arguments:
c        igrd             grid index
c        ifile            index of file in list for this grid
c        ncol             number of columns
c        nrow             number of rows
c        nlay             number of layers in grid
c        nlay_ems         number of layers in emissions foles
c        iunit            file unit for area emissions file
c        numspcs          number of model species
c
c     Output arguments:
c        aremis              area emissions rate (mole/s)
c
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
      include 'flags.inc'
      include 'vbs.inc'
      include 'netcdf.inc'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      integer igrd
      integer ifile 
      integer iunit
      integer ncol
      integer nrow
      integer nlay
      integer nlay_ems
      integer numspcs
      real    aremis(ncol,nrow,nlay_ems,numspcs)
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
      character*10  this_var
      integer       this_tstep, ierr, nlays_in, ispc, i, j, k, l, ibin
      integer       this_varid, data_start(4), data_count(4)
      integer       this_time_tflag, this_time_etflag, this_dimid
      integer       buffer_offset, itmp, jtmp
c
      real, allocatable, dimension(:,:,:) ::  argrid
      real, allocatable, dimension(:,:,:) ::  poa_gv_em
      real, allocatable, dimension(:,:,:) ::  poa_dv_em
      real, allocatable, dimension(:,:,:) ::  poa_mc_em
      real, allocatable, dimension(:,:,:) ::  poa_op_em
      real, allocatable, dimension(:,:,:) ::  poa_bb_em

c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c  --- set offset based on buffer or no buffer ---
c
      buffer_offset = buffer_offset_iarem(igrd,ifile)
c
      if( .NOT. is_netcdf_iarem(igrd,ifile) ) goto 9999
      write(action,'(2A,I3,A,I3)') 'Reading the NetCDF ',
     &                              'gridded emissions file. Grid: ',
     &                                           igrd,' File: ',ifile
c
c  ---- get the index for timestep containing this time ----
c
      this_tstep = ncf_get_tstep(iunit,action,date,time,
     &              this_time_tflag,this_time_etflag,le1day,.TRUE.)
c
c   ---- get the number of layers in this file ---
c
c  --- get the id for the layer dimension ---
c
      ierr = nf_inq_dimid(iunit, "LAY", this_dimid  )
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_inq_dimlen(iunit,this_dimid,nlays_in)
      if( ierr .NE. NF_NOERR ) goto 7000
c
c  ---- allocate the temporary array ---
c
      allocate( argrid(ncol-2*buffer_offset,nrow-2*buffer_offset,nlays_in) )
      allocate( poa_gv_em(ncol,nrow,nlays_in) )
      allocate( poa_dv_em(ncol,nrow,nlays_in) )
      allocate( poa_mc_em(ncol,nrow,nlays_in) )
      allocate( poa_op_em(ncol,nrow,nlays_in) )
      allocate( poa_bb_em(ncol,nrow,nlays_in) )
c
c   --- if number of layers is not 1 make sure it matches grid ---
c
      if( nlays_in .NE. 1 .AND. nlays_in .NE. nlay ) goto 7002
c
c   --- set the indexes for what to read ---
c
      data_start(1) = 1
      data_count(1) = ncol-2*buffer_offset
      data_start(2) = 1
      data_count(2) = nrow-2*buffer_offset
      data_start(3) = 1
      data_count(3) = nlays_in
      data_start(4) = this_tstep
      data_count(4) = 1
c
c  ---- loop over the list of emissions spcies ---
c
      do ispc=1,nemspc
c
c  --- load data into temporary arrray ---
c
          this_var = emspcname(ispc)
          ierr = nf_inq_varid(iunit, this_var, this_varid)
          if( ierr .NE. NF_NOERR) cycle
          ierr = nf_get_vara_real(iunit,this_varid,data_start,
     &                                               data_count,argrid)
          if( ierr .NE. NF_NOERR) goto 7001
c
c-----Put the data into the global array, if buffer cells are not
c     in file make edge cells same as neigboring cell ----
c
          do k=1,nlays_in
             do j=1+buffer_offset,nrow-buffer_offset
                jtmp = j-buffer_offset
                do i=1+buffer_offset,ncol-buffer_offset
                   itmp = i-buffer_offset
                   if (argrid(itmp,jtmp,k).lt.0.) then
                     write(iout,'(//,a)') 'ERROR in NCF_READAR:'
                     write(iout,'(A)') action(:istrln(action))
                     write(iout,'(a)') 'Negative emissions read.'
                     write(iout,'(a,3i3,a,i3)') 'Cell (i,j,k): ',
     &                           i+buffer_offset,j+buffer_offset,k,' Species: ',ispc
                     call camxerr()
                   endif
                   aremis(i,j,k,ispc) = aremis(i,j,k,ispc) + 
     &                                              argrid(itmp,jtmp,k)/(60.*dtems) 
                enddo
             enddo
          enddo
      enddo
c
c ---- handle VBS as special case ---
c
      if( lvbs .AND. LVBSPREPROC ) then
          this_var = 'POA_OP'
          ierr = nf_inq_varid(iunit, this_var, this_varid)
          if( ierr .EQ. NF_NOERR) then
             ierr = nf_get_vara_real(iunit,this_varid,data_start,
     &                                             data_count,argrid)
             if( ierr .NE. NF_NOERR) goto 7001
             do k=1,nlays_in
                do j=1+buffer_offset,nrow-buffer_offset
                   jtmp = j-buffer_offset
                   do i=1+buffer_offset,ncol-buffer_offset
                      itmp = i-buffer_offset
                      if (argrid(i,j,k).lt.0.) then
                        write(iout,'(//,a)') 'ERROR in NCF_READAR:'
                        write(iout,'(A)') action(:istrln(action))
                        write(iout,'(a)') 'Negative emissions read.'
                        write(iout,'(a,3i3,a,i3)') 'Cell (i,j,k): ',
     &                        i+buffer_offset,j+buffer_offset,k,' Species: ',ispc
                        call camxerr()
                      endif
                      poa_op_em(i,j,k) = argrid(itmp,jtmp,k)
                   enddo
                enddo
             enddo
          endif
c
          this_var = 'POA_GV'
          ierr = nf_inq_varid(iunit, this_var, this_varid)
          if( ierr .EQ. NF_NOERR) then
             ierr = nf_get_vara_real(iunit,this_varid,data_start,
     &                                             data_count,argrid)
             if( ierr .NE. NF_NOERR) goto 7001
             do k=1,nlays_in
                do j=1+buffer_offset,nrow-buffer_offset
                   jtmp = j-buffer_offset
                   do i=1+buffer_offset,ncol-buffer_offset
                      itmp = i-buffer_offset
                      if (argrid(i,j,k).lt.0.) then
                        write(iout,'(//,a)') 'ERROR in NCF_READAR:'
                        write(iout,'(A)') action(:istrln(action))
                        write(iout,'(a)') 'Negative emissions read.'
                        write(iout,'(a,3i3,a,i3)') 'Cell (i,j,k): ',
     &                          i+buffer_offset,j+buffer_offset,k,' Species: ',ispc
                        call camxerr()
                      endif
                      poa_gv_em(i,j,k) = argrid(itmp,jtmp,k)
                   enddo
                enddo
             enddo
          endif
c
          this_var = 'POA_DV'
          ierr = nf_inq_varid(iunit, this_var, this_varid)
          if( ierr .EQ. NF_NOERR) then
             ierr = nf_get_vara_real(iunit,this_varid,data_start,
     &                                             data_count,argrid)
             if( ierr .NE. NF_NOERR) goto 7001
             do k=1,nlays_in
                do j=1+buffer_offset,nrow-buffer_offset
                   jtmp = j-buffer_offset
                   do i=1+buffer_offset,ncol-buffer_offset
                      itmp = i-buffer_offset
                      if (argrid(i,j,k).lt.0.) then
                        write(iout,'(//,a)') 'ERROR in NCF_READAR:'
                        write(iout,'(A)') action(:istrln(action))
                        write(iout,'(a)') 'Negative emissions read.'
                        write(iout,'(a,3i3,a,i3)') 'Cell (i,j,k): ',
     &                      i+buffer_offset,j+buffer_offset,k,' Species: ',ispc
                        call camxerr()
                      endif
                      poa_dv_em(i,j,k) = argrid(itmp,jtmp,k)
                   enddo
                enddo
             enddo
          endif
c
          this_var = 'POA_MC'
          ierr = nf_inq_varid(iunit, this_var, this_varid)
          if( ierr .EQ. NF_NOERR) then
             ierr = nf_get_vara_real(iunit,this_varid,data_start,
     &                                             data_count,argrid)
             if( ierr .NE. NF_NOERR) goto 7001
             do k=1,nlays_in
                do j=1+buffer_offset,nrow-buffer_offset
                   jtmp = j-buffer_offset
                   do i=1+buffer_offset,ncol-buffer_offset
                      itmp = i-buffer_offset
                      if (argrid(i,j,k).lt.0.) then
                        write(iout,'(//,a)') 'ERROR in NCF_READAR:'
                        write(iout,'(A)') action(:istrln(action))
                        write(iout,'(a)') 'Negative emissions read.'
                        write(iout,'(a,3i3,a,i3)') 'Cell (i,j,k): ',
     &                       i+buffer_offset,j+buffer_offset,k,' Species: ',ispc
                        call camxerr()
                      endif
                      poa_mc_em(i,j,k) = argrid(itmp,jtmp,k)
                   enddo
                enddo
             enddo
          endif
c
          this_var = 'POA_BB'
          ierr = nf_inq_varid(iunit, this_var, this_varid)
          if( ierr .EQ. NF_NOERR) then
             ierr = nf_get_vara_real(iunit,this_varid,data_start,
     &                                             data_count,argrid)
             if( ierr .NE. NF_NOERR) goto 7001
             do k=1,nlays_in
                do j=1+buffer_offset,nrow-buffer_offset
                   jtmp = j-buffer_offset
                   do i=1+buffer_offset,ncol-buffer_offset
                      itmp = i-buffer_offset
                      if (argrid(i,j,k).lt.0.) then
                        write(iout,'(//,a)') 'ERROR in NCF_READAR:'
                        write(iout,'(A)') action(:istrln(action))
                        write(iout,'(a)') 'Negative emissions read.'
                        write(iout,'(a,3i3,a,i3)') 'Cell (i,j,k): ',
     &                           i+buffer_offset,j+buffer_offset,k,' Species: ',ispc
                        call camxerr()
                      endif
                      poa_bb_em(i,j,k) = argrid(itmp,jtmp,k)
                   enddo
                enddo
             enddo
          endif
c
          do ibin= 0,NVOLBIN
            do ispc=1,nemspc
               if( emspcname(ispc) .EQ. spname(kpap_c(ibin)) ) then
                  aremis(:,:,:,ispc) = aremis(:,:,:,ispc) +
     &                       ( poa_op_em(:,:,:) * poa_op_ef(ibin)
     &                       + poa_gv_em(:,:,:) * poa_gv_ef(ibin)
     &                       + poa_dv_em(:,:,:) * poa_dv_ef(ibin) ) /
     &                                              (60.*dtems)
              endif
           enddo
         enddo
         do ibin= 0,NVOLBIN
           do ispc=1,nemspc
              if( emspcname(ispc) .EQ. spname(kpcp_c(ibin)) ) then
                     aremis(:,:,:,ispc) = aremis(:,:,:,ispc) +
     &                        poa_mc_em(:,:,:) * poa_mc_ef(ibin) 
              endif
           enddo
         enddo
         do ibin= 0,NVOLBIN
           do ispc=1,nemspc
              if( emspcname(ispc) .EQ. spname(kpfp_c(ibin)) ) then
                aremis(:,:,:,ispc) = aremis(:,:,:,ispc) +
     &                        poa_bb_em(:,:,:) * poa_bb_ef(ibin)
              endif
           enddo
         enddo
c
      endif
c
c  ---- deallocate the temporary array ---
c
      deallocate( argrid )
      deallocate( poa_gv_em )
      deallocate( poa_dv_em )
      deallocate( poa_mc_em )
      deallocate( poa_op_em )
      deallocate( poa_bb_em )
c
c  ---- write message to out file ---
c
      write(iout,'(a35,i5,a,f7.0,i8.5)') 'Read area source file ',
     &                                         ifile,' at ',time,date
c
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,a)') 'ERROR in NCF_READAR:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(A)') 'Cannot find the dimension id for layer dimension.'
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7001 continue
      write(iout,'(//,a)') 'ERROR in NCF_READAR:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot read data for variable: ',
     &                                      this_var(:istrln(this_var))
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7002 continue
      write(iout,'(//,a)') 'ERROR in NCF_READAR:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(A)') 'Number of layers does not match: '
      write(iout,'(A,I4)') 'User supplied: ',nlay
      write(iout,'(A,I4)') 'Value in file: ',nlays_in
      call camxerr()
c
 9999 continue
      return
      end
