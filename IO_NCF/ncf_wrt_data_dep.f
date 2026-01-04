c**** NCF_WRT_DATA_DEP
c
      subroutine ncf_wrt_data_dep(action,iounit,igrd,
     &        num_cols,num_rows,num_lays_grid,
     &           nspcs_dep,nspcs_all,spec_name,vdep,depfld,height)
      use ncf_iomod
      use grid
      use chmstry
      use filunit
      implicit none
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routine writes the data for the timestep variables to the
c    NetCDF file
c
c      Copyright 1996 - 2022
c     Ramboll
c      Argument description:
c       Inputs:
c           action        C name of file to open
c           iounit        I NetCDF file ID of file
c           igrd          I grid number
c           num_cols      I number of columns in this grid
c           num_rows      I number of rows in this grid
c           num_lays_grid I number of layers in the array
c           nspcs_dep     I number of species in deposition array
c           nspcs_all     I number of species in deposition velocity array
c           spec_name     C names of the species 
c           vdep          R gridded array of deposition velocities
c           depfld        R gridded array of deposition velocities
c           height        R gridded array of layer heights
c       Outputs:
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
      include 'netcdf.inc'
      include 'flags.inc'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      character*(*) action
      integer       iounit
      integer       igrd
      integer       num_cols
      integer       num_rows
      integer       num_lays_grid
      integer       nspcs_dep
      integer       nspcs_all
      character*10  spec_name(*)
      real          vdep(num_cols,num_rows,nspcs_all)      
      real          depfld(num_cols,num_rows,nspcs_dep)      
      real          height(num_cols,num_rows,num_lays_grid)
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
      character*10 this_var
      integer      data_start(4), data_count(4), ispc, i, j, k
      integer      idxname, idxdata, itmp, jtmp, buffer_offset
      integer      data_start_z(4), data_count_z(4)
      integer      ierr, this_varid, this_species
      integer      ibegcel, iendcel, jbegcel, jendcel
      real,        allocatable, dimension(:,:,:) :: array_3d
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c  --- set the position in the NetCDF variable to write ---
c
      data_start(3) = 1
      data_count(3) = 1
      data_start(4) = ncf_cur_tstep
      data_count(4) = 1
      data_start_z(3) = 1
      data_count_z(3) = 1
      data_start_z(4) = ncf_cur_tstep
      data_count_z(4) = 1
c
c  --- allocate the array that will be used to write the data ---
c
      if( igrd .EQ. 1 ) then
         data_start(1) = 1
         data_count(1) = num_cols
         data_start(2) = 1
         data_count(2) = num_rows
         data_start_z(1) = 1
         data_count_z(1) = num_cols
         data_start_z(2) = 1
         data_count_z(2) = num_rows
         ibegcel = 1
         iendcel = num_cols
         jbegcel = 1
         jendcel = num_rows
         buffer_offset = 0
         allocate( array_3d(num_cols,num_rows,1) )
      else
         data_start(1) = 1
         data_count(1) = num_cols-2
         data_start(2) = 1
         data_count(2) = num_rows-2
         data_start_z(1) = 1
         data_count_z(1) = num_cols-2
         data_start_z(2) = 1
         data_count_z(2) = num_rows-2
         ibegcel = 2
         iendcel = num_cols-1
         jbegcel = 2
         jendcel = num_rows-1
         buffer_offset = 1
         allocate( array_3d(num_cols-2,num_rows-2,1) )
      endif
c
c  --- first do the z variable ---
c
      do j=jbegcel,jendcel
         jtmp = j-buffer_offset
         do i=ibegcel,iendcel
            itmp = i-buffer_offset
            array_3d(itmp,jtmp,1) = height(i,j,1)
         enddo
      enddo
c
c  --- get the id for the z variable and write the data ---
c
      this_var = 'z'
      ierr = nf_inq_varid(iounit,"z",this_varid)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_vara_real(iounit,this_varid,data_start_z,data_count_z,array_3d)
      if( ierr .NE. NF_NOERR ) goto 7001
c
c  --- loop over all variables in this file ---
c
      do ispc=1,ndepspc
         this_species = ldepmap(ispc)
c
c  --- get name for this variable, and get it's variable ID
c      in this file ---
c
         this_var = spec_name(ispc)
         ierr = nf_inq_varid(iounit,this_var,this_varid)
         if( ierr .NE. NF_NOERR ) goto 7000
c
c  --- load the data into the local array to write ---
c
         do j=jbegcel,jendcel
            jtmp = j-buffer_offset
            do i=ibegcel,iendcel
               itmp = i-buffer_offset
               array_3d(itmp,jtmp,1) = vdep(i,j,this_species)
            enddo
         enddo
c
         ierr = nf_put_vara_real(iounit,this_varid,data_start,data_count,array_3d)
         if( ierr .NE. NF_NOERR ) goto 7001
c
      enddo
c
c  --- loop over all variables in this file ---
c
      do ispc=1,3*ndepspc
         this_species = ispc + ndepspc
c
c  --- get name for this variable, and get it's variable ID
c      in this file ---
c
         this_var = spec_name(this_species)
         ierr = nf_inq_varid(iounit,this_var,this_varid)
         if( ierr .NE. NF_NOERR ) goto 7000
c
c  --- load the data into the local array to write ---
c
         do j=jbegcel,jendcel
            jtmp = j-buffer_offset
            do i=ibegcel,iendcel
               itmp = i-buffer_offset
               array_3d(itmp,jtmp,1) = depfld(i,j,ispc)
            enddo
         enddo
c
         ierr = nf_put_vara_real(iounit,this_varid,data_start,data_count,array_3d)
         if( ierr .NE. NF_NOERR ) goto 7001
c
      enddo
c
      if( lixemis ) then
         do ispc=1,2
            idxname = 4*ndepspc+ispc
            idxdata = 3*ndepspc+ispc
c
c  --- get name for this variable, and get it's variable ID
c      in this file ---
c
            this_var = spec_name(idxname)
            ierr = nf_inq_varid(iounit,this_var,this_varid)
            if( ierr .NE. NF_NOERR ) goto 7000
c
c  --- load the data into the local array to write ---
c
            do j=jbegcel,jendcel
               jtmp = j-buffer_offset
               do i=ibegcel,iendcel
                  itmp = i-buffer_offset
                  array_3d(itmp,jtmp,1) = depfld(i,j,idxdata)
               enddo
            enddo
c
            ierr = nf_put_vara_real(iounit,this_varid,data_start,data_count,array_3d)
            if( ierr .NE. NF_NOERR ) goto 7001
c
         enddo
      endif
c
c  ---- deallocate the array ---
c
      deallocate( array_3d )
c
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,a)') 'ERROR in NCF_WRT_DATA_DEP:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot find variable id for: ',
     &                                      this_var(:istrln(this_var))
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7001 continue
      write(iout,'(//,a)') 'ERROR in NCF_WRT_DATA_DEP:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot write data for the variable: ',
     &                                      this_var(:istrln(this_var))
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
 
