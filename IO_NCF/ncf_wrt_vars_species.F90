!**** NCF_WRT_VARS_SPECIES
!
      subroutine ncf_wrt_vars_species(action,iounit,numcols,numrows,numlays, &
                       nspcs,spnames,spunits,splong,spdesc,spcoord,num_dims)
      use camx_includes
      use chmstry
      use filunit
      use ncf_iomod
      implicit none
!
!----CAMx v7.20 220430
!
!-----------------------------------------------------------------------
!    Description:
!-----------------------------------------------------------------------
!
!   This routine writes the variable definitions and descriptions to 
!    the NetCDF file
!
!      Copyright 1996 - 2022
!     Ramboll
!      Argument description:
!       Inputs:
!           action   C  name of file to open
!           iounit   I  NetCDF file ID of file
!           numcols  I  number of columns in the grid
!           numrows  I  number of columns in the grid
!           numlays  I  number of columns in the grid
!           numcols  I  number of columns in the grid
!           nspcs    I  number of species in the file
!           spnames  C  names of each species
!           spunits  C  units for each species
!           splong   C  long name for each species
!           spdesc   C  description of each species
!           spcoord  C  description of each species
!           num_dims I  number of dimensions in variable
!       Outputs:
!
!-----------------------------------------------------------------------
!    LOG:
!-----------------------------------------------------------------------
!
!     02/20/17   --gwilson--    Original development
!
!-----------------------------------------------------------------------
!    Include files:
!-----------------------------------------------------------------------
!
      include 'netcdf.inc'
!
!-----------------------------------------------------------------------
!    Argument declarations:
!-----------------------------------------------------------------------
!
      character*(*) action
      integer       iounit
      integer       numcols
      integer       numrows
      integer       numlays
      integer       nspcs
      character*(*) spnames(nspcs)
      character*(*) spunits(nspcs)
      character*(*) splong(nspcs)
      character*(*) spdesc(nspcs)
      character*(*) spcoord(nspcs)
      integer       num_dims
!
!-----------------------------------------------------------------------
!    External functions:
!-----------------------------------------------------------------------
!
      integer istrln
!
!-----------------------------------------------------------------------
!    Local variables:
!-----------------------------------------------------------------------
!
      character*14 this_species
      integer      ierr, ispc, spec_dimid(4), spec_chunk(4), spec_varid
!
!-----------------------------------------------------------------------
!    Entry point:
!-----------------------------------------------------------------------
!
      spec_chunk(1) = INT(REAL(numcols)/NCF_CHUNK_SIZE_VAR_X)
      if( spec_chunk(1) .GT. numcols ) goto 7002
      spec_chunk(2) = INT(REAL(numrows)/NCF_CHUNK_SIZE_VAR_Y)
      if( spec_chunk(2) .GT. numrows ) goto 7003
      spec_chunk(3) = MAX(1,INT(numlays/2))
      spec_chunk(4) = 1
!
      spec_dimid(1) = ncf_col_dimid
      spec_dimid(2) = ncf_row_dimid
      spec_dimid(3) = ncf_lay_dimid
      spec_dimid(4) = ncf_tstep_dimid
!
      if( num_dims .EQ.  3) then
         spec_dimid(3) = ncf_tstep_dimid
         spec_chunk(3) = 1
      endif
      do ispc=1,nspcs
!
!  --- define everything for this species ---
!
        this_species = spnames(ispc)(:istrln(spnames(ispc)))
        ncf_species_units = spunits(ispc)
        ncf_species_long_name = splong(ispc)
        ncf_species_var_desc = spdesc(ispc)
        ncf_species_coordinates = spcoord(ispc)
!
!  --- define the variable ---
!
        ierr = nf_def_var(iounit,this_species(:istrln(this_species)),  &
                          NF_FLOAT, num_dims, spec_dimid, spec_varid)
        if( ierr .NE. NF_NOERR ) goto 7000
!
!  --- add the attributes ---
!
        ierr = nf_put_att_text(iounit,spec_varid,'long_name', &
                istrln(ncf_species_long_name),ncf_species_long_name)
        if( ierr .NE. NF_NOERR ) goto 7000

        ierr = nf_put_att_text(iounit,spec_varid,'units', &
                        istrln(ncf_species_units),ncf_species_units)
        if( ierr .NE. NF_NOERR ) goto 7000

        ierr = nf_put_att_text(iounit,spec_varid,'var_desc', &
                  istrln(ncf_species_var_desc),ncf_species_var_desc)
        if( ierr .NE. NF_NOERR ) goto 7000

        ierr = nf_put_att_text(iounit,spec_varid,'coordinates', &
                 istrln(ncf_species_coordinates),ncf_species_coordinates)
        if( ierr .NE. NF_NOERR ) goto 7000
!
!  --- chunkers go here ... Irie! ----
!
#ifdef CHUNK
        if( ncf_compress ) then
          ierr = nf_def_var_chunking(iounit, spec_varid, &
                                            NF_CHUNKED, spec_chunk)
          if( ierr .NE. NF_NOERR ) goto 7001
!
          ierr = nf_def_var_deflate(iounit, spec_varid, NCF_SHUFFLE,  &
                                      NCF_DEFLATE, NCF_DEFLATE_LEVEL )
          if( ierr .NE. NF_NOERR ) goto 7001
        endif
#endif
!
      enddo
!
      goto 9999
!
!-----------------------------------------------------------------------
!    Error messages:
!-----------------------------------------------------------------------
!
 7000 continue
      write(iout,'(//,A)') 'ERROR in NCF_WRT_VARS_SPECIES:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot create file variable for species: ', &
                                    this_species(:istrln(this_species))
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
!
 7001 continue
      write(iout,'(//,A)') 'ERROR in NCF_WRT_VARS_SPECIES:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot set chunk parameters for species: ', &
                                    this_species(:istrln(this_species))
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
!
 7002 continue
      write(iout,'(//,A)') 'ERROR in NCF_IODAT.INC:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(A)') 'Cannot set chunk parameters for this file.'
      write(iout,'(2A,/,A)') 'The NCF_CHUNK_SIZE_VAR_X parameter causes the ', &
           'chunk value to be larger','than the number of columns.'
      write(iout,'(A,I3,A,I4)') 'Number of columns: ',numcols
      call camxerr()
!
 7003 continue
      write(iout,'(//,A)') 'ERROR in NCF_IODAT.INC:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(A)') 'Cannot set chunk parameters for this file.'
      write(iout,'(2A,/,A)') 'The NCF_CHUNK_SIZE_VAR_Y parameter causes the ', &
           'chunk value to be larger','than the number of rows.'
      write(iout,'(A,I3,A,I4)') 'Number of rows: ',numrows
      call camxerr()
!
!-----------------------------------------------------------------------
!    Return point:
!-----------------------------------------------------------------------
!
 9999 continue
      return
      end
