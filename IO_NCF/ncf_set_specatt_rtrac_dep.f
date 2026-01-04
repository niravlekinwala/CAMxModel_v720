c**** NCF_SET_SPECATT_RTRAC_DEP
c
      subroutine ncf_set_specatt_rtrac_dep(spec_units,spec_long_name,
     &                                         spec_desc,spec_coords)
      use chmstry
      use tracer
      implicit none
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routine sets the species varaible attributes for the RTRAC
c   deposition output NetCDF file
c
c      Copyright 1996 - 2022
c     Ramboll
c      Argument description:
c       Inputs:
c       Outputs:
c            spec_units     C array of units for this each species
c            spec_long_name C array of "long names" for each each species
c            spec_desc      C array of desciption for this each species
c            spec_coords    C array of coordinates for this each species
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
      include 'ncf_iodat.inc'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
       character*20 spec_units(*)
       character*20 spec_long_name(*)
       character*60 spec_desc(*)
       character*60 spec_coords(*)
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
      integer nspcout, j
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
      nspcout = 0
      do j=1,ntotsp
          nspcout = nspcout + 1
          if( lsagas(j) ) then
             spec_units(nspcout) = "mol ha-1"
          else
             spec_units(nspcout) = "g ha-1"
          endif
          spec_long_name(nspcout) = "Dry dep mass"
          spec_desc(nspcout) = ptname(j)//" dry deposited mass"
          spec_coords(nspcout) = "latitude longitude"
      enddo
      do j=1,ntotsp
          nspcout = nspcout + 1
          if( lsagas(j) ) then
             spec_units(nspcout) = "mol ha-1"
          else
             spec_units(nspcout) = "g ha-1"
          endif
          spec_long_name(nspcout) = "Wet dep mass"
          spec_desc(nspcout) = ptname(j)//" wet deposited mass"
          spec_coords(nspcout) = "latitude longitude"
      enddo
c
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
 
