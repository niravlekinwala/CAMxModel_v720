c**** NCF_SET_SPECATT_SRFMOD
c
      subroutine ncf_set_specatt_srfmod(spec_units,spec_long_name,spec_desc,
     &                                          spec_coords,nspcs,specname)
      use chmstry
      implicit none
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routine sets the species varaible attributes for the surface
c   model output NetCDF file
c
c      Copyright 1996 - 2022
c     Ramboll
c      Argument description:
c       Inputs:
c            nspcs          I number of species in list
c            specname       C array of species names 
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
       integer      nspcs
       character*10 specname(*)
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
      integer l
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
      do l=1,nspcs
        spec_long_name(l) = "Mass on soilsnow"
        spec_units(l) = "mol ha-1"
        spec_desc(l) = specname(l)(2:istrln(specname(l)))//" mass on soil or snow"
        if( specname(l)(1:1) .EQ. 'V' ) then
           spec_long_name(l) = "Mass on veg"
           spec_desc(l) = specname(l)(2:istrln(specname(l)))//" mass on vegetation"
        endif
        spec_coords(l) = "latitude longitude"
      enddo
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
