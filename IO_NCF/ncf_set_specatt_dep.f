c**** NCF_SET_SPECATT_DEP
c
      subroutine ncf_set_specatt_dep(spec_units,spec_long_name,spec_desc,
     &                                                        spec_coords)
      use chmstry
      implicit none
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routine sets the species varaible attributes for the deposition
c   output NetCDF file
c
c      Copyright 1996 - 2022
c     Ramboll Environ
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
      include 'flags.inc'
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
      integer l
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
      do l = 1,ndepspc
        if( ldepmap(l) .GT. 0 ) then
c
           spec_units(l) = "m s-1"
           spec_long_name(l) = "Dry dep vel"
           spec_desc(l) = spname(ldepmap(l))(:istrln(spname(ldepmap(l))))
     &                              //" dry deposition velocity"
           spec_coords(l) = "latitude longitude"
c
           if( lgas(ldepmap(l)) ) then
             spec_units(ndepspc+l) = "mol ha-1"
           else
             spec_units(ndepspc+l) = "g ha-1"
           endif
           spec_long_name(ndepspc+l) = "Dry dep mass"
           spec_desc(ndepspc+l) = spname(ldepmap(l))(:istrln(spname(ldepmap(l))))
     &                              //" dry deposited mass"
           spec_coords(ndepspc+l) = "latitude longitude"
c
           if( lgas(ldepmap(l)) ) then
             spec_units(2*ndepspc+l) = "mol ha-1"
           else
             spec_units(2*ndepspc+l) = "g ha-1"
           endif
           spec_long_name(2*ndepspc+l) = "Wet dep mass"
           spec_desc(2*ndepspc+l) = spname(ldepmap(l))(:istrln(spname(ldepmap(l))))
     &                              //" wet deposited mass"
           spec_coords(2*ndepspc+l) = "latitude longitude"
c
           if( lgas(ldepmap(l)) ) then
             spec_units(3*ndepspc+l) = "mol l-1"
           else
             spec_units(3*ndepspc+l) = "g l-1"
           endif
           spec_long_name(3*ndepspc+l) = "Liquid conc"
           spec_desc(3*ndepspc+l) = spname(ldepmap(l))(:istrln(spname(ldepmap(l))))
     &                              //" liquid concentration"
           spec_coords(3*ndepspc+l) = "latitude longitude"
        endif
      enddo
c
      if (lixemis) then
        spec_units(4*ndepspc+1) = "mol m-2"
        spec_long_name(4*ndepspc+1) = "Oceanic halogen emissions"
        spec_desc(4*ndepspc+1) = "I2 oceanic emissions"
        spec_coords(4*ndepspc+1) = "latitude longitude"
        spec_units(4*ndepspc+2) = "mol m-2"
        spec_long_name(4*ndepspc+2) = "Oceanic halogen emissions"
        spec_desc(4*ndepspc+2) = "HOI oceanic emissions"
        spec_coords(4*ndepspc+2) = "latitude longitude"
      endif

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
 
