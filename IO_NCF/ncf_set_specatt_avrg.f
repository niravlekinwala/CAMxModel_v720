c**** NCF_SET_SPECATT_AVRG
c
      subroutine ncf_set_specatt_avrg(spec_units,spec_long_name,spec_desc,
     &                                                       spec_coords)
      use chmstry
      implicit none
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routine sets the species varaible attributes for the average
c   concentration output NetCDF file
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
      do l = 1,navspc
        if(lavmap(l) .GT. 0 ) then
          if( lgas(lavmap(l)) .AND. lgas_out_ppm ) then
             spec_units(l) = "ppmv"
          else
             spec_units(l) = "micrograms m-3"
          endif
          spec_long_name(l) = spavnam(l)
          spec_desc(l) = spavnam(l)(:istrln(spavnam(l)))//" air concentration"
          spec_coords(l) = "latitude longitude"
        endif
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
 
