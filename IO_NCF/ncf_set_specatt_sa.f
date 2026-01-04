c**** NCF_SET_SPECATT_SA
c
      subroutine ncf_set_specatt_sa(spec_units,spec_long_name,
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
c   This routine sets the species varaible attributes for the SA
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
      integer nspcout, j
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
      nspcout = 0
      do j=1,ntotsp
        if( loutsa(j) ) then
          nspcout = nspcout + 1
          if( lsagas(j) .AND. lgas_out_ppm ) then
             spec_units(nspcout) = "ppmv"
          else
             spec_units(nspcout) = "micrograms m-3"
          endif
          spec_long_name(nspcout) = ptname(j)
          if( ptname(j)(7:9) .EQ. 'IC' ) then
              spec_desc(nspcout) = ptname(j)(:3)
     &                         //" air concentration - Initial Conditions"
          else if( ptname(j)(4:9) .EQ. '000BC' ) then
              spec_desc(nspcout) = ptname(j)(:3)
     &                         //" air concentration - All Boundary Conditions"
          else if( ptname(j)(4:9) .EQ. 'NTHBC' ) then
              spec_desc(nspcout) = ptname(j)(:3)
     &                         //" air concentration - Northern Boundary Conditions"
          else if( ptname(j)(4:9) .EQ. 'ESTBC' ) then
              spec_desc(nspcout) = ptname(j)(:3)
     &                         //" air concentration - Eastern Boundary Conditions"
          else if( ptname(j)(4:9) .EQ. 'STHBC' ) then
              spec_desc(nspcout) = ptname(j)(:3)
     &                         //" air concentration - Southern Boundary Conditions"
          else if( ptname(j)(4:9) .EQ. 'WSTBC' ) then
              spec_desc(nspcout) = ptname(j)(:3)
     &                         //" air concentration - Western Boundary Conditions"
          else if( ptname(j)(4:9) .EQ. 'TOPBC' ) then
              spec_desc(nspcout) = ptname(j)(:3)
     &                         //" air concentration - Top Boundary Conditions"
          else
              spec_desc(nspcout) = ptname(j)(:3)
     &                         //" air concentration - Group: "//ptname(j)(4:6)//
     &                                                    " Region: "//ptname(j)(7:9)
          endif
          spec_coords(nspcout) = "latitude longitude"
       endif
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
 
