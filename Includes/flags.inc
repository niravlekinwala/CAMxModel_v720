c----CAMx v7.20 220430
c 
c     FLAGS.COM contains all model flags except those related to chemistry
c                           
c      Copyright 1996 - 2022
c     Ramboll
c           
c     Modifications: 
c        4/26/99   Removed LSMO, LBOT flags and replaced with IADVCT flag 
c        7/05/02   Removed lpig and replaced it with ipigflg
c        4/2/03    Removed option for UAM-V type cloud adjustment
c       10/09/08   Added ACM2 flag
c        8/09/09   Added dry deposition flag
c        7/14/10   Added in-line TUV flag
c       03/29/11   Support in-line TUV with aerosol optical depth
c       04/02/12   Removed RADM cloud adjustment option, cloud/aerosol
c                  adjustments now always done with in-line TUV
c       05/07/12   Added flexi-nesting flag
c       04/30/13   Added surface model
c        4/07/14   Added top con file
c       05/13/16   Added in-line Ix emissions flag
c       02/28/17   Added NetCDF output flag
c       07/23/18   Added Bi-Di NH3 drydep flag
c       07/10/19   Added stratospheric ozone treatment flag
c       07/23/19   Staggered wind flag is now grid-specific
c       04/15/21   Added flag for output gas species in PPM
c       01/05/22   Added flag for PPM vertical advection
c
c-----------------------------------------------------------------------
c     Integer flags:
c
c     ipigflg -- PiG flag: GREASD or IRON
c     iadvct  -- Advection solver flag
c                (2 = Bott, 3 = Piecewise Parabolic Method)
c     idrydep -- Dry deposition routine flag
c                (1 = Wesely 1989, 2 = Zhang 2003)
c
c     Logical flags:
c
c     lrstrt       -- simulation restart flag
c     lflexi       -- flexi-nesting flag
c     lchem        -- chemistry flag
c     ldry         -- dry deposition flag
c     lwet         -- wet deposition flag
c     lutm         -- UTM projection flag
c     llatlon      -- lat/lon coordinates flag
c     lrpolar      -- rotated (RAMS) polar stereographic projection flag
c     lambrt       -- Lambert conformal projection flag
c     lpolar       -- polar stereographic projection flag
c     lmerc        -- mercator projection flag
c     lstagw       -- staggered input wind field flag
c     larsrc       -- area source flag
c     lptsrc       -- point source flag
c     le1day       -- 1-day emission input flag
c     lairqul      -- output average file flag
c     ldiag        -- Diagnostic check flag
c     lcdfout      -- NetCDF output file flag
c     ncf_compress -- true if using NetCDF chunking
c     lmpi         -- flag to determine if the model is being run with MPI
c     lsuper       -- Super stepping advection flag (not allowed in MPI)
c     lacm2        -- ACM2 flag
c     lsrfmod      -- Surface model flag
c     ltopcon      -- Top concentration flag (T=input file)
c     lixemis      -- In-line Ix emissions flag
c     lixbypass    -- Bypass in-line Ix emissions flag
c     lbidinh3     -- Bi-Di NH3 drydep flag
c     lstrato3     -- Stratospheric ozone profile flag
c     lgas_out_ppm -- Gas species output in PPM
c     lzppm        -- Vertical PPM solver flag
c
c-----------------------------------------------------------------------
c
      integer ipigflg
      integer iadvct
      integer idrydep
      logical lrstrt
      logical lflexi
      logical lchem
      logical ldry
      logical lwet
      logical lutm
      logical llatlon
      logical lrpolar
      logical lambrt
      logical lpolar
      logical lmerc
      logical lstagw(MXGRID)
      logical larsrc
      logical lptsrc
      logical le1day
      logical lairqul
      logical ldiag
      logical lcdfout
      logical ncf_compress
      logical lmpi
      logical lsuper
      logical lacm2
      logical lsrfmod
      logical ltopcon
      logical lixemis
      logical lixbypass
      logical lbidinh3
      logical lstrato3
      logical lgas_out_ppm
      logical lzppm
c
      common /flags/ lrstrt, lchem, ldry, lwet, lutm, llatlon, lrpolar,
     &               lambrt, lpolar, lmerc, lstagw, larsrc, lptsrc, le1day, 
     &               lairqul, ldiag, lcdfout, ncf_compress, iadvct, ipigflg, lmpi, 
     &               lsuper, lacm2, idrydep, lflexi, lsrfmod, ltopcon, lixemis, 
     &               lixbypass,lbidinh3,lstrato3,lgas_out_ppm,lzppm
