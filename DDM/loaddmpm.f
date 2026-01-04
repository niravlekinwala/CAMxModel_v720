c*** LOADDMPM
c
      subroutine loaddmpm(filflg,numcol,numrow,numlay,nddm,grsens,
     &                    icl,jcl,kcl,nfam,nsen,sddm,convfac)
      use tracer
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c   Description:
c     This routine loads the sensitivies which are stored in a
c     4-D array from/into a 2-D array in conclusion/preperation 
c     of the DDM-PM chemistry routine.  The 2-D array contains the 
c     family of sensitivities for one cell and all modeled species.
c     The flag "filflg" determines wether the values in the 4-D
c     array are loaded into the 2-D array or vice-versa.  Sensitivities
c     for gases are converted from umol/m3 to ppm for chemistry and back.
c
c      Copyright 1996 - 2022
c     Ramboll
c
c   Argument descriptions:
c       filflg  L  flag for determining which direction to fill
c                  .TRUE.  = put 2-D values into 4-D gridded array
c                  .FALSE. = put 4-D values into 2-D array
c       numcol  I  number of cells in X direction
c       numrow  I  number of cells in Y direction
c       numlay  I  number of layers 
c       nddm    I  number of total DDM species
c       grsens  R  4-D array of DDM sensitivities
c       icl     I  the X grid location of current cell
c       jcl     I  the Y grid location of current cell
c       kcl     I  the vertical grid location of current layer
c       nfam    I  numbder of DDM families
c       nsen    I  number of modeled species
c       sddm    R  2-D array for this cell and species
c       convfac R  conversion from ppm to umol/m3
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     08/02/05   --bkoo--       Modified from loaddm for DDM-PM
c     06/25/08   --bkoo--       Added SDDM initialization to zero
c     09/01/12   --bkoo--       Updated for v5
c
c-----------------------------------------------------------------------
c   Include files:
c-----------------------------------------------------------------------
c
cbk      include 'camx.prm'
cbk      include 'tracer.com'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      logical filflg
      integer numcol
      integer numrow
      integer numlay
      integer nddm
      real    grsens(numcol,numrow,numlay,nddm)
      integer icl
      integer jcl
      integer kcl
      integer nfam
      integer nsen
      real    sddm(nfam,nsen)
      real    convfac
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      integer   ifam, ispc, iddm
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c  --- task 1: load 4-D array
c
      if ( filflg ) then
c
c  --- loop over the gas species ---
c
        do ispc = 1, NGDDMPM
          if ( jdpmap(ispc).eq.0 ) CYCLE
          do ifam = 1, nfam
            iddm = jdpmap(ispc) + ifam - 1
            grsens(icl,jcl,kcl,iddm) = sddm(ifam,ispc) * convfac
          enddo
        enddo
c
c  --- loop over the PM species ---
c
        do ispc = NGDDMPM+1, nsen
          if ( jdpmap(ispc).eq.0 ) CYCLE
          do ifam = 1, nfam
            iddm = jdpmap(ispc) + ifam - 1
            grsens(icl,jcl,kcl,iddm) = sddm(ifam,ispc)
          enddo
        enddo
c
c  --- task 2: load 2-D array
c
      else
c
c  --- SDDM initialization ---
c
        sddm = 0.0
c
c  --- loop over the gas species ---
c
        do ispc = 1, NGDDMPM
          if ( jdpmap(ispc).eq.0 ) CYCLE
          do ifam = 1, nfam
            iddm = jdpmap(ispc) + ifam - 1
            sddm(ifam,ispc) = grsens(icl,jcl,kcl,iddm) / convfac
          enddo
        enddo
c
c  --- loop over the PM species ---
c
        do ispc = NGDDMPM+1, nsen
          if ( jdpmap(ispc).eq.0 ) CYCLE
          do ifam = 1, nfam
            iddm = jdpmap(ispc) + ifam - 1
            sddm(ifam,ispc) = grsens(icl,jcl,kcl,iddm)
          enddo
        enddo
c
c  --- done ---
c
      endif
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
      return
      end
