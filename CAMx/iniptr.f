      subroutine iniptr(numcols,numrows)
      use filunit
      use grid
      use chmstry
      use tracer
      use procan
      implicit none
c
c----CAMx v7.20 220430
c
c     This routine calculates the pointers into the concentration
c     and meterology arrays for each grid.  The arrays are vectors
c     so this routine calculates the pointer that stores the 
c     first element of the 2-D or 3-D or 4-D field.  Some checks
c     are also made to ensure there is no array overflow of the
c     vectors.
c
c      Copyright 1996 - 2022
c     Ramboll
c          
c     Modifications:
c        10/8/99   Fixed a bug in assigning pointers for OSAT arrays for
c                  grids >= 3
c         9/3/02   Removed IPTRCL
c        1/13/03   Added IPTRDP for deposition output fields
c       11/10/03   Added IPSMP for RTRAC/PiG sampling grid output fields
c        7/29/05   Removed pointers for PiG sampling grid fields to READNML
c       04/30/13   Added surface model
c       09/02/14   Added subgrid convective model
c
c     Input arguments:
c        numcols   -- number of columns for each grid
c        numrows   -- number of columns for each grid
c
c     Output arguments:
c        nglay   -- maximum number of columns for each grid
c
c     Subroutine called:
c        none
c
c     Called by:
c        READNML
c
c-----Include files
c
      include 'camx.prm'
      include 'flags.inc'
      include 'deposit.inc'
c
c-----Argument declarations
c
      integer numcols(*)
      integer numrows(*)
c
c-----Local variables
c
      integer i
c
c-----Entry point
c
      iptr2d(1) = 1
      iptr2d_full(1) = 1
      iptr3d(1) = 1
      iptr3d_full(1) = 1
      iptr4d(1) = 1
      iptrav(1) = 1
      iptrem(1) = 1
      iptr1lay(1) = 1
      iptrlu(1) = 1
      iptrdp(1) = 1
      iptrsm(1) = 1
      ipsa2d(1) = 1
      ipsa2d_avrg(1) = 1
      ipsadep(1) = 1
      ipsa3d(1) = 1
      ipsa3d_ems(1) = 1
      iptrcig(1) = 1
      do i=2,ngrid
         iptr2d(i) = iptr2d(i-1) + numcols(i-1)*numrows(i-1)
         iptr2d_full(i) = iptr2d(i-1) + numcols(i-1)*numrows(i-1)
         iptr3d(i) = iptr3d(i-1) + numcols(i-1)*numrows(i-1)*nlay(i-1)
         iptr3d_full(i) = iptr3d(i)
         iptrcig(i) = iptrcig(i-1) + 
     &                2*nlay(i-1)*nlay(i-1)*numcols(i-1)*numrows(i-1)
         iptr4d(i) = iptr4d(i-1) + 
     &                       numcols(i-1)*numrows(i-1)*nlay(i-1)*nspec
         iptrem(i) = iptrem(i-1) + numcols(i-1)*numrows(i-1)*nlayers_ems*nemspc
         iptr1lay(i) = iptr1lay(i-1) + numcols(i-1)*numrows(i-1)*nspec
         if( nlu .GT. 0 ) then
            iptrlu(i) = iptrlu(i-1) + numcols(i-1)*numrows(i-1)*nlu
         else
            iptrlu(i) = 1
         endif
         iptrdp(i) = iptrdp(i-1) + numcols(i-1)*numrows(i-1)*
     &                             (ndepspc*3 + 2)
         iptrsm(i) = iptrsm(i-1) + numcols(i-1)*numrows(i-1)*nsmspc
         if( .NOT. l3davg(i-1) ) then 
             iptrav(i) = iptrav(i-1) + numcols(i-1)*numrows(i-1)*navspc
         else
             iptrav(i) = iptrav(i-1) + 
     &                       numcols(i-1)*numrows(i-1)*nlay(i-1)*navspc
         endif
         if( ltrace .OR. lddm .OR. lhddm .OR. lirr ) then
             if( lptdepout ) then
                 ipsadep(i) = ipsadep(i-1) + 
     &                              numcols(i-1)*numrows(i-1)*notimespc
             else
                 ipsadep(i) = 1
             endif
             ipsa2d(i) = ipsa2d(i-1) + numcols(i-1)*numrows(i-1)*ntotsp
             if( lsa_3davrg ) then
                ipsa2d_avrg(i) = ipsa2d_avrg(i-1) + 
     &                       numcols(i-1)*numrows(i-1)*nlay(i-1)*ntotsp
             else
                ipsa2d_avrg(i) = ipsa2d(i)
             endif
             ipsa3d(i) = ipsa3d(i-1) + DBLE(numcols(i-1))*
     &                      DBLE(numrows(i-1))*DBLE(nlay(i-1))*DBLE(ntotsp)
             ipsa3d_ems(i) = ipsa3d_ems(i-1) + DBLE(numcols(i-1))*
     &                    DBLE(numrows(i-1))*DBLE(nlayers_ems)*DBLE(ntotsp)
         else
             ipsadep(i) = 1
             ipsa2d(i) = 1
             ipsa2d_avrg(i) = 1
             ipsa3d(i) = 1
             ipsa3d_ems(i) = 1
         endif
      enddo
c
c----Return point
c
      return
      end
