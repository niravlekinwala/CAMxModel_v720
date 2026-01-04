      subroutine iniptr_node()
c
c-----------------------------------------------------------------------
c    Modules used:
c-----------------------------------------------------------------------
c
      use node_mod
      use grid
      use chmstry
      use filunit
      use tracer
      use procan
c
      implicit none
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c    Argument descriptions:
c     Input:  
c     Output:  
c        ngcol   -- maximum number of columns for each grid
c        ngrow   -- maximum number of columns for each grid
c        nglay   -- maximum number of columns for each grid
c
c    Called by: 
c       NODES_ALLOC
c
c      Copyright 1996 - 2022
c     Ramboll
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'deposit.inc'
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      integer :: i
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c  ---- set the pointers ---
c
      iptr2d(1) = 1
      iptr3d(1) = 1
      iptr4d(1) = 1
      iptrav(1) = 1
      iptrem(1) = 1
      iptr1lay(1) = 1
      iptrlu(1) = 1
      iptrdp(1) = 1
      ipsadep(1) = 1
      ipsa2d(1) = 1
      ipsa2d_avrg(1) = 1
      ipsa3d(1) = 1
      ipsa3d_ems(1) = 1
      iptrcig(1) = 1
c
      do i=2,ngrid
         iptr2d(i) = iptr2d(i-1) + mmxp(i-1)*mmyp(i-1)
         iptr3d(i) = iptr3d(i-1) + mmxp(i-1)*mmyp(i-1)*mmzp(i-1)
         iptrcig(i) = iptrcig(i-1) + 
     &                2*mmzp(i-1)*mmzp(i-1)*mmxp(i-1)*mmyp(i-1)
         iptr4d(i) = iptr4d(i-1) + mmxp(i-1)*mmyp(i-1)*mmzp(i-1)*nspec
         iptrem(i) = iptrem(i-1) + mmxp(i-1)*mmyp(i-1)*nlayers_ems*nemspc
         iptr1lay(i) = iptr1lay(i-1) + mmxp(i-1)*mmyp(i-1)*nspec
         if( nlu .GT. 0 ) then
            iptrlu(i) = iptrlu(i-1) + mmxp(i-1)*mmyp(i-1)*nlu
         else
            iptrlu(i) = 1
         endif
         iptrdp(i) = iptrdp(i-1) + mmxp(i-1)*mmyp(i-1)*(ndepspc*3 + 2)
         if( .NOT. l3davg(i-1) ) then
            iptrav(i) = iptrav(i-1) + mmxp(i-1)*mmyp(i-1)*navspc
         else
            iptrav(i) = iptrav(i-1) + mmxp(i-1)*mmyp(i-1)*mmzp(i-1)*navspc
         endif
         if( ltrace .OR. lddm .OR. lhddm .OR. lirr ) then
            if( lptdepout ) then
               ipsadep(i) = ipsadep(i-1) + mmxp(i-1)*mmyp(i-1)*notimespc
            else
               ipsadep(i) = 1
            endif
            ipsa2d(i) = ipsa2d(i-1) + mmxp(i-1)*mmyp(i-1)*ntotsp
            if( lsa_3davrg ) then
              ipsa2d_avrg(i) = ipsa2d_avrg(i-1) + 
     &                      mmxp(i-1)*mmyp(i-1)*mmzp(i-1)*ntotsp
            else
              ipsa2d_avrg(i) = ipsa2d(i)
            endif
            ipsa3d(i) = ipsa3d(i-1) + DBLE(mmxp(i-1))*
     &                         DBLE(mmyp(i-1))*DBLE(mmzp(i-1))*DBLE(ntotsp)
            ipsa3d_ems(i) = ipsa3d_ems(i-1) + DBLE(mmxp(i-1))*
     &                       DBLE(mmyp(i-1))*DBLE(nlayers_ems)*DBLE(ntotsp)
         else
            ipsadep(i) = 1
            ipsa2d(i) = 1
            ipsa2d_avrg(i) = 1
            ipsa3d(i) = 1
            ipsa3d_ems(i) = 1
         endif
      enddo
c
      return
      end
