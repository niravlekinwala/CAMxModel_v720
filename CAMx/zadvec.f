      subroutine zadvec(m1,m2,m3,i0,j0,ia,iz,ja,jz,ibcon,
     &                  losat,igrd,ncol,nrow,nlay,nrads,nspc,
     &                  nsen,delt,dx,
     &                  dy,idfin,depth,entrn,dilut,tempk,press,species,
     &                  conc,rhoscl,ctop,lvupsolv,ptrtop,fluxes,fluxtmp,
     &                  isaptr,ipa_xy,ipa_lay,iproc_id)
      use filunit
      use chmstry
      use bndary
      use procan
      use rtracchm
      use tracer
c
c----CAMx v7.20 220430
c
c     ZADVEC drives vertical transport of concentrations.  The operation
c     is performed on the regular model species (losat = False) or source
c     apportionment concentrations (losat = True).
c     This version also performs vertical transport of sensitivities 
c     if DDM is enabled.
c                          
c      Copyright 1996 - 2022
c     Ramboll
c          
c     Modifications:
c        12/3/99   Top boundary condition is converted from ppm to umol/m3
c                  using extrapolated density in an overlying layer
c                  that is assumed to be the same thickness as the top
c                  layer of the model
c       12/07/01   added instructions for OMP
c       01/30/02   Added code for RTRAC probing tool
c        3/06/06   Revised top BC approach; removed top dummy layer
c       05/01/07   Fixed bug that was causing PM PSAT species to
c                  be converted to micro-moles
c       07/16/07 -bkoo-     Added check for HDDM
c       07/16/08 -bkoo-     Added DDM turn-off flag
c       04/09/09   Improved definition of atmospheric density at top of model
c       11/04/09   Revised vertical advection solver technique, now
c                  employs zero-gradient top boundary condition
c       11/12/09   Added code to save fluxes and apply to the tracers
c        4/07/14   Added top con file
c       11/28/16 -bkoo-     Updated DDM for new top con
c        4/27/21   Added explicit PPM solver with time splitting 
c  
c     Input arguments:
c        igrd              grid index
c        ncol              number of columns
c        nrow              number of rows
c        nlay              number of layers
c        nrads             number of radicals
c        nspc              number of species
c        nsen              number of species X number of DDM parameters
c                          or number of tracer species
c                          or 1, whichever is larger
c        delt              time step (s)
c        dx                cell size in x-direction (m)
c        dy                cell size in y-direction (m)
c        idfin             map of nested grids in this grid
c        depth             layer depth (m)
c        entrn             entrainment rate (m/s)
c        dilut             dilution rate (m/s)
c        tempk             temperature (K)
c        press             pressure (mb)
c        species           species names
c        conc              species concentrations (umol/m3)
c        rhoscl            density scaling above top of model (unitless)
c        ctop              initial concentrations in top layer (umol/m3)
c        lvupsolv          3-D flag to use upstream donor solver
c        ptrtop            initial sensitivities in top layer
c        fluxtmp           temporary array for fluxes
c        ipa_xy            2-D gridded array to identify if cell is
c                          in a IPRM sub-domain
c        ipa_lay           3-D gridded array to identify IPRM sub-domain
c                          each layer is in
c        iproc_id          ID for this processor in MPI rank
c
c     Output arguments:
c        conc              species concentrations (umol/m3)
c        fluxes            fluxes across the boundaries (umol)
c
c     Routines Called:
c        VRTSLV
c
c     Called by:
c        EMISTRNS
c
      implicit none
      include "camx.prm"
      include "camx.inc"
      include "flags.inc"
c
      integer :: m1,m2,m3,i0,j0,ia,iz,ja,jz,ibcon
      integer :: m_zi1,m_zi2,m_zj1,m_zj2
c
      integer      ncol
      integer      nrow
      integer      nlay
      integer      nrads
      integer      nspc
      integer      igrd
      integer      nsen
      character*10 species(nspc)
      real         conc(m1,m2,m3,nspc)
      real         fluxtmp(m1,m2,m3,nspc)
      real         ctop(m1,m2,nspc)
      real         rhoscl(m1,m2)
      logical      lvupsolv(m1,m2,m3)
      real         ptrtop(m1,m2,nsen)
      integer*8    isaptr
      real, dimension(m1,m2,m3) :: entrn,dilut,tempk,press
      real, dimension(m1,m2,m3) :: depth
      real         dx(nrow)
      integer      idfin(m1,m2)
      logical      losat
      real*8       fluxes(nspc,11)
      real         delt
      real         dy
      real         entrnfac
      integer      iproc_id
      integer ispc,j,i,k,numspcs
      real convfac
      real c1d(MXLAYER+1)
      real d1d(MXLAYER)
      real ent1d(MXLAYER)
      real dil1d(MXLAYER)
      real fluxlay(MXLAYER)
      logical lvus(MXLAYER)
c
c===========================DDM Begin=================================
c
      integer isen,ioff
      integer*8 idx
      real sen1d((MXLAYER+1),MXFDDM)
c
c===========================DDM End===================================
c
c=================== Process Analysis Begin ==========================
c
      integer ipa_xy(ncol,nrow)
      integer ipa_lay(ncol,nrow,nlay)
      integer ipa_idx
      logical ldoipts
      real fc1(MXLAYER)
      real fc2(MXLAYER)
      real fc3(MXLAYER)
c
c==================== Process Analysis End ===========================
c
c
c-----Entry point
c
! Set grid point limits
      m_zi1= ia-1
      m_zi2= iz+1
      m_zj1= ja-1
      m_zj2= jz+1

      if(btest(ibcon,0)) m_zi1 = ia
      if(btest(ibcon,1)) m_zi2 = iz
      if(btest(ibcon,2)) m_zj1 = ja
      if(btest(ibcon,3)) m_zj2 = jz
c
c  ---- Initialize the temp flux array to zero ---
c
      if( .NOT. losat ) call zeros(fluxtmp,m1*m2*m3*nspc)
c
c-----Vertical advection, entrainment, and dilution
c
      entrnfac = 1.
      if (lzppm) entrnfac = -1.
      numspcs = nspc
      if( losat ) numspcs = nsaspc

      if( iproc_id .LE. 1 ) then
        if( .NOT. losat ) then
          write(*,'(a20,$)') 'z advection ......'
        else
          write(*,'(a20,$)') 'SA z advection ...'
        endif
      endif
      if( .NOT. losat ) then
         write(iout,'(a20,$)') 'z advection ......'
      else
         write(iout,'(a20,$)') 'SA z advection ...'
      endif
c
c$omp parallel default(shared)
c$omp&  private(ispc,i,j,k,ent1d,dil1d,d1d,c1d,lvus,fluxlay,
c$omp&          convfac,isen,ioff,idx,sen1d,
c$omp&          ldoipts,ipa_idx,fc1,fc2,fc3)
c
c$omp do schedule(guided)
c
      do 40 ispc = nrads+1,numspcs
           do 60 j = m_zj1, m_zj2   !2,nrow-1
             do 50 i = m_zi1, m_zi2  !2,ncol-1
c
c-----Skip cells occupied by child grids
c
            if (idfin(i,j).gt.igrd) goto 50
            do k = 1,nlay
              ent1d(k) = entrnfac*entrn(i,j,k)
              dil1d(k) = dilut(i,j,k)
              d1d(k)   = depth(i,j,k)
              c1d(k)   = conc(i,j,k,ispc)
              lvus(k)  = lvupsolv(i,j,k)
            enddo
            c1d(nlay+1) = rhoscl(i,j)*ctop(i,j,ispc)
            if (ltopcon) then
              convfac = 1.
              if (ispc .LE. ngas)
     &          convfac = densfac*(273./tempk(i,j,nlay))*
     &                            (press(i,j,nlay)/1013.)
              c1d(nlay+1) = convfac*c1d(nlay+1)
            endif
c
c================================DDM Begin============================
c
            if( (lddm .OR. lhddm) .AND. lddmcalc(igrd) ) then
              do isen = 1,nddmsp
                ioff = iptddm(ispc)+isen-1
                do k = 1,nlay
                  idx =  DBLE(isaptr-1) + DBLE(i) + 
     &               DBLE(m1)*DBLE(j-1) + DBLE(m1)*DBLE(m2)*DBLE(k-1) + 
     &                             DBLE(m1)*DBLE(m2)*DBLE(m3)*DBLE(ioff-1)
                  sen1d(k,isen) = ptconc(idx)
                enddo
                sen1d(k+1,isen) = rhoscl(i,j)*ptrtop(i,j,ioff)
                if (ltopcon) then
                  convfac = 1.
                  if (ispc .LE. ngas)
     &              convfac = densfac*(273./tempk(i,j,nlay))*
     &                                (press(i,j,nlay)/1013.)
                  sen1d(k+1,isen) = convfac*sen1d(k+1,isen)
                endif
              enddo
            endif
c
c================================DDM End==============================
c
c
c======================== Process Analysis Begin =====================
c
            ldoipts = .FALSE.
            if ( .NOT. ltrace .AND. lipr ) then
                if( i .GE. ia .AND. i .LE. iz .AND. j .GE. ja
     &                                         .AND. j .LE. jz ) then
                   if( ipa_xy(i+i0,j+j0) .GT. 0 ) ldoipts = .TRUE.
                endif
            endif
c
c
c========================= Process Analysis End ======================
c
c-----Solve vertical mass adjustments using Crank-Nicholson solver
c
            if (lzppm) then
              call vadvppm(nlay,igrd,delt,d1d,c1d,ent1d,dil1d,sen1d,fluxlay,
     &                     fc1,fc2,fc3,ldoipts)
            else
              call vrtslv(nlay,i,j,igrd,delt,ent1d,dil1d,d1d,c1d,
     &                    lvus,fluxlay(nlay),fluxlay,sen1d,
     &                    species(ispc),fc1,fc2,fc3,ldoipts)
            endif
c
c======================== Process Analysis Begin =====================
c
            if( ldoipts ) then
               do k=1,nlay
                  if( ipa_lay(i+i0,j+j0,k) .GT. 0 ) then
                    ipa_idx = ipa_lay(i+i0,j+j0,k)
c
c-----Concentration change from vertical advection
c
                    cipr(IPR_BADV, ipa_idx, ispc) =
     &                            cipr(IPR_BADV, ipa_idx, ispc) + fc1(k)
c
                    cipr(IPR_TADV, ipa_idx, ispc) =
     &                           cipr(IPR_TADV, ipa_idx, ispc) + fc2(k)
c
                    cipr(IPR_DADV, ipa_idx, ispc) =
     &                           cipr(IPR_DADV, ipa_idx, ispc) + fc3(k)
c
                  endif
               enddo
            endif
c
c========================= Process Analysis End ======================
c
            do k = 1,nlay
              conc(i,j,k,ispc) = c1d(k)
            enddo

            if( .NOT. losat ) then
               if(i .GE. ia .AND. i .LE. iz .AND. j .GE. ja 
     &                                     .AND. j .LE. jz ) then
                 if (fluxlay(nlay).lt.0) then
                   fluxtmp(i,j,1,ispc) = -fluxlay(nlay)*dx(j+j0)*dy*delt
                 else
                   fluxtmp(i,j,2,ispc) = -fluxlay(nlay)*dx(j+j0)*dy*delt
                 endif
               endif
            endif
c
c================================DDM Begin============================
c
            if( (lddm .OR. lhddm) .AND. lddmcalc(igrd) ) then
              do isen = 1, nddmsp
                ioff = iptddm(ispc)+isen-1
                do k = 1, nlay
                  idx =  DBLE(isaptr-1) + DBLE(i) + 
     &               DBLE(m1)*DBLE(j-1) + DBLE(m1)*DBLE(m2)*DBLE(k-1) + 
     &                           DBLE(m1)*DBLE(m2)*DBLE(m3)*DBLE(ioff-1)
                  ptconc(idx) = sen1d(k,isen)
                enddo
              enddo
            endif
c
c================================DDM End==============================
c

  50       continue
  60     continue
  40  continue
c
c$omp end parallel
c
c-----Put fluxes in global array
c
      if( .NOT. losat ) then
         do j = 2,m2-1
           do i = 2,m1-1
             do ispc = 1,nspc
               fluxes(ispc,9) = fluxes(ispc,9) + DBLE(fluxtmp(i,j,1,ispc))
               fluxes(ispc,10) = fluxes(ispc,10) + DBLE(fluxtmp(i,j,2,ispc))
             enddo
           enddo
         enddo
      endif
c
      if( iproc_id .LE. 1 ) then
        write(*,'(a)') '   Done'
        call flush(6)
      endif
      write(iout,'(a)') '   Done'
      call flush(iout)
c
      return
      end
