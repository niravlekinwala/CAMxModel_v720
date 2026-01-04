      subroutine strato3(m1,m2,nlayer,ia,iz,ja,jz,ibcon,igrd,nspc,
     &                   imapnest,height,tempk,press,conc,topo,ctop)
      use chmstry
      use grid
      use tracer
c
c----CAMx v7.20 220430
c
c     STRATO3 replaces O3 concentrations in all layers above the locally-
c     diagnosed tropopause based on a linear O3 gradient between the top
c     BC (model top) and a fixed value of 100 ppb at the top of the layer
c     containing the tropopause.
c
c      Copyright 1996 - 2022
c     Ramboll
c
c     Modifications:
c        05/11/20   Modified the ozone interpolation over altitude
c                   from linear to an exponential power of 2.5,
c                   convert altitude from AGL to MSL
c        None
c
c     Input arguments:
c        m1                number of columns in slice/domain
c        m2                number of rows in slice/domain
c        nlayer            number of layers
c        ia,iz             row start/end index
c        ja,jz             column start/end index
c        ibcon             MPI slice boundary flag
c        igrd              grid index
c        nspc              number of species
c        imapnest          map of nested grids in this grid
c        height            layer interface heights (m)
c        tempk             temperature (K)
c        press             pressure (mb)
c        conc              species concentrations (umol/m3)
c        topo              terrain elevation (m)
c        ctop              top concentrations above model (ppm)
c
c     Output arguments:
c        conc              species concentrations (umol/m3)
c
c     Routines Called:
c        None   
c
c     Called by:
c        EMISTRNS
c
      implicit none
      include "camx.prm"
      include "camx.inc"
c
      integer      m1,m2,ia,iz,ja,jz,ibcon
      integer      nlayer
      integer      nspc
      integer      igrd
      real         conc(m1,m2,nlayer,nspc)
      real         ctop(m1,m2,nspc)
      real         topo(m1,m2)
      real         tempk(m1,m2,nlayer)
      real         press(m1,m2,nlayer)
      real         height(m1,m2,nlayer)
      integer      imapnest(m1,m2)

      integer j,i,k,ktrop
      integer m_zi1,m_zi2,m_zj1,m_zj2
      integer idxcel,idx,lsao3
      real convfac
      real tropo3, topo3
      real do3dz
      real ztop
      real c1d(MXLAYER)
      real z1d(MXLAYER)
      real t1d(MXLAYER)
      real p1d(MXLAYER)
c
c-----Entry point
c
c-----Set grid point limits
c
      m_zi1= ia-1
      m_zi2= iz+1
      m_zj1= ja-1
      m_zj2= jz+1

      if(btest(ibcon,0)) m_zi1 = ia
      if(btest(ibcon,1)) m_zi2 = iz
      if(btest(ibcon,2)) m_zj1 = ja
      if(btest(ibcon,3)) m_zj2 = jz
c
      do 60 j = m_zj1, m_zj2   !2,nrow-1
        do 50 i = m_zi1, m_zi2  !2,ncol-1
c
c-----Skip cells occupied by child grids
c
          if (imapnest(i,j).gt.igrd) goto 50

          do k = 1,nlayer
            z1d(k) = topo(i,j) + height(i,j,k)
            t1d(k) = tempk(i,j,k)
            p1d(k) = press(i,j,k)
            c1d(k) = conc(i,j,k,ko3)
          enddo

c-----Diagnose tropopause

          call findtrop(nlayer,t1d,z1d,ktrop)
          if (ktrop.eq.nlayer) cycle

c-----Interpolate TC ozone to the layers above the tropopause
          
          tropo3 = 0.100
          topo3  = ctop(i,j,ko3)
          ztop = z1d(nlayer)
          do3dz = (topo3 - tropo3)/(ztop - z1d(ktrop))**2.5
          do k = ktrop+1,nlayer
            convfac = densfac*(273./t1d(k))*(p1d(k)/1013.)
            conc(i,j,k,ko3) = tropo3 + 
     &                   do3dz*((z1d(k) + z1d(k-1))/2. - z1d(ktrop))**2.5
            conc(i,j,k,ko3) = conc(i,j,k,ko3)*convfac
c
c=================== Source Apportionment Begin  =====================
c
c  --- Re-allocate new stratospheric ozone to O3N+O3V top BC tracers,
c      set all other O3N+O3V tracers to zero (lower bound) ---
c
            if( ltrace .AND. lozone ) then

              idxcel =  ipsa3d(igrd)-1 + 
     &                  DBLE(i) + 
     &                  DBLE(m1)*DBLE(j-1) + 
     &                  DBLE(m1)*DBLE(m2)*DBLE(k-1)
              do lsao3 = iptcls(idxipt(ITRO3N)),nptcls(idxipt(ITRO3V))
                idx = idxcel + 
     &                DBLE(m1)*DBLE(m2)*DBLE(nlayer)*DBLE(lsao3-1)
                if( (lbndry .AND. ptname(lsao3)(4:6) .EQ. 'TOP') .OR.
     &           (.NOT.lbndry .AND. ptname(lsao3)(7:8) .EQ. 'BC') )then
                  ptconc(idx) = 0.5*conc(i,j,k,ko3)
                else
                  ptconc(idx) = BNDLPT
                endif
              enddo

            endif
c
c=================== Source Apportionment End  =======================
c
          enddo
c
  50    continue
  60  continue
c
      return
      end
