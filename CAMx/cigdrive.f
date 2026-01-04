      subroutine cigdrive(m1,m2,m3,i0,j0,ia,iz,ja,jz,igrd,ncol,nrow,
     &                    nlay,nspcs,nspdep,idfin,deltat,densfac,
     &                    dtout,deltax,deltay,mapscl,ldark,depth,tempk,
     &                    press,water,cigfrc,cigwtr,cigph,cigpcp,cwc,
     &                    cph,pwr,pws,pwg,conc,cigmas,depfld,fluxes,
     &                    xmschem,iptrsa,ipa_cel,wetfld,fsurf,rtconc,
     &                    solmass,vegmass,m1rt,m2rt,m3rt,nrtsp)
      use chmstry
      use filunit
      use tracer
      use rtracchm
c
c----CAMx v7.20 220430
c 
c     CIGDRIVE is the transport and chemistry driver for the cloud-in-grid
c     algorithm.  It splits vertical column concentration profiles into in-cloud
c     and ambient, applies convective mass transport matrices to both, calls
c     aqueous chemistry and wet deposition for both, then re-combines final
c     profiles. 
c
c     If CiG is not invoked, this routine just performs wet deposition for
c     the regular gridded concentrations.
c 
c      Copyright 1996 - 2022
c     Ramboll
c           
c     Modifications:
c        10/13/17 -bkoo-     Added a dummy argument to AEROCHEM_CF call (aerosol pH)
c        10/16/17 -bkoo-     Added a dummy argument to AEROCHEM_CF call (time-weighted JNO2)
c                            Changed DUMARR2 size
c
c     Input arguments:
c        m1                  number of slice columns
c        m2                  number of slice rows
c        m3                  number of slice layers
c        i0                  slice i-offset
c        j0                  slice j-offset
c        ia                  slice computational starting i-index
c        iz                  slice computational ending i-index
c        ja                  slice computational starting j-index
c        jz                  slice computational ending j-index
c        igrd                grid index 
c        ncol                number of grid columns 
c        nrow                number of grid rows 
c        nlay                number of grid layers 
c        nspcs               number of total species
c        nspdep              number of deposition species
c        idfin               map of nested grids in this grid
c        deltat              time step (s)
c        densfac             factor to convert umol/m3 to ppm
c        dtout               output frequency (minutes)
c        deltax              cell size in x-direction (m)
c        deltay              cell size in y-direction (m)
c        mapscl              map scale factor
c        ldark               darkness flag (T=dark)
c        depth               cell depth (m)
c        tempk               temperature field (K)
c        press               pressure field (mb)
c        water               humidity field (ppm)
c        cigfrc              CiG column fractional coverage
c        cigwtr              CiG cloud water (g/m3)
c        cigph               CiG cloud pH
c        cigpcp              CiG precipitation water (g/m3)
c        cph                 cloud water pH
c        cwc                 cloud water content (g/m3)
c        pwr                 rain water content (g/m3)
c        pws                 snow water content (g/m3)
c        pwg                 graupel water content (g/m3)
c        conc                concentration field (umol/m3, ug/m3)
c        cigmas              CiG tracer transport matrix
c             
c     Output arguments: 
c        conc                concentration field (umol/m3, ug/m3)
c        depfld              2-D array of wet deposited mass (mol/ha, g/ha)
c                            and surface liquid concentrations (mol/l, g/l)
c        fluxes              boundary mass fluxes (umol, ug)
c        xmschem             change in grid mass from chemistry (umol)
c             
c     Routines called: 
c        SCAVRAT
c        HENRYFNC
c             
c     Called by: 
c        EMISTRNS
c 
      implicit none
      include 'camx.prm'
      include 'flags.inc'
      include 'deposit.inc'
      include 'soap.inc'
      include 'chmdbg.inc'
c
      integer m1,m2,m3,i0,j0,ia,iz,ja,jz
      integer igrd,ncol,nrow,nlay,nspcs,nspdep
      integer idfin(m1,m2)
      real deltat,densfac,dtout
      real deltax(nrow),deltay,mapscl(m1,m2)
      real depth(m1,m2,m3)
      real tempk(m1,m2,m3)
      real press(m1,m2,m3)
      real water(m1,m2,m3)
      real cigfrc(m1,m2)
      real cigwtr(m1,m2,m3)
      real cigph(m1,m2,m3)
      real cigpcp(m1,m2,m3)
      real cwc(m1,m2,m3)
      real cph(m1,m2,m3)
      real pwr(m1,m2,m3)
      real pws(m1,m2,m3)
      real pwg(m1,m2,m3)
      real conc(m1,m2,m3,nspcs)
      real depfld(m1,m2,3*nspdep)
      real cigmas(2,m3,m3,m1,m2)
      real*8 fluxes(nspcs,11),xmschem(nspcs)
      logical ldark(m1,m2)

      integer i,j,k,l,kk,ll,ks,ke,n,n1,n2
      integer aero_flag
      real dtchem,area,acig,aamb,tcell,pcell,wcell,
     &     convfac,cldph,cliq,volume,vol
      real tmassi,tmassf,terr
      real ftmp(MXSPEC),depfld1(MXSPEC),depfld2(MXSPEC)
      real con(MXLAYER,MXSPEC,2)
      real cnc(MXLAYER,MXSPEC)
      real dfac(MXLAYER)
      real cc(MXSPEC)
      real dtmp(MXSPEC)
      real dz(MXLAYER),tt(MXLAYER),pp(MXLAYER),cwtr(MXLAYER),
     &     ph(MXLAYER),rwtr(MXLAYER),swtr(MXLAYER),gwtr(MXLAYER)
      real dumvar
      real dumarr1(NNV),dumarr2(7,MXSPEC),dumarr3(MXFDDM*NNV),dumarr4(1)
      real fluxtmp(MXCELLS,MXCELLS,MXSPEC)
      real dchem(MXCELLS,MXCELLS,MXSPEC)
      common /comcigdrive/ fluxtmp, dchem
c
c======================== Probing Tool Begin ===========================
c
      integer*8 iptrsa
      integer nrtsp,m1rt,m2rt,m3rt
      integer ipa_cel(ncol,nrow,nlay)
      real wetfld(m1,m2,notimespc)
      real rtconc(m1rt,m2rt,m3rt,nrtsp)
      real fsurf(m1,m2,nlu)
      real solmass(m1rt,m2rt,nrtsp),vegmass(m1rt,m2rt,nrtsp)
      real cncrt(MXLAYER,MXTRSP)
c
c========================= Probing Tool End ============================
c
c-----Entry point
c
      do l = 1,nspec
        do j = 1,m2
          do i = 1,m1
            fluxtmp(i,j,l) = 0.
            dchem(i,j,l) = 0.
          enddo
        enddo
      enddo
      dtchem = deltat/3600.
      igrdchm = igrd
      aero_flag = 0
      if (aeropt.EQ.'CF') aero_flag = 1
c
c-----Loop over rows and columns
c
      do 10 j = ja,jz   !2,nrow-1
c
c$omp parallel default(shared)
c$omp&  private(i,l,ll,n,n1,n2,ftmp,depfld1,depfld2,k,kk,ks,ke,con,area,
c$omp&  dz,tt,pp,acig,cnc,cwtr,ph,rwtr,swtr,gwtr,aamb,cncrt,tcell,pcell,
c$omp&  wcell,convfac,cc,cldph,cliq,tmassi,tmassf,terr,dfac,volume,vol,
c$omp&  dtmp,dumvar,dumarr1,dumarr2,dumarr3,dumarr4)
c$omp&  copyin(/ijkgrd/)
c
c$omp do schedule(dynamic)
c
        do 20 i = ia,iz    !2,ncol-1
          if (idfin(i,j).gt.igrd) goto 20

          do n = 1,2
            do l = 1,nspec
              do k = 1,nlay
                con(k,l,n) = 0.
              enddo
            enddo
          enddo
          do k = 1,nlay
            dfac(k) = press(i,j,k)/tempk(i,j,k)
          enddo
c
c-----Apply convective transfer matrices to split the column profile
c
          if (lcig(igrd) .AND. cigfrc(i,j).gt.0.) then
            
            do l = nrad+1,nspec
              tmassi = 0.
              tmassf = 0.
              do ks = 1,nlay    !Starting layer
                tmassi = tmassi + conc(i,j,ks,l)*depth(i,j,ks)
                do ke = 1,nlay  !Ending layer
                  con(ke,l,1) = con(ke,l,1) + 
     &              conc(i,j,ks,l)*cigmas(1,ke,ks,i,j)*dfac(ke)/dfac(ks)
                  con(ke,l,2) = con(ke,l,2) + 
     &              conc(i,j,ks,l)*cigmas(2,ke,ks,i,j)*dfac(ke)/dfac(ks)
                enddo
              enddo
c
c-----Check for conc < 0 and mass conservation
c
              do k = 1,nlay
                if (con(k,l,1).lt.0. .or. con(k,l,2).lt.0.) then
                  write(iout,'(//,a)') 'ERROR in CIGDRIVE:'
                  write(iout,'(2a)') 'Negative concentration after',
     &                               ' convective transport'
                  write(iout,*) 'igrd, i, j, k = ', igrd,i,j,k
                  write(iout,*) 'Spec  Name   In-cloud  Ambient'
                  do kk = 1,nlay
                    write(iout,'(i3,2x,a7,2e10.3)') kk,spname(l),
     &                                          con(kk,l,1),con(kk,l,2)
                  enddo
                  call camxerr()
                endif
                tmassf = tmassf + depth(i,j,k)*(cigfrc(i,j)*con(k,l,1)
     &                                 + (1. - cigfrc(i,j))*con(k,l,2))
              enddo
              terr = abs(tmassf - tmassi)/tmassi
              if (terr.gt.0.002) then
                write(iout,'(//,a)')'ERROR in CIGDRIVE:'
                write(iout,*)'CiG transport matrix not conserving mass'
                write(iout,*)'Grid,I,J            : ',igrd,i,j
                write(iout,*)'Species             : ',spname(l)
                write(iout,*)'Starting/ending mass: ',tmassi,tmassf
                call camxerr()
              endif
            enddo
          else
            do l = nrad+1,nspec
              do k = 1,nlay
                con(k,l,2) = conc(i,j,k,l)
              enddo
            enddo
          endif
c
c-----Wet deposition
c
          if (lwet) then
            area = deltax(j+j0)*deltay/mapscl(i,j)**2
            do k = 1,nlay
              dz(k) = depth(i,j,k)
              tt(k) = tempk(i,j,k)
              pp(k) = press(i,j,k)
            enddo
c
c-----CIG profile first
c
            if (lcig(igrd) .AND. cigfrc(i,j).gt.0.) then
              acig = cigfrc(i,j)*area
              do l = nrad+1,nspec
                ftmp(l) = 0.
                depfld1(l) = 0.
                depfld2(l) = 0.
                do k = 1,nlay
                  cnc(k,l) = con(k,l,1)
                enddo
              enddo
              do k = 1,nlay
                cwtr(k) = cigwtr(i,j,k)
                ph(k)   = cigph(i,j,k)
                rwtr(k) = 0.
                swtr(k) = 0.
                gwtr(k) = 0.
                if (tt(k).ge.273.) then
                  rwtr(k) = cigpcp(i,j,k)
                else
                  gwtr(k) = cigpcp(i,j,k)
                endif
              enddo
              call wetdep(m1,m2,m3,i0,j0,i,j,igrd,ncol,nrow,nlay,nspec,
     &                    deltat,acig,dz,tt,pp,cwtr,rwtr,swtr,gwtr,ph,
     &                    densfac,dtout,ipa_cel,cnc,ftmp,
     &                    depfld1,depfld2,wetfld,iptrsa)
c
c-----Accumulate core model deposition fluxes
c
              do l = nrad+1,nspec
                do k = 1,nlay
                  con(k,l,1) = cnc(k,l)
                enddo
                fluxtmp(i,j,l) = fluxtmp(i,j,l) + ftmp(l)
                do ll = 1,ndepspc
                  if (l .eq. ldepmap(ll)) then
                    depfld(i,j,ndepspc+ll) = depfld(i,j,ndepspc+ll) +
     &                                       depfld1(l) 
                    depfld(i,j,2*ndepspc+ll) = depfld(i,j,2*ndepspc+ll)
     &                                         + depfld2(l)
                    goto 100
                  endif
                enddo
 100            continue
              enddo
c
            endif
c
c-----Ambient profile next
c
            aamb = area
            if (lcig(igrd) .AND. cigfrc(i,j).gt.0.)
     &         aamb = (1. - cigfrc(i,j))*area
            do l = nrad+1,nspec
              ftmp(l) = 0.
              depfld1(l) = 0.
              depfld2(l) = 0.
              do k = 1,nlay
                cnc(k,l) = con(k,l,2)
              enddo
            enddo
            do k = 1,nlay
              cwtr(k) = cwc(i,j,k)
              ph(k)   = cph(i,j,k)
              rwtr(k) = pwr(i,j,k)
              swtr(k) = pws(i,j,k)
              gwtr(k) = pwg(i,j,k)
            enddo
            call wetdep(m1,m2,m3,i0,j0,i,j,igrd,ncol,nrow,nlay,nspec,
     &                  deltat,aamb,dz,tt,pp,cwtr,rwtr,swtr,gwtr,ph,
     &                  densfac,dtout,ipa_cel,cnc,ftmp,
     &                  depfld1,depfld2,wetfld,iptrsa)
c
c-----Accumulate core model deposition fluxes
c
            do l = nrad+1,nspec
              do k = 1,nlay
                con(k,l,2) = cnc(k,l)
              enddo
              fluxtmp(i,j,l) = fluxtmp(i,j,l) + ftmp(l)
              do ll = 1,ndepspc
                if (l .eq. ldepmap(ll)) then
                  depfld(i,j,ndepspc+ll) = depfld(i,j,ndepspc+ll) +
     &                                     depfld1(l)
                  depfld(i,j,2*ndepspc+ll) = depfld(i,j,2*ndepspc+ll) +
     &                                       depfld2(l)
                  goto 101
                endif
              enddo
 101          continue
            enddo
c
c======================== Probing Tool Begin ===========================
c
c-----Wetdep for RTRAC
c
            if( ltrace .AND.
     &          (tectyp.EQ.RTRAC .OR. tectyp.EQ.RTCMC)) then
              do l = 1,nrtrac
                do k = 1,nlay
                  cncrt(k,l) = rtconc(i,j,k,l)
                enddo
              enddo
c 
              call wetdeprt(m1,m2,m1rt,m2rt,i,j,nlay,nrtrac,deltat,aamb,
     &                      dz,tt,pp,cwtr,rwtr,swtr,gwtr,ph,densfac,
     &                      cncrt,solmass,vegmass,fsurf)
c 
              do l = 1,nrtrac
                do k = 1,nlay
                  rtconc(i,j,k,l) = cncrt(k,l)
                enddo
              enddo
            endif
          endif
c
c========================= Probing Tool End ============================
c
c-----Aqueous chemistry
c
          if (lcig(igrd) .and. lchem .and. aeropt.EQ.'CF') then
            n1 = 2
            n2 = 2
            if (cigfrc(i,j).gt.0.) n1 = 1
            do k = 1,nlay
              ichm = i
              jchm = j
              kchm = k
c
              tcell = tempk(i,j,k)
              pcell = press(i,j,k)
              wcell = water(i,j,k)
              tchm = tcell
              wchm = wcell
              convfac = densfac*(273./tcell)*(pcell/1013.)
              volume = deltax(j+j0)*deltay*depth(i,j,k)/(mapscl(i,j)**2)
c
              do n = n1,n2   !Loop over in-cloud and ambient profiles
                if (n.eq.1) then
                  cldph = cigph(i,j,k)
                  cliq  = cigwtr(i,j,k)
                  vol   = cigfrc(i,j)*volume
                else
                  cldph = cph(i,j,k)
                  cliq  = cwc(i,j,k)
                  vol   = volume
                  if (cigfrc(i,j).gt.0.) vol = (1. - cigfrc(i,j))*volume
                endif
                if (tcell.lt.273.)
     &            cliq = amax1(0.,cliq*(tcell - tamin)/(273. - tamin))
c
c-----Sum mass and convert conc units before chemsitry
c
                do l = 1,nspec
                  dtmp(l) = con(k,l,n)*vol
                enddo
                do l = nrad+1,ngas
                  cc(l) = con(k,l,n)/convfac 
                enddo
                do l = ngas+1,nspec
                  cc(l) = con(k,l,n)
                enddo
                do l = nrad+1,nspec
                  if (cc(l).lt.0.) then
                    write(iout,'(//,a)') 'ERROR in CIGDRIVE:'
                    write(iout,*) 'Negative concentration before chem'
                    write(iout,*) 'igrd, i, j, k = ', igrd,i,j,k
                    do ll = nrad+1,nspec
                      write(iout,'(i3,2x,a7,e10.3)')ll,spname(ll),cc(ll)
                    enddo
                    call camxerr()
                  endif
                  cc(l) = amax1(bdnl(l),cc(l))
                enddo
c
                call aerochem_cf(wcell,tcell,pcell,cliq,cldph,
     &                           dumvar,dumvar,
     &                           dumarr1,cc,convfac,dtchem,dtchem,
     &                           dumarr2,dumarr3,dumarr4,1,1,.FALSE.,1,
     &                           aero_flag,.FALSE.,1,ldark(i,j))
c
c-----Update cloud pH, convert conc units, and sum mass after chemsitry
c
                if (n.eq.1) then
                  cigph(i,j,k) = cldph
                else
                  cph(i,j,k) = cldph
                endif
                do l = nrad+1,ngas
                  con(k,l,n) = amax1(cc(l)*convfac,bdnl(l))
                enddo
                do l = ngas+1,nspec
                  con(k,l,n) = amax1(cc(l),bdnl(l))
                enddo
                do l = 1,nspec
                  dtmp(l) = con(k,l,n)*vol - dtmp(l)
                  dchem(i,j,l) = dchem(i,j,l) + dtmp(l)
                enddo
c
              enddo
            enddo
          endif
c
c-----Update main concentration array
c
          if (lcig(igrd) .AND. cigfrc(i,j).gt.0.) then
            do l = nrad+1,nspec
              do k = 1,nlay
                conc(i,j,k,l) = cigfrc(i,j)*con(k,l,1) + 
     &                          (1. - cigfrc(i,j))*con(k,l,2)
              enddo
            enddo
          else
            do l = nrad+1,nspec
              do k = 1,nlay
                conc(i,j,k,l) = con(k,l,2) 
              enddo
            enddo
          endif
c
 20     continue
c
c  --- end of parallelized loop ---
c
c$omp end parallel
c
 10   continue
c
c  --- load fluxes into global array ---
c
      do l = 1,nspec
        do j = ja,jz
          do i = ia,iz    !2,ncol-1
            fluxes(l,11) = fluxes(l,11) + fluxtmp(i,j,l)
            xmschem(l) = xmschem(l) + dchem(i,j,l)
           enddo
        enddo
      enddo
c
      call flush(6)
c
      return
      end
