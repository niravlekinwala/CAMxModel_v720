      subroutine pigdrive(m1,m2,m3,ioff,joff,
     &                    igrd,iptr2d,ncol,nrow,nlay,nspcs,nspdep,
     &                    itzon,dt,dx,dy,delx,dely,meshfac,mapscl,
     &                    height,rkv,tempk,tsurf,press,water,windu,
     &                    windv,cldtrns,cwc,pwr,pws,pwg,cph,cellat,
     &                    cellon,topo,sfcz0,albedo,vdep,conc,pigdump,
     &                    pgmserr,fluxes,depfld,ipsa2d,ipsa3d,ipsadep,
     &                    ipa_cel,iproc_id)
      use camxcom
      use filunit
      use chmstry
      use o3colmap
      use pigsty
      use tracer
      use procan
      use rtracchm
      use rtcmcchm
      use node_mod
      implicit none
c
c----CAMx v7.20 220430
c
c     PIGDRIVE is the driver program for the PiG submodel.
c     It performs chemistry, puff growth, and mass dumping.
c
c      Copyright 1996 - 2022
c     Ramboll
c
c     Modifications:
c        02/02/06        Renamed PIGDRIVE from IRONDRIV
c        8/31/06         Added map scale factor
c        01/08/07        Added Stage 3 criteria for Mech6 (CB05)
c                        Added ETOH,MTBE,MBUT rxns to Stage 3 criteria for Mech5
c        06/22/07        Removed termination check for puffs extending to
c                        surface layer
c        07/11/07        Added code for RTCMC probing tool
c        02/21/08        Added coupling time for puff aerosol chemistry
c        03/10/08        Added avgcnc and acgrad to calls to IRONCHEM
c        03/31/08        Added average puff age when killed
c        05/14/08        Improved handling of puff PM increments
c        08/04/08        Added IPR code for PiG changes (lost starting at v4.40)
c        11/02/08        Removed z0 calculation to SRFRUF routine
c        04/15/09        Modifications to support OSAT/PSAT deposition output
c        06/17/09        Restructured puff sigmas to improve impacts of
c                        vertical directional shear on puff growth
c        05/28/10        Added mechanism 7 (CB6)
c        07/14/10        Added code for in-line TUV cloud adjustment
c        04/02/12        Removed RADM cloud adjustment option, cloud/aerosol
c                        adjustments now always done with in-line TUV; AHO
c                        file is now just ozone column
c        04/04/12        Introduced 4-dim interpolation of J values over
c                        zenith, height, albedo, terrain
c        04/12/12        Added T and P adjustments to photolysis rates
c        07/20/12        Added Mechanism 1 (CB06 with Iodine chemistry)
c        10/08/12        Added Mechanism 2 (CB6r1)
c        07/31/13        Restructured to improve GREASD PiG chemistry,
c                        wet dep, and dumping
c        10/14/13        Added pressure to KHETERO call
c        06/27/14        Added Mechanism 3 (CB6r2h)
c        08/11/14        Albedo now includes effects of snow
c        11/11/14        Added Mechanism 4 (CB6r3)
c        12/07/14        Modified for VBS
c        02/10/16        Updated for OSAT3
c        02/21/20        Added Mechanism 1 (CB6r5) -rlb-
c
c     Input arguments:
c        m1                  number of columns in this slice
c        m2                  number of rows in this slice
c        m3                  number of layers in this slice
c        ioff                X offset of this slice in complete domain
c        ioff                Y offset of this slice in complete domain
c        igrd                grid index 
c        iptr2d              pointers into vectors for 2-D fields
c        ncol                number of columns
c        nrow                number of rows
c        nlay                number of layers
c        nspcs               number of species
c        nspdep              number of species in deposition array
c        itzon               time zone 
c        dt                  time step (s)
c        dx                  cell size in x-direction (m)
c        dy                  cell size in y-direction (m)
c        delx                cell size in x-direction (km or deg)
c        dely                cell size in y-direction (km or deg)
c        meshfac             grid meshing factor relative to master grid
c        mapscl              map scale factor
c        height              gridded layer height (m)
c        rkv                 gridded vertical diffusivity (m2/s)
c        tempk               gridded temperature (K)
c        tsurf               gridded surface temperature (K)
c        press               gridded pressure (mb)
c        water               gridded water vapor concentration (ppm)
c        windu               gridded x-component windspeed (m/s)
c        windv               gridded y-component windspeed (m/s)
c        cldtrns             gridded energy transmission coefficient (fraction)
c        cwc                 gridded cloud water (g/m3)
c        pwr                 gridded precip rain (g/m3)
c        pws                 gridded precip snow (g/m3)
c        pwg                 gridded precip graupel (g/m3)
c        cph                 gridded cloud water pH
c        cellat              gridded cell centroid latitude (deg)
c        cellon              gridded cell centroid longitude (deg)
c        topo                gridded cell terrain altitude (m MSL)
c        sfcz0               gridded cell surface roughness (m)
c        albedo              gridded cell surface albedo
c        vdep                gridded species deposition velocity (m/s)
c        conc                species concentration in the grid (umol/m3,ug/m3)
c        ipa_cel             gridded array to identify if cell is
c                            in a IPRM sub-domain
c        iproc_id            process ID for this slice
c
c     Output arguments:
c        conc                species concentration in the grid (umol/m3,ug/m3)
c        pigdump             PiG dumped mass (umol)
c        pgmserr             PiG dumping mass error (umol)
c        fluxes              boundary mass fluxes (umol, ug)
c        depfld              2-D array of wet deposited mass (mol/ha, g/ha)
c                            and surface liquid concentrations (mol/l, g/l)
c
c     Subroutines Called:
c        PIGCOORD
c        VIRTDUMP
c        KTHERM
c        GETZNTH
c        KPHOTO
c        KHETERO
c        NOXCHEM
c        IRONCHEM
c        CALDATE
c        PIGGROW 
c
c     Called by:
c        PIGEVOL
c
      include "camx.prm"
      include "flags.inc"
      include "deposit.inc"
      include 'chmdbg.inc'
      include 'camx_aero.inc'
c
c
      integer m1
      integer m2
      integer m3
      integer ioff
      integer joff
      integer igrd
      integer iptr2d
      integer ncol
      integer nrow
      integer nlay
      integer nspcs
      integer nspdep
      integer itzon
      real    dt
      real    dx(nrow)
      real    dy
      real    delx
      real    dely
      integer meshfac
      real    mapscl(m1,m2)
      real    height(m1,m2,m3)
      real    rkv(m1,m2,m3)
      real    tempk(m1,m2,m3)
      real    tsurf(m1,m2)
      real    press(m1,m2,m3)
      real    water(m1,m2,m3)
      real    windu(m1,m2,m3)
      real    windv(m1,m2,m3)
      real    cldtrns(m1,m2,m3)
      real    cwc(m1,m2,m3)
      real    pwr(m1,m2,m3)
      real    pws(m1,m2,m3)
      real    pwg(m1,m2,m3)
      real    cph(m1,m2,m3)
      real    cellat(m1,m2)
      real    cellon(m1,m2)
      real    topo(m1,m2)
      real    sfcz0(m1,m2)
      real    albedo(m1,m2)
      real    vdep(m1,m2,nspcs)
      real    conc(m1,m2,m3,nspcs)
      real*8  pigdump(*)
      real*8  pgmserr(*)
      real*8  fluxes(nspcs,11)
      real    depfld(m1,m2,nspdep)
      integer ipsa2d
      integer*8 ipsa3d
      integer ipsadep
      integer iproc_id
c
c======================== Process Analysis Begin ===============================
c
      integer ipa_cel(ncol,nrow,nlay)
c
c========================= Process Analysis End ================================
c
      logical ldarkp,ldump,lvdmp(MXPIG),lkill(MXPIG),lstage3
      logical lchemp
      real*8  flxdry(MXSPEC), flxwet(MXSPEC), puff_flux(MXPIG,MXSPEC)
      real    puff_dep(MXPIG,2*MXSPEC)
      real    conpig(MXSPEC+1,MXRECTR), conback(MXSPEC+1)
      real    connew(MXSPEC+1)
      real    dumpmass(MXSPEC), pconc(MXLAYER)
      real    volvec(MXLAYER), wtfac(MXLAYER), pcvfac(MXLAYER)
      real    pcorrec(MXLAYER,MXSPEC), othrpuf(MXLAYER,MXSPEC)
      real    pbdnl(MXLAYER)
      real    depth(MXLAYER), rkpig(MXLAYER), vd(MXSPEC)
      real    dpfl(2*MXSPEC)
      real    delconc(MXLAYER,MXSPEC+1),avgcnc(MXSPEC+1)
      real    sendum(1,1)
c
      integer n,nj,nz,i,j,k,kk,kkk,l,ll,is,ierr,
     &        kpb,kpt,ij,iozon,nr,idum,
     &        idxcel,noxcnt,
     &        hinox(MXRECTR),ierr1,ierr2,aero_flag
      integer ipa_idx
      integer kohso2,kohco,kohno2
      real areamx,pufscl,delt,dtpig,tpuff,wpuff,ppuff,sumwt,
     &     deplyr,alb,trn,
     &     hght,ctrns,zenith,convfac,volpuff,uu,vv,wind,press0,
     &     delms,rctdmp,volfac,tmass,tvol,resid,negconc,
     &     z0,bmass,hshear,vshear,dudx,dudy,dudz,dvdx,
     &     dvdy,dvdz,avgup,avgum,avgvp,avgvm,dz,zzp,zzm,cwpuff,
     &     phpuff,dtchem,totgas
      real xpig,ypig,xdist,ydist,xlen,axisx,pi,test,atm,o2,ch4,h2
      real pufang,shrang
      real pufno(MXRECTR),pufno2(MXRECTR)
c
      common /pigdrivdat/ puff_flux, puff_dep
c
c======================== Source Apportion Begin =======================
c
      integer*8 idx
      real    rtcon(MXTRSP),rtdump(MXTRSP)
      real    dpflrt(MXTRSP),vdrt(MXTRSP)
      real*8  flxdryrt(MXTRSP)
c
c========================= Source Apportion End ========================
c
      pi = 3.1415927
      atm = 1.e6
      O2  = 2.095e5
      CH4 = 1.75
      H2  = 0.50
c
c-----Entry point
c
      do i=1,MXSPEC+1
        avgcnc(i) = 0.
      enddo
      do i=1,MXPIG
        do j=1,MXSPEC
          puff_flux(i,j) = 0.
          puff_dep(i,j) = 0.
          puff_dep(i,MXSPEC+j) = 0.
        enddo
      enddo 
c
c-----Is PiG to do chemistry?
c
      lchemp = .false.
      if (lchem .and. idmech.ne.10) lchemp = .true.
c
c-----Check for negative incoming CONC values and initialize puff overlap
c     array VPCONC
c
      do is = 1,nspec
        do k = 1,nlay
          do j = 1,m2
            do i = 1,m1
              vpconc(i,j,k,is) = 0.
              if (conc(i,j,k,is) .lt. 0.) then
                write(iout,'(/,a)') 
     &             'Start of PIGDRIVE ...negative value in grid conc'
                write(iout,*) 'species, grid, i, j, k: ',is,igrd,i,j,k
                call flush(iout)
                call camxerr()
              endif
            enddo
          enddo
        enddo
      enddo
      do n = 1,npig
        lkill(n) = .false.
      enddo
c
c-----Perform dry deposition
c       
      if (ldry) then
        do 10 n = 1,npig
        if( lmpi ) Lslice(igrd,iproc_id,n) = 0
c
c-----Skip puffs that are not in the current grid
c
          if (ingrd(n).ne.igrd) goto 10
c
c-----Perform deposition once for new puffs over their age
c
          if (lnewg(n)) then
            delt = (agepigf(n) + agepigb(n))/2.
c
c-----Perform deposition for old puffs over current grid's timestep
c
          else
            delt = dt
          endif
c
          xpig = (xpigf(n) + xpigb(n))/2.
          ypig = (ypigf(n) + ypigb(n))/2.
          call pigcoord(xpig,ypig,i,j,idum)
c
c--- adjust cell index and skip if PiG not in this slice ---
c
          i = i-ioff
          j = j-joff
          if( i .LT. ia .OR. i .GT. iz ) goto 10
          if( j .LT. ja .OR. j .GT. jz ) goto 10
          if( lmpi ) Lslice(igrd,iproc_id,n) = 1
c
          if(puffbot(n) .GE.
     &            amax1(10.,height(i,j,1)/2.)) goto 10
          dz = height(i,j,1) - puffbot(n)
          do l = 1,nspec
            vd(l) = vdep(i,j,l)
            flxdry(l) = 0.
            dpfl(l) = 0.
          enddo
c
          call irondry(nrad,nspec,nreactr,dx(j+joff),dy,dz,mapscl(i,j),
     &                 delt,puffmass(1,1,n),axisz(n),vd,dpfl,flxdry)
c
          do 11 l = 1,nspec
            fluxes(l,11) = fluxes(l,11) + flxdry(l)
            do ll = 1,ndepspc
              if (l .eq. ldepmap(ll)) then
                depfld(i,j,ll) = depfld(i,j,ll) + dpfl(l)
                goto 11
              endif
            enddo
 11       continue
c
c======================== Source Apportion Begin =======================
c
          if (ltrace ) then
             if ( tectyp .EQ. RTRAC .OR. tectyp .EQ. RTCMC) then
                do l = 1,nrtrac
                  idx = DBLE(ipsa2d-1) + 
     &                DBLE(i) + DBLE(m1)*DBLE(j-1) + 
     &                            DBLE(m1)*DBLE(m2)*DBLE(l-1)
                  vdrt(l) = vdeprt(idx)
                enddo
                call irondry(0,nrtrac,nreactr,dx(j+joff),dy,dz,
     &                       mapscl(i,j),delt,puffrt(1,1,n),axisz(n),
     &                       vdrt,dpflrt,flxdryrt)
             else
                if( lptdepout ) then
                   call pigdrysa(m1,m2,m3,nrad,nspec,notimespc,i,j,n,
     &                           dpfl,ptconc(ipsa3d),ptdryfld(ipsadep))
                endif
             endif
          endif
c
c========================= Source Apportion End ========================
c
 10     continue
      endif
c
c-----Perform chemistry and wet deposition
c
      if (.not.lchemp .AND. .not.lwet) goto 101
c
c-----Aerosol timing parameters:
c     time_aero : time to call aerosol routine (HHMM)
c     date_aero : date to call aerosol routine (YYJJJ)
c
      aero_flag = 0
      if( lchemp .AND. naero.GT.0 .AND. aeropt.EQ.'CF' ) then
         aero_flag = 1
      endif
c
c-----Treat overlapping puffs via a virtual dump of "large" puffs
c 
      do 20 n = 1,npig
        lvdmp(n) = .false.
        if (ingrd(n).ne.igrd .or. .not.loverlap) goto 20
        xpig = (xpigf(n) + xpigb(n))/2.
        ypig = (ypigf(n) + ypigb(n))/2.
        call pigcoord(xpig,ypig,i,j,idum)
c
c--- adjust cell index and skip if PiG not in this slice ---
c
        i = i-ioff
        j = j-joff
        if( i .LT. ia .OR. i .GT. iz ) goto 20
        if( j .LT. ja .OR. j .GT. jz ) goto 20
        if( lmpi ) Lslice(igrd,iproc_id,n) = 1
c
        xdist = (xpigf(n) - xpigb(n))/delx*dx(j+joff)*meshfac
        ydist = (ypigf(n) - ypigb(n))/dely*dy*meshfac
        xlen  = sqrt(xdist**2 + ydist**2)
        axisx = xlen + 3.*sigx(n)
        pufscl = pi*(axisx/2.)*(1.5*sigy(n))
c
        if (DXYMAX .gt. 0.) then
          areamx = DXYMAX*DXYMAX
        elseif (DXYMAX .lt. 0.) then
          areamx = MIN(DXYMAX*DXYMAX,dx(j+joff)*dy)
        else
          areamx = dx(j+joff)*dy
        endif
        if (pufscl .ge. FLEAK*areamx) then
          lvdmp(n) = .true.
          call virtdump(m1,m2,m3,i0,j0,n,ncol,nrow,nlay,nspec,i,j,
     &                  dx(j+joff),dy,mapscl,height,tempk,press,pcorrec)
          do is = nrad+1,nspec
            do kk = 1,nlay
              vpconc(i,j,kk,is) = vpconc(i,j,kk,is) + pcorrec(kk,is)
            enddo
          enddo
        endif
 20   continue
c
c-----CHEMISTRY + WET DEP loop
c
c$omp parallel default(shared)
c$omp&  private(n,delt,dtpig,xpig,ypig,i,j,idum,k,xdist,ydist,xlen,
c$omp&          axisx,volpuff,kpb,kpt,kk,pcorrec,is,othrpuf,
c$omp&          wtfac,sumwt,deplyr,tpuff,ppuff,wpuff,cwpuff,
c$omp&          phpuff,conback,convfac,
c$omp&          connew,nr,conpig,l,dpfl,flxwet,
c$omp&          noxcnt,hinox,ij,iozon,trn,alb,hght,
c$omp&          ctrns,zenith,ldarkp,lstage3,test,
c$omp&          ierr1,ierr2,rtcon,avgcnc,nj,
c$omp&          dtchem,nz,ierr,totgas,pufno,pufno2,kohso2,kohco,kohno2)
c
c$omp do schedule(dynamic)
c
      do 100 n = 1,npig
c
c-----Skip puffs that are not in the current grid
c
        if (ingrd(n) .ne. igrd) goto 100
c
c-----Perform chemistry once for new puffs over their age
c
        if (lnewg(n)) then
          delt = (agepigf(n) + agepigb(n))/2.
c
c-----Perform chemistry for old puffs over current grid's timestep
c
        else
          delt = dt 
        endif 
        dtpig = min( (agepigf(n) + agepigb(n))/(2.*3600.),
     &                aero_dt(igrd)+(dt/3600.) )
c
c-----Locate the puff in the grid
c
        xpig = (xpigf(n) + xpigb(n))/2.
        ypig = (ypigf(n) + ypigb(n))/2.
        call pigcoord(xpig,ypig,i,j,idum)
c
c--- adjust cell index and skip if PiG not in this slice ---
c
        i = i-ioff
        j = j-joff
        if( i .LT. ia .OR. i .GT. iz ) goto 100
        if( j .LT. ja .OR. j .GT. jz ) goto 100
        if( lmpi ) Lslice(igrd,iproc_id,n) = 1
c
        do k = 1,nlay
          if( height(i,j,k ) .GT. zpig(n)) goto 25
        enddo
        k = nlay
c
  25    xdist = (xpigf(n) - xpigb(n))/delx*dx(j+joff)*meshfac
        ydist = (ypigf(n) - ypigb(n))/dely*dy*meshfac
        xlen  = sqrt(xdist**2 + ydist**2)
        axisx = xlen + 3.*sigx(n)
        volpuff = pi*axisx*axisy(n)*axisz(n)/6.
c
c-----Find layers containing top and bottom of puff
c
        kpb = 1
        kpt = nlay
        do kk = 1,nlay
          if( height(i,j,kk) .GE. pufftop(n) ) then
            kpt = kk
            goto 26
          endif
        enddo
  26    continue
        do kk = nlay,1,-1
          if( height(i,j,kk) .LE. puffbot(n) ) then
            kpb = kk + 1
            goto 27
          endif
        enddo
  27    continue
c
c-----Get overlap vector for this puff if it contributed to the 
c     virtual dump above (subtract out its contribution to VPCONC).
c     Otherwise, do not apply any puff overlap
c
        if (lvdmp(n)) then
          call virtdump(m1,m2,m3,i0,j0,n,ncol,nrow,nlay,nspec,i,j,
     &                  dx(j+joff),dy,mapscl,height,tempk,press,pcorrec)
          do is = nrad+1,nspec
            do kk = kpb,kpt
              othrpuf(kk,is) = vpconc(i,j,kk,is) - pcorrec(kk,is)
            enddo
          enddo
        else
          do is = nrad+1,nspec
            do kk = kpb,kpt
              othrpuf(kk,is) = 0.
            enddo
          enddo
        endif
c 
c-----Average temp, press, water, and background concs (=overlap+gridded)
c     over layers contained in this puff; define layer weighting as the 
c     depth of puff coverage weighted by density
c
        if (kpb .eq. kpt) then
          wtfac(kpb) = 1.
        else
          sumwt = 0.
          deplyr = pufftop(n) - height(i,j,kpt-1)
          wtfac(kpt) = deplyr*press(i,j,kpt)/tempk(i,j,kpt)
          sumwt = sumwt + wtfac(kpt)
c
          deplyr = height(i,j,kpb) - puffbot(n)
          wtfac(kpb) = deplyr*press(i,j,kpb)/tempk(i,j,kpb)
          sumwt = sumwt + wtfac(kpb)
c
          do kk = kpb+1,kpt-1
            deplyr = height(i,j,kk) - height(i,j,kk-1)
            wtfac(kk) = deplyr*press(i,j,kk)/tempk(i,j,kk)
            sumwt = sumwt + wtfac(kk)
          enddo
c
          do kk = kpb,kpt
            wtfac(kk) = wtfac(kk)/sumwt
          enddo
        endif
c
        tpuff = 0.
        ppuff = 0.
        wpuff = 0.
        cwpuff = 0.
        phpuff = 0.
        do is = nrad+1,nspec+1
          conback(is) = 0.
        enddo
        do kk = kpb,kpt
          tpuff = tpuff + wtfac(kk)*tempk(i,j,kk)
          ppuff = ppuff + wtfac(kk)*press(i,j,kk)
          wpuff = wpuff + wtfac(kk)*water(i,j,kk)
          cwpuff = cwpuff + wtfac(kk)*cwc(i,j,kk)
          phpuff = phpuff + wtfac(kk)*cph(i,j,kk)
          do is = nrad+1,nspec
            conback(is) = conback(is) + wtfac(kk)*
     &                    (othrpuf(kk,is) + conc(i,j,kk,is))
          enddo
        enddo
        convfac = densfac*(273./tpuff)*(ppuff/1013.)
        do is = nrad+1,nspec
          if (is.le.ngas) conback(is) = conback(is)/convfac
          conback(is) = amax1(bdnl(is),conback(is))
          connew(is) = conback(is)
        enddo
        connew(nspec+1) = 0.
c
c-----Perform wet deposition
c
        if (lwet .AND. (pwr(i,j,1).ge.cwmin .or.
     &                  pws(i,j,1).ge.cwmin .or.
     &                  pwg(i,j,1).ge.cwmin)) then
c
          do nr = 1,nreactr
            do is = nrad+1,nspec
              conpig(is,nr) = puffmass(is,nr,n)/(volpuff/float(nreactr))
              if (is.le.ngas) conpig(is,nr) = conpig(is,nr)/convfac
            enddo
          enddo
c
          call ironwet(nspec,n,dx(j+joff),dy,mapscl(i,j),cwpuff,
     &                 pwr(i,j,1),pws(i,j,1),pwg(i,j,1),phpuff,tpuff,
     &                 ppuff,volpuff,conpig,conback,axisz(n),delt,dpfl,
     &                 flxwet)
c
          do l = nrad+1,nspec
            puff_flux(n,l) = puff_flux(n,l) + flxwet(l)
            puff_dep(n,l) = puff_dep(n,l) + dpfl(l)
            puff_dep(n,MXSPEC+l) = puff_dep(n,MXSPEC+l)
     &                                             + dpfl(nspec+l)
          enddo
c
c======================== Source Apportion Begin =======================
c
c  --- if doing tracer depostion, call the routine
c      to aadd to tracer arrays ---
c
          if( lptdepout ) then
             call pigwetsa(m1,m2,m3,nrad,nspec,notimespc,i,j,kpb,kpt,n,
     &                           dpfl,ptconc(ipsa3d),ptwetfld(ipsadep))
          endif
c
c======================== Source Apportion Begin =======================
c
        endif
c
c-----Perform chemistry
c
        if (.not.lchemp) goto 100
c
c-----Add background concentrations (CONBACK) to puff increments in
c     each reactor to get total concentration (CONPIG).
c     Background radicals taken from nominal pig center (no averaging)
c     Puff radicals taken from last time step to use as initial guess
c
        do nr = 1,nreactr
          totgas = 0.
          do is = nrad+1,nspec
            conpig(is,nr) = puffmass(is,nr,n)/(volpuff/float(nreactr))
            if (is.le.ngas) then
              conpig(is,nr) = conpig(is,nr)/convfac
              totgas = totgas + abs(conpig(is,nr))
            endif
            conpig(is,nr) = max(bdnl(is),conback(is) + conpig(is,nr))
          enddo
          if (totgas.lt.0.001) goto 100
          pufno(nr)  = conpig(kno,nr)  - conback(kno)
          pufno2(nr) = conpig(kno2,nr) - conback(kno2)
        enddo
        do nr = 1,nreactr
          conpig(nspec+1,nr) = 0.
        enddo
        if (lchemp) then
          do nr = 1,nreactr
            do l = 1,nrad
              conback(l) = conc(i,j,k,l)
              connew(l) = conback(l)
              conpig(l,nr) = puffmass(l,nr,n) ! puffmass(radical) in ppm
            enddo
          enddo
        endif
c
c-----Determine if any puff reactors have high NOx values
c
        noxcnt = 0
        do nr = 1,nreactr
          hinox(nr) = 0
          if (conpig(kno,nr) + conpig(kno2,nr) .gt. 0.2) then
            noxcnt = noxcnt + 1
            hinox(nr) = 1
          endif
        enddo
c
c-----Determine thermal rate constants
c
        call ktherm(tpuff,ppuff)
c
c-----Load local values of ozone, albedo, terrain ht, zenith angle, altitude
c
        call getznth(cellat(i,j),cellon(i,j),timec(igrd),datec(igrd),
     &               itzon,zenith,ldarkp)
c
        alb   = albedo(i,j)
        ij = i + (j - 1)*m1
        iozon = icdozn(iptr2d - 1 + ij)
c
        trn = topo(i,j)/1000.
        hght = height(i,j,k)/2000.
        if (k.gt.1) hght = (height(i,j,k) + height(i,j,k-1))/2000.
        ctrns = 0.
        if (.not.ldarkp) ctrns = cldtrns(i,j,k)
c
c-----Determine photolysis rates through interpolation of look-up table
c
        call kphoto(iozon,alb,trn,hght,zenith,ctrns,ldarkp,tpuff,ppuff)
c
c-----Determine if any GREASD puff reactors are in Stage 3 chemistry
c
        if (ipigflg.EQ.GRESPIG) then
          lstage3 = .false.
          do nr = 1,nreactr
            if (idmech.eq.5) then
              kohso2 = 44
              kohco  = 29
              kohno2 = 25
            elseif (idmech.eq.6) then
              kohso2 = 63
              kohco  = 66
              kohno2 = 28
            elseif (idmech.eq.3) then
              kohso2 = 52
              kohco  = 123
              kohno2 = 45
            elseif (idmech.eq.7) then
              kohso2 = 51
              kohco  = 50
              kohno2 = 41
            elseif (idmech.eq.1 .or. idmech.eq.4) then
              kohso2 = 52
              kohco  = 120
              kohno2 = 45
            else
              write(iout,'(/,a)')
     &         'In PIGDRIVE, selected mechanism not implemented for PiG'
              write(iout,'(a,i3)') 'Mechanism: ',idmech
              call flush(iout)
              call camxerr()
            endif
            test = 0.
            if( pufno2(nr) .NE. 0. ) 
     &            test = (rk(kohso2)*conpig(kSO2,nr) +
     &              rk(kohco)*conpig(kCO,nr))  /
     &             (rk(kohno2)*pufno2(nr))
            if (.not.ldarkp .AND. test.gt.1. .AND.
     &          pufno2(nr).gt.pufno(nr)) lstage3 = .true.
            if (pufno2(nr).lt.conback(kNO2)) lstage3 = .true.
c
            if (conpig(kNO2,nr).le.0. .OR. conpig(kNO,nr).le.0. .OR.
     &                                              lstage3) then
              lkill(n) = .true.
              nkill(5) = nkill(5) + 1
              goto 100
            endif
c
          enddo
        endif
c
c-----For IRONCHEM, perform chemistry for CONBACK if total NOx
c     (puff+background) is less than 200 ppb.
c     All reactors use same background, save results in CONNEW
c
        if (ipigflg .EQ. IRONPIG) then
          igrdchm = -n
          ichm    = i + ioff
          jchm    = j + joff
          kchm    = 0
          ierr1   = 0
          ierr2   = 0
          tchm = tpuff
          wchm = wpuff
          ldchm = ldarkp
c
c-----Calculate the rate constant for heterogeneous hydrolysis of N2O5
c
          call khetero(tpuff,ppuff,wpuff,conback,.FALSE.,1,sendum)
c
          if (noxcnt .ne. nreactr) then
            call ironchem(delt,dtpig,aero_flag,wpuff,tpuff,ppuff,cwpuff,
     &                    phpuff,ldarkp,connew,conback,avgcnc,convfac,
     &                    ierr1,ierr2)
            if (ierr1 .ne. 0) then
              lkill(n) = .true.
              nkill(3) = nkill(3) + 1
c             write(idiag,*)'Dumping puff ',n,': bad background LSODE'
              goto 100
            endif
          endif
        endif
c
c-----Loop over reactors, perform chemistry for CONPIG.
c     Check for NOx concentrations > 200 ppb; if found, use the 
c     reduced chemistry module
c
        do nr = 1,nreactr
          igrdchm = -n
          ichm    = i + ioff
          jchm    = j + joff
          kchm    = -nr
          ierr1   = 0
          ierr2   = 0
          tchm = tpuff
          wchm = wpuff
          ldchm = ldarkp
c
c-----Calculate the rate constant for heterogeneous hydrolysis of N2O5
c
          call khetero(tpuff,ppuff,wpuff,conpig(1,nr),.FALSE.,1,sendum)
c
c-----For GREASD PiG, remove background except for oxidants (O3,H2O2,FORM,CO)
c
          if (ipigflg.EQ.GRESPIG) then
            do is = nrad+1,nspec
              conpig(is,nr) = conpig(is,nr) - conback(is)
            enddo
            conpig(ko3,nr)   = MAX(bdnl(ko3),
     &                             conpig(ko3,nr) + conback(ko3))
            conpig(kh2o2,nr) = MAX(bdnl(kh2o2),
     &                             conpig(kh2o2,nr) + conback(kh2o2))
            conpig(kco,nr)   = MAX(bdnl(kco),
     &                             conpig(kco,nr) + conback(kco))
            if (idmech.eq.5) then
              conpig(khcho,nr)   = MAX(bdnl(khcho),
     &                             conpig(khcho,nr) + conback(khcho))
            else
              conpig(kform,nr)   = MAX(bdnl(kform),
     &                             conpig(kform,nr) + conback(kform))
            endif
          endif
c
          if (hinox(nr) .eq. 1) then
            call noxchem(delt,bdnl(ko3),bdnl(kno),conpig(kno,nr),
     &                   conpig(kno2,nr),conpig(ko3,nr),rk(ipigrxn))
          else
            call ironchem(delt,dtpig,aero_flag,wpuff,tpuff,ppuff,cwpuff,
     &                    phpuff,ldarkp,conpig(1,nr),conback,avgcnc,
     &                    convfac,ierr1,ierr2)
            if (ierr1 .ne. 0) then
              lkill(n) = .true.
              nkill(4) = nkill(4) + 1
c             write(idiag,*)'Dumping puff ',n,': bad puff LSODE'
              goto 100
            endif
            if (ierr2 .ne. 0) then
              lkill(n) = .true.
              nkill(9) = nkill(9) + 1
c             write(idiag,*)'Dumping puff ',n,': bad puff PM'
              goto 100
            endif
          endif
c
c======================== Source Apportion Begin =======================
c
c-----Original RTRAC
c
          if( ltrace .AND. tectyp .EQ. RTRAC .AND.
     &        (nrtherm .GT. 0 .OR. nrtphot .GT. 0) ) then
             do is = 1,nrtgas
               rtcon(is) = puffrt(is,nr,n)/(volpuff/float(nreactr))
               rtcon(is) = rtcon(is)/convfac
             enddo
             if (hinox(nr) .eq. 1) then
               avgcnc(koh)  = conpig(koh,nr)
               avgcnc(kno3) = conpig(kno3,nr)
               avgcnc(ko3)  = conpig(ko3,nr)
             endif
             dtchem = delt/3600.
             call pigchmrt(ppuff,tpuff,avgcnc(koh),avgcnc(ko3),
     &                     avgcnc(kno3),dtchem,rtcon)
             do is = 1,nrtgas
               rtcon(is) = rtcon(is)*convfac
               puffrt(is,nr,n) = rtcon(is)*volpuff/float(nreactr)
             enddo
          endif
c
c-----RTCMC
c
          if( ltrace .AND. tectyp .EQ. RTCMC ) then
             if( ktype .EQ. 1 ) then
                call krtc(tpuff,ppuff)
             else
                call ksci(tpuff,ppuff)
             endif
             if( njschm .GT. 0 ) then
                 do nj = 1,njschm
                    srkrtc(ijschm(nj)) = 0.
                    do nz = nzschm,1,-1
                       if( zenith .LE. zenschm(nz) )
     &                    srkrtc(ijschm(nj)) = rjschm(nz,nj)
                    enddo
                    rkrtc(ijschm(nj))  = DBLE( srkrtc(ijschm(nj)) )
                 enddo
             endif
             do is = 1,nrtgas
               rtcon(is) = puffrt(is,nr,n)/(volpuff/float(nreactr))
               rtcon(is) = rtcon(is)/convfac
             enddo
             if (hinox(nr) .eq. 1) then
               do is = 1,nspec
                 avgcnc(is)  = conpig(is,nr)
               enddo
             endif
             dtchem = delt/3600.
             if( isolv .EQ. 1) then
                call dlsdrv(dtchem,wpuff,atm,O2,CH4,H2,avgcnc,
     &                      rtcon,ierr,.true.)
             elseif( isolv .EQ. 2) then
                call slsdrv(dtchem,wpuff,atm,O2,CH4,H2,avgcnc,
     &                      rtcon,ierr,.true.)
             elseif( isolv .EQ. 3) then
                call rbkdrv(dtchem,wpuff,atm,O2,CH4,H2,avgcnc,
     &                      rtcon,ierr,.true.)
             endif
             do is = 1,nrtgas
               rtcon(is) = rtcon(is)*convfac
               puffrt(is,nr,n) = rtcon(is)*volpuff/float(nreactr)
             enddo
          endif
c
c========================= Source Apportion End ========================
c
c
c-----Subtract off the new background to get new increment.
c     CONPIG and PUFFMASS can be negative again
c
          if (ipigflg.EQ.IRONPIG) then
            do is = nrad+1,nspec
              if (hinox(nr).eq.1) then
                conpig(is,nr) = conpig(is,nr) - conback(is)
              else
                conpig(is,nr) = conpig(is,nr) - connew(is)
              endif
            enddo
          else
            conpig(kO3,nr)   = conpig(kO3,nr)   - connew(kO3)
            conpig(kh2o2,nr) = conpig(kh2o2,nr) - connew(kh2o2)
            conpig(kco,nr)   = conpig(kco,nr)   - connew(kco)
            if (idmech.eq.5) then
              conpig(khcho,nr) = conpig(khcho,nr) - connew(khcho)
            else
              conpig(kform,nr) = conpig(kform,nr) - connew(kform)
            endif
          endif
          do is = nrad+1,nspec
            if (is.le.ngas) conpig(is,nr) = conpig(is,nr)*convfac
            puffmass(is,nr,n) = conpig(is,nr)*volpuff/float(nreactr)
          enddo
          do l = 1,nrad
            puffmass(l,nr,n) = conpig(l,nr) ! puffmass(radical) in ppm
          enddo
        enddo
 100  continue
c
c  --- end of parallelized loop ---
c
c$omp end parallel
c
c  --- put the fluxes and wet deposition into the totals array ---
c
      do 110 n = 1,npig
        if( ingrd(n) .EQ. igrd ) goto 110
c
c-----Locate the puff in the grid
c
        xpig = (xpigf(n) + xpigb(n))/2.
        ypig = (ypigf(n) + ypigb(n))/2.
        call pigcoord(xpig,ypig,i,j,idum)
c
c--- adjust cell index and skip if PiG not in this slice ---
c
        i = i-ioff
        j = j-joff
        if( i .LT. ia .OR. i .GT. iz ) goto 110
        if( j .LT. ja .OR. j .GT. jz ) goto 110
c
c---- accumulate in totals array ---
c
        do 111 l=1,nspec
          fluxes(l,11) = fluxes(l,11) + puff_flux(n,l)
          do ll=1,ndepspc
             if( l .EQ. ldepmap(ll) ) then
                depfld(i,j,ndepspc+ll) = depfld(i,j,ndepspc+ll) + 
     &                                                    puff_dep(n,l)
                depfld(i,j,2*ndepspc+ll) = depfld(i,j,2*ndepspc+ll) +
     &                                              puff_dep(n,MXSPEC+l)
     &                                    
                goto 111
             endif
          enddo
  111   continue
  110 continue
c
c-----End CHEMISTRY loop
c
c-----Start GROWTH and DUMPING loop
c
 101  continue
      do 200 n = 1,npig
c
c-----Skip puffs that are not in the current grid
c
        if (ingrd(n) .ne. igrd) goto 200
c
c-----Perform growth once for new puffs over their age
c
        if (lnewg(n)) then
          delt = (agepigf(n) + agepigb(n))/2.
c
c-----Perform growth for old puffs over current grid's timestep
c
        else
          delt = dt 
        endif 
c
c-----Locate the pig in the grid
c
        xpig = (xpigf(n) + xpigb(n))/2.
        ypig = (ypigf(n) + ypigb(n))/2.
        call pigcoord(xpig,ypig,i,j,idum)
c
c--- adjust cell index and skip if PiG not in this slice ---
c
        i = i-ioff
        j = j-joff
        if( i .LT. ia .OR. i .GT. iz ) goto 200
        if( j .LT. ja .OR. j .GT. jz ) goto 200
        if( lmpi ) Lslice(igrd,iproc_id,n) = 1

        idxcel =  i + m1*(j-1)
        do k = 1,nlay
          if( height(i,j,k) .gt. zpig(n) ) goto 15
        enddo
        k = nlay
  15    xdist = (xpigf(n) - xpigb(n))/delx*dx(j+joff)*meshfac
        ydist = (ypigf(n) - ypigb(n))/dely*dy*meshfac
        xlen  = sqrt(xdist**2 + ydist**2)
        pufang = atan(ydist/(xdist+1.e-10))
c
c-----Find layers containing top and bottom of puff
c
        kpb = 1
        kpt = nlay
        do kk = 1,nlay
          if( height(i,j,kk) .GE. pufftop(n) ) then
            kpt = kk
            goto 151
          endif
        enddo
  151   continue
        do kk = nlay,1,-1
          if( height(i,j,kk) .LE. puffbot(n) ) then
            kpb = kk + 1
            goto 152
          endif
        enddo
  152   continue
c
c-----Get the drivers on puff growth
c
c-----Resolved shear metrics: use deformation for net horizontal shear
c
        dudx = (windu(i,j,k) - windu(i-1,j,k))/dx(j+joff)
        dvdy = (windv(i,j,k) - windv(i,j-1,k))/dy
        avgup = (windu(i,j,k) + windu(i,j+1,k) + windu(i-1,j,k) + 
     &           windu(i-1,j+1,k))/4.
        avgum = (windu(i,j,k) + windu(i,j-1,k) + windu(i-1,j,k) + 
     &           windu(i-1,j-1,k))/4.
        dudy = (avgup - avgum)/dy
        avgvp = (windv(i,j,k) + windv(i+1,j,k) + windv(i,j-1,k) + 
     &           windv(i+1,j-1,k))/4.
        avgvm = (windv(i,j,k) + windv(i-1,j,k) + windv(i,j-1,k) + 
     &           windv(i-1,j-1,k))/4.
        dvdx = (avgvp - avgvm)/dx(j+joff)
        hshear = sqrt((dudy + dvdx)**2 + (dudx - dvdy)**2)
c
        dudz = 0.
        dvdz = 0.
        dz = height(i,j,k)
        if (k.gt.1) dz = height(i,j,k) - height(i,j,k-1)
        if( axisz(n) .gt. dz .and. kpt-kpb .ge. 1 ) then
          zzp = (height(i,j,kpt) + height(i,j,kpt-1))/2.
          zzm = height(i,j,kpb)/2.
          if (kpb.gt.1) zzm = (height(i,j,kpb) + height(i,j,kpb-1))/2.
          avgup = (windu(i,j,kpt) + windu(i-1,j,kpt))/2.
          avgum = (windu(i,j,kpb) + windu(i-1,j,kpb))/2.
          dudz = (avgup - avgum)/(zzp - zzm)
          avgvp = (windv(i,j,kpt) + windv(i,j-1,kpt))/2.
          avgvm = (windv(i,j,kpb) + windv(i,j-1,kpb))/2.
          dvdz = (avgvp - avgvm)/(zzp - zzm)
        endif
        vshear = sqrt(dudz**2 + dvdz**2)
        shrang = atan(dvdz/(dudz+1.e-10))
c
c-----Surface layer met
c
        uu = 0.5*(windu(i,j,1) + windu(i-1,j,1))
        vv = 0.5*(windv(i,j,1) + windv(i,j-1,1))
        wind = sqrt(uu*uu + vv*vv)
        do kk = 1,nlay
          depth(kk) = height(i,j,kk)
          if (kk.gt.1) depth(kk) = height(i,j,kk) - height(i,j,kk-1)
          rkpig(kk) = rkv(i,j,kk)
        enddo
        press0 = press(i,j,1) -  depth(1)*(press(i,j,2) - 
     &                                   press(i,j,1))/height(i,j,2)
        z0 = sfcz0(i,j)
c
c-----Calculate new puff size and fraction of mass that has leaked (DELMS)
c
        if (DXYMAX .gt. 0.) then
          areamx = DXYMAX*DXYMAX
        elseif (DXYMAX .lt. 0.) then
          areamx = MIN(DXYMAX*DXYMAX,dx(j+joff)*dy)
        else
          areamx = dx(j+joff)*dy
        endif
        ldump = .false.
        call piggrow(n,nlay,k,delt,xlen,areamx,
     &               height(i,j,nlay),wind,
     &               tempk(i,j,1),tsurf(i,j),press(i,j,1),press0,rkpig,
     &               depth,z0,hshear,vshear,pufang,shrang,delms,ldump)
c
c-----Skip dumping if it's a new puff
c
        if (lnewg(n) .and. .not.lkill(n)) goto 40
c
c-----Find layers containing top and bottom of puff
c
        kpb = 1       
        kpt = nlay
        do kk = 1,nlay
          if (height(i,j,kk) .GE. pufftop(n)) then
            kpt = kk     
            goto 16
          endif       
        enddo
  16    continue
        do kk = nlay,1,-1
          if (height(i,j,kk) .LE. puffbot(n)) then
            kpb = kk + 1
            goto 17
          endif
        enddo
  17    continue
c
c-----Revise vertical averaging to account for growth; define layer weighting 
c     as the depth of puff coverage weighted by density
c
        if (kpb .eq. kpt) then
          wtfac(kpb) = 1.
        else
          sumwt = 0.
          deplyr = pufftop(n) - height(i,j,kpt-1)
          wtfac(kpt) = deplyr*press(i,j,kpt)/tempk(i,j,kpt)
          sumwt = sumwt + wtfac(kpt)

          deplyr = height(i,j,kpb) - puffbot(n)
          wtfac(kpb)= deplyr*press(i,j,kpb)/tempk(i,j,kpb)
          sumwt = sumwt + wtfac(kpb)
c
          do kk = kpb+1,kpt-1
            deplyr = height(i,j,kk) - height(i,j,kk-1)
            wtfac(kk) = deplyr*press(i,j,kk)/tempk(i,j,kk)
            sumwt = sumwt + wtfac(kk)
          enddo

          do kk = kpb,kpt
            wtfac(kk) = wtfac(kk)/sumwt
          enddo
        endif
c
        do kk = 1,nlay
          deplyr = height(i,j,1)
          if (kk.gt.1) deplyr = height(i,j,kk) - height(i,j,kk-1)
          volvec(kk) = deplyr*dx(j+joff)*dy/mapscl(i,j)**2
          pcvfac(kk) = densfac*(273./tempk(i,j,kk))*
     &                         (press(i,j,kk)/1013.)
        enddo
c
c-----First check for a full dump
c     Slaughter the puff if its horizontal size exceeds AREAMX
c     Slaughter the puff if remaining mass fraction < 0.1
c     Slaughter the puff according to chemistry flags
c
        axisx = xlen + 3.*sigx(n)
        pufscl = pi*(axisx/2.)*(1.5*sigy(n))
c
        if (lkill(n) .OR. pufscl.GT.areamx .OR. fmspig(n).lt.0.1) then
          if (.not.lkill(n)) then
            if (pufscl.gt.areamx) then
              nkill(1) = nkill(1) + 1
            elseif (fmspig(n).lt.0.1) then
              nkill(2) = nkill(2) + 1
            endif
          endif
          ingrd(n) = 0
          ldump = .true.
c
          do is = nrad+1,nspec
            dumpmass(is) = 0.
            do nr = 1,nreactr
              dumpmass(is) = dumpmass(is) + puffmass(is,nr,n)
              puffmass(is,nr,n) = 0.
            enddo
          enddo
          pigage(igrd) = pigage(igrd) + agepigf(n)
          nage(igrd) = nage(igrd) + 1
c
c======================== Source Apportion Begin =======================
c
          if( ltrace .AND. (tectyp .EQ. RTRAC .OR.
     &                      tectyp .EQ. RTCMC) ) then
             do is = 1,nrtrac
               rtdump(is) = 0.
               do nr = 1,nreactr
                 rtdump(is) = rtdump(is) + puffrt(is,nr,n)
                 puffrt(is,nr,n) = 0.
               enddo
             enddo
          endif
c
c========================= Source Apportion End ========================
c
          goto 60
        endif
c
c-----If no full dump, check for leakage
c     (LDUMP and DELMS were set in PIGGROW)
c
        if (ldump) then
          do is = nrad+1,nspec
            dumpmass(is) = 0.
            do nr = 1,nreactr
              rctdmp = puffmass(is,nr,n)*delms
              dumpmass(is) = dumpmass(is) + rctdmp
              puffmass(is,nr,n) = puffmass(is,nr,n) - rctdmp
            enddo
          enddo
c
c======================== Source Apportion Begin =======================
c
          if( ltrace .AND. (tectyp .EQ. RTRAC .OR.
     &                      tectyp .EQ. RTCMC) ) then
             do is = 1,nrtrac
               rtdump(is) = 0.
               do nr = 1,nreactr
                 rctdmp = puffrt(is,nr,n)*delms
                 rtdump(is) = rtdump(is) + rctdmp
                 puffrt(is,nr,n) = puffrt(is,nr,n) - rctdmp
               enddo
             enddo
          endif
c
c========================= Source Apportion End ========================
c
        endif
c
c-----Update total dumped mass 
c
 60     if (ldump) then
          do is = nrad+1,nspec
            pigdump(is) = pigdump(is) + DBLE(dumpmass(is))
          enddo
          delconc = 0.0 ! must be initialized
c
c-----Dump DUMPMASS into grid; allocate mass in proportion to layer densities
c     among layers contained within this puff
c
          do 300 is = nrad+1,nspec
            if (dumpmass(is).eq.0.) goto 300
            do kk = kpb,kpt
              volfac = wtfac(kk)/volvec(kk)
              pconc(kk) = conc(i,j,kk,is) + dumpmass(is)*volfac
              if (is.le.ngas) then
                pbdnl(kk) = bdnl(is)*pcvfac(kk)
              else
                pbdnl(kk) = bdnl(is)
              endif
c
c-----Woops! DUMPMASS leads to negative grid concentration
c
              if (pconc(kk) .lt. pbdnl(kk)) then
                tmass = 0.
                bmass = 0.
                tvol = 0.
                do kkk = kpb,kpt
                  tmass = tmass + conc(i,j,kkk,is)*volvec(kkk)
                  tvol  = tvol + volvec(kkk)
                  if (is.le.ngas) then
                    pbdnl(kkk) = bdnl(is)*pcvfac(kkk)
                  else
                    pbdnl(kkk) = bdnl(is)
                  endif
                  bmass = bmass + pbdnl(kkk)*volvec(kkk)
                enddo
                resid = tmass + dumpmass(is)
c
c-----Try spreading any positive residual mass evenly over layers contained 
c     within puff depth (if there is enough grid mass)
c
                if (resid .ge. bmass) then
                  do kkk = kpb,kpt
                    pconc(kkk) = resid/tvol
                  enddo
                  goto 301
c
c-----Not enough grid mass for a postive residual. Set grid conc to lower
c     bound and add negative mass increment to tracking array
c
                else
                  negconc = pconc(kk) - pbdnl(kk)
                  pgmserr(is) = pgmserr(is) + DBLE(negconc*volvec(kk))
                  pconc(kk) = pbdnl(kk)
                endif
              endif
            enddo
 301        do kk = kpb,kpt
              delconc(kk,is) = pconc(kk) - conc(i,j,kk,is)
              conc(i,j,kk,is) = pconc(kk)
c
c======================== Process Analysis Begin ====================================
c
              if( lipr ) then
                if( i .GE. ia .AND. i .LE. iz .AND. j .GE. ja
     &                                         .AND. j .LE. jz ) then
                     if( ipa_cel(i+i0,j+j0,kk) .GT. 0 ) then
                       ipa_idx = ipa_cel(i+i0,j+j0,kk)
                       cipr(IPR_PIGEMIS, ipa_idx, is) =
     &                   cipr(IPR_PIGEMIS, ipa_idx, is) + delconc(kk,is)
                     endif
                endif
              endif
c
c========================= Process Analysis End ================================
c
            enddo
 300      continue
c
c======================== Source Apportion Begin =======================
c
          if( ltrace ) then
c
c-----Update RTRAC tracer concentration in the grid ---
c
            if( tectyp .EQ. RTRAC .OR. tectyp .EQ. RTCMC ) then
              do is = 1,nrtrac
                do kk = kpb,kpt
                  volfac = wtfac(kk)/volvec(kk)
                  idx = DBLE(ipsa3d-1) + DBLE(i) + 
     &               DBLE(m1)*DBLE(j-1) + DBLE(m1)*DBLE(m2)*DBLE(kk-1) + 
     &                            DBLE(m1)*DBLE(m2)*DBLE(nlay)*DBLE(is-1)
                  ptconc(idx) = ptconc(idx) + rtdump(is)*volfac
                enddo
              enddo
c
c-----Call routine to update SA tracer concentration in the grid ---
c
            else
              call pigdumpsa(m1,m2,m3,nrad,nspec,ntotsp,
     &                       i,j,kpb,kpt,n,delconc,ptconc(ipsa3d) )
            endif
          endif
c
c========================= Source Apportion End ========================
c
        endif
c
c-----Reset ldump and set LNEWG to false
c
 40     continue
        ldump = .false.
        lnewg(n) = .false.
c
c-----End GROWTH and DUMPING loop
c
 200  continue
c
      return
      end
