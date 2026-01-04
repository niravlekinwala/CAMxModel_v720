      subroutine chemdriv(m1,m2,m3,i0,j0,ia,iz,ja,jz,igrd,
     &                    ncol,nrow,nlay,dt,itzon,idfin,dx,dy,
     &                    cldtrns,water,tempk,tsurf,press,height,uwind,
     &                    vwind,cwc,pwr,cph,aph,ajno2,
     &                    cigwtr,conc,cellat,cellon,
     &                    snow,snowrat,snowfrc,snowalb,topo,fsurf,albedo,
     &                    lrdlai,lai,ldark,solmas,vegmas,reemis,l3dflag,
     &                    iptr2d,iptrsa,ipa_cel,iproc_id)
      use filunit
      use chmstry
      use o3colmap
      use bndary
      use camxcom
      use procan
      use tracer
      use rtracchm
      use rtcmcchm
c
c----CAMx v7.20 220430
c
c     CHEMDRIV performs chemistry on the current grid for one time step.
c     It calls one of two ODE solvers: EBISOLV or LSCALL.
c       EBISOLV is the driver for the EBI solver with Hertel's solutions
c       LSCALL  is the driver for the LSODE solver
c
c     Mechanism specific subroutines are passed to the driver routines.
c
c     Local array element con(nspec+1) is used to store the concentrations
c     of non-used species, i.e., species that are in the chem solver but
c     not on the species list for this run
c
c      Copyright 1996 - 2022
c     Ramboll 
c
c     Modifications:
c        1/9/02        A minor bug fix related to fast aero routine
c                      (units conversion)
c        1/15/02       Added code to handle RTRAC
c        01/30/02      Added code for RTRAC probing tool
c        10/18/02      Updated AEROCHEM routine; added CWC for aqueous
c                      PM chemistry
c        03/14/03      Implemented Hg chemistry (Prakash Karamchandani, AER)
c        03/21/03      Removed the OSAT technology type OPPAT
c        04/2/03       Removed option for UAM-V type cloud adjustment
c        07/21/03      Added optional snow effect on UV albedo
c        04/21/04      Added CMU sectional aerosol chemistry
c        10/14/04      Added Mechanism 10
c        10/05/05      Added PA calls for Mech 4
c        10/06/05      Removed Mech 2
c        10/07/05      Cloud water content below 273 K is partitioned to
c                      liquid fraction like in wet dep
c        10/27/05      Added cloud water pH as a 3-D argument
c        01/12/06      Aqueous PM now called each timestep
c        11/28/06      Fixed bug in assigning ldchm flag
c        07/04/07      Added heterogeneous hydrolysis of N2O5
c        07/11/07      Added code for RTCMC probing tool
c        12/26/06 -bkoo-     Removed AEROCHEM_AQ (merged into AEROCHEM_CF/CMU)
c        12/29/06 -bkoo-     Added SOACHEM
c        01/08/07 -bkoo-     Added Mechanism 6 (CB05)
c                            Fixed titration rxn # for SAPRC99 (PA)
c        07/04/07 -bkoo-     Added heterogeneous hydrolysis of N2O5
c        07/16/07 -bkoo-     Revised for HDDM
c        04/24/08 -gyarwood- Added EBI solver
c        06/01/08 -bkoo-     Modified call to EBI for DDM
c        06/11/08 -bkoo-     Revised for rate constant sensitivity
c        07/16/08 -bkoo-     Added DDM turn-off flag
c        01/30/09 -bkoo-     Removed the CMC fast solver
c        12/08/09 -jjung-    Added new RTRAC surface model and
c                            partitioning algorithm
c        05/27/09 -jjung-    Added Mechanism 7 (CB6)
c        7/14/10  -cemery-   Added code for in-line TUV cloud adjustment
c        06/02/11 -pkaramchandani- Updated Hg chemistry code
c        04/02/12 -cemery-   Removed RADM cloud adjustment option, cloud/aerosol
c                            adjustments now always done with in-line TUV; AHO
c                            file is now just ozone column; replaced haze
c                            dimension with terrain height in photo file
c        04/04/12 -cemery-   Introduced 4-dim interpolation of J values over
c                            zenith, height, albedo, terrain
c        4/12/12  -cemery-   Added T and P adjustments to photolysis rates
c        07/23/12 -unopmongcol- Added Mechanism 1 (CB6-I)
c        10/08/12 -jjung-    Added Mechanism 2 (CB6r1)
c        04/29/13 -cemery-   Added surface model
c        05/24/13 -cemery-   RTCMC photolysis can point to core reactions
c        06/28/13 -bkoo-     Changed the last argument to LSCALL
c        08/23/13 -bkoo-     Revised for DDM-PM
c        09/21/13 -bkoo-     Updated DDM-PM for KHETERO
c        10/14/13 -bkoo-     Added pressure to KHETERO call
c        01/03/14 -bkoo-     Increased SDDM array size for organic nitrate 
c                            hydrolysis
c        03/18/14 -bkoo-     Revised for benzene SOA
c        06/27/14 -bkoo-     Added Mechanism 3 (CB6r2h)
c        08/11/14 -cemery-   Albedo now includes effects of snow
c        08/22/14 -cemery-   Extended Zhang LU to surface model, added snow
c        10/09/14 -cemery-   Added subgrid convective model
c        11/11/14 -bkoo-     Added Mechanism 4 (CB6r3)
c        12/07/14 -bkoo-     Revised for VBS
c        09/11/15 -bkoo-     Revised for SA v3
c        01/08/16 -bkoo-     Updated for DDMRATE
c                                        Revised VBS
c                                        Removal of IEH solver
c        02/29/16 -bkoo-     Revised SDDM array for hetero N2O5 + HCL rxn
c        07/05/16 -gyarwood- Added HR_XO
c        07/20/16 -bkoo-     Revised SDDM array for hetero INTR hydrolysis
c        08/25/16 -bkoo-     Updated for new SOAP
c        09/02/16 -bkoo-     Set SDDM array size with # of het rxns
c        09/16/16 -pk-       Added LDARK to AEROCHEM_CF & AEROCHEM_CMU calls
c        11/28/16 -bkoo-     Revised VBS precursors
c        10/13/17 -bkoo-     Added APH & AERPH to track aerosol pH
c        10/16/17 -bkoo-     Added AJNO2 to track time-weighted NO2 photolysis rate
c                            Expanded DELCON to store SOA-specific changes
c        01/12/18 -bkoo-     Removed BNZA/TOLA/XYLA/ISP/TRP
c        08/09/18 -bkoo-     Revised for rate term sensitivity
c        02/21/20 -rlb-      Added Mechanism 1 (CB6r5)
c        03/24/21 -gy-       Added CB7 (Mech 7)
c        04/05/21 -cemery-   RTCMC uses scaling factor for photolysis reactions
c                            that point to core mech reactions
c        04/23/21 -cemery-   Udpated CH4 to 1.85 ppm
c
c     Input arguments:
c        igrd                grid index
c        ncol                number of columns
c        nrow                number of rows
c        nlay                number of layers
c        dt                  timestep (s)
c        itzon               time zone
c        idfin               map of nested grids in this grid
c        dx                  cell size in x-direction (m)
c        dy                  cell size in y-direction (m)
c        cldtrns             energy transmission coefficient (fraction)
c        water               water vapor (ppm)
c        tempk               temperature (K)
c        tsurf               surface temperature (K)
c        press               pressure (mb)
c        height              layer interface height (m)
c        uwind               u-component wind (m/s)
c        vwind               v-component wind (m/s)
c        cwc                 cloud water content (g/m3)
c        pwr                 precipitation water content (g/m3)
c        cph                 cloud water pH
c        aph                 aerosol pH
c        ajno2               time-weighted accumulated NO2 photolysis rate ([1/hr]*[hr])
c        cigwtr              subgrid cloud water (g/m3)
c        conc                species concentration (umol/m3,ug/m3)
c        cellat              cell centroid latitude (deg)
c        cellon              cell centroid longitude (deg)
c        snow                snow cover water equivalent (m)
c        snowrat             snow water accumulation rate (m/hr)
c        snowfrc             snow cover fraction
c        snowalb             snow albedo
c        topo                cell terrain altitude (m MSL)
c        fsurf               landuse fraction
c        albedo              surface albedo
c        lrdlai              LAI was read
c        lai                 gridded LAI
c        ldark               darkness flag (T=dark)
c        solmas              surface model soil mass (mol/ha)
c        vegmas              surface model vegetation mass (mol/ha)
c        l3dflag             save 3-D average concentrations
c        iptr2d              pointers into vectors for 2-D fields
c        iptrsa              pointers into vectors for tracer conc
c        ipa_cel             gridded array to identify if cell is
c                            in a IPRM sub-domain
c        iproc_id            process ID for this slice
c
c     Output arguments:
c        conc                species concentration (umol/m3,ug/m3)
c        solmas              surface model soil mass (mol/ha)
c        vegmas              surface model vegetation mass (mol/ha)
c        reemiss             surface model reemission rate (mol/s)
c
c     Routines called:
c        CHEM10
c        KTHERM
c        KHETERO
c        GETZNTH
c        KPHOTO
c        RECALIB
c        LOADDM
c        HDDMADJ
c        EBISOLV
c        LSCALL
c        HDDMCHEM
c        VBSCHEM
c        SOACHEM
c        LOADDMPM
c        AEROCHEM_CF
c        AEROCHEM_CMU
c        CALDATE
c        SRFMOD
c        HGGASCHEM
c        HGAQSCHEM
c        HGADSORB
c        O3PRDSA
c        OSATSA
c        PSATSA
c        REPARTSA
c        CHEMRT
c        PARTITIONRT
c        SRFMODRT
c        KRTC
c        KSCI
c        DLSDRV
c        SLSDRV
c        RBKDRV
c        CPAMECH
c        CPARAD
c        UPTIME
c        CAMXERR
c
c     Called by:
c        CHEMRXN
c
      implicit none
      include "camx.prm"
      include 'camx_aero.inc'
      include 'soap.inc'
      include 'vbs.inc'
      include 'deposit.inc'
      include 'flags.inc'
      include 'chmdbg.inc'
      include 'rtracsrf.inc'
c
c----Argument declarations
c
      integer :: m1,m2,m3,i0,j0,ia,iz,ja,jz

      integer igrd
      integer ncol
      integer nrow
      integer nlay
      real    dt
      integer itzon
      integer idfin(m1,m2)
      integer iproc_id
c
      real, dimension(m1,m2,m3) :: cldtrns,cph,aph,ajno2,
     &             water,tempk,press,cwc,pwr,uwind,vwind,cigwtr

      real, dimension(m1,m2,m3) :: height
c
      real    conc(m1,m2,m3,nspec)
c
      real    cellat(m1,m2)
      real    cellon(m1,m2)
      real    snow(m1,m2)
      real    snowrat(m1,m2)
      real    snowfrc(m1,m2)
      real    snowalb(m1,m2)
      real    topo(m1,m2)
      real    tsurf(m1,m2)
      real    fsurf(m1,m2,NLU)
      real    albedo(m1,m2)
      real    lai(m1,m2)
      real    dx(nrow)
      real    dy
c
      real solmas(m1,m2,nsmspc)
      real vegmas(m1,m2,nsmspc)
      real reemis(m1,m2,nsmspc)
c
      logical   ldark(m1,m2)
      logical   lrdlai
      logical   l3dflag
      integer   iptr2d
      integer*8 iptrsa
      integer   ipa_cel(ncol,nrow,nlay)
c
c======================== Process Analysis Begin =======================
c
      logical ldoirr, ldoipr, ltrapirr
      real    titrt
      real    chg0, chg2, thg20
c
      real rrxn_irr(MXREACT)
      real patmp(MXTRSP)
c
c======================== Process Analysis End =========================
c
c========================= Hg Chemistry Begin ==========================
c
      real cl2c              ! Cl2 concentration (ppm)
      real chgiip, chgiipc   ! Adsorbed Hg concentrations (ug/m3)
      real brc, broc         ! Br and BrO concentrations (ppm)
      real henso2,tfso2      ! Henry's Law for SO2 (M/atm)
      real heno3,tfo3        ! Henry's Law for O3 (M/atm)
      real pm10              ! Total PM10 concentrations (ug/m3)
      integer iocean         ! Ocean flag: 1 = ocean, 0 = other
c
      real thg2                 ! Total Hg(2)
      real convhgfac            ! conv factor for hg (ppm to ug/m3)
c
c-----variables for Hg adsorption on primary particles
      real, parameter :: MWNA = 23.0 ! Mol. wt. of sodium
      real, parameter :: MWCL = 35.5 ! Mol. wt. of chlorine atom
      real, parameter :: MWNACL = MWNA + MWCL ! Mol. wt. of seasalt
      real, parameter :: RHOSEAS = 2.0E6 ! Density of seasalt (g/m3)

      real dfine, dcoarse       ! diameter of fine and coarse particles (um)
c
c-----surface area concentration (m2/m3) and volume concentrations (cm3/m3) of
c-----fine and coarse particles
      real safine,sacoarse,vfine,vcoarse
c
      real cna,cpcl,vseas,fseas ! Na,Cl concs; vol. conc of seasalt, fraction
                                ! of fine PM that is seasalt
c
c========================== Hg Chemistry End ===========================
c
c-----variables for SOACHEM/VBSCHEM
c
      real avgox(NOXID)    ! Average oxidant conc. (ppm)
      real con_prec(NPREC) ! Precursor conc. (ppm)
      real con_CG(NCG)     ! CG conc. (ppm)
      real con_NV(NNV)     ! Non-volatile product conc. (ppm)

      real yfac(3,NPREC_S) ! 3 CG products per each precursor
c
c======================== DDM Begin =======================
c
      real sddm(MXTRSP+NHETRXN*MXFDDM)
      real avgsen(MXFDDM*NOXID)   ! Average sensitivities of oxidants (ppm)
      real sen_prec(MXFDDM*NPREC) ! Sensitivities of precursors (ppm)
      real sen_CG(MXFDDM*NCG)     ! Sensitivities of CGs (ppm)
      real sen_NV(MXFDDM*NNV)     ! Sensitivities of NVs (ppm)
c
c======================== DDM End =======================
c
c-----subroutine names to be called by EBISOLV and LSCALL
c
      external ebirxn1, ebirate1, hr_hox1, hr_nox1, hr_nxy1, hr_pan1
      external ebirxn2, ebirate2, hr_hox2, hr_nox2, hr_nxy2, hr_pan2
      external ebirxn3, ebirate3, hr_hox3, hr_nox3, hr_nxy3, hr_pan3
      external ebirxn4, ebirate4, hr_hox4, hr_nox4, hr_nxy4, hr_pan4
      external ebirxn5, ebirate5, hr_hox5, hr_nox5, hr_nxy5, hr_pan5
      external ebirxn6, ebirate6, hr_hox6, hr_nox6, hr_nxy6, hr_pan6
      external ebirxn7, ebirate7, hr_hox7, hr_nox7, hr_nxy7, hr_pan7
      external hr_xo1, hr_xo3, hr_xo4, hr_xo7, hr_xo_dum
      external lsrate1
      external lsrate2
      external lsrate3
      external lsrate4
      external lsrate5
      external lsrate6
      external lsrate7
      external ddmrate1, hddmjac1
      external ddmrate2, hddmjac2
      external ddmrate3, hddmjac3
      external ddmrate4, hddmjac4
      external ddmrate5, hddmjac5
      external ddmrate6, hddmjac6
      external ddmrate7, hddmjac7
c
c-----Local variables
c
      integer   i,j,k,is,l,ij,irxn,ihg,iht,nj,nz,nn
      integer*8 idx
      integer   ierr
      integer   aero_flag
      integer   iozon,irt_cel
      integer   ipa_idx
      integer   lsm
      integer   nitsum(MXLAYER)
      integer   icig
c
      real dtchem,dtfact
      real tcell,pcell,convfac,cldph,aerph,ctrns,zenith,cliq
      real o2,ch4,atm,h2
      real con(MXSPEC+1)
      real avgcnc(MXSPEC+1)
      real hght,alb,trn
      real desto3,prdo3n,prdo3v
      real depth,press0,ucomp,vcomp,wind
      real avit(MXCELLS,MXCELLS)
      real sumit(MXLAYER)
c
c --- the delcon array stores the delta concentrations at
c     various stages of the chemistry process.
c
c    delcon(1,*)  -- change in concentration after just gas phase
c                    chemistry (no aerosol chemistry applied)
c    delcon(2,*)  -- change in concentration after gas phase
c                    and aqueous phase chemistry
c    delcon(3,*)  -- change in concentration after all chemistry
c                    is complete
c    delcon(4,*)  -- special instance to capture the change in
c                    concentration due to gas phase and aqueous 
c                    phase but also includes the transfer of
c                    some nitrate species during the aqueous phase
c    delcon(5,*)  -- change in concentration due to SOA polymerization
c    delcon(6,*)  -- change in concentration due to particle-phase SOA photolysis
c    delcon(7,*)  -- change in concentration due to in-cloud SOA formation
c
      integer ispc
      real    o3old, o3new
c
      real delcon(7,MXSPEC+1)
      real cold(MXSPEC)
      real cnew(MXSPEC)
      real rtcon(MXTRSP)
      real soil(MXSMSPC)
      real vegy(MXSMSPC)
      real reemit(MXSMSPC)
      real flu(NLU)
c
c-----Entry point
c
      call flush(6)
      call flush(iout)
c
      dtchem = dt/3600.
      con(nspec+1) = 0.
      ltrapirr = .FALSE.
c
c-----Aerosol timing parameters:
c     time_aero : time to call aerosol routine (HHMM)
c     date_aero : date to call aerosol routine (YYJJJ)
c     dt_aero   : time interval between calls to aerosol routines (min)
c     aero_dt   : actual time interval for each grid (hr) (accumulated dtchem)
c
      aero_flag = 0
      if( naero.GT.0 .AND. (aeropt.EQ.'CF' .OR. aeropt.EQ.'CMU')) then
         aero_flag = 1
         aero_dt(igrd) = aero_dt(igrd) + dtchem
         if ( (timec(igrd)-time_aero(igrd)) .ge. -0.01 .and.
     &        (datec(igrd) .eq. date_aero(igrd)) ) then
           aero_flag = 2
         endif
      endif
c
      igrdchm = igrd
      do 91 k = 1,nlay
c
c---Update list of variables private to each thread for updated
c---Hg treatment, PK, Ramboll 06/11
c
        nitsum(k)=1
        sumit(k)=1.
c
c$omp parallel default(shared)
c$omp&  private(i,j,l,is,ispc,con,
c$omp&  ij,iozon,trn,alb,hght,ctrns,
c$omp&  zenith,delcon,tcell,pcell,cliq,cldph,aerph,atm,O2,H2,CH4,
c$omp&  sddm,ldoipr,ldoirr,ipa_idx,titrt,irxn,rrxn_irr,
c$omp&  ierr,convfac,convhgfac,patmp,idx,nn,cold,cnew,avgcnc,o3old,
c$omp&  o3new,irt_cel,iocean,brc,broc,cl2c,henso2,tfso2,heno3,tfo3,
c$omp&  pm10,iht,ihg,ltrapirr,chg0,chg2,chgiip,chgiipc,thg2,
c$omp&  thg20,cna,cpcl,vseas,fseas,dfine,dcoarse,vfine,vcoarse,
c$omp&  safine,sacoarse,avgox,con_prec,con_CG,con_NV,yfac,
c$omp&  avgsen,sen_prec,sen_CG,sen_NV,
c$omp&  nj,nz,rtcon,desto3,prdo3n,prdo3v,
c$omp&  depth,press0,vcomp,ucomp,wind,lsm,
c$omp&  flu,soil,vegy,reemit)
c$omp&  copyin(/ijkgrd/)
c
c$omp do schedule(guided) collapse(2)
c
        do 90 j =  ja,jz   !2,nrow-1
          do 89 i = ia,iz   !2,ncol-1
            ltrapirr = .FALSE.
            ichm = i
            jchm = j
            kchm = k
            avit(i,j) = 0.
c
c-----skip chemistry if fine grid exists in this cell
c
            if (idfin(i,j).gt.igrd) goto 89
c
c-----For the gas phase species (numbered 1 to ngas)
c     Pass concentration to CON, converting from umol/m3 to ppm
c
            tcell = tempk(i,j,k)
            pcell = press(i,j,k)
            convfac = densfac*(273./tcell)*(pcell/1013.)
            cldph = cph(i,j,k)
            aerph = aph(i,j,k)
            tchm = tcell
            wchm = water(i,j,k)
            do is = 1,nrad
              con(is) = conc(i,j,k,is)
              if (con(is).lt.0.) then
                write(iout,'(//,a)') 'ERROR in CHEMDRIV:'
                write(iout,*) 'Negative concentration before chem'
                write(iout,*) 'igrd, i, j, k = ', igrd,i,j,k
                do l = 1,nspec
                  write(iout,'(i3,2x,a7,e10.3)') l,spname(l),con(l)
                enddo
                call camxerr()
              endif
              con(is) = amax1(bdnl(is),con(is))
            enddo
            do is = nrad+1,ngas
              con(is) = conc(i,j,k,is)/convfac
              if (con(is).lt.0.) then
                write(iout,'(//,a)') 'ERROR in CHEMDRIV:'
                write(iout,*) 'Negative concentration before chem'
                write(iout,*) 'igrd, i, j, k = ', igrd,i,j,k
                do l = 1,nspec
                  write(iout,'(i3,2x,a7,e10.3)') l,spname(l),con(l)
                enddo
                call camxerr()
              endif
              con(is) = amax1(bdnl(is),con(is))
            enddo
c
c-----Load any aerosols
c
            if (ngas.lt.nspec) then
              do is=ngas+1,nspec
                con(is) = conc(i,j,k,is)
                if (con(is).lt.0.) then
                  write(iout,'(//,a)') 'ERROR in CHEMDRIV:'
                  write(iout,*) 'Negative concentration before chem'
                  write(iout,*) 'igrd, i, j, k = ', igrd,i,j,k
                  do l = 1,nspec
                    write(iout,'(i3,2x,a7,e10.3)') l,spname(l),con(l)
                  enddo
                  call camxerr()
                endif
                con(is) = amax1(bdnl(is),con(is))
              enddo
            endif
c
c-----Mechanism 10 is done here
c
            if (idmech.eq.10) then
               call chem10(dtchem,tcell,pcell,water(i,j,k),con)
               goto 88
            endif
c
c======================== DDM Begin ====================================
c
c  ---- load the sensitivities for this cell into 2-D array ---
c
            if ( (lddm .OR. lhddm) .AND. lddmcalc(igrd) ) then
              call loaddm(.FALSE.,m1,m2,m3,ntotsp,ptconc(iptrsa),
     &                    i,j,k,nddmsp,nspec,ngas,nrad,NHETRXN,sddm,
     &                    convfac)
              if ( lhddm )
     &              call hddmadj(.FALSE.,nddmsp,con,sddm)
            endif
c
c======================== DDM End =======================================
c
c-----Determine thermal rate constants
c
            call ktherm(tcell,pcell)
c
c-----Calculate the rate constant for heterogeneous hydrolysis of N2O5
c
            call khetero(tcell,pcell,water(i,j,k),con,
     &                           (lddm.AND.lddmcalc(igrd)),nddmsp,sddm)
c
c-----Load local values of ozone, terrain ht, albedo, zenith angle, altitude
c
            call getznth(cellat(i,j),cellon(i,j),timec(igrd),
     &                   datec(igrd),itzon,zenith,ldark(i,j))
            ldchm = ldark(i,j)
c
            alb   = albedo(i,j)
            ij = i + (j-1)*m1
            iozon = icdozn(iptr2d-1+ij)
            iocean = icdocn(iptr2d-1+ij)
c
            trn = topo(i,j)/1000.
            hght = height(i,j,k)/2000.
            if (k.gt.1) hght = (height(i,j,k) + height(i,j,k-1))/2000.
            ctrns = 0.
            if (.not.ldark(i,j)) ctrns = cldtrns(i,j,k)
c
c-----Determine photolysis rates through interpolation of look-up table
c
            call kphoto(iozon,alb,trn,hght,zenith,ctrns,ldark(i,j),
     &                  tcell,pcell)
c
c-----Accumulate JNO2 ([1/hr]*[hr])
c
            ajno2(i,j,k) = ajno2(i,j,k) + rk(jno2rxn) * dtchem
c
c======================== Source Apportion Begin =======================
c
c  --- store the current value of the species needed for tracer "chemistry" ----
c
            if( ltrace ) then
               o3old = con(ko3)
               do ispc=1,nrad
                  delcon(1,ispc) = -con(ispc)*convfac
                  delcon(2,ispc) = -con(ispc)*convfac
                  delcon(3,ispc) = -con(ispc)*convfac
                  delcon(4,ispc) = -con(ispc)*convfac
                  cold(ispc) = con(ispc)
               enddo
               do ispc=nrad+1,ngas
                  delcon(1,ispc) = -con(ispc)*convfac
                  delcon(2,ispc) = -con(ispc)*convfac
                  delcon(3,ispc) = -con(ispc)*convfac
                  delcon(4,ispc) = -con(ispc)*convfac
                  cold(ispc) = con(ispc)*convfac
               enddo
               do ispc=ngas+1,nspec
                  delcon(1,ispc) = -con(ispc)
                  delcon(2,ispc) = -con(ispc)
                  delcon(3,ispc) = -con(ispc)
                  delcon(4,ispc) = -con(ispc)
                  cold(ispc) = con(ispc)
               enddo
               delcon(5:7,:) = 0.0 ! reset to zero
               delcon(1:7,nspec+1) = 0.0 ! species not in mechanism
c
c  --- call routine to recalibrate so that tracer species stays
c      on track with model species: needed because the numerics
c      cause tracer species to drift
c
               if( tectyp .NE. RTRAC .AND. tectyp .NE. RTCMC )
     &                call recalib(m1,m2,m3,ntotsp,
     &                                      ptconc(iptrsa),i,j,k,cold)
c
c   --- if ozone or nitrate tracers are requested,
c       set the flag for process analysis -- this will save needed
c       information from the solver routines ---
c
c  --- set variables for OSAT3 ----
c
               if( lozone .OR. lnitrate ) then
                 ltrapirr = .TRUE.
                 do irxn=1,nreact
                     rrxn_irr(irxn) = 0.0
                 enddo
               endif
            endif
c
c========================= Source Apportion End ========================
c
c-----Assign SOA precursor concentrations before chemistry integration
c
            if( aero_flag .GT. 0 ) then
              if ( lvbs ) then
                if (idmech.eq.1 .or. idmech.eq.3 .or. idmech.eq.4) then ! cb6
                  con_prec(1) = con(kBENZ)
                  con_prec(2) = con(kTOL )
                  con_prec(3) = con(kXYL )
                  con_prec(4) = con(kISOP)
                  con_prec(5) = con(kTERP)
                  con_prec(6) = con(kSQT )
                elseif (idmech.eq.5) then ! saprc07tc
                  con_prec(1) = con(kBENZ)
                  con_prec(2) = con(kTOLU) + con(kARO1)
                  con_prec(3) = con(kOXYL) + con(kMXYL) + con(kPXYL)
     &                                     + con(kB124) + con(kARO2)
                  con_prec(4) = con(kISOP)
                  con_prec(5) = con(kAPIN) + con(kTERP)
                  con_prec(6) = con(kSESQ)
                elseif (idmech.eq.6) then ! cb05
                  con_prec(1) = 0.0
                  con_prec(2) = con(kTOL )
                  con_prec(3) = con(kXYL )
                  con_prec(4) = con(kISOP)
                  con_prec(5) = con(kTERP)
                  con_prec(6) = con(kSQT )
                elseif (idmech.eq.7) then ! cb7
                  con_prec(1) = con(kBENZ)
                  con_prec(2) = con(kTOL )
                  con_prec(3) = con(kXYL )
                  con_prec(4) = con(kISOP)
                  con_prec(5) = con(kAPIN) + con(kTERP)
                  con_prec(6) = con(kSQT )
                endif
                con_prec( 7) = con(kIVOG)
                con_prec( 8) = con(kIVOD)
                con_prec( 9) = con(kIVOA)
                con_prec(10) = con(kIVOB)
              else ! soap
                if (idmech.eq.1 .or. idmech.eq.3 .or. idmech.eq.4) then ! cb6
                  con_prec(1) = con(kBENZ)
                  con_prec(2) = con(kTOL )
                  con_prec(3) = con(kXYL )
                  con_prec(5) = con(kISOP)
                  con_prec(6) = con(kTERP)
                  con_prec(7) = con(kSQT )
                elseif (idmech.eq.5) then ! saprc07tc
                  con_prec(1) = con(kBENZ)
                  con_prec(2) = con(kTOLU) + con(kARO1)
                  con_prec(3) = con(kOXYL) + con(kMXYL) + con(kPXYL)
     &                                     + con(kB124) + con(kARO2)
                  con_prec(5) = con(kISOP)
                  con_prec(6) = con(kAPIN) + con(kTERP)
                  con_prec(7) = con(kSESQ)
                elseif (idmech.eq.6) then ! cb05
                  con_prec(1) = 0.0
                  con_prec(2) = con(kTOL )
                  con_prec(3) = con(kXYL )
                  con_prec(5) = con(kISOP)
                  con_prec(6) = con(kTERP)
                  con_prec(7) = con(kSQT )
                elseif (idmech.eq.7) then ! cb7
                  con_prec(1) = con(kBENZ)
                  con_prec(2) = con(kTOL )
                  con_prec(3) = con(kXYL )
                  con_prec(5) = con(kISOP)
                  con_prec(6) = con(kAPIN) + con(kTERP)
                  con_prec(7) = con(kSQT )
                endif
                con_prec(4) = con(kIVOA)
c
c======================== DDM Begin =======================
c
                if ( lddm ) then
                  if (idmech.eq.1 .or. idmech.eq.3 .or. idmech.eq.4) then ! cb6
                    do is = 1, nddmsp
                      sen_prec(         is) = sddm(iptddm(kBENZ)-1+is)
                      sen_prec(  nddmsp+is) = sddm(iptddm(kTOL )-1+is)
                      sen_prec(2*nddmsp+is) = sddm(iptddm(kXYL )-1+is)
                      sen_prec(4*nddmsp+is) = sddm(iptddm(kISOP)-1+is)
                      sen_prec(5*nddmsp+is) = sddm(iptddm(kTERP)-1+is)
                      sen_prec(6*nddmsp+is) = sddm(iptddm(kSQT )-1+is)
                    enddo
                  elseif (idmech.eq.5) then ! saprc07tc
                    do is = 1, nddmsp
                      sen_prec(         is) = sddm(iptddm(kBENZ)-1+is)
                      sen_prec(  nddmsp+is) = sddm(iptddm(kTOLU)-1+is)
     &                                      + sddm(iptddm(kARO1)-1+is)
                      sen_prec(2*nddmsp+is) = sddm(iptddm(kOXYL)-1+is)
     &                                      + sddm(iptddm(kMXYL)-1+is)
     &                                      + sddm(iptddm(kPXYL)-1+is)
     &                                      + sddm(iptddm(kB124)-1+is)
     &                                      + sddm(iptddm(kARO2)-1+is)
                      sen_prec(4*nddmsp+is) = sddm(iptddm(kISOP)-1+is)
                      sen_prec(5*nddmsp+is) = sddm(iptddm(kAPIN)-1+is)
     &                                      + sddm(iptddm(kTERP)-1+is)
                      sen_prec(6*nddmsp+is) = sddm(iptddm(kSESQ)-1+is)
                    enddo
                  elseif (idmech.eq.6) then ! cb05
                    do is = 1, nddmsp
                      sen_prec(         is) = 0.0
                      sen_prec(  nddmsp+is) = sddm(iptddm(kTOL )-1+is)
                      sen_prec(2*nddmsp+is) = sddm(iptddm(kXYL )-1+is)
                      sen_prec(4*nddmsp+is) = sddm(iptddm(kISOP)-1+is)
                      sen_prec(5*nddmsp+is) = sddm(iptddm(kTERP)-1+is)
                      sen_prec(6*nddmsp+is) = sddm(iptddm(kSQT )-1+is)
                    enddo
                  elseif (idmech.eq.7) then ! cb7
                    do is = 1, nddmsp
                      sen_prec(         is) = sddm(iptddm(kBENZ)-1+is)
                      sen_prec(  nddmsp+is) = sddm(iptddm(kTOL )-1+is)
                      sen_prec(2*nddmsp+is) = sddm(iptddm(kXYL )-1+is)
                      sen_prec(4*nddmsp+is) = sddm(iptddm(kISOP)-1+is)
                      sen_prec(5*nddmsp+is) = sddm(iptddm(kAPIN)-1+is)
     &                                      + sddm(iptddm(kTERP)-1+is)
                      sen_prec(6*nddmsp+is) = sddm(iptddm(kSQT )-1+is)
                    enddo
                  endif
                  do is = 1, nddmsp
                    sen_prec(3*nddmsp+is) = sddm(iptddm(kIVOA)-1+is)
                  enddo
                endif
c
c======================== DDM End =======================
c
              endif
            endif
c
c-----Chemistry integration, pass subroutines for mechanism used
c
            atm = 1.e6
            O2  = 2.095e5
            CH4 = 1.85
            H2  = 0.60
c
c======================== Process Analysis Begin =======================
c
            ldoipr = .FALSE.
            ldoirr = .FALSE.
            if( lproca ) then
               if( lipr ) then
                  if( ipa_cel(i+i0,j+j0,k) .GT. 0 )  then
                     ipa_idx = ipa_cel(i+i0,j+j0,k)
                     ldoipr = .TRUE.
                  endif
               endif
               if( lirr ) then
                  ltrapirr = .TRUE.
                  if( ipa_cel(i+i0,j+j0,k) .GT. 0 ) then
                     ipa_idx = ipa_cel(i+i0,j+j0,k)
                     ldoirr = .TRUE.
                  endif
                  titrt = 0.0
                  do irxn=1,nreact
                     rrxn_irr(irxn) = 0.0
                  enddo
               endif
            endif
c
c======================== Process Analysis End =========================
c
            if ( idsolv .EQ. IDEBI ) then
               if (idmech.eq.1) then
                   call ebisolv(ebirxn1,ebirate1,hr_hox1,hr_nox1,
     &                  hr_nxy1,hr_pan1,hr_xo1,
     &                  dtchem,water(i,j,k),atm,O2,
     &                  CH4,H2,tcell,ldark(i,j),con,avgcnc,avit(i,j),
     &                  hddmjac1,nddmsp,sddm,NOXID,avgsen,
     &                  (lddm.AND.lddmcalc(igrd)),
     &                  nirrrxn,rrxn_irr,ltrapirr)
               elseif (idmech.eq.3) then
                   call ebisolv(ebirxn3,ebirate3,hr_hox3,hr_nox3,
     &                  hr_nxy3,hr_pan3,hr_xo3,
     &                  dtchem,water(i,j,k),atm,O2,
     &                  CH4,H2,tcell,ldark(i,j),con,avgcnc,avit(i,j),
     &                  hddmjac3,nddmsp,sddm,NOXID,avgsen,
     &                  (lddm.AND.lddmcalc(igrd)),
     &                  nirrrxn,rrxn_irr,ltrapirr)
               elseif (idmech.eq.4) then
                   call ebisolv(ebirxn4,ebirate4,hr_hox4,hr_nox4,
     &                  hr_nxy4,hr_pan4,hr_xo4,
     &                  dtchem,water(i,j,k),atm,O2,
     &                  CH4,H2,tcell,ldark(i,j),con,avgcnc,avit(i,j),
     &                  hddmjac4,nddmsp,sddm,NOXID,avgsen,
     &                  (lddm.AND.lddmcalc(igrd)),
     &                  nirrrxn,rrxn_irr,ltrapirr)
               elseif (idmech.eq.5) then
                   call ebisolv(ebirxn5,ebirate5,hr_hox5,hr_nox5,
     &                  hr_nxy5,hr_pan5,hr_xo_dum,
     &                  dtchem,water(i,j,k),atm,O2,
     &                  CH4,H2,tcell,ldark(i,j),con,avgcnc,avit(i,j),
     &                  hddmjac5,nddmsp,sddm,NOXID,avgsen,
     &                  (lddm.AND.lddmcalc(igrd)),
     &                  nirrrxn,rrxn_irr,ltrapirr)
               elseif (idmech.eq.6) then
                   call ebisolv(ebirxn6,ebirate6,hr_hox6,hr_nox6,
     &                  hr_nxy6,hr_pan6,hr_xo_dum,
     &                  dtchem,water(i,j,k),atm,O2,
     &                  CH4,H2,tcell,ldark(i,j),con,avgcnc,avit(i,j),
     &                  hddmjac6,nddmsp,sddm,NOXID,avgsen,
     &                  (lddm.AND.lddmcalc(igrd)),
     &                  nirrrxn,rrxn_irr,ltrapirr)
               elseif (idmech.eq.7) then
                   call ebisolv(ebirxn7,ebirate7,hr_hox7,hr_nox7,
     &                  hr_nxy7,hr_pan7,hr_xo7,
     &                  dtchem,water(i,j,k),atm,O2,
     &                  CH4,H2,tcell,ldark(i,j),con,avgcnc,avit(i,j),
     &                  hddmjac7,nddmsp,sddm,NOXID,avgsen,
     &                  (lddm.AND.lddmcalc(igrd)),
     &                  nirrrxn,rrxn_irr,ltrapirr)
               endif
c
            elseif ( idsolv .EQ. IDLSOD ) then
               if (idmech.eq.1) then
                  call lscall(lsrate1,dtchem,water(i,j,k),atm,
     &                        O2,CH4,H2,con,avgcnc,nirrrxn,
     &                        rrxn_irr,ierr,ltrapirr,0)
               elseif (idmech.eq.3) then
                  call lscall(lsrate3,dtchem,water(i,j,k),atm,
     &                        O2,CH4,H2,con,avgcnc,nirrrxn,
     &                        rrxn_irr,ierr,ltrapirr,0)
               elseif (idmech.eq.4) then
                  call lscall(lsrate4,dtchem,water(i,j,k),atm,
     &                        O2,CH4,H2,con,avgcnc,nirrrxn,
     &                        rrxn_irr,ierr,ltrapirr,0)
               elseif (idmech.eq.5) then
                  call lscall(lsrate5,dtchem,water(i,j,k),atm,
     &                        O2,CH4,H2,con,avgcnc,nirrrxn,
     &                        rrxn_irr,ierr,ltrapirr,0)
               elseif (idmech.eq.6) then
                  call lscall(lsrate6,dtchem,water(i,j,k),atm,
     &                        O2,CH4,H2,con,avgcnc,nirrrxn,
     &                        rrxn_irr,ierr,ltrapirr,0)
               elseif (idmech.eq.7) then
                  call lscall(lsrate7,dtchem,water(i,j,k),atm,
     &                        O2,CH4,H2,con,avgcnc,nirrrxn,
     &                        rrxn_irr,ierr,ltrapirr,0)
               endif
            endif
c
c======================== HDDM Begin ====================================
c
            if ( lhddm .AND. lddmcalc(igrd) ) then
              if (idmech.eq.1) then
                call hddmchem(hddmjac1,ebirxn1,ebirate1,ddmrate1,
     &                        nreact,ngas,ngas,
     &                        nddmsp,nrateddm,ntermddm,nhddm,
     &                        iprate,ipterm,wfterm,iphddm,
     &                        dtchem,water(i,j,k),atm,O2,CH4,H2,
     &                        rk,avgcnc,sddm,ierr)
              elseif (idmech.eq.3) then
                call hddmchem(hddmjac3,ebirxn3,ebirate3,ddmrate3,
     &                        nreact,ngas,ngas,
     &                        nddmsp,nrateddm,ntermddm,nhddm,
     &                        iprate,ipterm,wfterm,iphddm,
     &                        dtchem,water(i,j,k),atm,O2,CH4,H2,
     &                        rk,avgcnc,sddm,ierr)
              elseif (idmech.eq.4) then
                call hddmchem(hddmjac4,ebirxn4,ebirate4,ddmrate4,
     &                        nreact,ngas,ngas,
     &                        nddmsp,nrateddm,ntermddm,nhddm,
     &                        iprate,ipterm,wfterm,iphddm,
     &                        dtchem,water(i,j,k),atm,O2,CH4,H2,
     &                        rk,avgcnc,sddm,ierr)
              elseif (idmech.eq.5) then
                call hddmchem(hddmjac5,ebirxn5,ebirate5,ddmrate5,
     &                        nreact,ngas,ngas,
     &                        nddmsp,nrateddm,ntermddm,nhddm,
     &                        iprate,ipterm,wfterm,iphddm,
     &                        dtchem,water(i,j,k),atm,O2,CH4,H2,
     &                        rk,avgcnc,sddm,ierr)
              elseif (idmech.eq.6) then
                call hddmchem(hddmjac6,ebirxn6,ebirate6,ddmrate6,
     &                        nreact,ngas,ngas,
     &                        nddmsp,nrateddm,ntermddm,nhddm,
     &                        iprate,ipterm,wfterm,iphddm,
     &                        dtchem,water(i,j,k),atm,O2,CH4,H2,
     &                        rk,avgcnc,sddm,ierr)
              elseif (idmech.eq.7) then
                call hddmchem(hddmjac7,ebirxn7,ebirate7,ddmrate7,
     &                        nreact,ngas,ngas,
     &                        nddmsp,nrateddm,ntermddm,nhddm,
     &                        iprate,ipterm,wfterm,iphddm,
     &                        dtchem,water(i,j,k),atm,O2,CH4,H2,
     &                        rk,avgcnc,sddm,ierr)
              endif
              if (ierr.ne.0) then
                write(iout,'(//,a)') 'ERROR in CHEMDRIV:'
                write(iout,*)
     &                'Zero determinant in SGEFA in HDDMCHEM at ',ierr
                write(iout,*) 'igrd, i, j, k = ', igrd,i,j,k
                call camxerr()
              endif
            endif
c
c======================== HDDM End =======================================
c
c-----Perform gas-phase chemistry for SOA precursors
c
            if( aero_flag .GT. 0 ) then

              avgox(1) = amax1(avgcnc(kO)  , bdnl(kO)  )
              avgox(2) = amax1(avgcnc(kOH) , bdnl(kOH) )
              avgox(3) = amax1(avgcnc(kO3) , bdnl(kO3) )
              avgox(4) = amax1(avgcnc(kNO3), bdnl(kNO3))
              avgox(5) = amax1(avgcnc(kNO) , bdnl(kNO) )
              avgox(6) = amax1(avgcnc(kHO2), bdnl(kHO2))

              if ( lvbs ) then

                con_CG( 1) = con(kVAS1)
                con_CG( 2) = con(kVAS2)
                con_CG( 3) = con(kVAS3)
                con_CG( 4) = con(kVAS4)
                con_CG( 5) = con(kVBS1)
                con_CG( 6) = con(kVBS2)
                con_CG( 7) = con(kVBS3)
                con_CG( 8) = con(kVBS4)
                con_CG( 9) = con(kVAP1)
                con_CG(10) = con(kVAP2)
                con_CG(11) = con(kVAP3)
                con_CG(12) = con(kVAP4)
                con_CG(13) = con(kVCP1)
                con_CG(14) = con(kVCP2)
                con_CG(15) = con(kVCP3)
                con_CG(16) = con(kVCP4)
                con_CG(17) = con(kVFP1)
                con_CG(18) = con(kVFP2)
                con_CG(19) = con(kVFP3)
                con_CG(20) = con(kVFP4)
                call vbschem(iout,tcell,pcell,dtchem,
     &                       avgox,con_prec(:NPREC_V),
     &                       con_CG(:NVBS*NSET),con_NV(:NSET))
                if (.NOT.((idmech.eq.5 .and. kSESQ.lt.nspec+1) .OR.
     &                    (idmech.eq.7 .and. kSQT.lt.nspec+1)))
     &            con(kSQT) = amax1( con_prec( 6), bdnl(kSQT) )
                con(kIVOG) = amax1( con_prec( 7), bdnl(kIVOG) )
                con(kIVOD) = amax1( con_prec( 8), bdnl(kIVOD) )
                con(kIVOA) = amax1( con_prec( 9), bdnl(kIVOA) )
                con(kIVOB) = amax1( con_prec(10), bdnl(kIVOB) )

                con(kVAS1) = con_CG( 1)
                con(kVAS2) = con_CG( 2)
                con(kVAS3) = con_CG( 3)
                con(kVAS4) = con_CG( 4)
                con(kVBS1) = con_CG( 5)
                con(kVBS2) = con_CG( 6)
                con(kVBS3) = con_CG( 7)
                con(kVBS4) = con_CG( 8)
                con(kVAP1) = con_CG( 9)
                con(kVAP2) = con_CG(10)
                con(kVAP3) = con_CG(11)
                con(kVAP4) = con_CG(12)
                con(kVCP1) = con_CG(13)
                con(kVCP2) = con_CG(14)
                con(kVCP3) = con_CG(15)
                con(kVCP4) = con_CG(16)
                con(kVFP1) = con_CG(17)
                con(kVFP2) = con_CG(18)
                con(kVFP3) = con_CG(19)
                con(kVFP4) = con_CG(20)

              else

                con_CG(1) = con(kCG1)
                con_CG(2) = con(kCG2)
                con_CG(3) = con(kCG3)
                con_CG(4) = con(kCG4)
c
c======================== DDM Begin =======================
c
                if ( lddm ) then
                  do is = 1, nddmsp
                    sen_CG(         is) = sddm(iptddm(kCG1)-1+is)
                    sen_CG(  nddmsp+is) = sddm(iptddm(kCG2)-1+is)
                    sen_CG(2*nddmsp+is) = sddm(iptddm(kCG3)-1+is)
                    sen_CG(3*nddmsp+is) = sddm(iptddm(kCG4)-1+is)
                  enddo
                endif
c
c======================== DDM End =======================
c
                call soachem(iout,tcell,pcell,dtchem,
     &                       avgox,con_prec(:NPREC_S),
     &                       con_CG(:NSOAP),con_NV(:NNVOL),yfac,
     &                       nddmsp,avgsen,sen_prec,sen_CG,sen_NV,lddm)
c
c======================== DDM Begin =======================
c
                if ( lddm ) then
                  if (.NOT.((idmech.eq.5 .and. kSESQ.lt.nspec+1) .OR.
     &                      (idmech.eq.7 .and. kSQT.lt.nspec+1))) then
                    do is = 1, nddmsp
                      sddm(iptddm(kSQT)-1+is) = sen_prec(6*nddmsp+is)
                    enddo
                  endif
                  do is = 1, nddmsp
                    sddm(iptddm(kIVOA)-1+is) = sen_prec(3*nddmsp+is)

                    sddm(iptddm(kCG1)-1+is) = sen_CG(         is)
                    sddm(iptddm(kCG2)-1+is) = sen_CG(  nddmsp+is)
                    sddm(iptddm(kCG3)-1+is) = sen_CG(2*nddmsp+is)
                    sddm(iptddm(kCG4)-1+is) = sen_CG(3*nddmsp+is)
                  enddo
                endif
c
c======================== DDM End =======================
c
                if (.NOT.((idmech.eq.5 .and. kSESQ.lt.nspec+1) .OR.
     &                    (idmech.eq.7 .and. kSQT.lt.nspec+1)))
     &            con(kSQT) = amax1( con_prec(7), bdnl(kSQT) )
                con(kIVOA) = amax1( con_prec(4), bdnl(kIVOA) )

                con(kCG1) = con_CG(1)
                con(kCG2) = con_CG(2)
                con(kCG3) = con_CG(3)
                con(kCG4) = con_CG(4)
              endif

            endif
c
c======================== Source Apportion Begin =======================
c
c  --- subtract the current values of the species from the before values
c      to get the delta values ---
c
            if( ltrace ) then
               do ispc=1,nrad
                 delcon(1,ispc) = delcon(1,ispc) +
     &                     AMAX1(bdnl(ispc),con(ispc)) * convfac
                 delcon(4,ispc) = delcon(4,ispc) +
     &                     AMAX1(bdnl(ispc),con(ispc)) * convfac
                 if( aeropt .NE. 'CF' .OR. aero_flag .EQ. 0 ) then
                     delcon(2,ispc) = delcon(2,ispc) +
     &                     AMAX1(bdnl(ispc),con(ispc)) * convfac
                 endif
               enddo
               do ispc=nrad+1,ngas
                 delcon(1,ispc) = delcon(1,ispc) +
     &                     AMAX1(bdnl(ispc),con(ispc)) * convfac
                 delcon(4,ispc) = delcon(4,ispc) +
     &                     AMAX1(bdnl(ispc),con(ispc)) * convfac
                 if( aeropt .NE. 'CF' .OR. aero_flag .EQ. 0 ) then
                     delcon(2,ispc) = delcon(2,ispc) +
     &                     AMAX1(bdnl(ispc),con(ispc)) * convfac
                 endif
               enddo
               do ispc=ngas+1,nspec
                 delcon(1,ispc) = delcon(1,ispc) +
     &                               AMAX1(bdnl(ispc),con(ispc))
                 delcon(4,ispc) = delcon(4,ispc) +
     &                               AMAX1(bdnl(ispc),con(ispc))
                 if( aeropt .NE. 'CF' .OR. aero_flag .EQ. 0 ) then
                     delcon(2,ispc) = delcon(2,ispc) +
     &                               AMAX1(bdnl(ispc),con(ispc))
                 endif
               enddo
            endif
c
c======================== Source Apportion End =======================
c
c======================= Process Analysis Begin ========================
c
            if( ldoipr ) then
c
c --- Chemistry change in umol/m3 units (except for HG chemistry) ----
c
               do is=1,nrad
                  cipr(IPR_CHEM, ipa_idx, is) =
     &               cipr(IPR_CHEM, ipa_idx, is) +
     &                           con(is)-conc(i,j,k,is)
               enddo
               do is=nrad+1,ngas
                  cipr(IPR_CHEM, ipa_idx, is) =
     &               cipr(IPR_CHEM, ipa_idx, is) +
     &                  con(is)*convfac-conc(i,j,k,is)
               enddo
            endif
c
c======================== Process Analysis End =========================
c
c
c========================= DDM Begin ===================================
c
c  --- put the changed values back into the gridded array ---
c
            if( (lddm .OR. lhddm) .AND. lddmcalc(igrd) ) then
              if( lhddm )
     &              call hddmadj(.TRUE.,nddmsp,con,sddm)
              call loaddm(.TRUE.,m1,m2,m3,ntotsp,ptconc(iptrsa),
     &                    i,j,k,nddmsp,nspec,ngas,nrad,NHETRXN,sddm,
     &                    convfac)
            endif
c
c========================= DDM End =====================================
c
c
c-----Perform aerosol chemistry if it's time
c
            if ( aero_flag.gt.0 ) then
              cliq = cwc(i,j,k)
              if (lcig(igrd)) cliq = amax1(cliq,cigwtr(i,j,k))
              if (tcell.lt.273.)
     &          cliq = amax1(0.,cliq*(tcell - tamin)/(273. - tamin))
              if( aeropt.eq.'CF' ) then
c
c======================== DDM Begin ====================================
c
c  ---- load the sensitivities for this cell into 2-D array ---
c
                if ( lddm ) then
                  call loaddmpm(.FALSE.,
     &                          m1,m2,m3,ntotsp,ptconc(iptrsa),
     &                          i,j,k,nddmsp,NSDDMPM,sddm,convfac)
                endif
c
c======================== DDM End =======================================
c
                icig = 0
                if (lcig(igrd)) icig = 2
                call aerochem_cf(water(i,j,k),tcell,pcell,cliq,cldph,
     &                           aerph,ajno2(i,j,k),con_NV,
     &                           con,convfac,dtchem,aero_dt(igrd),
     &                           delcon,sen_NV,sddm,nddmsp,NSDDMPM,
     &                           ldoipr,ipa_idx,aero_flag,.FALSE.,icig,
     &                           ldark(i,j))
c
c========================= DDM Begin ===================================
c
c  --- put the changed values back into the gridded array ---
c
                if ( lddm ) then
                  call loaddmpm(.TRUE.,
     &                          m1,m2,m3,ntotsp,ptconc(iptrsa),
     &                          i,j,k,nddmsp,NSDDMPM,sddm,convfac)
                endif
c
c========================= DDM End =====================================
c
              elseif( aeropt.eq.'CMU' ) then
                call aerochem_cmu(water(i,j,k),tcell,pcell,cliq,cldph,
     &                            con_NV,
     &                            MXSPEC,con,bdnl,convfac,
     &                            timec(igrd),dtchem,aero_dt(igrd),
     &                            aero_flag,ldark(i,j))
              endif
            endif
c
c-----Prepare for surface model
c
            if (k.eq.1 .and. lsrfmod) then
              ij = i + (j-1)*m1
              depth = height(i,j,1)
              press0 = press(i,j,1) - depth*(press(i,j,2) - 
     &                        press(i,j,1))/height(i,j,2)
              ucomp = (uwind(i,j,1) + uwind(i-1,j,1))/2.
              vcomp = (vwind(i,j,1) + vwind(i,j-1,1))/2.
              wind = sqrt(ucomp**2 + vcomp**2)
              wind = amax1(0.1,wind)
              do l = 1,NLU
                flu(l) = fsurf(i,j,l)
              enddo
              do lsm = 1,nsmspc
                soil(lsm) = solmas(i,j,lsm)
                vegy(lsm) = vegmas(i,j,lsm)
              enddo

              call srfmod(datec(igrd),dtchem,dx(j+j0),dy,cellat(i,j),
     &                    cellon(i,j),snow(i,j),snowrat(i,j),
     &                    snowfrc(i,j),snowalb(i,j),zenith,ctrns,
     &                    ldark(i,j),pwr(i,j,1),tempk(i,j,1),tsurf(i,j),
     &                    press(i,j,1),press0,wind,depth,nsmspc,
     &                    flu,lrdlai,lai(i,j),soil,vegy,reemit)

              do lsm = 1,nsmspc
                solmas(i,j,lsm) = soil(lsm)
                vegmas(i,j,lsm) = vegy(lsm)
                reemis(i,j,lsm) = reemit(lsm)
              enddo
            endif
c
c========================= Hg Chemistry Begin ==========================
c
c-----Call Hg chemistry modules
c
            if (khg0.lt.nspec+1) then
c
c-----Assign Cl2, Br and BrO concentrations based on height(km), ocean
c-----and day/night
c
              do ihg = 1,NHTHAL
                if (hght .le. hthal(ihg)) then
                  iht = ihg
                  go to 115
                endif
              enddo
              iht = NHTHAL
 115          continue
              brc = 0.
              broc = 0.
              cl2c =  0.0

              if (.not. ldark(i,j)) then
                if (iocean .eq. 1) then
                  brc  = brwprof(iht)
                  broc = browprof(iht)
                  cl2c = cl2day(iht)
                else
                  brc  = brlprof(iht)
                  broc = brolprof(iht)
                endif
              else
                cl2c = cl2nite(iht)
              endif
c
c-----Gas-phase chemistry
c
c
c-----Put adsorbed Hg(II) into gas-phase Hg(2) for chemistry calculations
c-----convert adsorbed Hg(II) to ppm before adding
c
              convhgfac = MWHG*convfac
              thg2 = con(khg2) + (con(khgiip)+con(khgiipc))/convhgfac
              thg20 = thg2 ! For process analysis 
              call hggaschem(con(khg0),thg2,con(ko3),
     &                       con(kh2o2),con(koh),brc,broc,
     &                       tcell,pcell,dt)
              con(khg0) = max(bdnl(khg0),con(khg0))
              thg2 = max(bdnl(khg2),thg2)
c
c======================= Process Analysis Begin ========================
c
              if( ldoipr ) then
                chg0 = con(khg0)*convfac ! store HG0 in umol/m3
                chg2 = thg2*convfac ! store total HG2 in umol/m3
                cipr(IPR_CHEM, ipa_idx, khg0) =
     &            cipr(IPR_CHEM, ipa_idx, khg0) +
     &              chg0 - conc(i,j,k,khg0)
                cipr(IPR_CHEM, ipa_idx, khg2) =
     &            cipr(IPR_CHEM, ipa_idx, khg2) +
     &              chg2 - thg20*convfac
              endif
c
c======================== Process Analysis End =========================
c
c-----Aqueous-phase chemistry (conditional)
c
              cliq = cwc(i,j,k)
              if (tcell.lt.273.) then
                cliq = amax1(0.,cliq*(tcell - tamin)/(273. - tamin))
              endif
              if (cliq.ge.cwmin .and. tcell.ge.tamin ) then
c
c-----Calculate total PM concentrations
c
                henso2 = henry0(kso2)
                tfso2  = tfact(kso2)
                heno3  = henry0(ko3)
                tfo3   = tfact(ko3)
                pm10   = 0.
                do is = ngas+1,nspec
                  if (is.ne.kph2o) pm10 = pm10 + con(is)
                enddo
                call hgaqschem(con(khg0),thg2,con(ko3),
     &                         con(kso2),con(koh),con(kho2),
     &                         con(khcl),cl2c,pm10,cldph,cliq,
     &                         tcell,pcell,dt,henso2,tfso2,
     &                         heno3,tfo3)
                con(khg0) = max(bdnl(khg0),con(khg0))
                thg2 = max(bdnl(khg2),thg2)
c
c======================= Process Analysis Begin ========================
c
                if( ldoipr ) then
                  cipr(IPR_AQCHEM, ipa_idx, khg0) =
     &              cipr(IPR_AQCHEM, ipa_idx, khg0) +
     &                con(khg0)*convfac - chg0
                  cipr(IPR_AQCHEM, ipa_idx, khg2) =
     &              cipr(IPR_AQCHEM, ipa_idx, khg2) +
     &                thg2*convfac - chg2
                endif
c
c======================== Process Analysis End =========================
c
              endif
c
c-----Adsorption of Hg(2) to particles
c-----Calculate volume concentrations of fine & coarse primary particles
              vfine   = 0.
              if(kpec  < nspec+1)vfine = vfine + con(kpec)/roprt(kpec)
              if(kpoa  < nspec+1)vfine = vfine + con(kpoa)/roprt(kpoa)
              if(kfprm < nspec+1)vfine = vfine + con(kfprm)/roprt(kfprm)

c-----Add sea-salt to vfine

              if(kna == nspec+1)then
                vseas = nacl/RHOSEAS
              else  ! If Na is a modeled species, then PCl is too
                cna = con(kna)/MWNA
                cpcl = con(kpcl)/MWCL
                if(cna < cpcl) then
                  vseas = cna*MWNACL/roprt(kna)
                else
                  vseas = cpcl*MWNACL/roprt(kpcl)
                end if
              end if
              vfine = vfine + vseas
              dfine = sqrt(dcut(kph2o,1)*dcut(kph2o,2)) !H2O is a reqd. species
c
              vcoarse = 0.
              if (kccrs < nspec+1) then
                vcoarse = vcoarse + con(kccrs)/roprt(kccrs)
                dcoarse = sqrt(dcut(kccrs,1)*dcut(kccrs,2))
              end if
              if (kcprm < nspec+1) then
                vcoarse = vcoarse + con(kcprm)/roprt(kcprm)
                dcoarse = sqrt(dcut(kcprm,1)*dcut(kcprm,2))
              end if
c
c-----calculate number conc, surface area conc of fine and coarse primary PM
              if (vfine > 0.)then
                fseas = vseas/vfine
                fseas = MAX(0.0,MIN(1.0,fseas))
                safine = 6.0 * vfine / dfine  ! Units: m2/m3
              else
                fseas = 1.0
                safine = 0.0
              end if
              if (vcoarse > 0.)then
                sacoarse = 6.0 * vcoarse / dcoarse  ! Units: m2/m3
              else
                sacoarse = 0.0
              end if
c
              if (safine > 0. .or. sacoarse > 0.) then
                chgiip  = 0.
                chgiipc = 0.
                chg2    = thg2 * convhgfac ! in ug/m3
                call hgadsorb(tcell,chg2,chgiip,chgiipc,fseas,safine,
     &                        sacoarse)
c
c-----assign total hg2 back to individual species
                con(khg2)    = chg2/convhgfac  ! back to ppm
                con(khg2)    = max(bdnl(khg2),con(khg2))
                con(khgiip)  = max(bdnl(khgiip),chgiip)
                con(khgiipc) = max(bdnl(khgiipc),chgiipc)
              else
                con(khg2)    = thg2
                con(khgiip)  = bdnl(khgiip)
                con(khgiipc) = bdnl(khgiipc)
              endif

            endif
c
c========================== Hg Chemistry End ===========================
c
c
c======================== Source Apportion Begin =======================
c
c  --- subtract the current values of the species from the before values
c      to get the delta values ---
c
            if( ltrace ) then
               do ispc=1,nrad
                  cnew(ispc) = con(ispc)*convfac
                  delcon(3,ispc) = delcon(3,ispc) +
     &                     AMAX1(bdnl(ispc),con(ispc)) * convfac
               enddo
               do ispc=nrad+1,ngas
                  cnew(ispc) = con(ispc)*convfac
                  delcon(3,ispc) = delcon(3,ispc) +
     &                     AMAX1(bdnl(ispc),con(ispc)) * convfac
               enddo
               do ispc=ngas+1,nspec
                  cnew(ispc) = con(ispc)
                  delcon(3,ispc) = delcon(3,ispc) +
     &                               AMAX1(bdnl(ispc),con(ispc))
               enddo
               o3new = AMAX1(bdnl(ko3),con(ko3))
c
               if( ldark(i,j) ) then
                 desto3 = delcon(2,ko3)
                 prdo3n = 0.0
                 prdo3v = 0.0
               else
                 if( lozone ) then
                    call o3prdsa(rrxn_irr,convfac,delcon(2,ko3),
     &                                         prdo3n,prdo3v,desto3)
                    if( delcon(2,ko3) .lt. desto3 ) then
                       desto3 = delcon(2,ko3)
                       prdo3n = 0.0
                       prdo3v = 0.0
                    endif
                 endif
               endif
c
c  --- rrxn_irr: ppm -> umol/m3
c
               if( ltrapirr ) then
                  do irxn=1,nreact
                     rrxn_irr(irxn) = rrxn_irr(irxn) * convfac
                  enddo
               endif
c
c  --- SA: call routine to do the tracer species "chemistry", just makes
c      adjustments for production or decay of the regular model species ----
c
               if( lozone .OR. lnitrate ) then
                  call osatsa(m1,m2,m3,igrd,i,j,k,prdo3n,prdo3v,desto3,
     &                        nspec,delcon,convfac,cold,rrxn_irr,dtchem)
               endif
c
               if( lsulfate .OR. lnitrate .OR. lsoa .OR. lprimary .OR.
     &             lmercury ) then
                  call psatsa(m1,m2,m3,igrd,i,j,k,nspec,delcon,
     &                                                   NPREC_S,yfac)
                  call repartsa(m1,m2,m3,igrd,i,j,k,cnew,delcon)
               endif
            endif
c
c  --- RTRAC chemistry
c
            if( tectyp .EQ. RTRAC ) then
               if (nrtherm .GT. 0 .OR. nrtphot .GT. 0 ) then
                 if( ipa_cel(i+i0,j+j0,k) .GT. 0 )  then
                    irt_cel = ipa_cel(i+i0,j+j0,k)
                 else
                    irt_cel = -9
                 endif
                 call chemrt(m1,m2,m3,igrd,i,j,k,i0,j0,pcell,tcell,cold,cnew,
     &                     avgcnc(koh),avgcnc(ko3),avgcnc(kno3),dtchem,irt_cel)
               endif
               if( lparttn )
     &            call partitionrt(m1,m2,m3,nspec,igrd,i,j,k,cnew)
               if( k.eq.1 .and. lsrfmodrt )
     &            call srfmodrt(m1,m2,igrd,i,j,dtchem,zenith,ctrns,
     &                          ldark(i,j))
c
c  --- RTCMC chemistry
c
            else if( tectyp .EQ. RTCMC ) then
c
c  --- set rate constants
c
               if( ktype .EQ. 1 ) then
                  call krtc(tcell,pcell)
               else
                  call ksci(tcell,pcell)
               endif
               if( njschm .GT. 0 ) then
                  do nj = 1,njschm
                     srkrtc(ijschm(nj)) = 0.
                     if (ijpntr(nj).gt.0) then
                       srkrtc(ijschm(nj)) = rk(ijpntr(nj))*jscale(nj)
                     else
                       do nz = nzschm,1,-1
                          if( zenith .LE. zenschm(nz) )
     &                       srkrtc(ijschm(nj)) = rjschm(nz,nj)*ctrns
                       enddo
                     endif
                     rkrtc(ijschm(nj))  = DBLE( srkrtc(ijschm(nj)) )
                  enddo
               endif
c
c  --- do chemistry
c
               do l = 1,ntotsp
                  idx = DBLE(iptrsa-1) + DBLE(i) + 
     &                 DBLE(m1)*DBLE(j-1) + DBLE(m1)*DBLE(m2)*DBLE(k-1) + 
     &                                 DBLE(m1)*DBLE(m2)*DBLE(m3)*DBLE(l-1)
                  rtcon(l) = ptconc(idx)/convfac
               enddo
               if( isolv .EQ. 1) then
                  call dlsdrv(dtchem,water(i,j,k),atm,O2,CH4,H2,
     &                           avgcnc,rtcon,ierr,.false.)
               elseif( isolv .EQ. 2) then
                  call slsdrv(dtchem,water(i,j,k),atm,O2,CH4,H2,
     &                           avgcnc,rtcon,ierr,.false.)
               elseif( isolv .EQ. 3) then
                  call rbkdrv(dtchem,water(i,j,k),atm,O2,CH4,H2,
     &                           avgcnc,rtcon,ierr,.false.)
               endif
               do l = 1,ntotsp
                  idx = DBLE(iptrsa-1) + DBLE(i) + 
     &                DBLE(m1)*DBLE(j-1) + DBLE(m1)*DBLE(m2)*DBLE(k-1) + 
     &                               DBLE(m1)*DBLE(m2)*DBLE(m3)*DBLE(l-1)
                  ptconc(idx) = rtcon(l)*convfac
               enddo
            endif
c
c========================= Source Apportion End ========================
c
c
c======================= Process Analysis Begin ========================
c
c --- Chemical reaction tracking in ppm units
c     Account for titration reaction in TRAP at night ---
c
            if( ldoirr ) then
              do irxn=1,nirrrxn
                cirr(ipa_idx, irxn) =
     &                   cirr(ipa_idx, irxn) + rrxn_irr(irxn)
              enddo
              if ( titrt .GT. 0.0 ) then
                if ( idmech.EQ.5 ) then
                  cirr(ipa_idx, 7) =
     &                     cirr(ipa_idx, 7) + titrt
                else
                  cirr(ipa_idx, 3) =
     &                     cirr(ipa_idx, 3) + titrt
                endif
              endif
            endif
c
c   --- calculate the chemical process analysis (CPA) variables
c       time weight CPA variables that are not cumulative by dtfact ---
c
            if( lirr .AND. (l3dflag .OR. k .EQ. 1) ) then
               nn = 0
               dtfact = dt/dtout/60.
               if( idmech.EQ.1 ) then
                 call cpamech1(rrxn_irr,rk,dtfact,MXREACT,patmp,ntotsp,
     &                                                    nn,ldark(i,j))
               elseif( idmech.EQ.3 ) then
                 call cpamech3(rrxn_irr,rk,dtfact,MXREACT,patmp,ntotsp,
     &                                                    nn,ldark(i,j))
               elseif( idmech.EQ.4 ) then
                 call cpamech4(rrxn_irr,rk,dtfact,MXREACT,patmp,ntotsp,
     &                                                    nn,ldark(i,j))
               elseif( idmech.EQ.5 ) then
                 call cpamech5(rrxn_irr,rk,dtfact,MXREACT,patmp,ntotsp,
     &                                                    nn,ldark(i,j))
               elseif( idmech.EQ.6 ) then
                 call cpamech6(rrxn_irr,rk,dtfact,MXREACT,patmp,ntotsp,
     &                                                    nn,ldark(i,j))
               elseif( idmech.EQ.7 ) then
                 call cpamech7(rrxn_irr,rk,dtfact,MXREACT,patmp,ntotsp,
     &                                                    nn,ldark(i,j))
               else
                 write(iout,'(//,a)') 'ERROR in CHEMDRIV:'
                 write(iout,'(2A,I2)')'Chemical Process Analysis (CPA)',
     &                   ' does not support chemical mechanism: ',idmech
                 write(iout,*)
     &                  'Please choose another mechanism and try again.'
                 call camxerr()
               endif
c
c  --- save the cloud adjustment to photolysis rates ---
c
               nn = nn + 1
               ptname(nn)  = 'J_CLDADJ'
               patmp(nn) = ctrns * dtfact
c
c  --- save the aerosol pH ---
c
               nn = nn + 1
               ptname(nn)  = 'AER_PH'
               patmp(nn) = aerph * dtfact
c    
c  --- save radical concentrations unless accumulating throught run ---
c
               if( .NOT. lcpacum)
     &           call cparad(con, ngas, patmp, ntotsp, nn, dtfact)
c
c  --- add CPA values to gridded array in native units ---
c
               do l=1,ntotsp
                 idx = DBLE(iptrsa-1) + DBLE(i) + 
     &                 DBLE(m1)*DBLE(j-1) + DBLE(m1)*DBLE(m2)*DBLE(k-1) + 
     &                              DBLE(m1)*DBLE(m2)*DBLE(nlay)*DBLE(l-1)
                 ptconc(idx) = ptconc(idx) + patmp(l)
               enddo
            endif
c
c======================== Process Analysis End =========================
c
c-----Pass CON back to concentration, convert gases from ppm to umol/m3
c
 88         continue
            cph(i,j,k) = cldph
            aph(i,j,k) = aerph
            do is=1,nrad
              conc(i,j,k,is) = con(is)
            enddo
            do is=nrad+1,ngas
              conc(i,j,k,is) = con(is)*convfac
            enddo
            if (ngas.lt.nspec) then
              do is=ngas+1,nspec
                conc(i,j,k,is) = amax1(con(is),bdnl(is))
              enddo
            endif
c
  89      continue
  90    continue
c
c$omp end parallel
c       sumit(k) = 1.
c       nitsum(k) = 1
c       do i = ia,iz
c         do j =  ja,jz
c           sumit(k) = sumit(k) + avit(i,j)
c           if (idfin(i,j).le.igrd ) nitsum(k) = nitsum(k) + 1
c         enddo
c       enddo
c
  91  continue
c     write(idiag,'(a,2i5,2f8.1)')
c    &   ('EBI ',igrd,k,timec(igrd),sumit(k)/nitsum(k),k=1,nlay)
c
c-----If aerosol chemistry was called
c     - reset aero_dt & ajno2
c     - increase time_aero (& date_aero) by dt_aero (or multiple of dt_aero)
c
      if (aero_flag.EQ.2) then
        aero_dt(igrd) = 0.0
        ajno2 = 0.0
  92    continue
        call uptime(time_aero(igrd),date_aero(igrd),60.*dt_aero)
        if ( time_aero(igrd) .le. timec(igrd) .and.
     &       date_aero(igrd) .eq. datec(igrd) ) goto 92
      endif
c
      return
      end
