      subroutine aerochem_cf(h2o,tempk,press,cwc,cph,aph,ajno2,
     &                       cNV,con,convfac,
     &                       dtaq,dtaer,delcon,scNV,sddm,nfam,nsen,
     &                       ldoipr,ipa_cel,aero_flag,lpig,icig,ldark)
      use camxcom
      use chmstry
      use filunit
      use tracer
      use procan
c
c----CAMx v7.20 220430
c
c     AEROCHEM_CF drives the CF 2-section aerosol model in CAMx.
c     It calculates the following chemical transformations:
c        1. Condensible organic gasses to organic aerosol (SOAP or VBS)
c        2. Gaseous sulfate to aerosol sulfate (from gas-phase chem)
c        3. Inorganic aqueous reactions and aqueous formation of SOA
c           (CMAQ approach)
c        4. Inorganic gas-aerosol equilibrium partitioning for the
c           ammonium/nitrate/sulfate/sodium/chloride system (ISORROPIA)
c     The gas species are in ppm
c     The aerosol species are in (ug/m3)
c
c      Copyright 1996 - 2022
c     Ramboll 
c
c     Modifications: 
c        1/9/02         Minor bug fixes (units conversions)
c        11/19/02       Incorporated RADM aqueous chemistry and ISORROPIA
c        12/9/02        Incorporated SOAP routine
c        10/24/05       Cloud pH is now passed in/out; Henry's Law values
c                       for SO2 and O3 no longer passed out
c        8/23/06        Removed metastable option in ISOROPIA
c        12/26/06       AQ-CHEM now called each chemistry timestep
c        12/29/06       Revised for the updated SOA scheme
c        01/10/07       Added species mapping for CB05/SAPRC99
c        05/04/07       Added array (wt_save) to keep original total
c                       mass because ISOROPIA overwrites with lower
c                       bound value.
c        05/30/07       Skip call to RAQCHEM under high acid conditions.
c        05/13/08       Modified to skip ISOROPIA for certain PiG calls.
c        01/21/10       Added formation of mineral dust nitrate
c        07/31/13       Modified for new GREASD PiG chemistry 
c        08/23/13       Revised for DDM-PM; replaced hard-wired MW's with wtmol
c        03/18/14       Revised for benzene SOA
c        10/09/14       Added subgrid convective model
c        12/07/14       Revised for VBS
c        01/08/16       Updated for Revised VBS
c        08/25/16       Updated for new SOAP
c        09/16/16       Aqueous chemistry updates
c        09/14/17       Reestablished metastable option in ISOROPIA
c        10/13/17       Modified to get aerosol pH from ISOROPIA
c        10/16/17       Added particle-phase SOA photolysis to SOAP
c        11/22/17       Added PFE/PMN/PK/PCA/PMG
c                       Updated for revised WTMOL
c                       Added EQSAM option
c        10/26/18       Updated DDM for SOA photolysis & revised RADM
c
c     Input arguments:
c        h2o                 cell water vapor (ppm)
c        tempk               cell temperature (K)
c        press               cell pressure (mb)
c        cwc                 cloud water content (g/m3)
c        cph                 cloud water pH 
c        aph                 aerosol pH
c        ajno2               time-weighted accumulated NO2 photolysis rate ([1/hr]*[hr])
c        cNV                 NV product concentrations (ppm)
c        con                 species concentrations (ppm, ug/m3)
c        convfac             conversion factor: umol/m3 = ppm * convfac
c        dtaq                aqueous chemistry time step (hours)
c        dtaer               general PM chemistry time step (hours)
c        delcon              delta concentrations for source apportionment
c        scNV                NV product sensitivities
c        sddm                DDM-PM sensitivity array
c        nfam                number of DDM parameters
c        nsen                number of model species used in DDM-PM
c        ldoipr              flag to determine whether to do IPR calculation
c        ipa_cel             pointer into PA sub-domain array
c        aero_flag           1 - do only RAQCHEM (with ISOROPIA);
c                            2 - do SOAP, RAQCHEM, ISOROPIA
c        lpig                T - PiG call (do only RAQCHEM)
c        icig                0 - CiG not activated (default to aero_flag)
c                            1 - CiG call (do only RAQCHEM)
c                            2 - Grid call, CiG called prior (skip RAQCHEM)
c        ldark               darkness flag
c
c     Output arguments:
c        con                 species concentrations (ppm, ug/m3)
c        cph                 cloud water pH 
c        aph                 aerosol pH
c        delcon              delta concentrations for source apportionment
c        sddm                DDM-PM sensitivity array
c
c     Routines called: 
c        SOAP
c        RAQCHEM
c        ISOROPIA
c 
c     Called by: 
c        IRONCHEM
c        CIGDRIVE
c        CHEMDRIV 
c
      implicit none
      include 'camx.prm' 
      include 'camx_aero.inc'
      include 'soap.inc'
      include 'vbs.inc'
      include 'chmdbg.inc'
c
c========================== Process Analysis Begin =======================
c
      logical ldoipr
      integer ipa_cel
      real csulf,chno3,cnh3,cpso4,cpno3,cpnh4
c
c========================== Process Analysis End =======================
c
c========================= DDM Begin ===================================
c
      integer nfam, nsen, ierr, ip
      real scNV(nfam,NNV)
      real sddm(nfam,nsen)
c
c========================== DDM End ====================================
c
      real    h2o
      real    tempk
      real    press
      real    cwc
      real    cph
      real    aph
      real    ajno2
      real    cNV(NNV)
      real    con(*)
      real    convfac
      real    dtaq
      real    dtaer
      real    delcon(7,*)
      real    so2
      real    so4
      real    no3

      integer aero_flag
      integer icig

      logical lpig
      logical lnvoc
      logical lsvoc
      logical lradm
      logical lisor
      logical ldark ! darkness flag added for day-night difference in Fe_III
      logical ldosa

      integer ispc,is
      real qwatr,ev,es,rh
      real pres_pa
      real dt_sec
      real cold(MXSPEC), deln2o5
      real cw_kgm3
c
c-----Arrays for SOAP
c
      real soacon(NCG),cg(NCG),csatT(NCG)
      real preOM ! pre-existing OA mass (mole)
      real apoly(NSOAP), aphot(NSOAP+NNVOL)
      real fpoly, fphot
c
c-----Arrays for RADM aqueous chemistry
c
      real r_gas(15),r_aer(10),r_old(25)
      real r_1,r_2,r_3
      real csoac
c
c-----Variables for EQSAM4clim
c
      integer, parameter :: dp = SELECTED_REAL_KIND(12,307)
      real(dp) :: xTT,xAW,xWH2O,xYPa(5),xYMa(4)
      real(dp) :: xsPM,xaPM,xPMS,xPMT,xRHO,xVOL,xPH,xHp,xGF,
     &            xYPs(5),xYMs(4),xYG(4)
      real(dp), parameter :: Dd = 1.e-6_dp ! dry droplet diameter for the Kelvin effect (not used)
      integer, parameter :: jmeq = 1
c
c-----Variables for ISOROPIA need double precision
c
      integer imask(1,1)
      real*8 wi(5),wt(5),wt_save(5)
      real*8 rhi,tempi,cntrl(2)
      real*8 gasis(3),aerliq(12),aersld(9),other(6)
      character*15 scasi
c
c-----Variables for calcium nitrate
c
      real*8 dno3 ! calcium nitrate (mol/m3)
      logical lno3_isr ! flag to indicate non-zero nitrate for isorropia
c
c-----Variables for iron and manganese fractions in dust and other primary PM
c
      real fe_iii, mn_ii, fe_sol, mn_sol, fe_mmr, mn_mmr
c
c-----Speciation factors for elemental species (negative to fall back to background values)
c
      real, parameter :: f_FE_PRM = -1.0   , f_FE_CRS = -1.0    ! mass fraction of Fe in PRM or CRS
      real, parameter :: f_MN_PRM = -1.0   , f_MN_CRS = -1.0    ! mass fraction of Mn in PRM or CRS
      real, parameter :: f_CA_PRM = 0.0260 , f_CA_CRS = 0.0260  ! mass fraction of Ca in PRM or CRS
      real, parameter :: f_MG_PRM = -1.0   , f_MG_CRS = -1.0    ! mass fraction of Mn in PRM or CRS
      real, parameter :: f_K_PRM  = -1.0   , f_K_CRS  = -1.0    ! mass fraction of K in PRM or CRS
c 
c-----Entry point 
c
      lnvoc = .false.
      lsvoc = .false.
      lradm = .false.
      lisor = .false.
      ldosa = .false.

      if (.not.lpig .and. icig.ne.1) then
        lnvoc = .true.
        if (aero_flag.eq.2) then
          lsvoc = .true.
          lisor = .true.
        endif
        ldosa = ltrace
      endif
      if (cwc.ge.cwmin .and. tempk.gt.tamin) then
        so2 = con(kso2)*1.e-6
        so4 = (con(kpso4)/wtmol(1)/convfac)*1.e-6
        no3 = (con(kpno3)/wtmol(3)/convfac)*1.e-6
        if (so4.lt.5.0e-7 .and. no3.lt.1.5e-6 .and. so2.lt.1.0e-6) then
          lradm = .true.
          lisor = .true.
        endif
        if (lpig) lisor = .false.
        if (icig.eq.1) lisor = .false.
        if (icig.eq.2) lradm = .false.
      endif

      do ispc = 1,nspec
         cold(ispc) = con(ispc)
      enddo
      con(nspec+1) = 0.
      cph = 5.0
c
c-----Calculate relative humidity
c
      qwatr = 1.e-6*h2o*18./28.8 
      ev = qwatr*press/(qwatr + eps) 
      es = e0*exp((lv/rv)*(1./273. - 1./tempk)) 
      rh = 100.0*amin1(0.99,ev/es)
c
c-----Partitioning of non-volatile organic compounds.
c     Unit conversion is done here: [ppm] -> [ug/m3]
c
      if ( lnvoc ) then

        if ( lvbs ) then

          con(kPAS0) = con(kPAS0) + cNV(1) * convfac * mwPAS0
          con(kPBS0) = con(kPBS0) + cNV(2) * convfac * mwPBS0
          con(kPAP0) = con(kPAP0) + cNV(3) * convfac * mwPAP0
          con(kPCP0) = con(kPCP0) + cNV(4) * convfac * mwPCP0
          con(kPFP0) = con(kPFP0) + cNV(5) * convfac * mwPFP0
          if ( ldoipr ) then
            cipr(IPR_OAERO,ipa_cel,kPAS0) = cipr(IPR_OAERO,ipa_cel,kPAS0)
     &                                           + con(kPAS0)-cold(kPAS0)
            cipr(IPR_OAERO,ipa_cel,kPBS0) = cipr(IPR_OAERO,ipa_cel,kPBS0)
     &                                           + con(kPBS0)-cold(kPBS0)
            cipr(IPR_OAERO,ipa_cel,kPAP0) = cipr(IPR_OAERO,ipa_cel,kPAP0)
     &                                           + con(kPAP0)-cold(kPAP0)
            cipr(IPR_OAERO,ipa_cel,kPCP0) = cipr(IPR_OAERO,ipa_cel,kPCP0)
     &                                           + con(kPCP0)-cold(kPCP0)
            cipr(IPR_OAERO,ipa_cel,kPFP0) = cipr(IPR_OAERO,ipa_cel,kPFP0)
     &                                           + con(kPFP0)-cold(kPFP0)
          endif

        else ! SOAP

          con(ksopa) = con(ksopa) + cNV(1) * convfac * mwsopa
          con(ksopb) = con(ksopb) + cNV(2) * convfac * mwsopb
          if ( ldoipr ) then
            cipr(IPR_OAERO,ipa_cel,ksopa) = cipr(IPR_OAERO,ipa_cel,ksopa)
     &                                         + con(ksopa) - cold(ksopa)
            cipr(IPR_OAERO,ipa_cel,ksopb) = cipr(IPR_OAERO,ipa_cel,ksopb)
     &                                         + con(ksopb) - cold(ksopb)
          endif
c
c========================= DDM Begin ===================================
c
          if (lddm) then
            do ip = 1, nfam
              sddm(ip,jdpSOPA) = sddm(ip,jdpSOPA) + scNV(ip,1) * convfac * mwsopa
              sddm(ip,jdpSOPB) = sddm(ip,jdpSOPB) + scNV(ip,2) * convfac * mwsopb
            enddo
          endif
c
c========================== DDM End ====================================
c
        endif ! VBS or SOAP?

      endif ! lnvoc?
c
c-----Partitioning of semi-volatile organic compounds between the gas and
c     aerosol phases using the SOAP or VBS scheme.
c     CGs are in ppm; SOAs are in ug/m3; the unit conversion is done here.
c
      if ( lsvoc ) then

        if ( lvbs ) then

          cg(1)  = con(kVAS1) * convfac * mwVBS(1)
          cg(2)  = con(kVAS2) * convfac * mwVBS(2)
          cg(3)  = con(kVAS3) * convfac * mwVBS(3)
          cg(4)  = con(kVAS4) * convfac * mwVBS(4)
          cg(5)  = con(kVBS1) * convfac * mwVBS(5)
          cg(6)  = con(kVBS2) * convfac * mwVBS(6)
          cg(7)  = con(kVBS3) * convfac * mwVBS(7)
          cg(8)  = con(kVBS4) * convfac * mwVBS(8)
          cg(9)  = con(kVAP1) * convfac * mwVBS(9)
          cg(10) = con(kVAP2) * convfac * mwVBS(10)
          cg(11) = con(kVAP3) * convfac * mwVBS(11)
          cg(12) = con(kVAP4) * convfac * mwVBS(12)
          cg(13) = con(kVCP1) * convfac * mwVBS(13)
          cg(14) = con(kVCP2) * convfac * mwVBS(14)
          cg(15) = con(kVCP3) * convfac * mwVBS(15)
          cg(16) = con(kVCP4) * convfac * mwVBS(16)
          cg(17) = con(kVFP1) * convfac * mwVBS(17)
          cg(18) = con(kVFP2) * convfac * mwVBS(18)
          cg(19) = con(kVFP3) * convfac * mwVBS(19)
          cg(20) = con(kVFP4) * convfac * mwVBS(20)
          soacon(1)  = con(kPAS1)
          soacon(2)  = con(kPAS2)
          soacon(3)  = con(kPAS3)
          soacon(4)  = con(kPAS4)
          soacon(5)  = con(kPBS1)
          soacon(6)  = con(kPBS2)
          soacon(7)  = con(kPBS3)
          soacon(8)  = con(kPBS4)
          soacon(9)  = con(kPAP1)
          soacon(10) = con(kPAP2)
          soacon(11) = con(kPAP3)
          soacon(12) = con(kPAP4)
          soacon(13) = con(kPCP1)
          soacon(14) = con(kPCP2)
          soacon(15) = con(kPCP3)
          soacon(16) = con(kPCP4)
          soacon(17) = con(kPFP1)
          soacon(18) = con(kPFP2)
          soacon(19) = con(kPFP3)
          soacon(20) = con(kPFP4)
          preOM = con(kPAS0) / mwPAS0
     &          + con(kPBS0) / mwPBS0
     &          + con(kPAP0) / mwPAP0
     &          + con(kPCP0) / mwPCP0
     &          + con(kPFP0) / mwPFP0

          do ispc = 1, NVBS*NSET
            csatT(ispc) = csatV(ispc) * (cstempV(ispc)/tempk)
     &                                * exp((deltahV(ispc)/8.314)
     &                                * (1/cstempV(ispc)-1/tempk))
          enddo
          call soap(NSET*NVBS,mwVBS,soacon(:NSET*NVBS),cg(:NSET*NVBS),
     &              csatT(:NSET*NVBS),preOM,
     &              iout,igrdchm,ichm,jchm,kchm)

          con(kVAS1) = amax1(cg(1) /convfac/mwVBS(1) ,bdnl(kVAS1))
          con(kVAS2) = amax1(cg(2) /convfac/mwVBS(2) ,bdnl(kVAS2))
          con(kVAS3) = amax1(cg(3) /convfac/mwVBS(3) ,bdnl(kVAS3))
          con(kVAS4) = amax1(cg(4) /convfac/mwVBS(4) ,bdnl(kVAS4))
          con(kVBS1) = amax1(cg(5) /convfac/mwVBS(5) ,bdnl(kVBS1))
          con(kVBS2) = amax1(cg(6) /convfac/mwVBS(6) ,bdnl(kVBS2))
          con(kVBS3) = amax1(cg(7) /convfac/mwVBS(7) ,bdnl(kVBS3))
          con(kVBS4) = amax1(cg(8) /convfac/mwVBS(8) ,bdnl(kVBS4))
          con(kVAP1) = amax1(cg(9) /convfac/mwVBS(9) ,bdnl(kVAP1))
          con(kVAP2) = amax1(cg(10)/convfac/mwVBS(10),bdnl(kVAP2))
          con(kVAP3) = amax1(cg(11)/convfac/mwVBS(11),bdnl(kVAP3))
          con(kVAP4) = amax1(cg(12)/convfac/mwVBS(12),bdnl(kVAP4))
          con(kVCP1) = amax1(cg(13)/convfac/mwVBS(13),bdnl(kVCP1))
          con(kVCP2) = amax1(cg(14)/convfac/mwVBS(14),bdnl(kVCP2))
          con(kVCP3) = amax1(cg(15)/convfac/mwVBS(15),bdnl(kVCP3))
          con(kVCP4) = amax1(cg(16)/convfac/mwVBS(16),bdnl(kVCP4))
          con(kVFP1) = amax1(cg(17)/convfac/mwVBS(17),bdnl(kVFP1))
          con(kVFP2) = amax1(cg(18)/convfac/mwVBS(18),bdnl(kVFP2))
          con(kVFP3) = amax1(cg(19)/convfac/mwVBS(19),bdnl(kVFP3))
          con(kVFP4) = amax1(cg(20)/convfac/mwVBS(20),bdnl(kVFP4))
          con(kPAS1) = amax1(soacon(1) ,bdnl(kPAS1))
          con(kPAS2) = amax1(soacon(2) ,bdnl(kPAS2))
          con(kPAS3) = amax1(soacon(3) ,bdnl(kPAS3))
          con(kPAS4) = amax1(soacon(4) ,bdnl(kPAS4))
          con(kPBS1) = amax1(soacon(5) ,bdnl(kPBS1))
          con(kPBS2) = amax1(soacon(6) ,bdnl(kPBS2))
          con(kPBS3) = amax1(soacon(7) ,bdnl(kPBS3))
          con(kPBS4) = amax1(soacon(8) ,bdnl(kPBS4))
          con(kPAP1) = amax1(soacon(9) ,bdnl(kPAP1))
          con(kPAP2) = amax1(soacon(10),bdnl(kPAP2))
          con(kPAP3) = amax1(soacon(11),bdnl(kPAP3))
          con(kPAP4) = amax1(soacon(12),bdnl(kPAP4))
          con(kPCP1) = amax1(soacon(13),bdnl(kPCP1))
          con(kPCP2) = amax1(soacon(14),bdnl(kPCP2))
          con(kPCP3) = amax1(soacon(15),bdnl(kPCP3))
          con(kPCP4) = amax1(soacon(16),bdnl(kPCP4))
          con(kPFP1) = amax1(soacon(17),bdnl(kPFP1))
          con(kPFP2) = amax1(soacon(18),bdnl(kPFP2))
          con(kPFP3) = amax1(soacon(19),bdnl(kPFP3))
          con(kPFP4) = amax1(soacon(20),bdnl(kPFP4))
c
c========================= Process Analysis Begin ======================
c
          if ( ldoipr ) then
            cipr(IPR_OAERO,ipa_cel,kVAS1)=cipr(IPR_OAERO,ipa_cel,kVAS1)
     &                               + (con(kVAS1)-cold(kVAS1))*convfac
            cipr(IPR_OAERO,ipa_cel,kVAS2)=cipr(IPR_OAERO,ipa_cel,kVAS2)
     &                               + (con(kVAS2)-cold(kVAS2))*convfac
            cipr(IPR_OAERO,ipa_cel,kVAS3)=cipr(IPR_OAERO,ipa_cel,kVAS3)
     &                               + (con(kVAS3)-cold(kVAS3))*convfac
            cipr(IPR_OAERO,ipa_cel,kVAS4)=cipr(IPR_OAERO,ipa_cel,kVAS4)
     &                               + (con(kVAS4)-cold(kVAS4))*convfac
            cipr(IPR_OAERO,ipa_cel,kVBS1)=cipr(IPR_OAERO,ipa_cel,kVBS1)
     &                               + (con(kVBS1)-cold(kVBS1))*convfac
            cipr(IPR_OAERO,ipa_cel,kVBS2)=cipr(IPR_OAERO,ipa_cel,kVBS2)
     &                               + (con(kVBS2)-cold(kVBS2))*convfac
            cipr(IPR_OAERO,ipa_cel,kVBS3)=cipr(IPR_OAERO,ipa_cel,kVBS3)
     &                               + (con(kVBS3)-cold(kVBS3))*convfac
            cipr(IPR_OAERO,ipa_cel,kVBS4)=cipr(IPR_OAERO,ipa_cel,kVBS4)
     &                               + (con(kVBS4)-cold(kVBS4))*convfac
            cipr(IPR_OAERO,ipa_cel,kVAP1)=cipr(IPR_OAERO,ipa_cel,kVAP1)
     &                               + (con(kVAP1)-cold(kVAP1))*convfac
            cipr(IPR_OAERO,ipa_cel,kVAP2)=cipr(IPR_OAERO,ipa_cel,kVAP2)
     &                               + (con(kVAP2)-cold(kVAP2))*convfac
            cipr(IPR_OAERO,ipa_cel,kVAP3)=cipr(IPR_OAERO,ipa_cel,kVAP3)
     &                               + (con(kVAP3)-cold(kVAP3))*convfac
            cipr(IPR_OAERO,ipa_cel,kVAP4)=cipr(IPR_OAERO,ipa_cel,kVAP4)
     &                               + (con(kVAP4)-cold(kVAP4))*convfac
            cipr(IPR_OAERO,ipa_cel,kVCP1)=cipr(IPR_OAERO,ipa_cel,kVCP1)
     &                               + (con(kVCP1)-cold(kVCP1))*convfac
            cipr(IPR_OAERO,ipa_cel,kVCP2)=cipr(IPR_OAERO,ipa_cel,kVCP2)
     &                               + (con(kVCP2)-cold(kVCP2))*convfac
            cipr(IPR_OAERO,ipa_cel,kVCP3)=cipr(IPR_OAERO,ipa_cel,kVCP3)
     &                               + (con(kVCP3)-cold(kVCP3))*convfac
            cipr(IPR_OAERO,ipa_cel,kVCP4)=cipr(IPR_OAERO,ipa_cel,kVCP4)
     &                               + (con(kVCP4)-cold(kVCP4))*convfac
            cipr(IPR_OAERO,ipa_cel,kVFP1)=cipr(IPR_OAERO,ipa_cel,kVFP1)
     &                               + (con(kVFP1)-cold(kVFP1))*convfac
            cipr(IPR_OAERO,ipa_cel,kVFP2)=cipr(IPR_OAERO,ipa_cel,kVFP2)
     &                               + (con(kVFP2)-cold(kVFP2))*convfac
            cipr(IPR_OAERO,ipa_cel,kVFP3)=cipr(IPR_OAERO,ipa_cel,kVFP3)
     &                               + (con(kVFP3)-cold(kVFP3))*convfac
            cipr(IPR_OAERO,ipa_cel,kVFP4)=cipr(IPR_OAERO,ipa_cel,kVFP4)
     &                               + (con(kVFP4)-cold(kVFP4))*convfac
            cipr(IPR_OAERO,ipa_cel,kPAS1)=cipr(IPR_OAERO,ipa_cel,kPAS1)
     &                                  + con(kPAS1)-cold(kPAS1)
            cipr(IPR_OAERO,ipa_cel,kPAS2)=cipr(IPR_OAERO,ipa_cel,kPAS2)
     &                                  + con(kPAS2)-cold(kPAS2)
            cipr(IPR_OAERO,ipa_cel,kPAS3)=cipr(IPR_OAERO,ipa_cel,kPAS3)
     &                                  + con(kPAS3)-cold(kPAS3)
            cipr(IPR_OAERO,ipa_cel,kPAS4)=cipr(IPR_OAERO,ipa_cel,kPAS4)
     &                                  + con(kPAS4)-cold(kPAS4)
            cipr(IPR_OAERO,ipa_cel,kPBS1)=cipr(IPR_OAERO,ipa_cel,kPBS1)
     &                                  + con(kPBS1)-cold(kPBS1)
            cipr(IPR_OAERO,ipa_cel,kPBS2)=cipr(IPR_OAERO,ipa_cel,kPBS2)
     &                                  + con(kPBS2)-cold(kPBS2)
            cipr(IPR_OAERO,ipa_cel,kPBS3)=cipr(IPR_OAERO,ipa_cel,kPBS3)
     &                                  + con(kPBS3)-cold(kPBS3)
            cipr(IPR_OAERO,ipa_cel,kPBS4)=cipr(IPR_OAERO,ipa_cel,kPBS4)
     &                                  + con(kPBS4)-cold(kPBS4)
            cipr(IPR_OAERO,ipa_cel,kPAP1)=cipr(IPR_OAERO,ipa_cel,kPAP1)
     &                                  + con(kPAP1)-cold(kPAP1)
            cipr(IPR_OAERO,ipa_cel,kPAP2)=cipr(IPR_OAERO,ipa_cel,kPAP2)
     &                                  + con(kPAP2)-cold(kPAP2)
            cipr(IPR_OAERO,ipa_cel,kPAP3)=cipr(IPR_OAERO,ipa_cel,kPAP3)
     &                                  + con(kPAP3)-cold(kPAP3)
            cipr(IPR_OAERO,ipa_cel,kPAP4)=cipr(IPR_OAERO,ipa_cel,kPAP4)
     &                                  + con(kPAP4)-cold(kPAP4)
            cipr(IPR_OAERO,ipa_cel,kPCP1)=cipr(IPR_OAERO,ipa_cel,kPCP1)
     &                                  + con(kPCP1)-cold(kPCP1)
            cipr(IPR_OAERO,ipa_cel,kPCP2)=cipr(IPR_OAERO,ipa_cel,kPCP2)
     &                                  + con(kPCP2)-cold(kPCP2)
            cipr(IPR_OAERO,ipa_cel,kPCP3)=cipr(IPR_OAERO,ipa_cel,kPCP3)
     &                                  + con(kPCP3)-cold(kPCP3)
            cipr(IPR_OAERO,ipa_cel,kPCP4)=cipr(IPR_OAERO,ipa_cel,kPCP4)
     &                                  + con(kPCP4)-cold(kPCP4)
            cipr(IPR_OAERO,ipa_cel,kPFP1)=cipr(IPR_OAERO,ipa_cel,kPFP1)
     &                                  + con(kPFP1)-cold(kPFP1)
            cipr(IPR_OAERO,ipa_cel,kPFP2)=cipr(IPR_OAERO,ipa_cel,kPFP2)
     &                                  + con(kPFP2)-cold(kPFP2)
            cipr(IPR_OAERO,ipa_cel,kPFP3)=cipr(IPR_OAERO,ipa_cel,kPFP3)
     &                                  + con(kPFP3)-cold(kPFP3)
            cipr(IPR_OAERO,ipa_cel,kPFP4)=cipr(IPR_OAERO,ipa_cel,kPFP4)
     &                                  + con(kPFP4)-cold(kPFP4)
          endif
c
c========================== Process Analysis End =======================
c

        else ! SOAP

          cg(1)  = con(kcg1) * convfac * mwsoap(1)
          cg(2)  = con(kcg2) * convfac * mwsoap(2)
          cg(3)  = con(kcg3) * convfac * mwsoap(3)
          cg(4)  = con(kcg4) * convfac * mwsoap(4)
          soacon(1) = con(ksoa1)
          soacon(2) = con(ksoa2)
          soacon(3) = con(ksoa3)
          soacon(4) = con(ksoa4)
          preOM = con(kpoa) / mwpoa
     &          + con(ksopa) / mwsopa
     &          + con(ksopb) / mwsopb
c
c========================= DDM Begin ===================================
c
          if (lddm) then
            call DDMSOAINI (convfac,sddm,nfam,nsen)
          endif
c
c========================== DDM End ====================================
c
          do ispc = 1, NSOAP
            csatT(ispc) = csatS(ispc) * (cstempS(ispc)/tempk)
     &                                * exp((deltahS(ispc)/8.314)
     &                                * (1/cstempS(ispc)-1/tempk))
          enddo

          call soap(NSOAP,mwsoap,soacon(:NSOAP),cg(:NSOAP),
     &              csatT(:NSOAP),preOM,
     &              iout,igrdchm,ichm,jchm,kchm)
c
c-----Calculate SOA removal by particle-phase photolysis
c
          fphot = 1. - exp( -sfacJSOA * ajno2 ) ! fraction of removed SOA by photolysis
          do ispc = 1, NSOAP
            aphot(ispc) = soacon(ispc) * fphot
            soacon(ispc) = soacon(ispc) - aphot(ispc)
          enddo
          aphot(NSOAP+1) = con(ksopa) * fphot
          con(ksopa) = con(ksopa) - aphot(NSOAP+1)
          aphot(NSOAP+2) = con(ksopb) * fphot
          con(ksopb) = con(ksopb) - aphot(NSOAP+2)
c
c-----Calculate fraction of polymerized SOA
c
          fpoly = 1. - exp( -kpolya * dtaer ) ! polymerized fraction for anthro SOA
          do ispc = 1, NSOAA
            apoly(ispc) = soacon(ispc) * fpoly
            soacon(ispc) = soacon(ispc) - apoly(ispc)
            con(ksopa) = con(ksopa) + apoly(ispc)
          enddo
          fpoly = 1. - exp( -kpolyb * dtaer ) ! polymerized fraction for bio SOA
          do ispc = NSOAA + 1, NSOAP
            apoly(ispc) = soacon(ispc) * fpoly
            soacon(ispc) = soacon(ispc) - apoly(ispc)
            con(ksopb) = con(ksopb) + apoly(ispc)
          enddo
c
c========================= DDM Begin ===================================
c
          if (lddm) then
            call DDMSOAP (dtaer,ajno2,convfac,sddm,nfam,nsen,ierr,iout)
            if (ierr.lt.0) then
              write(iout,'(a,4i4)') ' igrd,i,j,k = ',
     &                                          igrdchm,ichm,jchm,kchm
              write(iout,*) ' T,DT = ',tempk,dtaer
              write(iout,*) ' CONVFAC = ',convfac
              write(iout,*) ' CG = ',con(kcg1),con(kcg2),
     &                               con(kcg3),con(kcg4)
              write(iout,*) ' SOA = ',con(ksoa1),con(ksoa2),
     &                                con(ksoa3),con(ksoa4),
     &                                con(ksopa),con(ksopb)
              write(iout,*) ' APHOT = ',(aphot(ispc),ispc=1,NSOAP+2)
              write(iout,*) ' APOLY = ',(apoly(ispc),ispc=1,NSOAP)
              write(iout,*) ' POA = ',con(kpoa)
              do ip = 1, nfam
                write(iout,*) ' SI(',ip,') = ',
     &                               (sddm(ip,jdpCG1+is-1),is=1,NSOAP),
     &                              (sddm(ip,jdpSOA1+is-1),is=1,NSOAP),
     &                               sddm(ip,jdpSOPA),sddm(ip,jdpSOPB),
     &                                                 sddm(ip,jdpPOA)
              enddo
              call camxerr()
            endif
          endif
c
c========================== DDM End ====================================
c
          con(kcg1)  = amax1(cg(1)/convfac/mwsoap(1),bdnl(kcg1))
          con(kcg2)  = amax1(cg(2)/convfac/mwsoap(2),bdnl(kcg2))
          con(kcg3)  = amax1(cg(3)/convfac/mwsoap(3),bdnl(kcg3))
          con(kcg4)  = amax1(cg(4)/convfac/mwsoap(4),bdnl(kcg4))
          con(ksoa1) = amax1(soacon(1),bdnl(ksoa1))
          con(ksoa2) = amax1(soacon(2),bdnl(ksoa2))
          con(ksoa3) = amax1(soacon(3),bdnl(ksoa3))
          con(ksoa4) = amax1(soacon(4),bdnl(ksoa4))
          con(ksopa) = amax1(con(ksopa),bdnl(ksopa))
          con(ksopb) = amax1(con(ksopb),bdnl(ksopb))
c
c========================= Process Analysis Begin ======================
c
          if ( ldoipr ) then
            cipr(IPR_OAERO,ipa_cel,kcg1) = cipr(IPR_OAERO,ipa_cel,kcg1)
     &                                 + (con(kcg1)-cold(kcg1))*convfac
            cipr(IPR_OAERO,ipa_cel,kcg2) = cipr(IPR_OAERO,ipa_cel,kcg2)
     &                                 + (con(kcg2)-cold(kcg2))*convfac
            cipr(IPR_OAERO,ipa_cel,kcg3) = cipr(IPR_OAERO,ipa_cel,kcg3)
     &                                 + (con(kcg3)-cold(kcg3))*convfac
            cipr(IPR_OAERO,ipa_cel,kcg4) = cipr(IPR_OAERO,ipa_cel,kcg4)
     &                                 + (con(kcg4)-cold(kcg4))*convfac
            cipr(IPR_OAERO,ipa_cel,ksoa1)=cipr(IPR_OAERO,ipa_cel,ksoa1)
     &                                 + con(ksoa1)-cold(ksoa1)
            cipr(IPR_OAERO,ipa_cel,ksoa2)=cipr(IPR_OAERO,ipa_cel,ksoa2)
     &                                 + con(ksoa2)-cold(ksoa2)
            cipr(IPR_OAERO,ipa_cel,ksoa3)=cipr(IPR_OAERO,ipa_cel,ksoa3)
     &                                 + con(ksoa3)-cold(ksoa3)
            cipr(IPR_OAERO,ipa_cel,ksoa4)=cipr(IPR_OAERO,ipa_cel,ksoa4)
     &                                 + con(ksoa4)-cold(ksoa4)
            cipr(IPR_OAERO,ipa_cel,ksopa)=cipr(IPR_OAERO,ipa_cel,ksopa)
     &                                 + con(ksopa)-cold(ksopa)
            cipr(IPR_OAERO,ipa_cel,ksopb)=cipr(IPR_OAERO,ipa_cel,ksopb)
     &                                 + con(ksopb)-cold(ksopb)
          endif
c
c========================== Process Analysis End =======================
c

        endif ! VBS or SOAP?

      endif ! lsvoc?
c
c-----Do RADM aqueous chemistry if CWC is above threshold
c     All conc units must be mol/mol (mixing ratio)
c
      if (lradm) then
        pres_pa = 100.*press
        dt_sec = dtaq*3600.
        cw_kgm3 = cwc/1000.
c
        r_gas(1)  = con(kso2)*1.e-6
        r_gas(2)  = con(khno3)*1.e-6
        r_gas(3)  = con(kn2o5)*1.e-6
        r_gas(4)  = co2*1.e-6
        r_gas(5)  = con(knh3)*1.e-6
        r_gas(6)  = con(khpo_c)*1.e-6
        r_gas(7)  = con(ko3)*1.e-6
        if (kfoa_c.lt.nspec+1) then
          r_gas(8) = con(kfoa_c)*1.e-6
        else
          r_gas(8) = foa*1.e-6
        endif
        if (kmhp_c.lt.nspec+1 .or. kohp_c.lt.nspec+1) then
          r_1 = con(kmhp_c) + con(kohp_c)
          r_gas(9) = r_1*1.e-6
          if( r_1 .NE. 0. ) r_1 = con(kmhp_c) / r_1
        else
          r_gas(9) = mhp*1.e-6
          r_1 = 0.0
        endif
        if (kpaa_c.lt.nspec+1 .or. kopa_c.lt.nspec+1) then
          r_2 = con(kpaa_c) + con(kopa_c)
          r_gas(10) = r_2*1.e-6
          if( r_2 .NE. 0. ) r_2 = con(kpaa_c) / r_2
        else
          r_gas(10) = paa*1.e-6
          r_2 = 0.0
        endif
        r_gas(11) = con(ksulf)*1.e-6
c
        if (kgly.lt.nspec+1) then
          r_gas(12) = con(kgly)*1.e-6
        else
          r_gas(12) = gly*1.e-6
        endif
        if (kmgly.lt.nspec+1) then
          r_gas(13) = con(kmgly)*1.e-6
        else
          r_gas(13) = mgly*1.e-6
        endif
        if (kglyd.lt.nspec+1) then
          r_gas(14) = con(kglyd)*1.e-6
        else
          r_gas(14) = glyd*1.e-6
        endif
        r_gas(15)  = con(koh)*1.e-6
c
        r_aer(1)  = (con(kpso4)/wtmol(1)/convfac)*1.e-6
        r_aer(2)  = (con(kpnh4)/wtmol(2)/convfac)*1.e-6
        r_aer(3)  = (con(kpno3)/wtmol(3)/convfac)*1.e-6

        if (f_CA_PRM.ge.0. .and. f_CA_CRS.ge.0.) then
          r_aer(4) = ((f_CA_PRM*con(kFPRM) + f_CA_CRS*con(kFCRS))/
     &                wtmol(8)/convfac)*1.e-6
        else
          r_aer(4) = (caco3/(wtmol(8)+wtmol(7))/convfac)*1.e-6
        endif
        if (kPCA.lt.nspec+1) then
          r_aer(4) = r_aer(4) + (con(kPCA)/wtmol(8)/convfac)*1.e-6
        endif

        if (f_MG_PRM.ge.0. .and. f_MG_CRS.ge.0.) then
          r_aer(5) = ((f_MG_PRM*con(kFPRM) + f_MG_CRS*con(kFCRS))/
     &                wtmol(9)/convfac)*1.e-6
        else
          r_aer(5) = (mgco3/(wtmol(9)+wtmol(7))/convfac)*1.e-6
        endif
        if (kPMG.lt.nspec+1) then
          r_aer(5) = r_aer(5) + (con(kPMG)/wtmol(9)/convfac)*1.e-6
        endif

        if (kNA.lt.nspec+1 .and. kPCL.lt.nspec+1) then
          if (con(kNA).gt.con(kPCL)) then
            r_aer(6) = (con(kPCL)/wtmol(4)/convfac)*1.e-6
          else
            r_aer(6) = (con(kNA)/wtmol(5)/convfac)*1.e-6
          endif
        else
          r_aer(6) = (nacl/(wtmol(5)+wtmol(4))/convfac)*1.e-6
        endif

        if ( ldark ) then
          fe_iii = 0.9  ! nighttime fraction partitioning of FE(III), GS 01july2011
        else
          fe_iii = 0.1  ! daytime fraction partitioning of FE(III), GS 01july2011
        endif
        fe_sol = 0.1    ! solubility of Fe, GS 01July2011
        if (f_FE_PRM.ge.0. .and. f_FE_CRS.ge.0.) then
          fe_mmr = ((f_FE_PRM*con(kFPRM) + f_FE_CRS*con(kFCRS))/
     &               wtmol(10)/convfac)*1.e-6
          r_aer(7) = fe_iii*fe_sol*fe_mmr ! GS 01July2011
        else
          r_aer(7) = (a3fe/wtmol(10)/convfac)*1.e-6
        endif
        if (kPFE.lt.nspec+1) then
          r_aer(7) = max(r_aer(7),
     &               fe_iii*fe_sol*(con(kPFE)/wtmol(10)/convfac)*1.e-6)
        endif

        mn_ii  = 1.0    ! fraction partitioning of MN(II), GS 01July2011
        mn_sol = 0.5    ! solubility of Mn, GS 28July2011
        if (f_MN_PRM.ge.0. .and. f_MN_CRS.ge.0.) then
          mn_mmr = ((f_MN_PRM*con(kFPRM) + f_MN_CRS*con(kFCRS))/ 
     &               wtmol(11)/convfac)*1.e-6
          r_aer(8) = mn_ii*mn_sol*mn_mmr ! GS 01July2011
        else
          r_aer(8) = (b2mn/wtmol(11)/convfac)*1.e-6
        endif
        if (kPMN.lt.nspec+1) then
          r_aer(8) = max(r_aer(8),
     &               mn_ii*mn_sol*(con(kPMN)/wtmol(11)/convfac)*1.e-6)
        endif

        if (f_K_PRM.ge.0. .and. f_K_CRS.ge.0.) then
          r_aer(9) = ((f_K_PRM*con(kFPRM) + f_K_CRS*con(kFCRS))/
     &                wtmol(6)/convfac)*1.e-6
        else
          r_aer(9) = (potcl/(wtmol(6)+wtmol(4))/convfac)*1.e-6
        endif
        if (kPK.lt.nspec+1) then
          r_aer(9) = r_aer(9) + (con(kPK)/wtmol(6)/convfac)*1.e-6
        endif
c
c-----Biogenic SOA tracer
c
        csoac = con(ksoac_c)
        r_aer(10) = (csoac/wtmol(12)/convfac)*1.e-6
c
c========================= DDM Begin ===================================
c
        if (lddm) then
          do is = 1, 15
            r_old(is) = r_gas(is)
          enddo
          do is = 1, 10
            r_old(is+15) = r_aer(is)
          enddo
          call DDMRADINI (con(kso2),con(kna),con(kpcl),
     &                    f_CA_PRM,f_CA_CRS,
     &                    convfac,wtmol,sddm,nfam,nsen,NWTMOL)
        endif
c
c========================== DDM End ====================================
c
        call raqchem(tempk,pres_pa,dt_sec,cw_kgm3,r_gas,r_aer,cph,
     &               idiag,iout,igrdchm,ichm,jchm,kchm,wtmol(12))
c
c========================= DDM Begin ===================================
c
        if (lddm) then
          call DDMRADM (r_1,r_gas(9)/r_old(9),r_2,r_gas(10)/r_old(10),
     &                  convfac,wtmol,sddm,nfam,nsen,NWTMOL,ierr,iout)
          if (ierr.lt.0) then
            write(iout,'(a,4i4)') ' igrd,i,j,k = ',
     &                                          igrdchm,ichm,jchm,kchm
            write(iout,*) ' T,P,CW,DT = ',tempk,pres_pa,cw_kgm3,dt_sec
            write(iout,*) ' CONVFAC = ',convfac
            write(iout,*) ' r1,r11,r2,r22 = ',r_1,r_gas(9)/r_old(9)
     &                                       ,r_2,r_gas(10)/r_old(10)
            write(iout,*) ' R_GAS = ',(r_old(is),is=1,15)
            write(iout,*) ' R_AER = ',(r_old(is+15),is=1,10)
            do ip = 1, nfam
              if (con(kna).gt.con(kpcl)) then
                r_3 = sddm(ip,jdpPCL )/wtmol(4)/convfac *1.e-6
              else
                r_3 = sddm(ip,jdpNA  )/wtmol(5)/convfac *1.e-6
              endif
              write(iout,*) ' SI(',ip,') = ',
     &      (sddm(ip,jdpNH3 )+sddm(ip,jdpPNH4)/wtmol(2)/convfac)*1.e-6,
     &      (sddm(ip,jdpHNO3)+sddm(ip,jdpN2O5)*2.0
     &                       +sddm(ip,jdpPNO3)/wtmol(3)/convfac)*1.e-6,
     &             sddm(ip,jdpFOA )*1.e-6, 0.0, sddm(ip,jdpSO2 )*1.e-6,
     &      (sddm(ip,jdpSULF)+sddm(ip,jdpPSO4)/wtmol(1)/convfac)*1.e-6,
     &                  sddm(ip,jdpH2O2)*1.e-6, sddm(ip,jdpO3  )*1.e-6,
     &                   ( sddm(ip,jdpMHP ) + sddm(ip,jdpOHP ) )*1.e-6,
     &                   ( sddm(ip,jdpPAA ) + sddm(ip,jdpOPA ) )*1.e-6,
     &          r_3,(f_CA_PRM*sddm(ip,jdpFPRM)
     &              +f_CA_CRS*sddm(ip,jdpFCRS))/wtmol(8)/convfac*1.e-6,
     &                    sddm(ip,jdpGLY)*1.e-6,sddm(ip,jdpMGLY)*1.e-6,
     &                     sddm(ip,jdpGLYD)*1.e-6,sddm(ip,jdpHO)*1.e-6,
     &                        sddm(ip,jdpSOPB)/wtmol(12)/convfac*1.e-6
            enddo
            call camxerr()
          endif
        endif
c
c========================== DDM End ====================================
c
        con(kso2)  = amax1(r_gas(1) *1.e6,bdnl(kso2))     ! SO2 (ppm)
        con(khno3) = amax1(r_gas(2) *1.e6,bdnl(khno3))    ! HNO3 (ppm)
        con(kn2o5) = amax1(r_gas(3) *1.e6,bdnl(kn2o5))    ! N2O5 gas (ppm)
        con(knh3)  = amax1(r_gas(5) *1.e6,bdnl(knh3))     ! NH3 (ppm)
        con(khpo_c)= amax1(r_gas(6) *1.e6,bdnl(khpo_c))   ! H2O2 (ppm)
        con(ko3)   = amax1(r_gas(7) *1.e6,bdnl(ko3))      ! O3 (ppm)
!bk        con(kfoa_c)= amax1(r_gas(8) *1.e6,bdnl(kfoa_c)) ! not changed by RADM
        con(kmhp_c)= amax1(r_gas(9) *1.e6*r_1,bdnl(kmhp_c))
        con(kohp_c)= amax1(r_gas(9) *1.e6*(1.-r_1),bdnl(kohp_c))
        con(kpaa_c)= amax1(r_gas(10)*1.e6*r_2,bdnl(kpaa_c))
        con(kopa_c)= amax1(r_gas(10)*1.e6*(1.-r_2),bdnl(kopa_c))
        con(ksulf) = amax1(r_gas(11)*1.e6,bdnl(ksulf))    ! H2SO4 (ppm)
c
        con(kgly)  = amax1(r_gas(12)*1.e6,bdnl(kgly))     ! GLY (ppm)
        con(kmgly) = amax1(r_gas(13)*1.e6,bdnl(kmgly))    ! MGLY (ppm)
        con(kglyd) = amax1(r_gas(14)*1.e6,bdnl(kglyd))    ! GLYD (ppm)
c
        con(kpso4) = amax1(r_aer(1)*convfac*wtmol(1)*1.e6,bdnl(kpso4)) ! PSO4 (ug/m3)
        con(kpnh4) = amax1(r_aer(2)*convfac*wtmol(2)*1.e6,bdnl(kpnh4)) ! PNH4 (ug/m3)
        con(kpno3) = amax1(r_aer(3)*convfac*wtmol(3)*1.e6,bdnl(kpno3)) ! PNO3 (ug/m3)
c
        con(ksoac_c) = amax1(r_aer(10)*convfac*wtmol(12)*1.e6,
     &                                                  bdnl(ksoac_c)) ! SOAC (ug/m3)
c
        con(nspec+1) = 0.
c
c========================= Process Analysis Begin ======================
c
        if ( ldoipr ) then
          cipr(IPR_AQCHEM,ipa_cel,kso2 )=cipr(IPR_AQCHEM,ipa_cel,kso2 )
     &                                +(con(kso2 )-cold(kso2 ))*convfac
          cipr(IPR_AQCHEM,ipa_cel,khno3)=cipr(IPR_AQCHEM,ipa_cel,khno3)
     &                                +(con(khno3)-cold(khno3))*convfac
          cipr(IPR_AQCHEM,ipa_cel,kn2o5)=cipr(IPR_AQCHEM,ipa_cel,kn2o5)
     &                                +(con(kn2o5)-cold(kn2o5))*convfac
          cipr(IPR_AQCHEM,ipa_cel,knh3 )=cipr(IPR_AQCHEM,ipa_cel,knh3 )
     &                                +(con(knh3 )-cold(knh3 ))*convfac
          cipr(IPR_AQCHEM,ipa_cel,khpo_c)
     &                                 =cipr(IPR_AQCHEM,ipa_cel,khpo_c)
     &                              +(con(khpo_c)-cold(khpo_c))*convfac
          cipr(IPR_AQCHEM,ipa_cel,ko3  )=cipr(IPR_AQCHEM,ipa_cel,ko3  )
     &                                +(con(ko3  )-cold(ko3  ))*convfac
!bk          if (kfoa_c.lt.nspec+1) cipr(IPR_AQCHEM,ipa_cel,kfoa_c) ! not changed by RADM
!bk     &                                 =cipr(IPR_AQCHEM,ipa_cel,kfoa_c)
!bk     &                              +(con(kfoa_c)-cold(kfoa_c))*convfac
          if (kmhp_c.lt.nspec+1) cipr(IPR_AQCHEM,ipa_cel,kmhp_c)
     &                                 =cipr(IPR_AQCHEM,ipa_cel,kmhp_c)
     &                              +(con(kmhp_c)-cold(kmhp_c))*convfac
          if (kohp_c.lt.nspec+1) cipr(IPR_AQCHEM,ipa_cel,kohp_c)
     &                                 =cipr(IPR_AQCHEM,ipa_cel,kohp_c)
     &                              +(con(kohp_c)-cold(kohp_c))*convfac
          if (kpaa_c.lt.nspec+1) cipr(IPR_AQCHEM,ipa_cel,kpaa_c)
     &                                 =cipr(IPR_AQCHEM,ipa_cel,kpaa_c)
     &                              +(con(kpaa_c)-cold(kpaa_c))*convfac
          if (kopa_c.lt.nspec+1) cipr(IPR_AQCHEM,ipa_cel,kopa_c)
     &                                 =cipr(IPR_AQCHEM,ipa_cel,kopa_c)
     &                              +(con(kopa_c)-cold(kopa_c))*convfac
          cipr(IPR_AQCHEM,ipa_cel,ksulf)=cipr(IPR_AQCHEM,ipa_cel,ksulf)
     &                                +(con(ksulf)-cold(ksulf))*convfac
          if (kgly.lt.nspec+1) cipr(IPR_AQCHEM,ipa_cel,kgly)
     &                                   =cipr(IPR_AQCHEM,ipa_cel,kgly)
     &                                  +(con(kgly)-cold(kgly))*convfac
          if (kmgly.lt.nspec+1) cipr(IPR_AQCHEM,ipa_cel,kmgly)
     &                                  =cipr(IPR_AQCHEM,ipa_cel,kmgly)
     &                                +(con(kmgly)-cold(kmgly))*convfac
          if (kglyd.lt.nspec+1) cipr(IPR_AQCHEM,ipa_cel,kglyd)
     &                                  =cipr(IPR_AQCHEM,ipa_cel,kglyd)
     &                                +(con(kglyd)-cold(kglyd))*convfac
          cipr(IPR_AQCHEM,ipa_cel,kpso4)=cipr(IPR_AQCHEM,ipa_cel,kpso4)
     &                                  +con(kpso4)-cold(kpso4)
          cipr(IPR_AQCHEM,ipa_cel,kpnh4)=cipr(IPR_AQCHEM,ipa_cel,kpnh4)
     &                                  +con(kpnh4)-cold(kpnh4)
          cipr(IPR_AQCHEM,ipa_cel,kpno3)=cipr(IPR_AQCHEM,ipa_cel,kpno3)
     &                                  +con(kpno3)-cold(kpno3)
          cipr(IPR_AQCHEM,ipa_cel,ksoac_c)
     &                                =cipr(IPR_AQCHEM,ipa_cel,ksoac_c)
     &                                  +con(ksoac_c)-csoac
        endif
c
c========================== Process Analysis End =======================
c
      endif ! lradm?
c
c======================== Source Apportion Begin =======================
c
      if( ldosa ) then
         do ispc=1,ngas
           delcon(2,ispc) = delcon(2,ispc) +
     &                     AMAX1(bdnl(ispc),con(ispc))*convfac
         enddo
         do ispc=ngas+1,nspec
           delcon(2,ispc) = delcon(2,ispc) +
     &                               AMAX1(bdnl(ispc),con(ispc))
         enddo
         deln2o5 = delcon(2,kn2o5) - delcon(1,kn2o5)
         delcon(4,khno3) = delcon(4,khno3) - deln2o5 * 2.0
         delcon(4,kn2o5) = delcon(4,kn2o5) + deln2o5
         if ( lsvoc .and. .not.lvbs ) then
           delcon(5,ksoa1) = -apoly(1)
           delcon(5,ksoa2) = -apoly(2)
           delcon(5,ksoa3) = -apoly(3)
           delcon(5,ksoa4) = -apoly(4)
           do ispc = 1, NSOAA
             delcon(5,ksopa) = delcon(5,ksopa) + apoly(ispc)
           enddo
           do ispc = NSOAA+1, NSOAP
             delcon(5,ksopb) = delcon(5,ksopb) + apoly(ispc)
           enddo
           delcon(6,ksoa1) = -aphot(1)
           delcon(6,ksoa2) = -aphot(2)
           delcon(6,ksoa3) = -aphot(3)
           delcon(6,ksoa4) = -aphot(4)
           delcon(6,ksopa) = -aphot(5)
           delcon(6,ksopb) = -aphot(6)
         endif
         if ( lradm ) then
           delcon(7,ksoac_c) = con(ksoac_c) - csoac
         endif
      endif
c
c======================== Source Apportion End =======================
c
c
c-----Inorganic aerosol equilibrium chemistry with ISOROPIA
c     Convert conc units to mol/m3 (double precision)
c
      if (lisor) then
c
c========================= Process Analysis Begin ======================
c
        if ( ldoipr ) then
          csulf = con(ksulf)
          chno3 = con(khno3)
          cnh3  = con(knh3 )
          cpso4 = con(kpso4)
          cpno3 = con(kpno3)
          cpnh4 = con(kpnh4)
        endif
c
c========================== Process Analysis End =======================
c
        if ( leqsam ) then ! EQSAM

        xTT = DBLE(tempk)
        xAW = DBLE(rh/100.)

        xYPa(1) = DBLE((con(knh3) *convfac + con(kpnh4)/wtmol(2))*1.e-6) ! NH4+  [mol/m3]
        if (kNA.lt.nspec+1 .and. kPCL.lt.nspec+1) then
          xYPa(2) = DBLE(                     con(kna )/wtmol(5) *1.e-6) ! Na+   [mol/m3]
          xYMa(4) = DBLE((con(khcl)*convfac + con(kpcl)/wtmol(4))*1.e-6) ! Cl-   [mol/m3]
        else
          xYPa(2) = DBLE(                nacl/(wtmol(5)+wtmol(4))*1.e-6)
          xYMa(4) = DBLE(con(khcl)*convfac*1.e-6) + xYPa(2)
        endif
        if (f_K_PRM.ge.0. .and. f_K_CRS.ge.0.) then
          xYPa(3) = DBLE((f_K_PRM*con(kFPRM) + f_K_CRS*con(kFCRS))
     &                                                 /wtmol(6) *1.e-6) ! K+    [mol/m3]
        else
          xYPa(3) = DBLE(               potcl/(wtmol(6)+wtmol(4))*1.e-6)
          xYMa(4) = xYMa(4) + xYPa(3)
        endif
        if (kPK.lt.nspec+1) then
          xYPa(3) = xYPa(3) + DBLE(           con(kpk )/wtmol(6) *1.e-6)
        endif
        if (f_CA_PRM.ge.0. .and. f_CA_CRS.ge.0.) then
          xYPa(4) = DBLE((f_CA_PRM*con(kFPRM) + f_CA_CRS*con(kFCRS))
     &                                                 /wtmol(8) *1.e-6) ! Ca++  [mol/m3]
        else
          xYPa(4) = DBLE(               caco3/(wtmol(8)+wtmol(7))*1.e-6)
        endif
        if (kPCA.lt.nspec+1) then
          xYPa(4) = xYPa(4) + DBLE(           con(kpca)/wtmol(8) *1.e-6)
        endif
        if (f_MG_PRM.ge.0. .and. f_MG_CRS.ge.0.) then
          xYPa(5) = DBLE((f_MG_PRM*con(kFPRM) + f_MG_CRS*con(kFCRS))
     &                                                 /wtmol(9) *1.e-6) ! Mg++  [mol/m3]
        else
          xYPa(5) = DBLE(               mgco3/(wtmol(9)+wtmol(7))*1.e-6)
        endif
        if (kPMG.lt.nspec+1) then
          xYPa(5) = xYPa(5) + DBLE(           con(kpmg)/wtmol(9) *1.e-6)
        endif
        xYMa(1) = DBLE((con(ksulf)*convfac + con(kpso4)/wtmol(1))*1.e-6) ! SO4-- [mol/m3]
        xYMa(2) = 0.0_dp                                                 ! HSO4- [mol/m3]
        xYMa(3) = DBLE((con(khno3)*convfac + con(kpno3)/wtmol(3))*1.e-6) ! NO3-  [mol/m3]

        xWH2O = 0.0_dp ! set to zero to make the EQSAM metastable flag effective

        wt_save(2) = xYMa(1)
        wt_save(3) = xYPa(1)
        wt_save(4) = xYMa(3)
        wt_save(5) = xYMa(4)

        call EQSAM4clim(xTT,xAW,xsPM,xaPM,xPMS,xPMT,xRHO,xVOL,xPH,xHp,
     &                  xGF,xWH2O,xYPa,xYPs,xYMa,xYMs,xYG,
     &                  imask,1,1,jmeq,Dd)

        con(ksulf) = bdnl(ksulf)
        con(kpso4) = wt_save(2) - DBLE(con(ksulf)*1.e-6*convfac)
        con(kpso4) = amax1(con(kpso4)*wtmol(1)*1.e6,bdnl(kpso4))

        con(khno3) = REAL(xYG(1))*1.e6/convfac
        con(khno3) = amax1(con(khno3),bdnl(khno3))
        con(kpno3) = wt_save(4) - DBLE(con(khno3)*1.e-6*convfac)
        con(kpno3) = amax1(con(kpno3)*wtmol(3)*1.e6,bdnl(kpno3))

        con(knh3 ) = REAL(xYG(3))*1.e6/convfac
        con(knh3 ) = amax1(con(knh3),bdnl(knh3))
        con(kpnh4) = wt_save(3) - DBLE(con(knh3)*1.e-6*convfac)
        con(kpnh4) = amax1(con(kpnh4)*wtmol(2)*1.e6,bdnl(kpnh4))

        if (khcl.lt.nspec+1) then
          con(khcl ) = REAL(xYG(2))*1.e6/convfac
          con(khcl ) = amax1(con(khcl),bdnl(khcl))
        endif
        if (kpcl.lt.nspec+1) then
          con(kpcl ) = wt_save(5) - DBLE(con(khcl)*1.e-6*convfac)
          con(kpcl ) = amax1(con(kpcl)*wtmol(4)*1.e6,bdnl(kpcl))
        endif

        con(kph2o) = amax1(REAL(xWH2O),bdnl(kph2o))

        aph = REAL(xPH)
        aph = MIN(MAX(aph,-2.0),10.0)

        else ! ISORROPIA

        rhi = DBLE(rh/100.)
        tempi = DBLE(tempk)
        cntrl(1) = 0.d0                   ! 0 = forward problem
        cntrl(2) = 1.d0                   ! 1 = metastable liquid state

        if (kNA.lt.nspec+1 .and. kPCL.lt.nspec+1) then
          wi(1) = DBLE(con(kna)/wtmol(5)*1.e-6)                        ! total sodium
          wi(5) = DBLE((con(khcl)*convfac + con(kpcl)/wtmol(4))*1.e-6) ! total chloride
        else
          wi(1) = DBLE(nacl/(wtmol(5)+wtmol(4))*1.e-6)                 ! total sodium
          wi(5) = DBLE(con(khcl)*convfac*1.e-6) + wi(1)                ! total chloride
        endif
        wi(2) = DBLE((con(ksulf)*convfac + con(kpso4)/wtmol(1))*1.e-6) ! total sulfate
        wi(3) = DBLE((con(knh3) *convfac + con(kpnh4)/wtmol(2))*1.e-6) ! total ammonium
        wi(4) = DBLE((con(khno3)*convfac + con(kpno3)/wtmol(3))*1.e-6) ! total nitrate

        if (f_CA_PRM.ge.0. .and. f_CA_CRS.ge.0.) then
          dno3 = DBLE((f_CA_PRM*con(kFPRM) + f_CA_CRS*con(kFCRS))
     &                                               /wtmol(8) *1.e-6)
        else
          dno3 = DBLE(caco3/(wtmol(8)+wtmol(7))*1.e-6)
        endif
        if (kPCA.lt.nspec+1) then
          dno3 = dno3 + DBLE(con(kPCA)/wtmol(8)*1.e-6) ! half of CaCO3 is assumed to be replaced by Ca(NO3)2          
        endif
        if (dno3.lt.wi(4)) then
          wi(4) = wi(4) - dno3 ! exclude calcium nitrate
          lno3_isr = .true.
        else ! all nitrate is calcium nitrate
          con(khno3) = bdnl(khno3)
          con(kpno3) = amax1( (REAL(wi(4))*1.e6 - con(khno3)*convfac)
     &                                  * wtmol(3), bdnl(kpno3) )
          wi(4) = 0.d0
          lno3_isr = .false.
        endif
c
        do is = 1, 5
          wt_save(is) = wi(is)
        enddo
c
c========================= DDM Begin ===================================
c
        if (lddm) then
          call DDMISRINI (lno3_isr,f_CA_PRM,f_CA_CRS,
     &                    convfac,wtmol,sddm,nfam,nsen,NWTMOL,iout)
        endif
c
c========================== DDM End ====================================
c
        call isoropia(wi,rhi,tempi,cntrl,wt,gasis,aerliq,aersld,
     &                scasi,other)
c
c========================= DDM Begin ===================================
c
        if (lddm) then
          call DDMISRPIA (lno3_isr,f_CA_PRM,f_CA_CRS,
     &                    convfac,wtmol,sddm,nfam,nsen,NWTMOL,ierr,iout)
          if (ierr.lt.0) then
            write(iout,'(a,4i4)') ' igrd,i,j,k = ',
     &                                          igrdchm,ichm,jchm,kchm
            write(iout,*) ' RH,T = ',rhi,tempi
            write(iout,*) ' CONVFAC = ',convfac
            write(iout,*) ' WI = ',(wt_save(is),is=1,5)
            write(iout,*) ' LNO3_ISR = ',lno3_isr,f_CA_PRM,f_CA_CRS,
     &                     (sddm(ip,jdpFPRM),sddm(ip,jdpFCRS),ip=1,nfam)
            do ip = 1, nfam
              if (lno3_isr) then
                r_1 = (          sddm(ip,jdpHNO3)*convfac
     &                +          sddm(ip,jdpPNO3)/wtmol(3)
     &                -(f_CA_PRM*sddm(ip,jdpFPRM)
     &                 +f_CA_CRS*sddm(ip,jdpFCRS))/wtmol(8))*1.e-6
              else
                r_1 = 0.0
              endif
              write(iout,*) ' SI(',ip,') = ',
     &                                sddm(ip,jdpNA  )/wtmol(5) *1.e-6,
     &      (sddm(ip,jdpSULF)*convfac+sddm(ip,jdpPSO4)/wtmol(1))*1.e-6,
     &      (sddm(ip,jdpNH3 )*convfac+sddm(ip,jdpPNH4)/wtmol(2))*1.e-6,
     &                                                             r_1,
     &      (sddm(ip,jdpHCL )*convfac+sddm(ip,jdpPCL )/wtmol(4))*1.e-6
            enddo
            call camxerr()
          endif
        endif
c
c========================== DDM End ====================================
c 
c-----Load results back to local CON array (ppm for gas, ug/m3 for aerosol)
c 
        con(ksulf) = bdnl(ksulf)                                  ! sulfuric acid gas
        con(kpso4) = wt_save(2) - DBLE(con(ksulf)*1.e-6*convfac)  ! sulfate
        con(kpso4) = amax1(con(kpso4)*wtmol(1)*1.e6,bdnl(kpso4))
c
        if (lno3_isr) then
          con(khno3) = REAL(gasis(2))*1.e6/convfac                ! nitric acid gas
          con(khno3) = amax1(con(khno3),bdnl(khno3))    
          con(kpno3) = wt_save(4) - DBLE(con(khno3)*1.e-6*convfac)! nitrate
          con(kpno3) = con(kpno3) + dno3                          ! add calcium nitrate back
          con(kpno3) = amax1(con(kpno3)*wtmol(3)*1.e6,bdnl(kpno3))
        endif
c
        con(knh3 ) = REAL(gasis(1))*1.e6/convfac                  ! ammonia gas
        con(knh3 ) = amax1(con(knh3),bdnl(knh3))    
        con(kpnh4) = wt_save(3) - DBLE(con(knh3)*1.e-6*convfac)   ! ammonium
        con(kpnh4) = amax1(con(kpnh4)*wtmol(2)*1.e6,bdnl(kpnh4))
c
        if (khcl.lt.nspec+1) then
          con(khcl) = REAL(gasis(3))*1.e6/convfac                 ! HCl gas
          con(khcl) = amax1(con(khcl),bdnl(khcl))
        endif
        if (kpcl.lt.nspec+1) then
          con(kpcl) = wt_save(5) - DBLE(con(khcl)*1.e-6*convfac)  ! chloride
          con(kpcl) = amax1(con(kpcl)*wtmol(4)*1.e6,bdnl(kpcl))
        endif
c
        con(kph2o) = REAL(aerliq(8))                              ! aerosol water
        con(kph2o) = amax1(con(kph2o)*18.*1.e6,bdnl(kph2o))
c
        aph = LOG10(REAL(aerliq(8)*18.D-3/MAX(aerliq(1),1.D-20))) ! aerosol pH
        aph = MIN(MAX(aph,-2.0),10.0)

        endif ! leqsam?
c
c========================= Process Analysis Begin ======================
c
        if ( ldoipr ) then
          cipr(IPR_IAERO,ipa_cel,ksulf) = cipr(IPR_IAERO,ipa_cel,ksulf)
     &                                  + (con(ksulf)-csulf)*convfac
          cipr(IPR_IAERO,ipa_cel,khno3) = cipr(IPR_IAERO,ipa_cel,khno3)
     &                                  + (con(khno3)-chno3)*convfac
          cipr(IPR_IAERO,ipa_cel,knh3 ) = cipr(IPR_IAERO,ipa_cel,knh3 )
     &                                  + (con(knh3 )-cnh3 )*convfac
          cipr(IPR_IAERO,ipa_cel,kpso4) = cipr(IPR_IAERO,ipa_cel,kpso4)
     &                                  + con(kpso4)-cpso4
          cipr(IPR_IAERO,ipa_cel,kpno3) = cipr(IPR_IAERO,ipa_cel,kpno3)
     &                                  + con(kpno3)-cpno3
          cipr(IPR_IAERO,ipa_cel,kpnh4) = cipr(IPR_IAERO,ipa_cel,kpnh4)
     &                                  + con(kpnh4)-cpnh4
          if (khcl.lt.nspec+1)
     &      cipr(IPR_IAERO,ipa_cel,khcl) = cipr(IPR_IAERO,ipa_cel,khcl)
     &                                 + (con(khcl)-cold(khcl))*convfac
          if (kpcl.lt.nspec+1)
     &      cipr(IPR_IAERO,ipa_cel,kpcl) = cipr(IPR_IAERO,ipa_cel,kpcl)
     &                                   + con(kpcl)-cold(kpcl)
          cipr(IPR_IAERO,ipa_cel,kph2o) = cipr(IPR_IAERO,ipa_cel,kph2o)
     &                                  + con(kph2o)-cold(kph2o)
        endif
c
c========================== Process Analysis End =======================
c
      endif ! lisor?
c
      return
c
      end
