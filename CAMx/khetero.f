      subroutine khetero(tcell,pcell,H2O,con,lddm,nfam,sddm)
      use chmstry
      use filunit
      implicit none
c
c----CAMx v7.20 220430
c
c     KHETERO estimates the rate constants for heterogeneous reactions of N2O5:
c       N2O5 + H2O -> 2 HNO3
c       N2O5 + HCL -> CLN2 + HNO3
c     KHETERO also sets the rate constant for hydrolysis of organic nitrate:
c       NTR2       -> HNO3
c       INTR       -> HNO3
c     KHETERO also sets the rate constant for heterogeneous reaction of SO2:
c       SO2        -> SULF
c
c     Parameterization of N2O5 on aqueous particles:
c     Bertram and Thornton (2009) Atmos. Chem. Phys., 9, 8351-8363.
c
c     Mean molecular speed formula - Seinfeld and Pandis (1998)
c     Atmospheric Chemistry and Physics, p453.
c
c     NTR2 hydrolysis rate:
c     Liu et al. (2012) Aerosol Sci. Technol., 46, 1359-1369.
c
c     INTR hydrolysis rate:
c     Fisher et al. (2016) Atmos. Chem. Phys., 16, 5969-5991.
c
c     Organic nitrate partitioning ratio:
c     Rollins et al. (2013) J. Geophys. Res., 118, 6651-6662.
c
c     Heterogeneous SO2 reaction rate:
c     Zheng et al. (2015) Atmos. Chem. Phys., 15, 2031-2049.
c
c     Molecular diffusion coefficient:
c     Tang et al. (2014) Atmos. Chem. Phys., 14, 9233-9247.
c     (https://sites.google.com/site/mingjintang/home/diffusion)
c
c      Copyright 1996 - 2022
c     Ramboll
c
c     Modifications:
c        01/30/2006  - bkoo -     Original development
c        07/04/2007  - bkoo -     Modified to use pointer to N2O5 hydrolysis
c                                 Extended for CMU_AERO
c        03/29/2011  -cemery-     Revised to allow for inert PM with gas-phase
c                                 chemistry
c        03/28/2012  -cemery-     Revised to set rate of 0.25/hr in the
c                                 absence of PM chemistry
c        09/21/2013  - bkoo -     Revised for DDM-PM
c        10/14/2013  - bkoo -     Added hydrolysis of organic nitrate
c        01/03/2014  - bkoo -     Updated DDM-PM for organic nitrate hydrolysis
c        03/18/2014  - bkoo -     Updated for benzene SOA
c        12/07/2014  - bkoo -     Modified for VBS
c                                 Cap khetN2O5 sens
c        01/08/2016  - bkoo -     Updated for Revised VBS
c        02/29/2016  - bkoo -     Updated heterogeneous rxn of N2O5
c        07/20/2016  - bkoo -     Added heterogeneous hydrolysis of INTR
c        08/25/2016  - bkoo -     Updated for new SOAP
c        09/02/2016  - bkoo -     Added het SO2 -> SULF; revised het rxn pointers
c        11/28/2016  - bkoo -     Revised N2O5 uptake rxn formula
c
c     Input arguments:
c        tcell          temperature (K)
c        pcell          pressure (mb)
c        H2O            water vapor concentration (ppm)
c        con            species concentrations (ppm or ug/m3)
c        lddm           logical flag for DDM sensitivities
c        nfam           number of sensitivity families
c        sddm           sensitivity matrix (ppm; ug/m3; 1/hr or 1/ppm-hr)
c                          sddm(:,nspec+1:nspec+2) [1/ppm-hr]
c                          sddm(:,nspec+3:nspec+5) [1/hr]
c
c     Output arguments:
c        sddm
c
c     Routines called:
c        NONE
c
c     Called by:
c        CHEMDRIV
c
      include 'camx.prm'
      include 'camx_aero.inc'
      include 'vbs.inc'
c
c-----Arguments
c
      real  tcell, pcell, H2O
      real  con(*)
c
c-----Local variables:
c
      real, parameter :: khon0  = 0.16667       ! [1/hr]; lifetime = 6 hr (particle-phase NTR2 hydrolysis)
      real, parameter :: khon20 = 1.0           ! [1/hr]; lifetime = 1 hr (particle-phase INTR hydrolysis)
      real, parameter :: alfa1  = 0.34,         ! partioning parameters for total ON
     &                   alfa2  = 0.66,
     &                   cstar1 = 0.73,
     &                   cstar2 = 1000.
      real, parameter :: fntr2  = 0.71          ! NTR2 / total ON
      real, parameter :: fintr  = 0.08          ! INTR / total ON

      real :: fcond                             ! condensable ON / total ON
      real :: khon                              ! psuedo gas-phase NTR2 hydrolysis rate constant [1/hr]
      real :: khon2                             ! psuedo gas-phase INTR hydrolysis rate constant [1/hr]
      real :: fp                                ! fraction of NTR2 or INTR in the particle phase
                                                ! (assuming the same fp for NTR2 and INTR)
      real :: ctoa                              ! total OA concentration [ug/m3]

      real :: qwatr, ev, es, rh                 ! relative humidity calculation

      real, parameter :: Rgc    = 8.314         ! gas constant [J/mol-K]
      real, parameter :: Pi     = 3.1415927     ! Pi
      real, parameter :: mwN2O5 = 0.108         ! MW of N2O5 [kg/mol]
      real, parameter :: mwSO2  = 0.064         ! MW of SO2 [kg/mol]
      real, parameter :: rtau   = 1.66667e-3    ! reciprocal of minimum lifetime of N2O5 (10 min = 600 sec) [1/s]

      real, parameter :: Afac   = 3.2e-8        ! empirical pre-factor for gamma [s]
      real, parameter :: beta   = 1.15e6        ! parameter for k2f' [1/s]
      real, parameter :: delta  = 1.3e-1        ! parameter for k2f' [1/M]
      real, parameter :: k3k2b  = 6.0e-2        ! k3/k2b
      real, parameter :: k4k2b  = 29.0          ! k4/k2b
      real, parameter :: molarPH2O = 55.        ! molar concentration of aerosol water [M]; assume diluted aqueous phase

      real, parameter :: diff_N2O5 = 65.        ! Diffusivity of N2O5 at 296 K [Torr-cm2/s]
      real, parameter :: diff_SO2 = 94.         ! Diffusivity of SO2 at 296 K [Torr-cm2/s]
      real, parameter :: g1_SO2 = 2.e-5         ! lower limit of SO2 uptake coefficient
      real, parameter :: g2_SO2 = 5.e-5         ! upper limit of SO2 uptake coefficient

      real, parameter :: h2s = 3600.0           ! hour to second

      real :: vt                                ! total dry volume of fine particles
      real :: diam                              ! diameter of a fine particle [um]
      real :: nden                              ! number density of fine particles [#/cm3_air]
      real :: saden                             ! surface area density [um2/cm3_air]
      real :: cbarN2O5                          ! mean molecular speed of N2O5 [m/s]
      real :: cbarSO2                           ! mean molecular speed of SO2 [m/s]
      real :: k2fprime                          ! k2f' = k2f * [H2O] [1/s]
      real :: dc_N2O5                           ! diffusion coefficient of N2O5 [m2/s]
      real :: dc_SO2                            ! diffusion coefficient of SO2 [m2/s]
      real :: gamma_N2O5                        ! uptake coefficient of N2O5
      real :: gamma_SO2                         ! uptake coefficient of SO2
      real :: y_clno2                           ! product yield for ClNO2

      real*8 :: r1, r2, rdenom
      real :: rscl
      real :: khet0, khet1, khet2, khet5

      integer :: l, isempty
      integer :: kwtr, isec, iaero
c
c======================== DDM Begin =======================
c
      logical :: lddm
      integer :: nfam
      integer :: ip
      real :: sddm(nfam,MXSPEC+NHETRXN)
      real :: stoa,sfp
      real :: const0,sad0
      real :: dvt,dsad,dr1r
c
c======================== DDM End =======================
c
c-----No PM chemistry:
c       N2O5 + H2O -> 2 HNO3      : reset to the IUPAC 2006 rate
c       N2O5 + HCL -> CLN2 + HNO3 : use the rate set in the CHEMPARAM input file
c       NTR2       -> HNO3        : use the rate set in the CHEMPARAM input file
c       INTR       -> HNO3        : use the rate set in the CHEMPARAM input file
c       SO2        -> SULF        : use the rate set in the CHEMPARAM input file
c
      if (aeropt.eq.'NONE' .or. aeropt.eq.'INERT') then
        rk(ihetrxn(1)) = 2.5e-22 * 8.87e16    ! [cm3/molecule-sec] -> [1/ppm-hr]
        return
      endif
c
c-----Calculate RH
c
      qwatr = 1.e-6*H2O*18./28.8
      ev = qwatr*pcell/(qwatr + eps)
      es = e0*exp((lv/rv)*(1./273. - 1./tcell))
      rh = amin1(1.,ev/es)
c
c-----Hydrolysis of organic nitrate in the organic solution phase
c
      if (ihetrxn(3).gt.0) then
        if (aeropt.eq.'CF') then
          ctoa = con(kpoa)  + con(ksoa1) + con(ksoa2)
     &                      + con(ksoa3) + con(ksoa4)
     &                      + con(ksopa) + con(ksopb)
          if (lvbs) then
            ctoa = con(kpas0)+con(kpas1)+con(kpas2)+con(kpas3)+con(kpas4)
     &            +con(kpbs0)+con(kpbs1)+con(kpbs2)+con(kpbs3)+con(kpbs4)
     &            +con(kpap0)+con(kpap1)+con(kpap2)+con(kpap3)+con(kpap4)
     &            +con(kpcp0)+con(kpcp1)+con(kpcp2)+con(kpcp3)+con(kpcp4)
     &            +con(kpfp0)+con(kpfp1)+con(kpfp2)+con(kpfp3)+con(kpfp4)
          endif
        elseif (aeropt.eq.'CMU') then
          ctoa = 0.0
          do isec = 1, nbin
            ctoa = ctoa + con(kpoa_1+isec-1)
     &                  + con(ksoa1_1+isec-1) + con(ksoa2_1+isec-1)
     &                  + con(ksoa3_1+isec-1) + con(ksoa4_1+isec-1)
     &                  + con(ksopa_1+isec-1) + con(ksopb_1+isec-1)
          enddo
        else
          write(iout,'(//,a)') 'ERROR in KHETERO: invalid AEROPT'
          call camxerr()
        endif
        fp = 0.0
        if (ctoa.ge.0.0001) then
          fp = alfa1 / (1. + cstar1/ctoa) + alfa2 / (1. + cstar2/ctoa)
          fcond = fntr2
          if (ihetrxn(4).gt.0) fcond = fcond + fintr
          fp = fp / fcond
        endif
        if (rh.lt.0.2) then
          rk(ihetrxn(3)) = 0.0
          if (ihetrxn(4).gt.0) rk(ihetrxn(4)) = 0.0
        elseif (rh.ge.0.4) then
          rk(ihetrxn(3)) = fp * khon0
          if (ihetrxn(4).gt.0) rk(ihetrxn(4)) = fp * khon20
        else
          rk(ihetrxn(3)) = fp * khon0 * (rh - 0.2) / 0.2
          if (ihetrxn(4).gt.0) rk(ihetrxn(4)) = fp * khon20
     &                                * (rh - 0.2) / 0.2
        endif
c
c======================== DDM Begin =======================
c
        if ( lddm .and. aeropt.eq.'CF' ) then
          do ip = 1, nfam
            stoa = sddm(ip,kpoa)
     &           + sddm(ip,ksoa1) + sddm(ip,ksoa2)
     &           + sddm(ip,ksoa3) + sddm(ip,ksoa4)
     &           + sddm(ip,ksopa) + sddm(ip,ksopb)
            sfp = 0.0
            if (ctoa.ge.0.0001) then
              sfp = ( fp - ( 
     &                alfa1 / ((1. + cstar1/ctoa)*(1. + cstar1/ctoa))
     &              + alfa2 / ((1. + cstar2/ctoa)*(1. + cstar2/ctoa))
     &                      ) / fcond ) * stoa / ctoa
            endif
            if (rh.lt.0.2) then
              sddm(ip,nspec+3) = 0.0
              if (ihetrxn(4).gt.0) sddm(ip,nspec+4) = 0.0
            elseif (rh.ge.0.4) then
              sddm(ip,nspec+3) = sfp * khon0
              if (ihetrxn(4).gt.0) sddm(ip,nspec+4) = sfp * khon20
            else
              sddm(ip,nspec+3) = sfp * khon0 * (rh - 0.2) / 0.2
              if (ihetrxn(4).gt.0) sddm(ip,nspec+4) = sfp * khon20
     &                                       * (rh - 0.2) / 0.2
            endif
          enddo
        endif
c
c======================== DDM End =======================
c
      endif
c
c-----Heterogeneous rxn of N2O5 / SO2 on aerosol surface area
c
      khet1 = 0.0
      khet2 = 0.0
      khet5 = 0.0

      k2fprime = beta * ( 1. - exp(-delta*molarPH2O) )      ! [1/s]
      cbarN2O5 = sqrt( 8.0*Rgc*tcell / (Pi*mwN2O5) )        ! [m/s]
      dc_N2O5 = (diff_N2O5/7.6e6) * (1013.25/pcell) * (tcell/296.)**1.75 ! [m2/s]

      cbarSO2 = sqrt( 8.0*Rgc*tcell / (Pi*mwSO2) )          ! [m/s]
      dc_SO2 = (diff_SO2/7.6e6) * (1013.25/pcell) * (tcell/296.)**1.75 ! [m2/s]
      if (rh.lt.0.5) then
        gamma_SO2 = g1_SO2
      else
        gamma_SO2 = g1_SO2 + (g2_SO2 - g1_SO2)*(rh - 0.5)/0.5
      endif
c
c-----CF PM chemistry
c
      if (aeropt.eq.'CF') then
        diam = sqrt(dcut(kph2o,1)*dcut(kph2o,2))            ! geometric mean dry diameter [um]
        !bk - limit hetero SO2 rxn only at dust particle surface
        saden = (6./diam) * ( con(kfprm)/roprt(kfprm)       ! dust particle surface area [m2/m3_air]
     &                      + con(kfcrs)/roprt(kfcrs) 
     &                      + con(kpfe) /roprt(kfcrs) 
     &                      + con(kpmn) /roprt(kfcrs) 
     &                      + con(kpca) /roprt(kfcrs) 
     &                      + con(kpmg) /roprt(kfcrs) 
     &                      + con(kpk)  /roprt(kfcrs) 
     &                      + con(kpal) /roprt(kfcrs) 
     &                      + con(kpsi) /roprt(kfcrs) 
     &                      + con(kpti) /roprt(kfcrs) )
        khet5 = saden / ( diam * 0.5e-6 / dc_SO2 +
     &                    4. / ( cbarSO2 * gamma_SO2 ) )    ! [1/s]
c
c-----Adjust dry diameter to wet diameter for N2O5 hydrolysis
c
        if (con(kph2o).le.bdnl(kph2o)) goto 800

        isempty = 1
        vt = 0.
        do l = ngas+1,nspec
          if (dcut(l,2).lt.(dcut(kph2o,2)+1.e-5)) then      ! fine
            if (l.ne.kph2o) then                            ! dry
              if (con(l).gt.bdnl(l)) isempty = 0
              vt = vt + con(l)/roprt(l)                     ! [ug/m3_air]/[g/m3]
            endif
          endif
        enddo
        if (isempty.eq.1) goto 800

        nden = vt / (diam*diam*diam*Pi/6.0)                 ! [ug/m3_air]/[um3-g/m3] = [#/m3_air] * 1.E-12
c
c======================== DDM Begin =======================
c
        if ( lddm ) then !bk - broken
        endif
c
c======================== DDM End =======================
c
        diam = diam * ( 1.0 +
     &         con(kph2o)/roprt(kph2o)/vt )**0.33333        ! wet diameter [um]
        saden = Pi*diam*diam*nden                           ! [m2/m3_air]
        r1 = k3k2b * (con(kph2o)/18.)/(con(kpno3)/62.)
        r2 = k4k2b * (con(kpcl)/35.)/(con(kpno3)/62.)
        rdenom = r1 + 1. + r2
        gamma_N2O5 = Afac * k2fprime * (1. - 1. / rdenom)
        y_clno2 = r2 / (r1 + r2)
        khet0 = saden / ( diam * 0.5e-6 / dc_N2O5 +
     &                    4. / ( cbarN2O5 * gamma_N2O5 ) )  ! [1/s]
        khet0 = MIN( khet0, rtau )                          ! check against upper limit [1/s]
        khet1 = (1. - y_clno2) * khet0
        khet2 = y_clno2 * khet0
c
c======================== DDM Begin =======================
c
        if ( lddm ) then !bk - broken
        endif
c
c======================== DDM End =======================
c
c-----CMU PM chemistry
c
      elseif (aeropt.eq.'CMU') then
        kwtr = (kph2o_1 - ngas)/nbin + 1
        if (nbin.eq.1) kwtr = kph2o_1 - ngas
        do isec = 1, nbin
          diam = sqrt(dcut(ngas+isec,1)*dcut(ngas+isec,2))
          !bk - limit hetero SO2 rxn only at dust particle surface
          saden = (6./diam) * con(kcrst_1+isec-1)/roprt(kcrst_1+isec-1)
          khet5 = khet5 + saden / ( diam * 0.5e-6 / dc_SO2 +
     &                              4. / ( cbarSO2 * gamma_SO2 ) )

          if (con(kph2o_1+isec-1).le.bdnl(kph2o_1+isec-1)) cycle
          isempty = 1
          vt = 0.
          do iaero = 1, naero
            if (iaero.ne.kwtr) then
              l = ngas + (iaero - 1)*nbin + isec
              if (con(l).gt.bdnl(l)) isempty = 0
              vt = vt + con(l)/roprt(l)
            endif
          enddo
          if (isempty.eq.1) cycle

          nden = vt / (diam*diam*diam*Pi/6.0)
          diam = diam * ( 1.0 +
     &           con(kph2o_1+isec-1)/roprt(kph2o_1+isec-1)/vt )**0.33333
          saden = Pi*diam*diam*nden

          r1 = k3k2b * (con(kph2o_1+isec-1)/18.)
     &               / (con(kpno3_1+isec-1)/62.)
          r2 = k4k2b * (con(kpcl_1+isec-1)/35.)
     &               / (con(kpno3_1+isec-1)/62.)
          rdenom = r1 + 1. + r2
          gamma_N2O5 = Afac * k2fprime * (1. - 1. / rdenom)
          y_clno2 = r2 / (r1 + r2)

          khet0 = saden / ( diam * 0.5e-6 / dc_N2O5 +
     &                      4. / ( cbarN2O5 * gamma_N2O5 ) )
          khet1 = khet1 + (1. - y_clno2) * khet0
          khet2 = khet2 + y_clno2 * khet0
        enddo
        rscl = MIN( 1., rtau / (khet1 + khet2) )
        khet1 = rscl * khet1
        khet2 = rscl * khet2
      else
        write(iout,'(//,a)') 'ERROR in KHETERO: invalid AEROPT'
        call camxerr()
      endif
c
c --- adjust rk(N2O5); rk(N2O5) will be reset by KTHERM every timestep
c
 800  continue
      rk(ihetrxn(1)) = rk(ihetrxn(1)) + khet1 * h2s / H2O   ! add gas-phase rate; [1/ppm-hr]
      if (ihetrxn(2).gt.0) rk(ihetrxn(2)) = khet2 * h2s / con(khcl)
      if (ihetrxn(5).gt.0) rk(ihetrxn(5)) = khet5 * h2s     ! [1/hr]

      return

      end
