c*** EBISOLV
c
      subroutine ebisolv(ebirxn,ebirate,hr_hox,hr_nox,hr_nxy,hr_pan,
     &                   hr_xo,
     &                   dtin,H2O,atm,O2,CH4,H2,tcell,ldark,conc,
     &                   avgcnc,avgitr,hddmjac,nfam,sddm,noxin,avgsen,
     &                   lddm,nirrrxn,rrxn_irr,ldoirr)
      use chmstry
      use filunit
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c   Description:
c
c     EBISOLV solves gas-phase chemistry for one time step using
c     the Euler Backward Iterative method with Hertel's explicit
c     solutions for four groups of species.  Subroutines are generated
c     by the Chemical Mechanism Compiler (CMC).
c
c      Copyright 1996 - 2022
c     Ramboll
c
c   Argument descriptions:
c     Inputs:
c        ebirxn            subroutine to calculate reaction rates
c        ebirate           subroutine to calculate species rates
c        hr_hox            subroutine to solve HOx group
c        hr_nox            subroutine to solve NOx group
c        hr_nxy            subroutine to solve NO3 and N2O5
c        hr_pan            subroutine to solve PAN and C2O3
c        hr_xo             subroutine to solve halogen oxides
c        dtin              time duration to be integrated (hr)
c        H2O               water vapor concentration (ppm)
c        atm               total gas concentration (M) in ppm
c        O2                oxygen concentration (ppm)
c        CH4               methane concentration (ppm)
c        H2                hydrogen concentration (ppm)
c        tcell             temperature (K)
c        ldark             logical flag for darkness
c        conc              state species concentrations (ppm)
c        hddmjac           subroutine to calculate Jacobian
c        nfam              number of sensitivity families
c        sddm              sensitivity matrix (ppm)
c        noxin             number of oxidants to store in avgsen
c        lddm              logical flag for DDM sensitivities
c        nirrrxn           number of reactions for IRR
c        ldoirr            logical flag to calculate IRR data
c
c    Outputs:
c        conc              state species concentrations (ppm)
c        avgcnc            average state species concentrations (ppm)
c        avgitr            average inetrations taken per step
c        sddm              sensitivity matrix (ppm)
c        avgsen            average sensitivities of oxidants (ppm)
c        rrxn_irr          integrated reaction rates for IRR
c
c    Subroutines called:
c        EBIRXN            calculate reaction rates
c        EBIRATE           calculate species rate of change
c        HR_HOX            solve HOx group (OH, HO2, HONO, PNA)
c        HR_NOX            solve NOx group (O1D, O, NO, NO2, O3)
c        HR_NXY            solve NO3 and N2O5
c        HR_PAN            solve PAN and C2O3
c        HR_XO             solve halogen oxides
c        DDMEBI
c        CAMXERR
c
c     Called by:
c        CHEMDRIV
c
c-----------------------------------------------------------------------
c   LOG:
c-----------------------------------------------------------------------
c
c     03/14/08   --gyarwood--    Original development
c     06/01/08   --bkoo--        Added code for DDM
c     08/23/13   --bkoo--        Added avgsen for DDM-PM
c     09/21/13   --bkoo--        Updated DDM-PM for KHETERO
c     01/03/14   --bkoo--        Added organic nitrate hydrolysis rate sensitivities
c     03/18/14   --bkoo--        Added more oxidants to avgsen for benzene SOA
c     07/06/15   --gyarwood--    Report average iterations taken per step
c     07/06/15   --gyarwood--    Revise convergence test
c     02/29/16   --bkoo--        Revised for hetero N2O5 + HCL rxn
c     07/05/16   --gyarwood--    Added HR_XO
c     07/20/16   --bkoo--        Revised for hetero INTR hydrolysis
c     09/02/16   --bkoo--        Set SDDM/KHETSEN array size with # of het rxns
c
c-----------------------------------------------------------------------
c   Include files:
c-----------------------------------------------------------------------
c
      implicit none
      include "camx.prm"
      include "ddmchm.inc"
      include "chmdbg.inc"
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      real     dtin
      real     H2O, atm, O2, CH4, H2, tcell, avgitr
      real     conc(MXSPEC+1), avgcnc(MXSPEC+1)
      logical  ldark
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      integer    MXY, MXITER, MXSTEP
      real       RTOL, ATOL, YMIN, YMAX, DTMIN
      parameter (MXY = MXSPEC+1)
      parameter (MXITER = 1000)
      parameter (MXSTEP = 1000)
      parameter (RTOL = 1.0e-3)
      parameter (ATOL = 1.0e-8)
      parameter (YMIN = 1.0e-18)
      parameter (YMAX = 1.0e3)
      parameter (DTMIN = 1.0e-5)
c
      integer  i, l, ny
      integer  istep, nstep, iter, niter, noconv
      real     y0(MXY), y1(MXY),  yh(MXY), yinit(MXY)
      real     loss(MXY), gain(MXY), rerr(MXY)
      real     rrxn(MXREACT)
      real     dtuse, dtmx, clock, wt
      logical  lavg
c
c-----------------------------------------------------------------------
c    Externals:
c-----------------------------------------------------------------------
c
      external ebirxn, ebirate, hr_hox, hr_nox, hr_nxy, hr_pan, hr_xo
c
c======================== DDM Begin =======================
c
      external hddmjac
      integer  nfam, noxin, ifam, ierr
      real     sddm(nfam,MXSPEC+NHETRXN)
      real     avgsen(nfam,noxin)
      real     ysen(nfam,MXSPEC+1)
      real     khetsen(nfam,NHETRXN)
      logical  lddm, lavgddm
c
c======================== DDM End =======================
c
c========================== Process Analysis Begin =====================
c
      integer  nirrrxn
      real     rrxn_irr(nirrrxn)
      logical  ldoirr
c
c=========================== Process Analysis End ======================
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c --- Load concentrations
c
      ny = ngas
      yh(ny+1) = 0.0 
      y1(ny+1) = 0.0 
      y0(ny+1) = 0.0 
      do i = 1, ngas
        y1(i) = conc(i)
        yinit(i) = conc(i)
        y0(i) = conc(i)
        gain(i) = 0.0
        loss(i) = 0.0
      enddo
c
c --- Set solver parameters
c
      lavg = .true.
      clock = 0.0
      dtmx = 5.01/60.0
      dtmx = dtin/MAX(2, (1+INT(dtin/dtmx)))
      dtuse = dtmx
c
c======================== DDM Begin =======================
c
      if( lddm ) then
        do i = 1, ngas
          do ifam = 1, nfam
            ysen(ifam,i) = sddm(ifam,i)
          enddo
        enddo
        do ifam = 1, nfam
          ysen(ifam,ngas+1) = 0.0
        enddo
        do i = 1, NHETRXN
          do ifam = 1, nfam
            khetsen(ifam,i) = sddm(ifam,nspec+i)
          enddo
        enddo
        lavgddm = .false.
      endif
c
c======================== DDM End =======================
c
c========================== Process Analysis Begin =====================
c
      if (ldoirr) then
        do i = 1, nirrrxn
          rrxn_irr(i) = 0.0
        enddo
      endif
c
c=========================== Process Analysis End ======================
c
c
c --- Integrate in steps until clock reaches dtin
c     istep counts both succesful and non-converged steps
c     abort if istep = MXSTEP or dtuse < DTMIN
c
      niter = 0
      do istep = 1, MXSTEP
        nstep = istep
        if( istep .EQ. MXSTEP) goto 900
c
c --- Iterate this step to convergence
c
        do iter = 1, MXITER
          niter = niter + 1
          do l = 1, ny
            yh(l) = y1(l)
            y1(l) = 0.0
          enddo
c
c --- Calculate species rates
c
          call ebirxn(ny,nreact,yh,H2O,atm,O2,CH4,H2,rk,rrxn)
          call ebirate(ny,nreact,rrxn,gain,loss)
c
c --- Solve all species using Backward Euler method
c
          do 100 l = 1, ny
            y1(l) = (y0(l) + gain(l)*dtuse) 
     &                              / (1.0 + loss(l)*dtuse/yh(l))
            y1(l) = MIN(YMAX, MAX(YMIN,y1(l)))
 100      continue
c
c --- Solve four groups of species using Hertel's equations
c     for the Backward Euler method
c         y0 is the initial conc for this step
c         yh is the last iteration
c         y1 is the next iteration
c
          call hr_hox(y0,y1,yh,H2O,atm,O2,CH4,H2,ny,rk,rrxn,
     &                 nreact,dtuse)
          call hr_nox(y0,y1,yh,H2O,atm,O2,CH4,H2,ny,rk,rrxn,
     &                 nreact,dtuse)
          call hr_nxy(y0,y1,yh,H2O,atm,O2,CH4,H2,ny,rk,rrxn,
     &                 nreact,dtuse)
          call hr_pan(y0,y1,yh,H2O,atm,O2,CH4,H2,ny,rk,rrxn,
     &                 nreact,dtuse)
          call hr_xo(y0,y1,yh,H2O,atm,O2,CH4,H2,ny,rk,rrxn,
     &                 nreact,dtuse)
c
c --- Check for convergence
c
          noconv = 0
          do l = 1, ny
             rerr(l) = MAX(YMIN,ABS(y1(l)-yh(l)))/MAX(ATOL,yh(l))
             noconv = MAX(noconv, INT(MIN(10.0, rerr(l)/RTOL)))
          enddo
          if( noconv .EQ. 0) goto 200
        enddo
 200    continue
c
c --- Update clock or reduce timestep and try again
c
        if( noconv .GT. 0) then
          dtuse = dtuse*0.5
          if( dtuse .LT. DTMIN) goto 900
          do l = 1, ny
            y1(l) = y0(l)
          enddo
        else
          clock = clock + dtuse
          if( lavg .AND. clock .GT. 0.5*dtin) then
            lavg = .false.
            wt = (clock - 0.5*dtin) / dtuse
            do i = 1, ngas
              avgcnc(i) = MAX(bdnl(i), (wt*y0(i) + (1.0-wt)*y1(i)) )
            enddo
c
c======================== DDM Begin =======================
c
            if ( lddm ) then
              do ifam = 1, nfam ! save sensitivities of previous step
                avgsen(ifam,1) = ysen(ifam,lO)
                avgsen(ifam,2) = ysen(ifam,lOH)
                avgsen(ifam,3) = ysen(ifam,lO3)
                avgsen(ifam,4) = ysen(ifam,lNO3)
                avgsen(ifam,5) = ysen(ifam,lNO)
                avgsen(ifam,6) = ysen(ifam,lHO2)
              enddo
              lavgddm = .true.
            endif
c
c======================== DDM End =======================
c
          endif
c
c======================== DDM Begin =======================
c
          if ( lddm ) then
            call ddmebi(hddmjac,nreact,nfam,ngas,NHETRXN,
     &                  dtuse,H2O,atm,O2,CH4,H2,
     &                  rk,y1,ysen,khetsen,
     &                  lN2O5,lHCL,lNTR2,lHNO3,lCLN2,lINTR,ierr)
            if (ierr.ne.0) then
              write(iout,'(//,a)') 'ERROR in EBISOLV:'
              write(iout,*) 
     &                 'Zero determinant in SGEFA in DDMEBI at ',ierr
              write(iout,*) 'igrd, i, j, k = ', igrdchm,ichm,jchm,kchm
              call camxerr()
            endif

            if ( lavgddm ) then ! calculate average sensitivities
              do ifam = 1, nfam
                avgsen(ifam,1) = ( wt * avgsen(ifam,1) +
     &                               (1.0-wt) * ysen(ifam,lO) )
                avgsen(ifam,2) = ( wt * avgsen(ifam,2) +
     &                               (1.0-wt) * ysen(ifam,lOH) )
                avgsen(ifam,3) = ( wt * avgsen(ifam,3) +
     &                               (1.0-wt) * ysen(ifam,lO3) )
                avgsen(ifam,4) = ( wt * avgsen(ifam,4) +
     &                               (1.0-wt) * ysen(ifam,lNO3) )
                avgsen(ifam,5) = ( wt * avgsen(ifam,5) +
     &                               (1.0-wt) * ysen(ifam,lNO) )
                avgsen(ifam,6) = ( wt * avgsen(ifam,6) +
     &                               (1.0-wt) * ysen(ifam,lHO2) )
              enddo
              lavgddm = .false.
            endif
          endif
c
c======================== DDM End =======================
c
c========================== Process Analysis Begin =====================
c
          if (ldoirr) then
            call ebirxn(ny,nreact,yh,H2O,atm,O2,CH4,H2,rk,rrxn)
            do i = 1, nirrrxn
              rrxn_irr(i) = rrxn_irr(i) + dtuse*rrxn(i)
            enddo
          endif
c
c=========================== Process Analysis End ======================
c
          if( clock/dtin .GT. 0.995) goto 300
          dtuse = MIN(dtuse*2.0, dtmx, dtin-clock)
          do l = 1, ny
            y0(l) = y1(l)
          enddo
        endif
c
c --- Next step
c
      enddo
c
c --- Save species concentrations in CAMx arrays
c
 300  continue
      avgitr = float(niter)/float(nstep)
      do i = 1, ngas
        conc(i) = MAX(bdnl(i), y1(i))
      enddo
c
c======================== DDM Begin =======================
c
      if ( lddm ) then
        do i = 1, ngas
          do ifam = 1, nfam
            sddm(ifam,i) = ysen(ifam,i)
          enddo
        enddo
      endif
c
c======================== DDM End =======================
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
      return
c
c-----------------------------------------------------------------------
c    Error handling:
c-----------------------------------------------------------------------
c
 900  continue
      write(iout,*)
      write(iout,*)' ERROR in the EBISOLV gas-phase chemistry solver'
      write(iout,*)
      write(iout,*) 'ldark, temp(K), water = ',
     &               ldark, tcell, H2O
      write(iout,*) 'M, O2, CH4, H2  = ',
     &               atm, O2, CH4, H2
      write(iout,'(a,i10,3i4)')' igrd, i, j, k =',
     &                           igrdchm,ichm,jchm,kchm
      write(iout,'(a,3f8.5,2i6)')' dtin, clock, dtuse, nstep, niter =',
     &                           dtin, clock, dtuse, nstep, niter
      write(iout,*) 'No  Name      New Conc    Init Conc   Rel Err'
      do l = 1,ngas
        write(iout,905) l, spname(l), y1(l), yinit(l), rerr(l)
      enddo
      write(iout,*) 'No   Rate Constant'
      do l = 1,nreact
        write(iout,910) l, rk(l)
      enddo
      write(iout,*)
c
      call camxerr()
c
 905  format(i3,2x,a7,1p5e12.3)
 910  format(i3,2x,1p3e12.3)
      end

