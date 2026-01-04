      subroutine ironchem(dt,dtpig,aero_flag,wpuff,tpuff,ppuff,cwpuff,
     &                    phpuff,ldark,conpin,cback,avgcnc,convfac,
     &                    ierr1,ierr2)
      use chmstry
      use filunit
      use tracer
      implicit none
c
c----CAMx v7.20 220430
c
c     IRONCHEM performs chemistry on the current puff reactor for one
c     time step.  It uses the LSODE solver exclusively for the gas-phase.
c
c     Local array element con(nspec+1) is used to store the concentrations
c     of non-used species, i.e., species that are in the chem solver but
c     not on the species list for this run.
c
c      Copyright 1996 - 2022
c     Ramboll
c
c     Modifications:
c        10/06/05        Removed mechanism 2
c        01/12/06        Added PM chemistry calls for AERO_CF
c        12/26/06        Removed AEROCHEM_AQ (merged into AEROCHEM_CF)
c        01/08/07        Added mechanism 6 (CB05)
c        02/21/08        Added coupling time for puff aerosol chemistry
c        03/31/08        Puff killed before PM chemistry for concs < lower bound
c        05/13/08        Modified to skip ISOROPIA for background chemistry.
c        05/28/10        Added mechanism 7 (CB6)
c        07/20/12        Added Mechanism 1 (CB06 with Iodine chemistry)
c        10/10/12        Added mechanism 2 (CB6r1)
c        06/28/13        Modified for new GREASD PiG chemistry (RADM-AQ only)
c        03/18/14        Revised for benzene SOA (AEROCHEM_CF interface change)
c        06/27/14        Added Mechanism 3 (CB6r2h)
c        11/11/14        Added Mechanism 4 (CB6r3)
c        10/13/17        Added a dummy argument to AEROCHEM_CF call (aerosol pH)
c        10/16/17        Added a dummy argument to AEROCHEM_CF call (time-weighted JNO2)
c                        Changed DUMARR2 size
c        02/21/20        Added Mechanism 1 (CB6r5) -rlb-
c
c     Input arguments:
c        dt                  gas/aqueous chemistry timestep (s)
c        dtpig               PiG coupling timestep (hr)
c        aero_flag           flag to determine which PM algorithm to call
c        wpuff               water vapor (ppm)
c        tpuff               temperature (K)
c        ppuff               pressure (mb)
c        cwpuff              cloud liquid water content (g/m3)
c        phpuff              cloud pH
c        ldark               darkness flag (T=dark)
c        conpin              puff species concentration (ppm,ug/m3)
c        cback               background species concentration (ppm,ug/m3)
c        convfac             conversion factor: umol/m3 = ppm * convfac
c
c     Output arguments:
c        conpin              puff species concentration (ppm,ug/m3)
c        avgcnc              time-step average puff concentrations
c        ierr1               error code (=1 if LSODE reports error)
c        ierr2               error code (=1 if neg concs for AEROCHEM)
c
c     Routines called:
c        LSCALL
c        AEROCHEM_CF
c
c     Called by:
c        PIGDRIVE
c
      include 'camx.prm'
      include 'flags.inc'
      include 'soap.inc'
c
c-----Subroutine names to be called by LSODE
c
      external lsrate1
      external lsrate2
      external lsrate3
      external lsrate4
      external lsrate5
      external lsrate6
      external lsrate7
      external lsrate_pig
c
      real conpin(*)
      real cback(*)
      real avgcnc(*)
c
      real rrxn(MXREACT)
      real dt,dtpig,wpuff,tpuff,ppuff,cwpuff,phpuff,convfac,dtchem
      real atm,o2,ch4,h2,cliq
      real dumvar
      real dumarr1(NNV),dumarr2(7,MXSPEC),dumarr3(MXFDDM*NNV),dumarr4(1)
      real cold_so2,cold_pso4,cold_hno3,cold_pno3,cold_nh3,cold_pnh4,
     &     cold_na, cold_pcl
      real coldb_so2,coldb_pso4
      real so2frc

      integer is,ierr1,ierr2,aero_flag
      logical ldark
c
      real con(MXSPEC+1)
c
c-----Entry point
c
      dtchem = dt/3600.
      con(nspec+1) = 0.
      do is = 1,nspec
        con(is) = conpin(is)
      enddo
c
c-----Chemistry integration, pass subroutines for mechanism used
c
      atm = 1.e6
      O2  = 2.095e5
      CH4 = 1.75
      H2  = 0.50
c
      if (ipigflg.eq.GRESPIG) then
        call lscall(lsrate_pig,dtchem,wpuff,atm,O2,CH4,H2,con,
     &              avgcnc,nreact,rrxn,ierr1,.FALSE.,ipigflg)
      elseif (idmech.eq.1) then
        call lscall(lsrate1,dtchem,wpuff,atm,O2,CH4,H2,con,
     &              avgcnc,nreact,rrxn,ierr1,.FALSE.,ipigflg)
      elseif (idmech.eq.3) then
        call lscall(lsrate3,dtchem,wpuff,atm,O2,CH4,H2,con,
     &              avgcnc,nreact,rrxn,ierr1,.FALSE.,ipigflg)
      elseif (idmech.eq.4) then
        call lscall(lsrate4,dtchem,wpuff,atm,O2,CH4,H2,con,
     &              avgcnc,nreact,rrxn,ierr1,.FALSE.,ipigflg)
      elseif (idmech.eq.5) then
        call lscall(lsrate5,dtchem,wpuff,atm,O2,CH4,H2,con,
     &              avgcnc,nreact,rrxn,ierr1,.FALSE.,ipigflg)
      elseif (idmech.eq.6) then
        call lscall(lsrate6,dtchem,wpuff,atm,O2,CH4,H2,con,
     &              avgcnc,nreact,rrxn,ierr1,.FALSE.,ipigflg)
      elseif (idmech.eq.7) then
        call lscall(lsrate7,dtchem,wpuff,atm,O2,CH4,H2,con,
     &              avgcnc,nreact,rrxn,ierr1,.FALSE.,ipigflg)
      endif
      if (ierr1 .ne. 0) return
c
c-----Perform aerosol chemistry if it's time
c
      if (aero_flag .gt. 0) then
        cliq = cwpuff
        if (tpuff.lt.273.)
     &    cliq = amax1(0.,cliq*(tpuff - tamin)/(273. - tamin))
        if (cliq.ge.cwmin .and. tpuff.gt.tamin) then
c
c-----Minimum conc check: kill the puff if below BDNL values
c
          ierr2 = 0
          if (con(kso2) .lt. 0. .or.
     &        con(khno3).lt. 0. .or.
     &        con(kn2o5).lt. 0. .or.
     &        con(knh3) .lt. 0. .or.
     &        con(kh2o2).lt. 0. .or.
     &        con(ko3)  .lt. 0. .or.
     &        con(ksulf).lt. 0.) ierr2 = 1
          do is = ngas+1,nspec
            if (is.ne.kph2o .and. con(is).lt.0.) ierr2 = 1
          enddo
          if (ierr2 .ne. 0) return
c
c-----Convert all SULF->PSO4, N2O5->HNO3 in cloud water before PM chemistry
c
          con(kpso4) = con(kpso4) + con(ksulf)*96.*convfac
          con(ksulf) = 0.
          con(khno3) = con(khno3) + 2.*con(kn2o5)
          con(kn2o5) = 0.
          cold_so2   = con(kso2)
          cold_pso4  = con(kpso4)
          cold_hno3  = con(khno3)
          cold_pno3  = con(kpno3)
          cold_nh3   = con(knh3)
          cold_pnh4  = con(kpnh4)
          cold_na    = con(kna)
          cold_pcl   = con(kpcl)
          con(kso2)  = con(kso2)  + cback(kso2)
          con(kpso4) = con(kpso4) + cback(kpso4)
          coldb_so2  = con(kso2)
          coldb_pso4 = con(kpso4)
          con(khno3) = con(khno3) + cback(khno3)
          con(kpno3) = con(kpno3) + cback(kpno3)
          con(knh3)  = con(knh3)  + cback(knh3)
          con(kpnh4) = con(kpnh4) + cback(kpnh4)
          con(kna)   = con(kna)   + cback(kna)
          con(kpcl)  = con(kpcl)  + cback(kpcl)
          so2frc = cold_so2/max(con(kso2),bdnl(kso2))
c
          call aerochem_cf(wpuff,tpuff,ppuff,cliq,phpuff,
     &                     dumvar,dumvar,
     &                     dumarr1,con,convfac,dtchem,dtpig,
     &                     dumarr2,dumarr3,dumarr4,1,1,
     &                     .FALSE.,1,aero_flag,.TRUE.,0,ldark)
c
          con(kso2)  = cold_so2  + (con(kso2)  - coldb_so2)*so2frc
          con(kpso4) = cold_pso4 + (con(kpso4) - coldb_pso4)*so2frc
          con(khno3) = cold_hno3
          con(kpno3) = cold_pno3
          con(knh3)  = cold_nh3
          con(kpnh4) = cold_pnh4
          con(kna)   = cold_na
          con(kpcl)  = cold_pcl
        endif
      endif
c
c-----Update input species arrays
c
      do is = 1,nrad
        con(is) = amax1(con(is),bdnl(is))
      enddo
      do is = 1,nspec
        conpin(is) = con(is)
      enddo
      conpin(nspec+1) = 0.
c
      return
      end
