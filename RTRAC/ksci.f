c**** KSCI
c
      subroutine ksci(temp, pres)
      use rtcmcchm
c 
c----CAMx v7.20 220430
c 
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c     KSCI calculates temperature/pressure-dependent rate constants
c     for the RTRAC CMC solver.  The rate expressions follow the
c     SCICHEM conventions
c
c     The SCICHEM expression types supported are:
c
c     ID_RAD  = 0 Radiation dependent reaction rate identifier
c     ID_CONS = 1 Constant reaction rate identifier
c     ID_TEMP = 2 Temperature dependent reaction rate identifier
c     ID_PRES = 3 Pressure dependent reaction rate identifier
c     ID_H2O  = 4 H2O dependent reaction rate identifier
c     ID_M    = 5 M dependent reaction rate identifier
c     ID_LWC  = 6 LWC (and drop diam) dependent reaction rate identifier
c     ID_PRES2= 7 Pressure (type 2) dependent reaction rate identifier
c     ID_EQM  = 8 reverse decomposition rate dependent on forward rate
c     ID_O2   = 9 O2 dependent reaction rate identifier
c     ID_N2   =10 N2 dependent reaction rate identifier
c     ID_FOFF1=11 Falloff (type 1) reaction rate identifier
c     ID_FOFF2=12 Falloff (type 2) reaction rate identifier
c     ID_FOFF3=13 Falloff (type 3) reaction rate identifier
c     ID_CH4  =14 CH4 dependent reaction rate identifier
c     ID_H2OB =15 H2O dependent (type 2) reaction rate identifier
c 
c      Copyright 1996 - 2022
c    Ramboll
c
c     Argument description:
c      Inputs:
c       temp         R     temperature (K)
c       pres         R     pressure (mbar)
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c    07/06/07   --gyarwood--    Original development
c    13/04/07   --gyarwood--    Rate expressions for SCICHEM2012
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c    Argument declaration:
c-----------------------------------------------------------------------
c
      real    temp, pres
c 
c-----------------------------------------------------------------------
c    External functions:
c-----------------------------------------------------------------------
c
      double precision darren
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      real*8     cf1
      parameter (cf1 = 2.462D13)
c
      integer    i, iord, ityp, iref
      real*8     cf, cfm, rktmp, ppmpres, factor, t2h
      real*8     dpres, dtemp, dtemp2, dtemp3
      real*8     ratio, rk0, rk1, rk2, rk3
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
      dpres = DBLE(pres)
      dtemp = DBLE(temp)
      dtemp2 = dtemp*dtemp
      dtemp3 = dtemp2*dtemp
c
c --- Conversions from molecule cm-3 to ppm and from hr to min to sec
c       cf1 from gas law: (n/V) = NP/RT
c       cf1 = (6.022e23 * 1.013e5 * 1e-12 ) / (8.314 * 298.0)
c     The cf1 conversion is at 298 K and 1013 mb and the
c     adjustment "factor" below generalises to other conditions
c
      if( itunit .EQ. 1 ) then
         t2h = 3600.0D0
      elseif( itunit .EQ. 2 ) then
         t2h = 60.0D0
      else
         t2h = 1.0D0
      endif
c
      do i = 1,nrxnrtc
        ityp = ityprtc(i)
        iord = max(1, nrct(i))
        if (icunit .EQ. 1) then
          cf  = cf1**DBLE((iord-1))
          cfm = cf1**DBLE((iord))
        else
          cf  = 1.0D0
          cfm = 1.0D0
        endif
        ppmpres = 1.0D6*dpres/1013.0D0
        factor = ( (298.D0/dtemp)*(dpres/1013.D0) )**DBLE((iord-1))
c
c --- Note that SCICHEM rate expression types are offset by +100
c     For reference: function darren(a,ea,b,tref,temp)
c                    darren = a*((temp/tref)**b)*dexp(-ea/temp)
c
c Intialize rate constants to zero
c
        rkrtc(i) = 0.0D0
c
c Type 0 is photolysis so nothing to do here
c
c Type 1 is a constant value; adjust for conc and time units
c
        if (ityp.EQ.101) then
          rkrtc(i) = DBLE(rkprmrtc(i,1))*cf*t2h
c
c Type 2 is k = A * B^T * exp(C/T)
c
        elseif (ityp.EQ.102) then
          rkrtc(i) = darren(DBLE(rkprmrtc(i,1)),DBLE(-rkprmrtc(i,2)),
     &                      DBLE(rkprmrtc(i,3)),300.0D0,dtemp)*cf*t2h
c
c Type 3 is Troe with k0 and kinf set = A * B^T and F = 0.6
c
        elseif (ityp.EQ.103) then
          rk0 = darren(DBLE(rkprmrtc(i,1)),0.0D0,DBLE(rkprmrtc(i,2)),
     &                               300.0D0,dtemp)*cfm*ppmpres
          ratio = rk0/(darren(DBLE(rkprmrtc(i,3)),0.0D0,
     &                         DBLE(rkprmrtc(i,4)),300.0D0,dtemp)*cf)
          rkrtc(i) = (rk0/(1.0D0+ratio))*0.6D0**(1.0D0/(1.0D0+
     &                                        LOG10(ratio)**2))*t2h
c
c Type 4 is type 2 with embedded [H2O] in ppm
c
        elseif (ityp.EQ.104) then
          rkrtc(i) = darren(DBLE(rkprmrtc(i,1)),DBLE(-rkprmrtc(i,2)),
     &                  DBLE(rkprmrtc(i,3)),300.0D0,dtemp)*t2h*cf/cf1
c
c Use of type 2 with embedded M, O2, N2, CH4, H2O, O2*M, H2
c
        elseif (ityp.EQ.105 .OR. ityp.EQ.109 .OR. ityp.EQ.110
     &             .OR. ityp.EQ.114 .OR. ityp.EQ.115
     &             .OR. ityp.EQ.117
     &             .OR. ityp.EQ.121 ) then
          rkrtc(i) = darren(DBLE(rkprmrtc(i,1)),DBLE(-rkprmrtc(i,2)),
     &                      DBLE(rkprmrtc(i,3)),300.0D0,dtemp)*cf*t2h
c
c Type 6 involves cloud water and is not supported
c
c
c k = k0 * (1 + 0.6*P) used for OH + CO
c
        elseif (ityp.EQ.107) then
          rkrtc(i) = DBLE(rkprmrtc(i,1))
     &                    * ( 1.0D0 + 0.6D0*ppmpres/1.0D6 )*t2h
c
c Type 8 sets a ratio to another rate constant k2 = k1 / ratio
c with ratio set using type 3 expression
c
        elseif (ityp.EQ.108) then
          factor = 1.0D0
          iref = NINT(rkprmrtc(i,1))
          if (icunit .EQ. 2) 
     &       cf  = cf1**(nrct(iref)-iord)
          rktmp = darren(DBLE(rkprmrtc(i,2)),
     &               DBLE(-rkprmrtc(i,3)),0.0D0,300.0D0,dtemp)*cf
          rkrtc(i) = rkrtc(iref) / rktmp 
c
c Troe with k0 and kinf set  = A * (B/300)^T
c
        elseif (ityp.EQ.111 .OR. ityp.EQ.112) then
          rk0 = darren(DBLE(rkprmrtc(i,1)),0.0D0,DBLE(rkprmrtc(i,2)),
     &                               300.0D0,dtemp)*cfm*ppmpres
          ratio = rk0/(darren(DBLE(rkprmrtc(i,3)),0.0D0,
     &                        DBLE(rkprmrtc(i,4)),300.0D0,dtemp)*cf)
          rkrtc(i) = (rk0/(1.0D0+ratio))*0.6D0**(1.0D0/(1.0D0+
     &                                        LOG10(ratio)**2))*t2h
c
c Type 13 is Lindeman-Hinshelwood used for OH + HNO3
c
        elseif (ityp.EQ.113) then
          rk1 = darren(DBLE(rkprmrtc(i,1)),DBLE(-rkprmrtc(i,2)),
     &                      0.0D0,300.0D0,dtemp)*cf
          rk2 = darren(DBLE(rkprmrtc(i,3)),DBLE(-rkprmrtc(i,4)),
     &                      0.0D0,300.0D0,dtemp)*cf
          rk3 = darren(DBLE(rkprmrtc(i,5)),DBLE(-rkprmrtc(i,6)),
     &                      0.0D0,300.0D0,dtemp)*cfm
          rkrtc(i) = ( rk1 + (rk3*ppmpres /
     &                    (1.0D0 + (rk3*ppmpres / rk2) ) ) )*t2h
c
c Type 22 is  k = A * (B/300)^T * exp(C/T); type 16 has M & O2 embedded
c (like Type 3 but with Tref = 300)
c
        elseif (ityp.EQ.116 .OR. ityp.EQ.122) then
          rkrtc(i) = darren(DBLE(rkprmrtc(i,1)),DBLE(-rkprmrtc(i,2)),
     &                      DBLE(rkprmrtc(i,3)),300.0D0,dtemp)*cf*t2h
c
c Type 18 is the full Troe with Fc and k = A * B^T * exp(C/T)
c
        elseif (ityp.EQ.118) then
          rk0 = darren(DBLE(rkprmrtc(i,1)),DBLE(-rkprmrtc(i,3)),
     &              DBLE(rkprmrtc(i,2)),300.0D0,dtemp)*cfm*ppmpres
          ratio = rk0/(darren(DBLE(rkprmrtc(i,4)),DBLE(-rkprmrtc(i,6)),
     &                        DBLE(rkprmrtc(i,5)),300.0D0,dtemp)*cf)
          rkrtc(i) = (rk0/(1.0D0+ratio))*DBLE(rkprmrtc(i,7))**(1.0D0/(1.0D0+
     &                                        LOG10(ratio)**2))*t2h
c
c Type 19 and 20 are pressure dependent k = k0 + k1 * M
c with 20 having H2O embedded
c
        elseif (ityp.EQ.119 .OR. ityp.EQ.120) then
          rk1 = darren(DBLE(rkprmrtc(i,1)),DBLE(-rkprmrtc(i,2)),
     &                      0.0D0,300.0D0,dtemp)*cf
          rk2 = darren(DBLE(rkprmrtc(i,3)),DBLE(-rkprmrtc(i,4)),
     &                      0.0D0,300.0D0,dtemp)*cfm
          rkrtc(i) = ( rk1 + (rk2*ppmpres) )*t2h
c
        endif
        rkrtc(i) =  rkrtc(i)*factor
        srkrtc(i) =  SNGL( rkrtc(i) )
      enddo
c
      return
      end
