      subroutine exptbl(iunits,rxntyp,rxnord,rxnpar)
      use chmstry
c 
c----CAMx v7.20 220430
c 
c     EXPTBL calculates a lookup table for temperature/pressure-dependent
c     rate constants using parameters from the chemistry parameters file.
c     The pressure and temperature dependance of the conversion factor from
c     cm3 molec-1 units to ppm units is accounted for.
c     The expression types supported are:
c       Type 1 Temperature independent
c       Type 2 UAM/OZIPM format Arrhenius
c       Type 3 Generalized temperature dependent
c       Type 4 Troe teperature/pressure dependent
c       Type 5 Equilibrium constant ratio to another reaction
c       Type 6 Lindemann-Hinshelwood used for OH + HNO3
c       Type 7 k = k1 + k2[M] used for OH + CO and HO2 self-reaction
c     During the model run rate constants are always interpolated from 
c     the lookup table.  The output rate constants in RKTBL have units
c     ppm-n hour-1.
c 
c      Copyright 1996 - 2022
c     Ramboll
c           
c     Modifications: 
c        09/15/04            Accept rate constants in cms units
c        10/20/05            Restructured the definition of the bins in
c                            the temperature and pressure lookup table
c        07/21/06            Fixed bug in unit conversion for rxn type 7
c        01/08/07            Fixed bug in pressure dependence term
c        12/07/14            Fixed bug in adjustment factor for rxn type 5
c        06/29/15            Explicitly cast variables in expressions
c 
c     Input arguments: 
c        iunits              input units for rate expressions 
c                              1 = ppm and minutes (ppm)
c                              2 = molec/cm3 and sec (cms)
c        rxntyp              rate constant expression type
c        rxnord              order of reaction
c        rxnpar              rate constant expression parameters
c 
c     Output arguments: 
c        none
c             
c     Routines Called:
c        none
c             
c     Called by: 
c        READCHM
c
      include "camx.prm"
c
c-----Conversions from molecule cm-3 to ppm and from hr to min to sec
c       CF1 from gas law: (n/V) = NP/RT
c       CF1 = (6.022e23 * 1.013e5 * 1e-6 ) / (8.314 * 298.0)
c     The CF1 conversion is at 298 K and 1013 mb and the
c     adjustment "factor" below generalises to other conditions
c
      double precision darren
      double precision CF1, M2H, S2M
      parameter (CF1 = 2.462d13)
      parameter (M2H = 60.d0)
      parameter (S2M = 60.d0)
c
      double precision rxnpar(MXREACT,12)
      integer rxntyp(MXREACT)
      integer rxnord(MXREACT)
c
      double precision cf, cfm, rktmp, factor, af, afm
      double precision dtmp, dprs
      double precision rk0, rk1, rk2, rk3, ratio
      integer i, j, k, iref
      real del_t, del_p
c
c-----Entry point
c
      del_t = (TEMPHI - TEMPLO) / float(NTEMPR - 1)
      del_p = (PRESHI - PRESLO) / float(NPRESR - 1)
      do j = 1,NTEMPR
        tempr(j) = TEMPLO + float(j-1)*del_t
      enddo
      do k = 1,NPRESR
        presr(k) = PRESLO + float(k-1)*del_p
      enddo
c
c-----Fill lookup table
c
      do i = 1,nreact
        if (iunits .EQ. 1) then
          cf  = 1.d0
          cfm = 1.d0
        else
          cf  = S2M * CF1**dble( (rxnord(i)-1) )
          cfm = S2M * CF1**dble( (rxnord(i)) )
        endif
c
        do j = 1,NTEMPR
          dtmp = dble( tempr(j) )
c
          do k = 1,NPRESR
            dprs = dble( presr(k) )
            factor = (298.d0/dtmp)*(dprs/1013.d0)
            af     = cf  * factor**dble( rxnord(i)-1 )
            afm    = cfm * factor**dble( rxnord(i) )
c
            if (ABS(rxntyp(i)).eq.1) then
              rktbl(i,j,k) = sngl( rxnpar(i,1)*af*M2H )
c
            elseif (ABS(rxntyp(i)).eq.2) then
              rktmp = rxnpar(i,2)*(1.d0/298.d0 - 1.d0/dtmp)
              rktbl(i,j,k) = sngl( rxnpar(i,1)*dexp(rktmp)*af*M2H )
c
            elseif (ABS(rxntyp(i)).eq.3) then
              rktbl(i,j,k) = sngl( darren(rxnpar(i,1),rxnpar(i,2),
     &                       rxnpar(i,3),rxnpar(i,4),dtmp)*af*M2H )
c
            elseif (ABS(rxntyp(i)).eq.4) then
              rk0 = darren(rxnpar(i,1),rxnpar(i,2),
     &              rxnpar(i,3),rxnpar(i,4),dtmp)*afm*1.d6
              ratio = rk0/(darren(rxnpar(i,5),rxnpar(i,6),
     &              rxnpar(i,7),rxnpar(i,8),dtmp)*af)
              rktbl(i,j,k) =  sngl(
     &          (rk0/(1.d0+ratio))*rxnpar(i,9)**(1.d0/(1.d0+
     &          (0.43429d0*dlog(ratio)/rxnpar(i,10))**2.d0))*M2H )
c
            elseif (ABS(rxntyp(i)).eq.5) then
              iref = nint(rxnpar(i,1))
              if (iunits .EQ. 1) then
                af = factor**dble( rxnord(iref)-rxnord(i) )
              else
                af = (CF1*factor)**dble( rxnord(iref)-rxnord(i) )
              endif
              rktmp = darren(rxnpar(i,2),rxnpar(i,3),
     &                rxnpar(i,4),rxnpar(i,5),dtmp)*af
              rktbl(i,j,k) = sngl( rktbl(iref,j,k) / rktmp )
c
            elseif (ABS(rxntyp(i)).eq.6) then
              rk1 = darren(rxnpar(i,1),rxnpar(i,2),
     &              rxnpar(i,3),rxnpar(i,4),dtmp)*af
              rk2 = darren(rxnpar(i,5),rxnpar(i,6),
     &              rxnpar(i,7),rxnpar(i,8),dtmp)*af
              rk3 = darren(rxnpar(i,9),rxnpar(i,10),
     &              rxnpar(i,11),rxnpar(i,12),dtmp)*afm*1.d6
              rktbl(i,j,k) = sngl( (rk1 + (rk3 / 
     &                             (1.d0+(rk3/rk2))))*M2H )
c
            elseif (ABS(rxntyp(i)).eq.7) then
              rk1 = darren(rxnpar(i,1),rxnpar(i,2),
     &              rxnpar(i,3),rxnpar(i,4),dtmp)*af
              rk2 = darren(rxnpar(i,5),rxnpar(i,6),
     &              rxnpar(i,7),rxnpar(i,8),dtmp)*afm*1.d6
              rktbl(i,j,k) = sngl( ( rk1 + rk2 )*M2H )
c
            endif
          enddo
        enddo
      enddo
c
      return
      end
c
      function darren(a,ea,b,tref,dtmp)
c
c --- Calculate a rate constant using k = A*(T/Tref)^B*exp(-Ea/T)
c
      double precision darren
      double precision a,ea,b,tref,dtmp
c
      darren = a*((dtmp/tref)**b)*dexp(-ea/dtmp)
c
      end
