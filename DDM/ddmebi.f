      subroutine ddmebi(hddmjac,nrxn,nfam,njac,nhet,
     &                  dt,H2O,atm,O2,CH4,H2,
     &                  rk,y1,ysen,khetsen,
     &                  lN2O5,lHCL,lNTR2,lHNO3,lCLN2,lINTR,ierr)
c
c----CAMx v7.20 220430
c
c     DDMEBI advances the 1st order sensitivty coefficients one time
c     step using Backward Euler method.
c
c      Copyright 1996 - 2022
c     Ramboll
c
c     History:
c        06/01/08   --bkoo--       Original development
c        09/21/13   --bkoo--       Updated for KHETERO sens
c        01/03/14   --bkoo--       Updated for organic nitrate hydrolysis rate sensitivities
c        12/07/14   --bkoo--       Revised the cutoff code (EPS)
c        02/29/16   --bkoo--       Revised for hetero N2O5 + HCL rxn
c        07/20/16   --bkoo--       Revised for hetero INTR hydrolysis
c        09/02/16   --bkoo--       Set KHETSEN array size with # of het rxns; revised for hetero SO2
c
c     Input arguments:
c        hddmjac           name of subroutine to calc full Jacobian
c        nrxn              number of reactions
c        nfam              number of sensitivity families
c        njac              number of gas and radical species (ngas)
c        nhet              number of heterogeneous rxns
c        dt                EBI solver time step (hr)
c        H2O               water vapor concentration (ppm)
c        atm               total gas concentration (M) in ppm
c        O2                oxygen concentration (ppm)
c        CH4               methane concentration (ppm)
c        H2                hydrogen concentration (ppm)
c        rk                reaction rate constant (ppm/hr)
c        y1                concentrations (radical+gas) at t+dt (ppm)
c        ysen              sensitivity (radical+gas) matrix (ppm)
c        khetsen           sensitivity to k_hetero (1/hr or 1/ppm-hr)
c        lN2O5             spc idx to N2O5
c        lHCL              spc idx to HCL
c        lNTR2             spc idx to NTR2
c        lHNO3             spc idx to HNO3
c        lCLN2             spc idx to CLN2
c        lINTR             spc idx to INTR
c
c     Output arguments:
c        ysen              sensitivity (radical+gas) matrix (ppm)
c        ierr              error flag from LU decomposition
c
c     Routines called:
c        HDDMJAC
c        SGEFA
c        SGESL
c
c     Called by:
c        EBISOLV
c
      implicit none

      integer  nrxn,nfam,njac,nhet,ierr
      integer  lN2O5,lHCL,lNTR2,lHNO3,lCLN2,lINTR
      real     dt,H2O,atm,O2,CH4,H2
      real     rk(nrxn),y1(njac+1),ysen(nfam,njac)
      real     khetsen(nfam,nhet)
      real     stmp1,stmp2,stmp3,stmp4

      external hddmjac

      real     r_eps
      parameter ( r_eps = 1.e25 ) ! 1/eps

      real     rjac(njac+1,njac+1),aa(njac,njac),bb(njac)
      real     flg

      integer  ipvt(njac)
      integer  i,j,ip
c
c---  Entry point
c
      ierr = 0
      y1(njac+1) = 0.0
c
c---  Get the Jacobian
c
      flg = 1.0   ! keep (effective) 1st-order rxns
      rjac = 0.0  ! initialize the Jacobian
      call hddmjac(njac,nrxn,H2O,atm,O2,CH4,H2,flg,rjac,y1,rk)
c
c---  Fill the coeff matrix
c
      do i = 1, njac
        do j = 1, njac
          aa(i,j) = -dt * rjac(i,j)
        enddo
        aa(i,i) = aa(i,i) + 1.0
      enddo
c
c---  Factor the coeff matrix
c
      call sgefa(aa,njac,njac,ipvt,ierr)
      if (ierr.ne.0) goto 990 ! zero determinant
c
c---  Solve for 1st order sensitivities
c
      do ip = 1, nfam

        stmp1 = y1(lN2O5)*H2O*dt*khetsen(ip,1)
        stmp2 = y1(lN2O5)*y1(lHCL)*dt*khetsen(ip,2)
        stmp3 = y1(lNTR2)*dt*khetsen(ip,3)
        stmp4 = y1(lINTR)*dt*khetsen(ip,4)
        do i = 1, njac
          bb(i) = ysen(ip,i)
        enddo
        if ( lN2O5.le.njac ) bb(lN2O5) = bb(lN2O5) - stmp1 - stmp2      ! add k_het sens term
        if ( lHCL .le.njac ) bb(lHCL)  = bb(lHCL)  - stmp2              ! add k_het sens term
        if ( lNTR2.le.njac ) bb(lNTR2) = bb(lNTR2) - stmp3              ! add k_het sens term
        if ( lINTR.le.njac ) bb(lINTR) = bb(lINTR) - stmp4              ! add k_het sens term
        if ( lHNO3.le.njac ) bb(lHNO3) = bb(lHNO3) + 2.*stmp1
     &                                     + stmp2 + stmp3 + stmp4      ! add k_het sens term
        if ( lCLN2.le.njac ) bb(lCLN2) = bb(lCLN2) + stmp2              ! add k_het sens term
        do i = 1, njac
          bb(i) = bb(i) * MIN(1.0, AINT( ABS(bb(i))*r_eps ))
        enddo

        call sgesl(aa,njac,njac,ipvt,bb,0)

        do i = 1, njac
          bb(i) = bb(i) * MIN(1.0, AINT( ABS(bb(i))*r_eps ))
          ysen(ip,i) = bb(i)
        enddo

      enddo
c
c---  Return point
c
990   return
      end
