      subroutine soap(ntot,mw,caer,cgas,csatT,cpx,
     &                iout,igrdchm,ichm,jchm,kchm)
      implicit none
c
c**********************************************************************
c                                                                     * 
c                                                                     * 
c  SOAP - A SUBROUTINE TO CALCULATE GAS-AEROSOL PARTITIONING OF       *
c                     SECONDARY ORGANIC SPECIES                       *
c                                                                     *
c                           DEVELOPED BY:                             *
c                                                                     * 
c                           ROSS STRADER                              *
c            DEPARTMENT OF CIVIL/ENVIRONMENTAL ENGINEERING            *
c             CARNEGIE MELLON UNIVERSITY, 5000 FORBES AVE             *
c                   PITTSBURGH, PENNSYLVANIA, 15213                   *
c                                                                     * 
c                               AND                                   *
c                                                                     * 
c                         SPYROS N. PANDIS                            *
c       DEPARTMENTS OF CHEMICAL ENGINEERING AND ENGINEERING AND       *
c                          PUBLIC POLICY                              *
c             CARNEGIE MELLON UNIVERSITY, 5000 FORBES AVE             *
c                   PITTSBURGH, PENNSYLVANIA, 15213                   *
c                                                                     * 
c                                                                     * 
c                                                                     * 
c                                                                     * 
c      ASSOCIATED SUBROUTINES:                                        *
c      CALC       - CALCULATES THE MOLE FRACTION OF EACH SOLUTION-    *
c                   FORMING COMPOUND IN THE PSEUDO-IDEAL SOLUTION     * 
c                   (AEROSOL PHASE)                                   *
c      FCN        - NONLINEAR SYSTEM TO BE SOLVED BY HYBRD            *
c      RNDNUM     - GENERATES RANDOM NUMBERS                          *
c      HYBRD      - SOLVES A SYSTEM OF NONLINEAR EQUATIONS            *
c                                                                     * 
c                                                                     * 
c**********************************************************************
c                                                                     * 
c                                                                     * 
c     THIS MODEL USES A PSEUDO-IDEAL SOLUTION THEORY TO PARTITION     *
c     SECONDARY ORGANIC SPECIES BETWEEN THE AEROSOL AND GAS PHASE.    *
c     THE THEORY IS OUTLINED IN:                                      * 
c                                                                     * 
c     STRADER, R. (1998), 'EVALUATION OF SECONDARY ORGANIC AEROSOL    *
c          FORMATION IN WINTER,' M.S. THESIS, DEPT. OF CIVIL/         *
c          ENVIRONMENTAL ENGINEERING, CARNEGIE MELLON UNIVERSITY      *
c                                                                     * 
c                                                                     * 
c**********************************************************************
c
c
c  MODIFICATIONS
c
c     Jan 2003 -gyarwood- Revised for CAMx
c     01/20/03 -bkoo-     Revised for CAMx
c     03/09/03 -bkoo-     Revised for new SOAP driver
c     05/27/05 -bkoo-     Revised for new bi-section solver
c     12/29/06 -bkoo-     Revised for the updated SOA scheme
c     04/13/10 -bkoo-     Approximate solutions if POA is dominant
c     08/23/13 -bkoo-     Revised for DDM-PM
c     03/18/14 -bkoo-     Revised for benzene SOA
c     08/25/16 -bkoo-     Updated for new SOAP
c     10/30/19 -gyarwood- Relax xtol and revise convergence tests at 99
c
c DEFINITION OF VARIABLES:
c
c  INPUTS
c
c     ntot    - total number of CG/SOA pairs
c     mw      - molecular weights of SOA species (g/mole)
c     caer    - aerosol-phase concentrations of SOA species (ug/m3)
c     cgas    - gas-phase concentrations of CG species (ug/m3)
c     csatT   - saturation concentrations of CG/SOA pairs at current T (ug/m3)
c     cpx     - pre-existing OA mass (mole)
c     iout    - standard output file unit
c     igrdchm - index for grid containing the grid cell
c     ichm    - i index of grid cell
c     jchm    - j index of grid cell
c     kchm    - k index of grid cell
c
c   OUTPUTS
c
c     caer    - aerosol-phase concentrations of SOA species (ug/m3)
c     cgas    - gas-phase concentrations of CG species (ug/m3)
c
c   VARIABLES USED WITHIN SOAP
c
c     i       - counter
c     icont   - counter
c     sum     - counter
c     nsol    - total number of solution-forming SOA species
c     smw     - molecular weights of solution-forming SOA species (g/mole)
c     scaer   - aerosol-phase concentrations of solution-forming SOA species (ug/m3)
c     scgas   - gas-phase concentrations of solution-forming SOA species (ug/m3)
c     scsat   - saturation concentrations of solution-forming SOA species (ug/m3)
c     sctot   - total concentrations of solution-forming SOA species (ug/m3)
c     znum    - counter for number of iterations
c     conmin  - use simple solution for species below this level (ug/m3)
c     cpxmin  - no pre-existing organics if cpx < cpxmin (mole) 
c     xtol    - error tolerance for bi-section method
c
c***********************************************************************
c
c VARIABLE DECLARATION
c
      real, parameter :: conmin  = 1.e-6
      real, parameter :: cpxmin  = 1.e-10
      real, parameter :: xtola   = 1.e-8
      real, parameter :: xtolr   = 1.e-5

      integer :: ntot
      real, dimension(ntot) :: mw, caer, cgas, ctot, csatT,
     &                         smw, scaer, scgas, sctot, scsat
      integer, dimension(ntot) :: idx

      real    :: cpx, sum
      integer :: iout, igrdchm, ichm, jchm, kchm
      integer :: i, icont, nsol, znum

      real    :: bb, cc, xend, fend, xmid, fmid, dx
      real    :: caer_in(ntot), cgas_in(ntot)

      logical, save :: lfirsttime = .true.
c
c***********************************************************************
c
c Entry point
c
      caer_in = caer
      cgas_in = cgas

!      if (lfirsttime) then
!        lfirsttime = .false.
!      endif

      do i=1,ntot
        ctot(i) = caer(i) + cgas(i)
      enddo
c
c CALCULATE AEROSOL-PHASE CONCENTRATION (CAER) AND GAS-PHASE 
c CONCENTRATION (CGAS) FOR NON-SOLUTION-FORMING COMPOUNDS
c COMPOUNDS THAT HAVE A CONCENTRATION OF LESS THAN conmin ARE IGN0RED
c MAP COMPOUNDS THAT FORM SOLUTIONS ONTO ARRAYS
c
      icont=0
      do i=1,ntot
         if (ctot(i).lt.conmin) then
            cgas(i) = ctot(i)
            caer(i) = 0.0
c
c========================= DDM Begin ===================================
c
            CALL DDMSOATRV (ctot(i),csatT(i),0.0,0.0,0.0,i,2,iout)
c
c========================== DDM End ====================================
c
         else
            icont=icont+1
            idx(icont) = i
            smw(icont)=mw(i)
            scsat(icont)=csatT(i)
            sctot(icont)=ctot(i)
            scaer(icont)=caer(i)
         endif
      enddo
      nsol=icont
c
c Check for a trivial solution
c
      if (nsol.eq.0) goto 1000
      if (nsol.eq.1) then
         if (cpx.lt.cpxmin) then
           scgas(1) = amin1(sctot(1), scsat(1))
           scaer(1) = sctot(1) - scgas(1)
c
c========================= DDM Begin ===================================
c
           CALL DDMSOATRV (sctot(1),scsat(1),0.0,0.0,0.0,idx(1),1,iout)
c
c========================== DDM End ====================================
c
         else ! This case has an analytical solution
           bb = scsat(1)-sctot(1)+cpx*smw(1)
           cc = -sctot(1)*cpx*smw(1)
           scaer(1) = amin1( sctot(1), .5*(-bb+SQRT(bb*bb-4.*cc)) )
           scgas(1) = sctot(1) - scaer(1)
c
c========================= DDM Begin ===================================
c
           CALL DDMSOATRV (sctot(1),scsat(1),
     &                     cpx,smw(1),scaer(1),idx(1),3,iout)
c
c========================== DDM End ====================================
c
         endif
         goto 900
      endif
      sum=0.0
      do i=1,nsol
         sum = sum + sctot(i)/scsat(i)
      enddo
      if (cpx.lt.cpxmin .and. sum.le.1.0) then
         do i=1,nsol
            scgas(i)=sctot(i)
            scaer(i)=0.0
c
c========================= DDM Begin ===================================
c
            CALL DDMSOATRV (sctot(i),scsat(i),0.0,0.0,0.0,idx(i),2,iout)
c
c========================== DDM End ====================================
c
         enddo
         goto 900
      endif
c
c Find the solution using a bi-section method (approach from max)
c
      xend = 0.0
      do i = 1, nsol
        xend = xend + sctot(i)/smw(i)
      enddo
      xend = xend + cpx
      call spfcn (nsol,sctot,scsat,scaer,smw,cpx,xend,fend)
      if ( abs(fend).le.xtola .OR.
     &    (xend-cpx).lt.cpx*xtolr ) goto 99       ! converged, no iteration needed
      if (fend.gt.cpx*xtolr) then
        write (iout,'(//,a)') ' ERROR in SOAP:'   ! first iteration will diverge
        write (iout,'(/,a)') ' ERROR: positive end point'
        write (iout,'(a,2e15.7)') 'xend,fend ',xend,fend
        goto 50
      endif
      dx = xend - cpx
      do znum = 1, 200
        dx = 0.5 * dx
        xmid = xend - dx
        call spfcn (nsol,sctot,scsat,scaer,smw,cpx,xmid,fmid)
        if (abs(fmid).le.xtola .or.
     &      dx.le.max(xtola,xtolr*xmid)) goto 100 ! converged after iteration
        if (fmid.lt.0.0) xend = xmid
      enddo
      write (iout,'(//,a)') ' ERROR in SOAP:'     ! failed to converge
      write (iout,'(/,a)') ' ERROR: max number of iterations reached'
      write (iout,'(a,2e15.7)') ' xmid,fmid ',xmid,fmid
 50   write (iout,'(a,i10,a,3i4)') ' Grid = ', igrdchm,
     &                 ' cell(i,j,k) = ', ichm, jchm, kchm
      write (iout,'(a5,5a)')
     &       ' spec','  gas in[ug/m3]','    aerosol in[ug/m3]',
     &           '     total[ug/m3]','    aerosol[ug/m3]','      c* [ug/m3]'
      write (iout,'(i5,5e18.7)')
     &       (idx(i),cgas_in(i),caer_in(i),sctot(i),scaer(i),scsat(i),i=1,nsol)
      write (iout,'(a,e15.7)') ' cpx ',cpx
      call camxerr()
c
c Converged
c
  99  xmid = xend
 100  continue
      do i=1,nsol
         scaer(i) = amin1( sctot(i), scaer(i) )
         scgas(i) = sctot(i) - scaer(i)
      enddo
c
c========================= DDM Begin ===================================
c
      CALL DDMSOAAP (sctot,scsat,smw,xmid,nsol,idx,iout)
c
c========================== DDM End ====================================
c
c REMAP COMPOUNDS THAT FORM SOLUTIONS BACK ONTO ORIGINAL ARRAYS
c      
 900  continue
      do i=1,nsol
         caer(idx(i))=scaer(i)
         cgas(idx(i))=scgas(i)
      enddo
c
c Return point
c
 1000 continue
c
      return
      end
