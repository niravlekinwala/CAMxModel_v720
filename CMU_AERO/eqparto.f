c=====================================================================
c     EQILIBRIUM PARTITIONING OF CONDENSING OR EVAPORATING ORGANICS
c     - Updated for new SOAP by bkoo (08/25/2016)
c     - revised for benzene SOA by bkoo (03/18/2014)
c     - modified for the updated SOA scheme by bkoo (12/29/2006)
c     - revised for CAMx by bkoo (03/09/03)
c     - created by bkoo (04/18/02)
c
c     Predefined pointers:
c     - ICG1  : pointer of the first condensable gas in the gas array
c     - KSOA1 : pointer of the first SOA species in the aerosol array
c     - KSOPA : pointer of polymerized anthro SOA
c     - KSOPB : pointer of polymerized biogenic SOA
c     - KPOM  : pointer of primary organic aerosol in the aerosol array
c=====================================================================
      subroutine eqparto(t,q)

      include 'dynamic.inc'
      include 'soap.inc'

      parameter (itermax = 1000) ! comparable to ITMAXEQ in equaer.inc

      real*8  q(ntotal),dq,accom,frq(nsec),frqtot,rlambda,frt
      real    ctot(nsoap),caer(nsoap),cgas(nsoap),tempk
      real    apoly(nsoap),fpoly
      real    preOM

      integer iter
      logical prevsp

      real*8 xsum(nsec),mnk,mxk ! bkoo (11/14/02)
      real   csatT(nsoap)       ! bkoo (11/14/02)

      logical, save :: lfirsttime = .true.
c
c     check consistency in MW
c
      if (lfirsttime) then
        lfirsttime = .false.
        if (nsoap.ne.ngas-4) then
          write(*,*) 'ERROR in EQPARTO: ',
     &               'nsoap .ne. ngas-4 (check dynamic.inc)'
          stop
        endif
        do i=1,nsoap
          if ( gmw(ICG1+i-1) .ne.mwsoap(i) .or. 
     &         emw(KSOA1+i-1).ne.mwsoap(i) ) then
            write(*,*) 'ERROR in EQPARTO: Inconsistent MW for CG/SOA',i
            write(*,*) '  BLOCK.F:  ',gmw(ICG1+i-1),emw(KSOA1+i-1)
            write(*,*) '  SOAP.INC: ',mwsoap(i)
            stop
          endif
        enddo
        if ( emw(KSOPA).ne.mwsopa ) then
          write(*,*) 'ERROR in EQPARTO: Inconsistent MW for SOPA'
          write(*,*) '  BLOCK.F:  ',emw(KSOPA)
          write(*,*) '  SOAP.INC: ',mwsopa
          stop
        endif
        if ( emw(KSOPB).ne.mwsopb ) then
          write(*,*) 'ERROR in EQPARTO: Inconsistent MW for SOPB'
          write(*,*) '  BLOCK.F:  ',emw(KSOPB)
          write(*,*) '  SOAP.INC: ',mwsopb
          stop
        endif
        if ( emw(KPOM).ne.mwpoa ) then
          write(*,*) 'ERROR in EQPARTO: Inconsistent MW for POA'
          write(*,*) '  BLOCK.F:  ',emw(KPOM)
          write(*,*) '  SOAP.INC: ',mwpoa
          stop
        endif
      endif

      prevsp = .false. ! currently SOAP do not utilize PREVSP

      ng = nsp*nsecx2                 ! gases
      prod=rgas*temp/(1.01325D5*pres) ! conversion from ppm to umoles/m3
      tempk = SNGL(temp)
      do i=1,nsoap
        cgas(i) = SNGL(q(ng+ICG1-1+i) * gmw(ICG1-1+i) / prod) ! ug/m3
        caer(i) = 0.0
        do j=1,nsecx2
          caer(i) = caer(i) + SNGL(q((j-1)*nsp+KSOA1-1+i))
        enddo
        ctot(i) = cgas(i) + caer(i)
        csatT(i) = csatS(i) * (cstempS(i)/tempk)
     &                      * exp((deltahS(i)/8.314)
     &                      * (1/cstempS(i)-1/tempk))
      enddo

      preOM = 0.0
      do j=1,nsecx2
        preOM = preOM + SNGL( q((j-1)*nsp+KSOPA)/emw(KSOPA) )
     &                + SNGL( q((j-1)*nsp+KSOPB)/emw(KSOPB) )
     &                + SNGL( q((j-1)*nsp+KPOM) /emw(KPOM) )
      enddo

      call get_param(igrdchm,ichm,jchm,kchm,iout,idiag)
      call soap(nsoap,mwsoap,caer,cgas,csatT,preOM,
     &          iout,igrdchm,ichm,jchm,kchm)

      rlambda=0.065d0 ! this should be consistent with that in eqpart.f

      do isec = 1,nsecx2                                   ! bkoo (11/14/02)
        indx=(isec-1)*nsp                                  ! bkoo (11/14/02)
        xsum(isec)=q(indx+KPOM)/emw(KPOM)                  ! bkoo (11/14/02)
        do isp = 1,nsoap                                   ! bkoo (11/14/02)
          xsum(isec)=xsum(isec)+q(indx+KSOA1-1+isp)/emw(KSOA1-1+isp) ! bkoo (11/14/02)
        enddo                                              ! bkoo (11/14/02)
        xsum(isec)=xsum(isec)+q(indx+KSOPA)/emw(KSOPA)
     &                       +q(indx+KSOPB)/emw(KSOPB)
        ! give any non-zero value to xsum if it is zero in the section
        if(xsum(isec).eq.0.0) xsum(isec) = tinys
      enddo                                                ! bkoo (11/14/02)

      do isp = 1,nsoap
        ! calculate DQ & ACCOM
        dq = q(ng+ICG1-1+isp)*gmw(ICG1-1+isp)/prod - DBLE(cgas(isp)) ! ug/m3
        accom = delta(ICG1-1+isp)
        ! calculate frq
        frqtot = 0.d0
        mnk = 0.d0 ! bkoo (11/14/02)
        mxk = 0.d0 ! bkoo (11/14/02)
        do isec = 1,nsecx2
          frq(isec) = qn(isec)*
     &        ( q(ng+ICG1-1+isp)*gmw(ICG1-1+isp)/prod                   ! bkoo (11/14/02)
     &        - q((isec-1)*nsp+KSOA1-1+isp)/emw(KSOA1-1+isp)/xsum(isec) ! bkoo (11/14/02)
     &        * csatT(isp) )*                                       ! bkoo (11/14/02)
     &        dsec(isec)/(1.d0+rlambda/(accom*dsec(isec)))
c          frqtot = frqtot + frq(isec)
          mnk = min(mnk,frq(isec)) ! bkoo (11/14/02)
          mxk = max(mxk,frq(isec)) ! bkoo (11/14/02)
        enddo

        if(dq.gt.0. .and. mnk.lt.0. .and. mxk.gt.0.) then     ! bkoo (11/14/02)
          do isec = 1,nsecx2                                  ! bkoo (11/14/02)
            frq(isec)=max(frq(isec)-mnk,0.d0)                 ! bkoo (11/14/02)
          enddo                                               ! bkoo (11/14/02)
        elseif(dq.lt.0. .and. mxk.gt.0. .and. mnk.lt.0.) then ! bkoo (11/14/02)
          do isec = 1,nsecx2                                  ! bkoo (11/14/02)
            frq(isec)=min(frq(isec)-mxk,0.d0)                 ! bkoo (11/14/02)
          enddo                                               ! bkoo (11/14/02)
        endif                                                 ! bkoo (11/14/02)
        do isec = 1,nsecx2                                    ! bkoo (11/14/02)
          frqtot = frqtot + frq(isec)                         ! bkoo (11/14/02)
        enddo                                                 ! bkoo (11/14/02)
        ! normalize frq
        do isec = 1,nsecx2
          frq(isec) = frq(isec) / frqtot
        enddo
        ! condense all condensing species
        if(dq.gt.0.d0) then
          do isec = 1,nsecx2
            indx=(isec-1)*nsp
            q(indx+KSOA1-1+isp) = q(indx+KSOA1-1+isp) + dq * frq(isec)
          enddo
          ! set the gas species concentration to the equilibrium value
          q(ng+ICG1-1+isp) = DBLE(cgas(isp)) * prod / gmw(ICG1-1+isp) ! ppm
        ! evaporate all evaporating species
        elseif(dq.lt.0.d0) then
          iter = 0
 100      frt = 1.d0
          do isec = 1,nsecx2
            indx=(isec-1)*nsp
            if(frq(isec).gt.0.d0) frt = MAX( MIN(q(indx+KSOA1-1+isp) ! bkoo (10/07/03)
     &                                  /(-dq*frq(isec)),frt), 0.d0)
          enddo
          frqtot = 0.d0
          do isec = 1,nsecx2
            indx=(isec-1)*nsp
            q(indx+KSOA1-1+isp) = 
     &               MAX(q(indx+KSOA1-1+isp)+frt*dq*frq(isec),0.d0)
            if(q(indx+KSOA1-1+isp).lt.tinys) frq(isec) = 0.d0
            frqtot = frqtot + frq(isec)
          enddo
          ! check if we should evaporate more
          dq = (1.d0 - frt) * dq
          if(dq.lt.-1.d-7) then
            if(frqtot.gt.tinys) then ! we have sections which are not empty
              if(iter.le.itermax) then ! check infinite loop
                iter = iter + 1
                do isec = 1,nsecx2
                  frq(isec) = frq(isec) / frqtot
                enddo
                goto 100
              endif
            endif
            ! we need to evaporate more to achieve equilibrium
            ! but we completely evaporate the species in all sections
            ! or exceed itermax
cbk            call errprt('EQUIO-F:DQ    ',t,KSOA1-1+isp,-1,dq,0)
c            if(dq.lt.-1.d-3) then                                        ! bkoo_dbg
c              write(6,'(A9,E15.5,4I4)')'EQUIO-F: ',dq,isp,ichm,jchm,kchm ! bkoo_dbg
c            endif                                                        ! bkoo_dbg
          endif
          ! now set the gas species concentration conservatively
          q(ng+ICG1-1+isp) = DBLE(ctot(isp)) ! ug/m3
          do isec = 1,nsecx2
            indx=(isec-1)*nsp
            q(ng+ICG1-1+isp) = q(ng+ICG1-1+isp) - q(indx+KSOA1-1+isp)
          enddo
          q(ng+ICG1-1+isp) = q(ng+ICG1-1+isp) * prod / gmw(ICG1-1+isp) ! ppm
        endif
      enddo
c
c Calculate fraction of polymerized SOA
c
! currently apoly is not sectionalized
! need to correct if PSAT is implemented in CMU_AERO
      fpoly = 1. - exp( -kpolya * SNGL(dt/3600.d0) ) ! polymerized fraction for anthro SOA
      do isec = 1,nsecx2
        indx=(isec-1)*nsp
        do isp = 1, nsoaa
          apoly(isp) = fpoly * q(indx+KSOA1-1+isp)
          q(indx+KSOA1-1+isp) = q(indx+KSOA1-1+isp) - apoly(isp)
          q(indx+KSOPA) = q(indx+KSOPA) + apoly(isp)
        enddo
      enddo
      fpoly = 1. - exp( -kpolyb * SNGL(dt/3600.d0) ) ! polymerized fraction for bio SOA
      do isec = 1,nsecx2
        indx=(isec-1)*nsp
        do isp = nsoaa + 1, nsoap
          apoly(isp) = fpoly * q(indx+KSOA1-1+isp)
          q(indx+KSOA1-1+isp) = q(indx+KSOA1-1+isp) - apoly(isp)
          q(indx+KSOPB) = q(indx+KSOPB) + apoly(isp)
        enddo
      enddo

      return
      end


c=====================================================================
c     This subroutine calculates the organic water absorption using
c     fits from UNIFAC for representative SOA species.
c     - bkoo (04/18/02)
c     Currently not used in PMCAMx
c=====================================================================
      subroutine orgwtr(orgtotal,orgwater)

      include 'dynamic.inc'

      real*8 orgtotal(7),orgwater ! bkoo (09/29/02)
      real*8        soawtr(7,70) ! bkoo (09/29/02)
      integer       iswtr ! bkoo (09/29/02)
      common / soawtr / soawtr,iswtr ! bkoo (09/29/02)
      real*8 wtr
      integer idrh

      if(iswtr.ne.1) then ! bkoo (09/29/02)
        orgwater = 0.d0
        return
      endif

c     modified - bkoo (09/29/02)
      idrh = IDINT(rh*100.d0)
      orgwater = 0.d0
      do i=1,7
        if(idrh.lt.30) then
          wtr = soawtr(i,30) * rh / 0.3d0
        elseif(idrh.ge.99) then
          wtr = soawtr(i,70)
        else
          wtr = soawtr(i,idrh-29)
     &        + (rh*100.d0 - DFLOAT(idrh))
     &        * (soawtr(i,idrh-28) - soawtr(i,idrh-29))
        endif
        orgwater = orgwater + wtr * orgtotal(i) * 1.d-9 ! kg/m3
      enddo

      return
      end


c=====================================================================
c     This subroutine calculates partial pressures of organics at the 
c     surface of the particles.
c     To calculate partial pressure of organics at equilibrium, we use the
c     pseudo-ideal solution theory by Bowman et al. (1997) with parameters
c     taken from Strader et al. (1999) for saturation concentrations.
c     We use the Clausius Clapeyron Equation here to calculate
c     the temperature dependence of the saturation concentration.
c     We divide the code into two subroutines so that the temperature
c     dependence calculations are done only once before DIFFUND.
c     - bkoo (04/18/02)
c     - revised for PMCAMx by bkoo (03/09/03)
c=====================================================================
      subroutine orgps1(q)

      include 'dynamic.inc'
      include 'soap.inc'

      real*8  q(ntotal),ctotp(nsoap),ctot,rsum
      real*8  csatT(nsoap)
      real*8  conmin,cpxmin,cpx
      integer icont
      integer nsol,idx(nsoap)
      logical ispre
      logical, save :: lfirsttime = .true.
      common / csatx / csatT,ispre,nsol,idx

      if (lfirsttime) then
        lfirsttime = .false.
        if(nsoap.ne.ngas-4) stop
     &      'ERROR in ORGPS1: nsoap .ne. ngas-4 (check dynamic.inc)'
      endif

      ispre = .TRUE.
      cpx = 0.0d0
      do k=1,nsecx
        cpx = cpx + q((k-1)*nsp+KPOM) /emw(KPOM)
     &            + q((k-1)*nsp+KSOPA)/emw(KSOPA)
     &            + q((k-1)*nsp+KSOPB)/emw(KSOPB)
      enddo

      do i=1,nsoap
        csatT(i) = csatS(i) * rgas * temp / mwsoap(i) * 1.d-6 ! Pa
        csatT(i) = csatT(i)
     &           * exp((deltahS(i)/8.314)*(1/cstempS(i)-1/temp))
      enddo

      ng = nsp*nsecx                  ! gases
      prod=rgas*temp/(1.01325D5*pres) ! conversion from ppm to umoles/m3
      conmin = 1.d-6                  ! should be consistent with that of SOAP
      cpxmin = 1.d-10                 ! should be consistent with that of SOAP
      icont = 0
      rsum  = 0.d0
      do i=1,nsoap
        ctot = q(ng+ICG1-1+i) * gmw(ICG1-1+i) / prod ! ug/m3
        do k=1,nsecx
          ctot = ctot + q((k-1)*nsp+KSOA1-1+i)
        enddo
        ctotp(i) = ctot * rgas * temp / mwsoap(i) * 1.d-6 ! Pa
        ! as in SOAP
        if(ctot.ge.conmin) then
          icont=icont+1
          idx(icont) = i
          rsum = rsum + ctotp(i) / csatT(i)
        else ! for ignored species flux should be zero
          do k=1,nsecx
            ps(k,ICG1-1+i) = q(ng+ICG1-1+i)*1.01325d-1*pres ! Pa
          enddo
        endif
      enddo
      nsol = icont
      if(nsol.eq.0) return
      if(rsum.le.1.d0 .and. cpx.lt.cpxmin) then
        do i=1,nsol
          do k=1,nsecx
            ps(k,ICG1-1+idx(i)) = ctotp(idx(i))
          enddo
        enddo
        nsol = 0
      endif

      return
      end

      subroutine orgps2(q)

      include 'dynamic.inc'

      real*8  q(ntotal),orgtot
      real*8  csatT(ngas-4)
      integer nsol,idx(ngas-4)
      logical ispre
      common / csatx / csatT,ispre,nsol,idx

      if(nsol.eq.0) return

      do i = 1,nsecx
        orgtot = 0.d0
        if(ispre) orgtot=q((i-1)*nsp+KPOM)/emw(KPOM) ! primary OC
        do k = 1,nsol
          orgtot=orgtot+q((i-1)*nsp+KSOA1-1+idx(k))/emw(KSOA1-1+idx(k))
        enddo
        orgtot=orgtot+q((i-1)*nsp+KSOPA)/emw(KSOPA)
     &               +q((i-1)*nsp+KSOPB)/emw(KSOPB)
        if(orgtot.lt.tinys) then
          do k = 1,nsol
            ps(i,ICG1-1+idx(k)) = 0.d0
          enddo
        else
          do k = 1,nsol
            ps(i,ICG1-1+idx(k)) = dmax1( ( q((i-1)*nsp+KSOA1-1+idx(k))
     &                / emw(KSOA1-1+idx(k)) ) / orgtot * csatT(idx(k)) ! Pa
     &                , 0.d0) ! bkoo (09/29/02)
          enddo
        endif
      enddo

      return
      end

