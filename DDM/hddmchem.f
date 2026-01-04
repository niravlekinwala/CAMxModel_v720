      subroutine hddmchem(hddmjac,ebirxn,ebirate,ddmrate,
     &                    nrxn,nspc,njac,
     &                    nfam,nrfam,ntfam,nhfam,
     &                    ipr,ipt,wft,iph,
     &                    dt,H2O,atm,O2,CH4,H2,
     &                    rk,avgcnc,sddm,ierr)
c
c----CAMx v7.20 220430
c
c     HDDMCHEM advances the 1st and 2nd order sensitivty coefficients
c     one chemistry time step. This follows the CMAQ DDM-3D approach.
c
c      Copyright 1996 - 2022
c     Ramboll
c
c     History:
c        07/16/07   --bkoo--       Original development
c        06/11/08   --bkoo--       Added rate constant sensitivity
c        12/07/14   --bkoo--       Revised the cutoff code (EPS)
c        01/08/16   --bkoo--       Updated for DDMRATE
c        08/09/18   --bkoo--       Added rate term sensitivity
c
c     Input arguments:
c        hddmjac           name of subroutine to calc full Jacobian
c        ebirxn            name of subroutine to calculate reaction rates
c        ebirate           name of subroutine to calculate species rates needed for the EBI solver
c        ddmrate           name of subroutine to calculate the remaining species rates
c        nrxn              number of reactions
c        nspc              number of state species
c        njac              nspc
c        nfam              number of sensitivity families
c        nrfam             number of rate constant sensitivity families
c        ntfam             number of rate term sensitivity families
c        nhfam             number of 2nd order sensitivity families
c        ipr               pointer array for rxns in rate constant groups;
c                            ipr(0,i) stores the number of rxns in the i-th group
c        ipt               pointer array for terms in rate term groups;
c                            ipt(0,j,i) = 1 if i-th group includes j-th rxn; otherwise, 0
c                            ipt(1,j,i) = # of reactants of j-th rxn included in the i-th group
c                            ipt(1+k,j,i) = pointer to k-th reactant of j-th rxn included ...
c                            ipt(2+ipt(1,j,i),j,i) = # of products of j-th rxn included ...
c                            ipt(2+ipt(1,j,i)+k,j,i) = pointer to k-th product of j-th rxn included ...
c        wft               weighting factors of reactant term sens in net product formation sens formula
c        iph               pointer array for 1st order sensitivity families
c                            to which 2nd order sensitivity is computed
c        dt                chemistry time step (hr)
c        H2O               water vapor concentration (ppm)
c        atm               total gas concentration (M) in ppm
c        O2                oxygen concentration (ppm)
c        CH4               methane concentration (ppm)
c        H2                hydrogen concentration (ppm)
c        rk                reaction rate constant (ppm/hr)
c        avgcnc            average state species concentrations (ppm)
c        sddm              sensitivity matrix sddm (ppm)
c
c     Output arguments:
c        sddm              updated sensitivity matrix
c        ierr              error flag from LU decomposition
c
c     Routines called:
c        HDDMJAC
c        EBIRXN
c        EBIRATE
c        DDMRATE
c        SGEFA
c        SGESL
c
c     Called by:
c        CHEMDRIV
c
      implicit none

      integer  nrxn,nspc,njac,nfam,nrfam,ntfam,nhfam,ierr
      integer  ipr(0:nrxn,nrfam),ipt(0:nspc,nrxn,ntfam),iph(2,nhfam)
      real     wft(nspc,nrxn,ntfam)
      real     dt,H2O,atm,O2,CH4,H2
      real     rk(nrxn),avgcnc(nspc+1)
      real     sddm(nfam,nspc)

      external hddmjac,ebirxn,ebirate,ddmrate

      real     r_eps
      parameter ( r_eps = 1.e18 ) ! 1/eps

      real     ymid(njac+1),rjac(njac+1,njac+1)
      real     aa(njac,njac),a1(njac,njac),bb(njac)
      real     sen1d(njac+1),smid(nfam,njac)
      real     rrxn(nrxn),rtmp(nrxn)
      real     gain(njac+1),loss(njac+1),prod(njac)
      real     flg

      integer  ipvt(njac)
      integer  i,j,n,ip,ipx,ipz1,ipz2,ip1,ip2,iptmp1,iptmp2
c
c---  Entry point
c
      ierr = 0
c
c---  Load concentrations
c
      do i = 1, nspc
        ymid(i) = avgcnc(i)
      enddo
      ymid(njac+1) = 0.0
c
c---  Get the Jacobian
c
      flg = 1.0   ! keep (effective) 1st-order rxns
      rjac = 0.0  ! initialize the Jacobian
      call hddmjac(njac,nrxn,H2O,atm,O2,CH4,H2,flg,rjac,ymid,rk)
c
c---  Fill the coeff matrix
c
      do i = 1, njac
        do j = 1, njac
          aa(i,j) = -0.5 * dt * rjac(i,j)
          a1(i,j) = -aa(i,j)
        enddo
        aa(i,i) = aa(i,i) + 1.0
        a1(i,i) = a1(i,i) + 1.0
      enddo
c
c---  Factor the coeff matrix
c
      call sgefa(aa,njac,njac,ipvt,ierr)
      if (ierr.ne.0) goto 990 ! zero determinant
c
c---  Fill the rxn rate vector if rate const (or term) sens to be computed
c
      if ( nrfam.gt.0 .or. ntfam.gt.0 ) then
        call ebirxn(njac,nrxn,ymid,H2O,atm,O2,CH4,H2,rk,rrxn)
      endif
c
c---  Solve for 1st order sensitivities
c
      ipz1 = nfam - nrfam - ntfam - nhfam
      ipz2 = nfam - ntfam - nhfam

      do ip = 1, nfam - nhfam

        do i = 1, nspc
          sen1d(i) = sddm(ip,i)
        enddo
        do i = 1, njac
          sen1d(i) = sen1d(i) * MIN(1.0, AINT( ABS(sen1d(i))*r_eps ))
        enddo

        prod = 0.0
        if ( ip .gt. ipz2 ) then ! rate terms sens
          ipx = ip - ipz2
          do j = 1, nrxn
            if ( ipt(0,j,ipx).eq.1 ) then
              rtmp = 0.0
              rtmp( j ) = rrxn( j )
              call ebirate(njac,nrxn,rtmp,gain,loss)
              call ddmrate(njac,nrxn,rtmp,gain,loss)
              do i = 1, ipt(1,j,ipx) ! loop over reactants selected
                prod(ipt(1+i,j,ipx)) = prod(ipt(1+i,j,ipx))
     &                               - loss(ipt(1+i,j,ipx)) * dt
     &                               * wft(1+i,j,ipx) ! weighting factors for net product formation sens
              enddo
              do i = 1, ipt(2+ipt(1,j,ipx),j,ipx) ! loop over products selected
                prod(ipt(2+ipt(1,j,ipx)+i,j,ipx))
     &                  = prod(ipt(2+ipt(1,j,ipx)+i,j,ipx))
     &                  + gain(ipt(2+ipt(1,j,ipx)+i,j,ipx)) * dt
              enddo
            endif
          enddo
        elseif ( ip .gt. ipz1 ) then ! rate constant sens
          rtmp = 0.0
          ipx = ip - ipz1
          do j = 1, ipr(0,ipx)
            rtmp( ipr(j,ipx) ) = rrxn( ipr(j,ipx) )
          enddo
          call ebirate(njac,nrxn,rtmp,gain,loss)
          call ddmrate(njac,nrxn,rtmp,gain,loss)
          do i = 1, njac
            prod(i) = ( gain(i) - loss(i) ) * dt
          enddo
        endif

        do i = 1, njac
          bb(i) = 0.0
          do j = 1, njac
            bb(i) = bb(i) + a1(i,j) * sen1d(j)
          enddo
          bb(i) = bb(i) + prod(i)
        enddo

        call sgesl(aa,njac,njac,ipvt,bb,0)

        do i = 1, njac
          bb(i) = bb(i) * MIN(1.0, AINT( ABS(bb(i))*r_eps ))
          smid(ip,i) = 0.5 * (sen1d(i) + bb(i))
        enddo
        do i = 1, nspc
          sddm(ip,i) = bb(i)
        enddo

      enddo
c
c---  Solve for 2nd order sensitivities
c
      do ip = nfam - nhfam + 1, nfam

        ip1 = iph(1, ip - nfam + nhfam)
        ip2 = iph(2, ip - nfam + nhfam)

        prod = 0.0
        iptmp1 = ip1
        iptmp2 = ip2
        do n = 1, 2 ! loop for two 1st-order sens parameters
          if ( iptmp1 .gt. ipz2 ) then ! rate terms sens
          ! NOT SUPPORTED YET (see RDOPTDDM.F)
          elseif ( iptmp1 .gt. ipz1 ) then ! rate constant sens
            rtmp = 0.0
            ipx = iptmp1 - ipz1
            do j = 1, ipr(0,ipx)
              rtmp( ipr(j,ipx) ) = rk( ipr(j,ipx) )
            enddo
            flg = 1.0   ! keep (effective) 1st-order rxns
            rjac = 0.0  ! initialize the Jacobian
            call hddmjac(njac,nrxn,H2O,atm,O2,CH4,H2,flg,rjac,ymid,rtmp)
            do i = 1, njac
              do j = 1, njac
                prod(i) = prod(i) + rjac(i,j) * smid(iptmp2,j)
              enddo
            enddo
            if ( iptmp1.eq.iptmp2 ) then
              prod = 2.0 * prod
              EXIT
            endif
          endif
          iptmp1 = ip2
          iptmp2 = ip1
        enddo

        do i = 1, njac
          sen1d(i) = smid(ip1,i)
        enddo
        sen1d(njac+1) = 0.0

        flg = 0.0   ! remove (effective) 1st-order rxns
        rjac = 0.0  ! initialize the Jacobian
        call hddmjac(njac,nrxn,H2O,atm,O2,CH4,H2,flg,rjac,sen1d,rk)

        do i = 1, nspc
          sen1d(i) = sddm(ip,i)
        enddo
        do i = 1, njac
          sen1d(i) = sen1d(i) * MIN(1.0, AINT( ABS(sen1d(i))*r_eps ))
        enddo

        do i = 1, njac
          bb(i) = 0.0
          do j = 1, njac
            bb(i) = bb(i) + a1(i,j) * sen1d(j)
            prod(i) = prod(i) + rjac(i,j) * smid(ip2,j)
          enddo
          bb(i) = bb(i) + dt * prod(i)
        enddo

        call sgesl(aa,njac,njac,ipvt,bb,0)

        do i = 1, njac
          bb(i) = bb(i) * MIN(1.0, AINT( ABS(bb(i))*r_eps ))
        enddo
        do i = 1, nspc
          sddm(ip,i) = bb(i)
        enddo

      enddo
c
c---  Return point
c
990   return
      end
