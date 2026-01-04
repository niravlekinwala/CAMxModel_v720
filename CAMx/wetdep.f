      subroutine wetdep(m1,m2,m3,i0,j0,ixcl,jycl,igrd,ncol,nrow,nlay,
     &                  nspcs,deltat,area,depth,tempk,press,
     &                  cwc,pwr,pws,pwg,cph,densfac,dtout,ipa_cel,conc,
     &                  fluxtmp,depfld1,depfld2,wetfld,iptrsa)
      use chmstry
      use procan
      use tracer
c
c----CAMx v7.20 220430
c 
c     WETDEP modifies vertical concentration profiles for a given grid via 
c     precipitation processes.  This subroutine has been completely rewritten
c     for CAMx v4.
c 
c      Copyright 1996 - 2022
c     Ramboll
c           
c     Modifications:
c        06/15/03   Fixed bug in scaling of tmass
c        02/11/04   Updated empirical rainfall relationships
c        04/21/04   Incorporated sectional PM
c        05/03/05   Revised gas scavenging calculation
c        05/04/05   Revised DDM adjustments
c        06/20/05   Revised to handle liquid/frozen cloud and precip water
c        10/21/05   Added dissociation to Henry's Law for NH3, HNO3, SO2
c        08/31/06   Added map scale factor
c        07/16/07   Fixed PA bug - adjust cipr if evap
c        07/16/07 -bkoo-     Added check for HDDM
c        05/12/08 -cemery-   Revised deposited liquid concentration to use 
c                            correct rainfall volume
c        07/16/08 -bkoo-     Added DDM turn-off flag
c        07/16/08 -gwilson-  Now accumulates mass flux summary outside OMP loop
c        12/22/09 -cemery-   Improved treatment of "c0" array
c        09/23/13 -bkoo-     Placed levap-setting code outside of the spc loop
c                            Added DDM-PM code
c        10/09/14 -cemery-   Restructured for subgrid convective model
c        12/07/14 -bkoo-     Modified for VBS
c        03/17/16 -cemery-   Fixed bug to dissallow tmass<0 (Ondrej Vlcek, CHMI)
c        06/06/16 -cemery-   Simplified gas scavenging approach, removed
c                            rain evap
c        08/04/16 -cemery-   Revised mean PM size in the fine mode (<2.5um) to
c                            account for more efficient PM scavenging at the
c                            high end of the size distribution
c        08/16/21 -cemery-   Fixed calculation errors in wet deposition
c                            liquid concentration
c
c     Input arguments:
c        m1                  number of slice columns
c        m2                  number of slice rows
c        m3                  number of slice layers
c        i0                  slice i-offset
c        j0                  slice j-offset
c        ixcl                i-cell index
c        jycl                j-cell index
c        igrd                grid index 
c        ncol                number of grid columns
c        nrow                number of grid rows
c        nlay                number of grid layers 
c        nspcs               number of total species
c        deltat              time step (s)
c        area                cell area (m2)
c        depth               cell depth (m)
c        tempk               temperature field (K)
c        press               pressure field (mb)
c        cwc                 cloud water content (g/m3)
c        pwr                 rain water content (g/m3)
c        pws                 snow water content (g/m3)
c        pwg                 graupel water content (g/m3)
c        cph                 cloud water pH
c        densfac             factor to convert to umol/m3
c        dtout               output frequency (minutes)
c        ipa_cel             identify if cell is in a IPRM sub-domain
c        conc                concentration field (umol/m3, ug/m3)
c             
c     Output arguments: 
c        conc                concentration field (umol/m3, ug/m3)
c        fluxtmp             temporary array for fluxes
c        depfld1             wet deposited mass (mol/ha, g/ha)
c        depfld2             surface liquid concentrations (mol/l, g/l)
c        wetfld              wet depostion array for tracer species
c             
c     Routines called: 
c        SCAVRAT
c        HENRYFNC
c             
c     Called by: 
c        CIGDRIVE
c 
      implicit none
      include 'camx.prm'
      include 'flags.inc'
c
      integer nlay,nspcs,igrd
      real deltat,area,densfac,dtout
      real, dimension(nlay) :: tempk,press,cwc,pwr,pws,pwg,cph,depth
      real, dimension(nspcs) :: depfld1,depfld2,fluxtmp
      real, dimension(MXLAYER,MXSPEC) :: conc
c
      logical lcloud,ltop,lgraupl
      integer kbot,ktop,k,kzcl,l,isemptyf,isemptyc,kwtr,
     &        isec,isempty,iaero
      real rd,rhoh2o,cellvol,rhoair,delc,delm,convfac,cmin,hlaw,
     &     dscav,gscav,ascav,qtf,qtc,vtf,vtc,psizec,ascavf,
     &     roprta,psize,ascavc,qt,vt,cwat,pwat,pwtr,rconst,volume,
     &     dtfall,drpvel,delc0,ceqc,ceqp,ceqpp,pliq
c
      real pp
      real rr
      real volrat
      real c0(MXSPEC)
      real tmass(MXSPEC)
      real delr(MXSPEC)
c
c======================== Probing Tool Begin ===========================
c
      integer m1,m2,m3,i0,j0,ixcl,jycl
      integer ncol,nrow
      integer icls
      integer*8 iptrsa
      real fc2r,fr2c,fc2rc0,fr2cc0,cyctr
      logical lzerc
c
      real c0trac(MXTRSP)
      real tmtrac(MXTRSP)
      real delcls(MXALCLS)
      real c0cls(MXALCLS)
      real concls(MXALCLS)
      real cytcls(MXALCLS)
c
      integer ipa_cel(ncol,nrow,nlay), ipa_idx
      real, dimension(m1,m2,notimespc) :: wetfld
c
c========================= Probing Tool End ============================
c
      data rd /287./         ! Dry air gas constant (J/K/kg)
      data rhoh2o /1.e6/     ! water density (g/m3)
      data rconst /8.206e-2/ ! gas constant (l.atm/mol.K)
c
c-----Entry point
c
c-----If precip is not getting to the ground, there is no wet dep
c
      if (pwr(1).lt.cwmin .and. pws(1).lt.cwmin .and.
     &    pwg(1).lt.cwmin) return
c
c-----Scan column for layers containing cloud top
c
      kbot = 1
      pp = pwr(1) + pws(1) + pwg(1)
      volrat = pp/rhoh2o                 ! drop volume/air volume
      rr = (volrat/1.0e-7)**1.27         ! rainfall rate (mm/hr)
c
      ktop = 0
      do k = 1,nlay
        if (cwc(k).lt.cwmin .and. ktop.ne.0) then
          goto 26
        endif
        if (cwc(k).ge.cwmin) ktop = k
      enddo
      ktop = nlay
  26  continue
c
      do l = 1,nspec
        c0(l)    = 0.
        tmass(l) = 0.
      enddo
c
c-----Loop over layers and species for precipitating columns
c
      ltop = .true.
      pliq = 0.
      do 30 k = ktop,kbot,-1
        kzcl = k
        lcloud = .false.
        lgraupl = .false.
        if (cwc(k).ge.cwmin) lcloud = .true.
        if (pwg(k).ge.cwmin) lgraupl = .true.
        cellvol = area*depth(k)
        rhoair = 100.*press(k)/(rd*tempk(k))
        cwat = cwc(k)
        pwtr = pp + pliq
        if (tempk(k).lt.273.)
     &    cwat = max(0.,cwat*(tempk(k) - tamin)/(273. - tamin))
c
c======================== Probing Tool Begin ===========================
c
        if( ((lddm .OR. lhddm) .AND. lddmcalc(igrd)) .OR.
     &      (ltrace .AND. tectyp.NE.RTRAC .AND. tectyp.NE.RTCMC) ) then
          if (ltop) then
            do l = 1,ntotsp
              tmtrac(l) = 0.
              c0trac(l) = 0.
            enddo
          endif
          do icls = 1,ntrcls
            delcls(icls) = 0.
            concls(icls) = 0.
            c0cls(icls)  = 0.
            cytcls(icls) = 0.
          enddo
        endif
c
c======================== Probing Tool End =============================
c
c-----Calculate scavenging for soluble gas species
c
        do 40 l = nrad+1,ngas
          if( henry0(l).LT.1.e-6 ) goto 40
          convfac = densfac*(273./tempk(k))*(press(k)/1013.)
          cmin = bdnl(l)*convfac
          conc(k,l) = max(cmin,conc(k,l))
          call henryfnc(l,henry0(l),tfact(l),tempk(k),cph(k),knh3,
     &                  khno3,kso2,hlaw)
          hlaw = hlaw*rconst*tempk(k)
          ceqc = conc(k,l)*(1. - 1./(1. + hlaw*cwat/rhoh2o))
          call henryfnc(l,henry0(l),tfact(l),tempk(k),5.,knh3,
     &                  khno3,kso2,hlaw)
          hlaw = hlaw*rconst*tempk(k)
         
          pwat = pwtr
          if (tempk(k).lt.273. .and. rscale(l).gt.0.) pwat = 0.
          ceqp  = conc(k,l)*(1. - 1./(1. + hlaw*pwat/rhoh2o))
          ceqpp = (conc(k,l) - ceqc)*(1. - 1./(1. + hlaw*pwat/rhoh2o))
c
          call scavrat(.false.,lcloud,lgraupl,tamin,rr,tempk(k),
     &                 rhoair,diffrat(l),0.,0.,dscav,gscav,ascav,drpvel)
c
          dtfall = depth(k)/drpvel
          if (c0(l).lt.ceqp .or. ceqp.eq.0.) then
            delc  = ceqc*(1. - exp(-dscav*deltat)) +
     &              ceqpp*(1. - exp(-gscav*deltat))
            delc0 = ceqc*(1. - exp(-dscav*min(dtfall,deltat))) +
     &              ceqpp*(1. - exp(-gscav*min(dtfall,deltat)))
            if (ceqp.eq.0.) then
              delc  = min(delc,conc(k,l)-cmin)
              delc0 = min(delc0,conc(k,l)-cmin)
            else
              delc  = min(delc,ceqp-c0(l),conc(k,l)-cmin)
              delc0 = min(delc0,ceqp-c0(l),conc(k,l)-cmin)
            endif
          else
            delc  = (ceqp - c0(l))*(1. - exp(-gscav*deltat))
            delc0 = (ceqp - c0(l))*(1. - exp(-gscav*min(dtfall,deltat)))
          endif
          delm = delc*cellvol
          if (delm.lt.0.) then
            delm = max(delm,-tmass(l))
            delc = delm/cellvol
          endif
c
c======================== Probing Tool Begin ===========================
c
c  --- Accumulate OSAT/PSAT tracer adjustments by tracer class ---
c
          if( ltrace .AND. tectyp.NE.RTRAC .AND. tectyp.NE.RTCMC ) then
            cyctr = 0.5*min(conc(k,l),c0(l))
            if (pwat.eq.0) then 
              cyctr = cyctr*(1. - exp(-dscav*deltat))
            else
              cyctr = cyctr*(1. - exp(-(dscav+gscav)*deltat))
            endif
            do icls = 1,ntrcls
              concls(icls) = concls(icls) + conc(k,l)*fluxmap(l,icls)
              c0cls(icls)  = c0cls(icls)  + c0(l)*fluxmap(l,icls)
              delcls(icls) = delcls(icls) + delc*fluxmap(l,icls)
              cytcls(icls) = cytcls(icls) + cyctr*fluxmap(l,icls)
            enddo
          endif
c
c  --- For DDM convert delc to flux ---
c
          if( (lddm .OR. lhddm) .AND. lddmcalc(igrd) ) then
            if( delc.GT.0. ) then
               fc2r = delc/conc(k,l)
               fr2c = 0.
               fc2rc0 = delc0/conc(k,l)
               fr2cc0 = 0.
            elseif( c0(l) .GT. 1.0e-20 ) then
               fc2r = 0.
               fr2c = -delc/c0(l)
               fc2rc0 = 0.
               fr2cc0 = -delc0/c0(l)
            else
               fc2r = 0.
               fr2c = 0.
               fc2rc0 = 0.
               fr2cc0 = 0.
            endif

            lzerc = .FALSE.
            if( conc(k,l)-delc .LE. 1.1*cmin ) lzerc = .TRUE.

            call adjddmc0(l,ixcl,jycl,kzcl,fc2r,fr2c,fc2rc0,fr2cc0,
     &                    c0trac,tmtrac,cellvol,lzerc,
     &                    m1,m2,m3,ntotsp,ptconc(iptrsa))
          endif
c
c-----PA change from wet deposition
c
          if( lipr ) then
            if( ipa_cel(ixcl+i0,jycl+j0,k) .GT. 0 ) then
              ipa_idx = ipa_cel(ixcl+i0,jycl+j0,k)
              cipr(IPR_WDEP,ipa_idx,l) = cipr(IPR_WDEP,ipa_idx,l) - delc
            endif
          endif
c
c========================= Probing Tool End ============================
c
c-----Update the cell concentration and rain mass
c
          conc(k,l) = conc(k,l) - delc
          c0(l) = max(c0(l)+delc0,0.)
          tmass(l) = max(tmass(l)+delm,0.)
 40     continue
        pliq = pliq + cwat*(1. - exp(-dscav*min(dtfall,deltat)))
c
c-----Calculate scavenging for particulate species
c
        if (naero .gt. 0) then
c
c-----Recalculate particle size for wet diameter
c     Adjust geometric mean particle diameters <2.5um by 10x to account for
c     mass scavenging over the whole size distribution (assuming geometric
c     standard deviation of 2.0)
c
c-----CF 2-bin scheme
c
          if (lchem .AND. aeropt.eq.'CF') then
            isemptyf = 1
            isemptyc = 1
            qtf = 0. ! fine dry total mass
            qtc = 0. ! coarse dry total mass
            vtf = 0. ! fine dry total volume
            vtc = 0. ! coarse dry total volume
            do l = ngas+1,nspec
              if (dcut(l,2).lt.(dcut(kph2o,2)+1.e-5)) then     ! fine
                if (l.ne.kph2o) then
                  if (conc(k,l).gt.bdnl(l)) isemptyf = 0
                  qtf = qtf + conc(k,l)
                  vtf = vtf + conc(k,l)/roprt(l)
                endif
              else                                             ! coarse
                if (conc(k,l).gt.bdnl(l)) isemptyc = 0
                qtc = qtc + conc(k,l)
                vtc = vtc + conc(k,l)/roprt(l)
                psizec = sqrt(dcut(l,1)*dcut(l,2))
              endif
            enddo
            ascavf = 0.
            if (isemptyf.eq.0) then
              roprta = (qtf + conc(k,kph2o)) /
     &                 (vtf + conc(k,kph2o)/roprt(kph2o))
              psize = sqrt(dcut(kph2o,1)*dcut(kph2o,2))
              psize = 1.e-6*psize*(1. + 
     &                conc(k,kph2o)/roprt(kph2o)/vtf)**0.33333
              if (psize.lt.2.5e-6) psize = psize*10.
              call scavrat(.true.,lcloud,lgraupl,tamin,rr,tempk(k),
     &                 rhoair,0.,psize,roprta,dscav,gscav,ascavf,drpvel)
            endif
            ascavc = 0.
            if (isemptyc.eq.0) then
              roprta = qtc/vtc
              psize = psizec*1.e-6
              call scavrat(.true.,lcloud,lgraupl,tamin,rr,tempk(k),
     &                 rhoair,0.,psize,roprta,dscav,gscav,ascavc,drpvel)
            endif
            do l = ngas+1,nspec
              if (dcut(l,2).lt.(dcut(kph2o,2)+1.e-5)) then
                delr(l) = 1. - exp(-ascavf*deltat)
              else
                delr(l) = 1. - exp(-ascavc*deltat)
              endif
            enddo
c
c-----CMU multi-section scheme
c
          elseif (lchem .AND. aeropt.eq.'CMU') then
            kwtr = (kph2o_1 - ngas)/nbin + 1
            if (nbin.eq.1) kwtr = kph2o_1 - ngas
            do isec = 1, nbin
              isempty = 1
              qt = 0. ! dry total mass
              vt = 0. ! dry total volume
              do iaero = 1,naero
                if (iaero.ne.kwtr) then
                  l = ngas + (iaero-1)*nbin + isec
                  if (conc(k,l).gt.bdnl(l)) isempty = 0
                  qt = qt + conc(k,l)
                  vt = vt + conc(k,l)/roprt(l)
                endif
              enddo
              ascav = 0.
              if (isempty.eq.0) then
                roprta = (qt + conc(k,kph2o_1-1+isec))/
     &                   (vt + conc(k,kph2o_1-1+isec)/
     &                   roprt(kph2o_1))
                psize = sqrt(dcut(ngas+isec,1)*dcut(ngas+isec,2))
                psize = 1.e-6*psize*(1. + conc(k,kph2o_1-1+isec)/
     &                  roprt(kph2o_1)/vt)**0.33333
                call scavrat(.true.,lcloud,lgraupl,tamin,rr,tempk(k),
     &                  rhoair,0.,psize,roprta,dscav,gscav,ascav,drpvel)
              endif
              do iaero = 1,naero
                l = ngas + (iaero-1)*nbin + isec
                delr(l) = 1. - exp(-ascav*deltat)
              enddo
            enddo
c
c-----Other scheme
c
          else
            do l = ngas+1,nspec
              psize = 1.e-6*sqrt(dcut(l,1)*dcut(l,2))
              if (psize.lt.2.5e-6) psize = psize*10.
              call scavrat(.true.,lcloud,lgraupl,tamin,rr,tempk(k),
     &                rhoair,0.,psize,roprt(l),dscav,gscav,ascav,drpvel)
              delr(l) = 1. - exp(-ascav*deltat)
            enddo
          endif
c
          do 50 l = ngas+1,nspec
            delc = 0.
            delm = 0.
            cmin = bdnl(l)
            conc(k,l) = max(cmin,conc(k,l))
            delc = conc(k,l)*delr(l)
            delc = min(delc,conc(k,l)-cmin)
c
            conc(k,l) = conc(k,l) - delc
            delm = delc*cellvol
            tmass(l) = tmass(l) + delm
c
c========================= Probing Tool Begin ============================
c
c-----PA change from wet deposition
c
            if( lipr ) then
              if( ipa_cel(ixcl+i0,jycl+j0,k) .GT. 0 ) then
                ipa_idx = ipa_cel(ixcl+i0,jycl+j0,k)
                cipr(IPR_WDEP,ipa_idx,l) = cipr(IPR_WDEP,ipa_idx,l) - 
     &                                     delc
              endif
            endif
c
            if( ltrace .AND. tectyp.NE.RTRAC .AND.
     &                       tectyp.NE.RTCMC ) then
              do icls = 1,ntrcls
                delcls(icls) = delcls(icls) + delc*fluxmap(l,icls)
              enddo
            endif
c
c  --- Adjust DDM sensitivities for PM species ---
c
            if ( lddm .AND. lddmcalc(igrd) ) then
              fc2r = delr(l)
              fr2c = 0.
              fc2rc0 = 0.
              fr2cc0 = 0.

              lzerc = .FALSE.
              if( conc(k,l) .LE. 1.1*cmin ) lzerc = .TRUE.

              call adjddmc0(l,ixcl,jycl,kzcl,fc2r,fr2c,fc2rc0,fr2cc0,
     &                      c0trac,tmtrac,cellvol,lzerc,
     &                      m1,m2,m3,ntotsp,ptconc(iptrsa))
            endif
c
c========================= Probing Tool End ============================
c
 50       continue
        endif
c
c======================== Source Apportion Begin =======================
c
c  --- Adjust OSAT/PSAT tracer by tracer class ---
c
        if( ltrace .AND. tectyp.NE.RTRAC .AND. tectyp.NE.RTCMC ) then
           call adjstc0(m1,m2,m3,igrd,ixcl,jycl,kzcl,delcls,concls,
     &                  cytcls,c0cls,tmtrac,cellvol)
        endif
c
c======================== Source Apportion End =========================
c
        ltop = .false.
 30   continue
c
c-----Increment deposition flux arrays
c
      volume = 1.e-3*rr*area               ! mm/hr * m2 -> m3 (in 1 hr)
      do l = nrad+1,nspec
        fluxtmp(l) = -tmass(l)
        depfld1(l) = 1.e-2*tmass(l)/area   ! ug(mol) / m2 -> g(mol)/ha
        depfld2(l) = 1.e-9*tmass(l)/volume ! ug(mol) / m3 -> g(mol)/l
      enddo
c
c======================== Source Apportion Begin =======================
c
      if( lptdepout) then
        do l = 1,notimespc
          wetfld(ixcl,jycl,l) = wetfld(ixcl,jycl,l) + 
     &                          1.e-2*tmtrac(l)/area
        enddo
      endif
c
c======================== Source Apportion End =======================
c
      return
      end
