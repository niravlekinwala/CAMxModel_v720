      subroutine wetdeprt(m1,m2,m1rt,m2rt,ixcl,jycl,nlay,nrtsp,deltat,area,
     &                    depth,tempk,press,cwc,pwr,pws,pwg,cph,densfac,
     &                    conc,solmass,vegmass,fsurf)
      use chmstry
      use rtracchm
c
c----CAMx v7.20 220430
c 
c     This version is for the RTRAC species.
c     WETDEP modifies vertical concentration profiles for a given grid via 
c     precipitation processes.  This subroutine has been completely rewritten
c     for CAMx v4.
c 
c      Copyright 1996 - 2022
c     Ramboll
c           
c     Modifications:
c        05/01/03   Taken from WETDEP for regular model
c        06/15/03   Fixed bug in scaling of tmass
c        02/11/04   Updated empirical rainfall relationships
c        05/03/05   Revised gas scavenging calculation
c        06/20/05   Revised to handle liquid/frozen cloud and precip water
c        10/21/05   Added dissociation to Henry's Law for NH3, HNO3, SO2
c        08/31/06   Added map scale factor
c        10/29/09   Added flux of RTRAC species to surface model
c        12/22/09   Improved treatment of "c0" array
c        10/09/14   Restructured for subgrid convective model
c        03/17/16   Fixed bug to dissallow tmass<0 (Ondrej Vlcek, CHMI)
c        06/06/16   Simplified gas scavenging approach, removed rain evap
c        08/04/16   Revised mean PM size in the fine mode (<2.5um) to
c                   account for more efficient PM scavenging at the
c                   high end of the size distribution
c        08/26/21   Surface file doubles as (default) dep output or (flagged)
c                   surface model output
c
c     Input arguments:
c        m1                  number of slice columns
c        m2                  number of slice rows
c        m1rt                number of slice columns for RTRAC
c        m2rt                number of slice rows for RTRAC
c        ixcl                i-cell index
c        jycl                j-cell index
c        nlay                number of layers 
c        nrtsp               number of species in RTRAC array
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
c        conc                concentration field (umol/m3, ug/m3)
c        solmass             surface soil mass field (mol/ha, g/ha)
c        vegmass             surface veg mass field (mol/ha, g/ha)
c        fsurf               fractional landuse coverage
c             
c     Output arguments: 
c        conc                concentration field (umol/m3, ug/m3)
c        solmass             surface soil mass field (mol/ha, g/ha)
c        vegmass             surface veg mass field (mol/ha, g/ha)
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
      include 'deposit.inc'
      include 'rtracsrf.inc'
c
      integer m1,m2,m1rt,m2rt,ixcl,jycl,nlay,nrtsp
      real deltat,area,densfac
      real, dimension(nlay) :: tempk,press,cwc,pwr,pws,pwg,cph,depth
      real, dimension(MXLAYER,MXTRSP) :: conc
      real, dimension(m1rt,m2rt,nrtsp) :: solmass,vegmass
      real, dimension(m1,m2,NLUW89) :: fsurf

      logical lcloud,lgraupl
      integer kbot,ktop,k,l,n
      integer iwater
      real rd,rhoh2o,cellvol,rhoair,delc,delm,convfac,cmin,hlaw,
     &     dscav,gscav,ascav,psize,cwat,pwat,pwtr,rconst,
     &     depmas,dtfall,drpvel,delc0,ceqc,ceqp,ceqpp,pliq
c
      real pp
      real rr
      real volrat
      real tmass(MXTRSP)
      real c0(MXTRSP)
c
      data rd /287./         ! Dry air gas constant (J/K/kg)
      data rhoh2o /1.e6/     ! water density (g/m3)
      data rconst /8.206e-2/ ! gas constant (l.atm/mol.K)
c
c-----Entry point
c
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
      do l = 1,nrtsp
        tmass(l) = 0.
        c0(l)    = 0.
      enddo
c
c-----Loop over layers and species for precipitating columns
c
      pliq = 0.
      do 30 k = ktop,kbot,-1
        lcloud = .false.
        lgraupl = .false.
        if (pwg(k).ge.cwmin) lgraupl = .true.
        if (cwc(k).ge.cwmin) lcloud = .true.
        cellvol = area*depth(k)
        rhoair = 100.*press(k)/(rd*tempk(k))
        cwat = cwc(k)
        pwtr = pp + pliq
        if (tempk(k).lt.273.)
     &    cwat = max(0.,cwat*(tempk(k) - tamin)/(273. - tamin))
c
c-----Calculate scavenging for soluble gas species
c
        do 40 l = 1,nrtgas
          if( rthlaw(l).LT.1.e-6 ) goto 40
          convfac = densfac*(273./tempk(k))*(press(k)/1013.)
          cmin = rtlbnd(l)*convfac
          conc(k,l) = max(cmin,conc(k,l))
          call henryfnc(0,rthlaw(l),rttfact(l),tempk(k),cph(k),knh3,
     &                  khno3,kso2,hlaw)
          hlaw = hlaw*rconst*tempk(k)
          ceqc = conc(k,l)*(1. - 1./(1. + hlaw*cwat/rhoh2o))
          call henryfnc(0,rthlaw(l),rttfact(l),tempk(k),5.,knh3,
     &                  khno3,kso2,hlaw)
          hlaw = hlaw*rconst*tempk(k)
          pwat = pwtr
          if (tempk(k).lt.273. .and. rscale(l).gt.0.) pwat = 0.
          ceqp = conc(k,l)*(1. - 1./(1. + hlaw*pwat/rhoh2o))
          ceqpp = (conc(k,l) - ceqc)*(1. - 1./(1. + hlaw*pwat/rhoh2o))
c
          call scavrat(.false.,lcloud,lgraupl,tamin,rr,tempk(k),
     &                 rhoair,rtdrate(l),0.,0.,dscav,gscav,ascav,drpvel)
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
c-----Update the cell concentration and rain mass
c
          conc(k,l) = conc(k,l) - delc
          tmass(l) = max(tmass(l)+delm,0.)
          c0(l) = max(c0(l)+delc0,0.)
 40     continue
        pliq = pliq + cwat*(1. - exp(-dscav*min(dtfall,deltat)))
c
c-----Calculate scavenging for particulate species
c     Adjust geometric mean particle diameters <2.5um by 10x to account for
c     mass scavenging over the whole size distribiution (assuming geometric
c     standard deviation of 2.0)
c
        if (nrtaero .gt. 0) then
          do 50 l = nrtgas+1,nrtrac
            delc = 0.
            delm = 0.
            cmin = rtlbnd(l)
            conc(k,l) = amax1(cmin,conc(k,l))
            psize = 1.e-6*sqrt(rtlcut(l)*rtucut(l))
            if (psize.lt.2.5e-6) psize = psize*10.
c
            call scavrat(.true.,lcloud,lgraupl,tamin,rr,tempk(k),
     &               rhoair,0.,psize,rtdens(l),dscav,gscav,ascav,drpvel)
c
            delc = conc(k,l)*(1. - exp(-ascav*deltat))
            delc = amin1(delc,conc(k,l)-cmin)
c
            conc(k,l) = conc(k,l) - delc
            delm = delc*cellvol
            tmass(l) = tmass(l) + delm
 50       continue
        endif
 30   continue
c
c-----Increment surface mass array
c
      do l = 1,nrtrac
        depmas = 1.e-2*tmass(l)/area
        if( lsrfmodrt ) then
          do n = 1,NLUW89
            iwater = 7
            if (n.ne.iwater) then
              solmass(ixcl,jycl,l) = solmass(ixcl,jycl,l) +
     &                             fsoil(n)*fsurf(ixcl,jycl,n)*depmas
              vegmass(ixcl,jycl,l) = vegmass(ixcl,jycl,l) +
     &                      (1. - fsoil(n))*fsurf(ixcl,jycl,n)*depmas
            endif
          enddo
        else
          vegmass(ixcl,jycl,l) = vegmass(ixcl,jycl,l) + depmas
        endif
      enddo
c
      return
      end
