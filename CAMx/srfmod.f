      subroutine srfmod(date,dtime,dx,dy,cellat,cellon,swe,snowrat,
     &                  snowfrc,snowalb,zenith,ctrns,ldrk,pwr,tempk,
     &                  temp0,press,press0,wind,depth,nspc,fsurf,
     &                  lrdlai,lai,solmas,vegmas,reemis)
c
      use chmstry
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c   Description:
c     This routine determines fraction of deposited material that is
c     sorbed to surface elements.  Unsorbed fraction is available for
c     re-emission.  Sorbed fraction undergoes chemistry on soil/vegetation
c     using information from the chemistry parameters input file.
c
c      Copyright 1996 - 2022
c     Ramboll
c
c   Argument descriptions:
c     Inputs:
c       date    I  model date (YYJJJ)
c       dtime   R  time step (hours)
c       dx,dy   R  cell size (m)
c       cellat  R  latitude (deg)
c       cellon  R  longitude (deg)
c       swe     R  snow cover water equivalent (m)
c       snowrat R  snow water accumulation rate (m/hr)
c       snowfrc R  snow cover fraction
c       snowalb R  snow albedo
c       zenith  R  zenith angle (deg)
c       ctrns   R  cloud transmissivity (unitless)
c       ldrk    L  darkness flag
c       pwr     R  rain water content (g/m3)
c       tempk   R  layer 1 temperature (K)
c       temp0   R  surface temperature (K)
c       press   R  layer 1 pressure (mb)
c       press0  R  surface pressure (mb)
c       wind    R  layer 1 wind (m/s)
c       depth   R  layer 1 thickness (m)
c       nspc    I  number of surface model species
c       fsurf   R  fraction LU
c       lrdlai  L  LAI inputs were read
c       lai     R  gridded LAI
c       solmas  R  Soil mass (mol/ha)
c       vegmas  R  Vegetation mass (mol/ha)
c     Outputs:
c       solmas  R  Soil mass (mol/ha)
c       vegmas  R  Vegetation mass (mol/ha)
c       reemis  R  Re-emission rate (mol/s)
c
c   Called by:
c       CHEMDRIV 
c
c   Routines called:
c       CALDATE
c       MICROMET
c       VE_GAS
c
c-----------------------------------------------------------------------
c   LOG:
c-----------------------------------------------------------------------
c
c     04/24/13   ---cemery---    Original development (from RTRAC)
c     08/22/14   ---cemery---    Extended Zhang LU to surface model
c     08/25/14   ---cemery---    Added snow effects to surface model
c
c-----------------------------------------------------------------------
c   Include files:
c-----------------------------------------------------------------------
c
      implicit none
      include 'camx.prm'
      include 'deposit.inc'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      integer date
      integer nspc
      real    dtime
      real    dx, dy
      real    cellat, cellon
      real    swe
      real    snowrat
      real    snowfrc
      real    snowalb
      real    zenith
      real    ctrns
      real    pwr
      real    fsurf(NLU)
      real    tempk, temp0
      real    press, press0
      real    wind
      real    depth
      real    lai
      real    solmas(nspc)
      real    vegmas(nspc)
      real    reemis(nspc)
      logical ldrk
      logical lrdlai
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      integer   l
      integer   isp
      integer   m
      integer   mbin, latbin
      integer   isesn
      integer   idate, iyear, month, iday
      integer   nday(12)
      real      fshads
      real      fshadv
      real      fzenith
      real      deg2rad
      real      zenang
      real      totland
      real      smkrat
      real      smjrat
      real      effrates
      real      effratev
      real      delcons
      real      delconv
      real      eps
      real      rhoh2o
      real      volrat
      real      rain
      real      snow
      real      pbl,ustar,el,psih,wstar
      real      ve
      real      z0
      real      ssorb,slch
      real      solemis,tmpsol
      real      vegemis,tmpveg
      real      smas,vmas
      real      vemis(nspc)
      real      ref_lai, lai_ref_intpl, rlai, lai_f
      logical lstable
      logical lsnow
c
      data deg2rad /0.01745329/
      data eps /1.0e-20/
      data rhoh2o /1.e6/     ! water density (g/m3)
      data nday/31,28,31,30,31,30,31,31,30,31,30,31/
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
      idate = date
      call caldate(idate)
      iyear = idate/10000
      month = (idate - 10000*iyear)/100
      if (mod(iyear,4).eq.0) nday(2)=29
      iday = idate - 100*(idate/100)

      lsnow = .false.
      if (swe.ge.0.01) lsnow = .true.  ! Snow effects occur when > 10 cm deep
c
c-----Prepare for LU-dependent surface roughness calculation
c
      mbin = month
      if (cellat.lt.0.) then
        mbin = mod(month+6,12)
        if (mbin.eq.0) mbin = 12
      endif
      latbin = 1
      if (abs(cellat).gt.20.) then
        latbin = 2
      elseif (abs(cellat).gt.35.) then
        latbin = 3
      elseif (abs(cellat).gt.50.) then
        latbin = 4
      elseif (abs(cellat).gt.75.) then
        latbin = 5
      endif
      if ((cellat.gt.50. .and. cellat.lt.75.) .and.
     &    (cellon.gt.-15. .and. cellon.lt.15.)) latbin = 3
      isesn = iseason(latbin,mbin)
c
c-----Use input snow cover to set season, if specified
c
      if (swe.ge.0.001) isesn = 4   ! Snow cover check > 1 cm
c
c-----Determine cell-average relative LAI
c
      if (lrdlai .and. NLU.eq.NLUZ03) then
        ref_lai = 0.
        do m = 1,NLU
          lai_ref_intpl = lai_ref(m,mbin) +
     &                    float(min(nday(mbin),iday))/float(nday(mbin))*
     &                    (lai_ref(m,mbin+1) - lai_ref(m,mbin))
          ref_lai = ref_lai + fsurf(m)*lai_ref_intpl
        enddo
        rlai = lai/(ref_lai + 1.e-10)
      endif
c
c-----Loop over land use; surface roughness for water is dependent on
c     wind speed
c
      do l = 1,nspc
        vemis(l) = 0.
      enddo
      do 10 m = 1,NLU
        if (fsurf(m).lt.0.01) goto 10
c
c-----Set surface roughness for Wesely (1989) scheme
c
        if (NLU.eq.NLUW89) then
          z0 = z0lu(m,isesn)
          if (m.eq.7) z0 = amax1(z0,2.0e-6*wind**2.5)
c
c-----Set surface roughness and LAI for Zhang (2003) scheme
c
        else
          lai_f = lai_ref(m,mbin) +
     &            float(min(nday(mbin),iday))/float(nday(mbin))*
     &            (lai_ref(m,mbin+1) - lai_ref(m,mbin))
          if (lrdlai) then
            lai_f = lai_f*rlai
            lai_f = amin1(lai_ref(m,15),lai_f)
            lai_f = amax1(lai_ref(m,14),lai_f)
          endif
          if (m.eq.1 .or. m.eq.3) then
            z0 = 2.0e-6*wind**2.5
          else
            if (z02(m).gt.z01(m)) then
              z0 = z01(m) + (lai_f - lai_ref(m,14))/
     &                      (lai_ref(m,15) - lai_ref(m,14))*
     &                      (z02(m) - z01(m))
            else
              z0 = z01(m)
            endif
          endif
        endif
c
c-----Get surface layer micrometeorological parameters for this cell and
c     landuse type
c
        call micromet(tempk,temp0,press,press0,depth/2.,wind,z0,pbl,
     &                ustar,el,psih,wstar,lstable)
c
c-----Loop over species, and calculate re-emission velocity for this cell,
c     landuse, and current species.  Aggregated into cell-average re-emission
c     velocity.
c
        do l = 1,nspc
          isp = idsmsp(l)
          if ((NLU.eq.NLUW89 .and. m.eq.7) .or. (NLU.eq.NLUZ03 .and.
     &                                    (m.eq.1 .or. m.eq.3)))then
            ve = 0.
          else
            call ve_gas(z0,depth/2.,psih,ustar,diffrat(isp),ve)
          endif
          vemis(l) = vemis(l) + ve*fsurf(m)
        enddo
 10   continue
c
c-----Calculate fraction of sorbed and un-sorbed material:
c     Un-sorbed fraction is available for re-emission
c     Sorbed fraction remains on surface and undergoes chemistry, leaching
c     and vegetative penetration
c
      do l = 1,nspc
        ssorb = smssrb(l)
        if (lsnow) ssorb = smisrb(l)
        solemis = (vemis(l)/depth)*solmas(l)/(1. + ssorb) ! (mol/ha/s)
        vegemis = (vemis(l)/depth)*vegmas(l)/(1. + smvsrb(l)) ! (mol/ha/s)
        tmpsol  = MAX(eps,solmas(l) - solemis*dtime*3600.) ! (mol/ha)
        tmpveg  = MAX(eps,vegmas(l) - vegemis*dtime*3600.) ! (mol/ha)
        solemis = (solmas(l) - tmpsol)/(dtime*3600.)          ! (mol/ha/s)
        vegemis = (vegmas(l) - tmpveg)/(dtime*3600.)          ! (mol/ha/s)
        solmas(l) = tmpsol                                    ! (mol/ha)
        vegmas(l) = tmpveg                                    ! (mol/ha)
        reemis(l) = (solemis + vegemis)*dx*dy*1.e-4           ! (mol/s)
      enddo
c
c --- Calculate photolysis factors (shade, zenith, clouds, snow)
c
      zenang = deg2rad*zenith
      fshads = 0.0
      fshadv = 0.0
      fzenith = 0.0
      totland = 0.0

      if (.not.ldrk) then
        fshadv = 0.5
        do m = 1,NLU
          if (NLU.eq.NLUW89) then
            fshads = fshads + fsoil(m)*fsurf(m)
          else
            fshads = fshads + fsoilz(m)*fsurf(m)
          endif
          totland = totland + fsurf(m)
        enddo
        fshads = fshads/totland
        if (lsnow) then
          fshads = fshads + snowfrc*snowalb
          fshadv = fshadv + snowfrc*snowalb
        endif
        fzenith = COS(zenang)
      endif
c
c-----Chemical decay of precursor species to product species
c
      do l = 1,nsmrxn
        isp = idsmpre(l)
        smkrat = smskrat(l)
        smjrat = smsjrat(l)
        if (lsnow) then
          smkrat = smikrat(l)
          smjrat = smijrat(l)
        endif
        effrates = smkrat + smjrat*fshads*fzenith*ctrns
        effratev = smvkrat(l) + smvjrat(l)*fshadv*fzenith*ctrns
        delcons = solmas(isp)*(1. - EXP(-effrates*dtime*60.))
        delconv = vegmas(isp)*(1. - EXP(-effratev*dtime*60.))
        smas = MAX(eps,solmas(isp) - delcons)
        vmas = MAX(eps,vegmas(isp) - delconv)
        delcons = solmas(isp) - smas
        delconv = vegmas(isp) - vmas
        solmas(isp) = smas
        vegmas(isp) = vmas
        if (idsmprd(l).ne.0) then
          isp = idsmprd(l)
          solmas(isp) = solmas(isp) + delcons
          vegmas(isp) = vegmas(isp) + delconv
        endif
      enddo
c
c --- Apply decay rates for leaching, snow melt, and plant penetration
c     Accelerate if raining (1 mm/hr = e-folding in 1 hr)
c     Accelerate if snowing (1 cm/hr = e-folding in 1 day)
c
      rain = (pwr/rhoh2o/1.0e-7)**1.27       !rain rate (mm/hr)
      if (rain.lt.0.2) rain = 0.
      rain = rain/60.                        !mm/hr -> 1/min
      snow = 0.
      if (snowrat.gt.0.) 
     &  snow = 100.*10.*snowrat/(24.*60.)    !m/hr SWE -> cm/hr depth -> 1/min
      do l = 1,nspc
        slch = smlch(l)
        if (lsnow) then
          slch = 0.
          if (temp0.gt.272.) slch = smmlt(l)
        endif
        delcons = EXP(-(slch + rain + snow)*dtime*60.)
        delconv = EXP(-(smpen(l) + rain + snow)*dtime*60.)
        solmas(l) = MAX(eps,solmas(l)*delcons)
        vegmas(l) = MAX(eps,vegmas(l)*delconv)
      enddo
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
      return
      end
