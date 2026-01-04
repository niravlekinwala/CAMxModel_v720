      subroutine ironwet(nspc,n,dx,dy,mapscl,cwc,pwr,pws,pwg,cph,
     &                   temp,press,puffvol,conpig,conbck,depth,dt,
     &                   depfld,flxwet)
      use camxcom
      use chmstry
      use pigsty
      implicit none
c
c----CAMx v7.20 220430
c
c     IRONWET calculates wet deposition for a given IRON PiG puff.
c     Material is depleted (+ puff increment) or added (- puff increment)
c     over the entire volume of the puff using a simplified scavenging
c     mechanism.
c
c      Copyright 1996 - 2022
c     Ramboll
c
c     Modifications:
c        8/31/06         Added map scale factor
c        5/12/08         Revised deposited liquid concentration to use
c                        correct rainfall volume
c       12/22/09         Corrected equilibrium gas concentration in rain
c        7/26/13         Improved for new GREASD chemistry
c       12/07/14         Modified for VBS
c       06/06/16         Scavenging modifications following WETDEP mods
c       08/16/21         Fixed calculation errors in wet deposition
c                        liquid concentration
c
c     Input arguments:
c        nspc                number of species
c        n                   puff index
c        dx                  grid cell size in x-direction (m)
c        dy                  grid cell size in y-direction (m)
c        mapscl              map scale factor
c        cwc                 cloud water content (g/m3)
c        pwr                 rain water content (g/m3)
c        pws                 snow water content (g/m3)
c        pwg                 graupel water content (g/m3)
c        cph                 cloud water pH
c        temp                temperature (K)
c        press               pressure (mb)
c        puffvol             puff volume (m3)
c        conpig              puff increment (umol/m3,ug/m3)
c        conbck              background concentration (umol/m3,ug/m3)
c        depth               puff depth (m)
c        dt                  time step (s)
c        depfld              wet deposited mass (mol/ha, g/ha)
c                            and surface liquid concentrations (mol/l, g/l)
c        flxwet              deposition mass flux (umol, ug)
c  
c     Output arguments:
c        depfld              wet deposited mass (mol/ha, g/ha)
c                            and surface liquid concentrations (mol/l, g/l)
c        flxwet              deposition mass flux (umol, ug)
c
c     Routines Called:
c        HENRYFNC
c        SCAVRAT
c
c     Called by:
c        PIGDRIVE
c
      include 'camx.prm'
      include 'flags.inc'
c
      integer nspc
      integer n
      real cwc,pwr,pws,pwg,cph,temp,press,puffvol,depth,
     &     dt,dx,dy,mapscl
      real depfld(*)
      real conpig(MXSPEC+1,MXRECTR)
      real conbck(MXSPEC+1)
      real*8 flxwet(*)
c
      integer nr,l,isemptyf,isemptyc
      real rd,rhoh2o,rr,volrat,rhoair,convfac,cmin,conc,
     &     hlaw,dscav,gscav,ascav,delm,qtf,qtc,vtf,vtc,psize,psizec,
     &     ascavf,ascavc,roprta,cph2o,pwc,rconst,cwat,pwat,pwtr,delc
      real area, volume
      real drpvel,ceqc,ceqp
      real fdep(MXSPEC)
      logical lcloud,lgraupl
c
      data rd /287./         ! Dry air gas constant (J/K/kg)
      data rhoh2o /1.e6/     ! water density (g/m3)
      data rconst /8.206e-2/ ! gas constant (l.atm/mol.K)
c
c-----Entry point
c
      do l = 1,nspec
        depfld(l) = 0.
        depfld(nspec+l) = 0.
        flxwet(l) = 0.
      enddo
c
      pwc = pwr + pws + pwg
      volrat = pwc/rhoh2o                 ! drop volume/air volume
      rr     = (volrat/1.0e-7)**1.27      ! rainfall rate (mm/hr)
c
      lcloud = .false.
      lgraupl = .false.
      if (cwc .ge. cwmin) lcloud = .true.
      if (pwg .ge. cwmin) lgraupl = .true.
      rhoair = 100.*press/(rd*temp)
      cwat = cwc
      pwtr = pwc
      if (temp.lt.273.)
     &  cwat = max(0.,cwc*(temp - tamin)/(273. - tamin))
c
c-----Loop over puff reactors
c
      do 10 nr = 1,nreactr
c
c-----Calculate scavenging for soluble gas species
c
        do 40 l = nrad+1,ngas
          if( henry0(l) .LT. 1.e-6 ) goto 40
          convfac = densfac*(273./temp)*(press/1013.)
          cmin = bdnl(l)*convfac
          conc = MAX(cmin,conpig(l,nr) + conbck(l))
          call henryfnc(l,henry0(l),tfact(l),temp,cph,knh3,khno3,
     &                  kso2,hlaw)
          hlaw = hlaw*rconst*temp
          ceqc = conc*(1. - 1./(1. + hlaw*cwat/rhoh2o))
          call henryfnc(l,henry0(l),tfact(l),temp,5.,knh3,khno3,
     &                  kso2,hlaw)
          hlaw = hlaw*rconst*temp
          pwat = pwtr
          if (temp.lt.273. .and. rscale(l).gt.0.) pwat = 0.
          ceqp = (conc - ceqc)*(1. - 1./(1. + hlaw*pwat/rhoh2o))
c
          call scavrat(.false.,lcloud,lgraupl,tamin,rr,temp,rhoair,
     &                 diffrat(l),0.,0.,dscav,gscav,ascav,drpvel)
c
          delc = ceqc*(1. - exp(-dscav*dt)) +
     &           ceqp*(1. - exp(-gscav*dt))
          delc = min(delc,conc-cmin)
          delc = MAX(0.,delc*MAX(0.,conpig(l,nr))/conc)
          delm = delc*puffvol/float(nreactr)
          puffmass(l,nr,n) = puffmass(l,nr,n) - delm
c
          flxwet(l) = flxwet(l) - DBLE(delm)
          area = dx*dy/mapscl**2
          volume = 1.e-3*rr*area
          depfld(l) = depfld(l) + 1.e-2*delm/area
          depfld(nspec+l) = depfld(nspec+l) +
     &                      1.e-9*(delm/volume)
 40     continue
c
c-----Calculate scavenging for particulate species
c
        if (naero .gt. 0) then
c
c-----Recalculate particle size for wet diameter (CF scheme)
c
          if (lchem .AND. aeropt.eq.'CF') then
            isemptyf = 1
            isemptyc = 1
            qtf = 0. ! fine dry total mass
            qtc = 0. ! coarse dry total mass
            vtf = 0. ! fine dry total volume
            vtc = 0. ! coarse dry total volume
            do l = ngas+1,nspec
              conc = MAX(bdnl(l),conpig(l,nr) + conbck(l))
              if (dcut(l,2) .lt. (dcut(kph2o,2)+1.e-5)) then   ! fine
                if (l .ne. kph2o) then
                  if (conc.gt.bdnl(l)) isemptyf = 0
                  qtf = qtf + conc
                  vtf = vtf + conc/roprt(l)
                endif
              else                                             ! coarse
                if (conc.gt.bdnl(l)) isemptyc = 0
                qtc = qtc + conc
                vtc = vtc + conc/roprt(l)
                psizec = sqrt(dcut(l,1)*dcut(l,2))
              endif
            enddo
c
            ascavf = 0.
            if (isemptyf.eq.0) then
              cph2o = amax1(bdnl(kph2o),conpig(kph2o,nr))
              roprta = (qtf + cph2o)/(vtf + cph2o/roprt(kph2o))
              psize = sqrt(dcut(kph2o,1)*dcut(kph2o,2))
              psize = 1.e-6*psize*
     &                (1. + cph2o/roprt(kph2o)/vtf)**0.33333
              if (psize.lt.2.5e-6) psize = psize*10.
              call scavrat(.true.,lcloud,lgraupl,tamin,rr,temp,
     &                 rhoair,0.,psize,roprta,dscav,gscav,ascavf,drpvel)
            endif
c
            ascavc = 0.
            if (isemptyc.eq.0) then
              roprta = qtc/vtc
              psize = psizec*1.e-6
              call scavrat(.true.,lcloud,lgraupl,tamin,rr,temp,
     &                 rhoair,0.,psize,roprta,dscav,gscav,ascavc,drpvel)
            endif
c
            do l = ngas+1,nspec
              if (dcut(l,2) .lt. (dcut(kph2o,2)+1.e-5)) then
                fdep(l) = 1. - exp(-ascavf*dt)
              else
                fdep(l) = 1. - exp(-ascavc*dt)
              endif
            enddo
c
c-----Other scheme
c
          else
            do l = ngas+1,nspec
               psize = 1.e-6*sqrt(dcut(l,1)*dcut(l,2))
               if (psize.lt.2.5e-6) psize = psize*10.
               call scavrat(.true.,lcloud,lgraupl,tamin,rr,temp,
     &                rhoair,0.,psize,roprt(l),dscav,gscav,ascav,drpvel)
               fdep(l) = 1. - exp(-ascav*dt)
             enddo
          endif

          do 50 l = ngas+1,nspec
            delm = MAX(0.,puffmass(l,nr,n)*fdep(l))
            puffmass(l,nr,n) = puffmass(l,nr,n) - delm
c
            flxwet(l) = flxwet(l) - DBLE(delm)
            area = dx*dy/mapscl**2
            volume = 1.e-3*rr*area
            depfld(l) = depfld(l) + 1.e-2*delm/area
            depfld(nspec+l) = depfld(nspec+l) +
     &                      1.e-9*delm/volume
 50       continue
c
        endif
c
 10   continue
c
      return
      end
