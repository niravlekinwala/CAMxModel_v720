      subroutine convtrans(igrid,dtinp,ncol,nrow,nlay,deltat,idfin,
     &                     depth,press,temp,cigfrc,cigtim,cigent,cigdet,
     &                     cigmas)
      use filunit
c
c----CAMx v7.20 220430
c
c     CONVTRANS integrates cloud-in-grid mass transport matrices based on input
c     K-F entrainment/detrainment fluxes.  Assuming a steady-state sub-grid
c     cloud environment, these CiG matrices are calculated for each convective
c     column over an entire grid for the duration of the column-specific
c     convective time scale and then used in CIGDRIVE for all timesteps between
c     the meteorological update times (usually 1 hour).
c
c      Copyright 1996 - 2022
c     Ramboll
c
c     Modifications:
c        None
c
c     Input arguments:
c        igrid               grid index
c        dtinp               met input frequency (min)
c        ncol                number of columns
c        nrow                number of rows
c        nlay                number of layers
c        deltat              time step (s)
c        idfin               map of nested grids in this grid
c        depth               layer depth field (m)
c        press               pressure field (mb)
c        temp                temperature field (K)
c        cigfrc              CiG cloud fraction field
c        cigtim              CiG cloud timescale (s)
c        cigent              CiG sub-grid cloud entrainment rate (kg/m2/s)
c        cigdet              CiG sub-grid cloud detrainment rate (kg/m2/s)
c
c     Output arguments:
c        cigmas              CiG transport matrix (unitless)
c
c     Routines called:
c        None
c
c     Called by:
c        TSTEP_INIT
c
      implicit none
      include "camx.prm"

      integer igrid,ncol,nrow,nlay
      integer idfin(ncol,nrow)
      real deltat,dtinp
      real depth(ncol,nrow,nlay)
      real press(ncol,nrow,nlay)
      real temp(ncol,nrow,nlay)
      real cigfrc(ncol,nrow)
      real cigtim(ncol,nrow)
      real cigent(ncol,nrow,nlay)
      real cigdet(ncol,nrow,nlay)
      real cigmas(2,nlay,nlay,ncol,nrow)

      integer i,j,k,kk,l,n,ntimes,ks,ke
      real timescale,dt,dcflx,daflx
      real fm,fp,crat
      real tmassi,tmassf,terr
      real dz(MXLAYER)
      real rho(MXLAYER)
      real entc(MXLAYER),detc(MXLAYER)
      real enta(MXLAYER),deta(MXLAYER)
      real cflx(MXLAYER),aflx(MXLAYER)
      real aatrac(MXLAYER,MXLAYER),actrac(MXLAYER,MXLAYER)
      real cctrac(MXLAYER,MXLAYER),catrac(MXLAYER,MXLAYER)
      real aanew(MXLAYER,MXLAYER),acnew(MXLAYER,MXLAYER)
      real ccnew(MXLAYER,MXLAYER),canew(MXLAYER,MXLAYER)
c
c-----Entry point
c
      do j = 1,nrow
        do i = 1,ncol
          do l = 1,nlay
            do k = 1,nlay
              do n = 1,2
                cigmas(n,k,l,i,j) = 0.
              enddo
            enddo
          enddo
        enddo
      enddo
c
c-----Loop over vertical grid columns, skipping any covered by a nest or 
c     without active sub-grid convection
c
      do 20 j = 2,nrow-1
        do 10 i = 2,ncol-1
          if (idfin(i,j).gt.igrid) goto 10
          if (cigfrc(i,j).eq.0. .or. cigtim(i,j).eq.0.) goto 10
c
c-----Set hydrostatic density profiles (kg/m3) and entrainment/detrainment
c     rates (kg/m2/s)
c
          timescale = min(cigtim(i,j),60.*dtinp)
          crat  = cigfrc(i,j)/(1. - cigfrc(i,j))
          do k = 1,nlay
            dz(k) = depth(i,j,k)
            rho(k) = 100.*press(i,j,k)/(287.*temp(i,j,k))
            entc(k) = cigent(i,j,k)
            detc(k) = cigdet(i,j,k)
            enta(k) = entc(k)*crat
            deta(k) = detc(k)*crat
          enddo
c
c-----Calculate in-cloud and ambient net vertical fluxes (kg/m2/s)
c
          cflx(1) = entc(1) - detc(1)
          aflx(1) = deta(1) - enta(1)
          do k = 2,nlay
            cflx(k) = cflx(k-1) + entc(k) - detc(k)
            aflx(k) = aflx(k-1) - enta(k) + deta(k)
          enddo
c
c-----Initialize tracer matrix
c
          do l = 1,nlay
            do k = 1,nlay
              cctrac(k,l) = 0.
              aatrac(k,l) = 0.
              catrac(k,l) = 0.
              actrac(k,l) = 0.
            enddo
            cctrac(l,l) = 1.
            aatrac(l,l) = 1.
          enddo
c
c-----Check time step, determine integration step for first-order upstream
c     transport technique
c
          dt = timescale
 100      continue
          do k = 1,nlay
            dcflx = cflx(k)
            daflx = aflx(k)
            if (k.gt.1) then
              dcflx = max(abs(cflx(k-1)),abs(cflx(k)))
              daflx = max(abs(aflx(k-1)),abs(aflx(k)))
            endif
            if ((dcflx + detc(k))*dt/dz(k).gt.rho(k)/2. .or.
     &          (daflx + enta(k))*dt/dz(k).gt.rho(k)/2.) then
              dt = dt/2.
              goto 100
            endif
          enddo
c
c-----Loop for the duration of the input timestep
c
          ntimes = nint(timescale/dt)
          dt = timescale/float(ntimes)
          do n = 1,ntimes
c
c-----Cloud tracer change
c
          do l = 1,nlay
            do k = 1,nlay
              if (k.eq.1) then
                dcflx = -cflx(k)*cctrac(k,l)
                if (cflx(k).lt.0.) dcflx = -cflx(k)*cctrac(k+1,l)
              else
                fm = cflx(k-1)*cctrac(k-1,l)
                if (cflx(k-1).lt.0.) fm = cflx(k-1)*cctrac(k,l)
                fp = cflx(k)*cctrac(k,l)
                if (cflx(k).lt.0.) then
                  fp = 0.
                  if (k.lt.nlay) fp = cflx(k)*cctrac(k+1,l)
                endif
                dcflx = fm - fp
              endif
              ccnew(k,l) = cctrac(k,l) + (dt/dz(k))*(1./rho(k))*
     &          (dcflx + entc(k)*catrac(k,l)*crat - detc(k)*cctrac(k,l))

              if (k.eq.1) then
                daflx = -aflx(k)*catrac(k,l)
                if (aflx(k).lt.0.) daflx = -aflx(k)*catrac(k+1,l)
              else
                fm = aflx(k-1)*catrac(k-1,l)
                if (aflx(k-1).lt.0.) fm = aflx(k-1)*catrac(k,l)
                fp = aflx(k)*catrac(k,l)
                if (aflx(k).lt.0.) then
                  fp = 0.
                  if (k.lt.nlay) fp = aflx(k)*catrac(k+1,l)
                endif
                daflx = fm - fp
              endif
              canew(k,l) = catrac(k,l) + (dt/dz(k))*(1./rho(k))*
     &          (daflx - enta(k)*catrac(k,l) + deta(k)*cctrac(k,l)/crat)
            enddo
            do k = 1,nlay
              cctrac(k,l) = ccnew(k,l)
              catrac(k,l) = canew(k,l)
            enddo
          enddo
c
c-----Ambient tracer change
c
          do l = 1,nlay
            do k = 1,nlay
              if (k.eq.1) then
                daflx = -aflx(k)*aatrac(k,l)
                if (aflx(k).lt.0.) daflx = -aflx(k)*aatrac(k+1,l)
              else
                fm = aflx(k-1)*aatrac(k-1,l)
                if (aflx(k-1).lt.0.) fm = aflx(k-1)*aatrac(k,l)
                fp = aflx(k)*aatrac(k,l)
                if (aflx(k).lt.0.) then
                  fp = 0.
                  if (k.lt.nlay) fp = aflx(k)*aatrac(k+1,l)
                endif
                daflx = fm - fp
              endif
              aanew(k,l) = aatrac(k,l) + (dt/dz(k))*(1./rho(k))*
     &         (daflx - enta(k)*aatrac(k,l) + deta(k)*actrac(k,l)/crat)

              if (k.eq.1) then
                dcflx = -cflx(k)*actrac(k,l)
                if (cflx(k).lt.0.) dcflx = -cflx(k)*actrac(k+1,l)
              else
                fm = cflx(k-1)*actrac(k-1,l)
                if (cflx(k-1).lt.0.) fm = cflx(k-1)*actrac(k,l)
                fp = cflx(k)*actrac(k,l)
                if (cflx(k).lt.0.) then
                  fp = 0.
                  if (k.lt.nlay) fp = cflx(k)*actrac(k+1,l)
                endif
                dcflx = fm - fp
              endif
              acnew(k,l) = actrac(k,l) + (dt/dz(k))*(1./rho(k))*
     &          (dcflx + entc(k)*aatrac(k,l)*crat - detc(k)*actrac(k,l))
            enddo
            do k = 1,nlay
              aatrac(k,l) = aanew(k,l)
              actrac(k,l) = acnew(k,l)
            enddo
          enddo

          enddo  ! End time loop
c
c-----Load CiG transport matrix array
c
          do ks = 1,nlay    !Starting layer
            do ke = 1,nlay  !Ending layer
              cigmas(1,ke,ks,i,j) = cctrac(ke,ks) + actrac(ke,ks)/crat
              cigmas(2,ke,ks,i,j) = aatrac(ke,ks) + catrac(ke,ks)*crat
            enddo
          enddo
c
c-----Check that transport matrices conserve mass
c
          do ks = 1,nlay
            tmassi = rho(ks)*dz(ks)
            tmassf = 0.
            do ke = 1,nlay
              if (cigmas(1,ke,ks,i,j).lt.0. .or.
     &            cigmas(2,ke,ks,i,j).lt.0.) then
                write(iout,'(//,a)') 'ERROR in CONVTRANS:'
                write(iout,'(2a)') 'Negative tracer after',
     &                             ' convective transport'
                write(iout,*) 'igrid, i, j, ks, ke = ',igrid,i,j,ks,ke
                write(iout,*) 'Spec  Name   In-cloud  Ambient'
                do kk = 1,nlay
                  write(iout,'(2i3,2e10.3)')
     &                   ks,kk,cigmas(1,kk,ks,i,j),cigmas(2,kk,ks,i,j)
                enddo
                call camxerr()
              endif
              tmassf = tmassf + rho(ke)*dz(ke)*
     &                 (cigmas(1,ke,ks,i,j)*cigfrc(i,j) +
     &                  cigmas(2,ke,ks,i,j)*(1.-cigfrc(i,j)))
            enddo
            terr = abs(tmassf - tmassi)/tmassi
            if (terr.gt.0.001) then
              write(iout,'(//,a)')'ERROR in CONVTRANS:'
              write(iout,*)'CiG transport matrix not conserving mass'
              write(iout,*)'Grid,I,J,K          : ',igrid,i,j,ks
              write(iout,*)'Starting/ending mass: ',tmassi,tmassf
              call camxerr()
            endif
          enddo
c
 10     continue
 20   continue

      return
      end
