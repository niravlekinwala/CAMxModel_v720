      subroutine vdiffacmx(ns,nddm,nz,deltat,temp0,press0,z0,delz,rkv,
     &                     qv,temp,press,wind,kpbl,fconv,con)
c
c----CAMx v7.20 220430
c
c     VDIFFACMX performs vertical diffusion of concentrations using
c     the Asymmetric Convective Model v2 (ACM2/ACM1; Pleim, 2007;2014). 
c
c     ACM2 portion of code is based on: 
c         Models-3/CMAQ vdiffacm2.F,v 1.13 2012/01/19 14:37:47 
c
c      Copyright 1996 - 2022
c     Ramboll
c
c     Modifications:
c        5/27/16   Updated to parallel ACM2 speed improvements in CMAQ v5.1,
c                  remove K-theory solution and dry deposition (both done using
c                  CAMx K-theory solver), and extend to H/DDM
c
c     Input arguments:
c        ns                number of species
c        nddm              number of parameters for which sensitivities are
c                          calculated - should be zero if DDM is not enabled
c        nz                number of layers
c        deltat            timestep (s)
c        temp0             Surface temperature (K)
c        press0            Surface pressure (mb)
c        z0                Surface roughness (m)
c        delz              Layer depth (m)
c        rkv               Kv profile (m2/s)
c        qv                Humidity profile (mixing ratio) 
c        temp              Temperature profile (K)
c        press             Pressure profile (mb)
c        wind              Wind speed profile
c        con               Species concentrations (umol/m3, ug/m3) followed by
c                          sensitivities (umol/m3/parameter unit).
c
c     Output arguments:
c        kpbl              Layer index at the top of the PBL
c        fconv             Convective fraction to scale Kv
c        con               Species concentrations (umol/m3, ug/m3) followed by
c                          sensitivities (umol/m3/parameter unit).
c
c     Routines Called:
c        MICROMET
c        MATRIX
c
c     Called by:
c        DIFFUS
c
      implicit none
      include "camx.prm"
c
      integer ns,nz,nddm,kpbl
      real deltat,temp0,press0,z0,fconv
      real delz(nz),rkv(nz),qv(nz),temp(nz),press(nz),
     &     wind(nz)
      real con(nz+nz*nddm,ns)
c
      real theta,thbar
      parameter (theta = 0.5)         !Crank-Nicholson = 0.5 (semi-implicit)
      parameter (thbar = 1 - theta)
      logical lconv,lstable
      integer k,kk,nstep,n,l,m
      real gamma,vk
      real tv,pbl,zf,critk,ustar,el,psih,wstar,hoverl,dt,mbar,meddy,
     &     delc,lfac1,lfac2,rhopbl,rhosum,roplus,rz,dtacm,
     &     dfacp,dfacq
      real zz(0:MXLAYER)
      real dz(MXLAYER),dzinv(MXLAYER),zm(MXLAYER),thetav(MXLAYER),
     &     rho(MXLAYER),dzm(MXLAYER),dzminv(MXLAYER),seddy(MXLAYER)
      real mbarks(MXLAYER),mdwn(MXLAYER),mfac(MXLAYER),dfsp(MXLAYER),
     &     dfsq(MXLAYER)
      real aa(MXLAYER),bb(MXLAYER),ee(MXLAYER),dd(MXLAYER),uu(MXLAYER)
      real cc(MXLAYER,1+MXTRSP)
c
      data gamma /0.286/, vk /0.4/
c
c-----Entry point
c
c-----Determine vertical grid structure variables and thermodynamic profiles.
c
      zz(0) = 0.
      do k = 1,nz
        dz(k) = delz(k)
        zz(k) = zz(k-1) + dz(k)
        zm(k)  = (zz(k) + zz(k-1))/2.
        tv = temp(k)*(1. + 0.608*qv(k))
        thetav(k) = tv*(1000./press(k))**gamma
        rho(k) = press(k)/tv
        seddy(k) = rkv(k)
      enddo
      do k = 1,nz-1
        dzm(k) = zm(k+1) - zm(k)
      enddo
      dzm(nz) = 0.
c
c-----Get the PBL parameters
c
      kpbl = 1
      pbl = dz(1)
      do k = 2,nz-1
        if (kpbl .ne. k-1) goto 100
        zf = dz(k)/dz(k-1)
        critk = 0.03*dz(k-1)*dzm(k)*(1. + zf)/200.
        if (rkv(k-1) .gt. critk) then
          kpbl = k
          pbl = pbl + dz(k)
        endif
      enddo
 100  continue

      call micromet(temp(1),temp0,press(1),press0,dz(1)/2.,wind(1),z0,
     &              pbl,ustar,el,psih,wstar,lstable)

      lconv = .false.
      hoverl = pbl/el
      if (((thetav(1) - thetav(2)) .gt. 1.e-8 ) .and.
     &    (hoverl .lt. -0.1) .and.
     &    (kpbl .gt. 3)) lconv = .true.
      if (.not.lconv) return
c
c-----Initialize ACM2 arrays
c
      do k = 1,nz - 1
        mbarks(k) = 0.
        mdwn(k) = 0.
      enddo
      mdwn(nz) = 0.
c
c-----Calculate ACM2 non-local mixing rates
c
      do k = 1,nz
        dz(k) = dz(k)*rho(k)
        dzinv(k) = 1./dz(k)
      enddo
      do k = 1,nz-1
        roplus = (rho(k)+rho(k+1))/2.
        dzm(k) = dzm(k)*roplus
        dzminv(k) = 1./dzm(k)
        seddy(k) = seddy(k)*roplus*roplus*dzminv(k)
      enddo
      dzminv(nz) = 0.

      mbar = 0.
      rhopbl = 0.
      do k = 2,kpbl
        rhopbl = rhopbl + dz(k)
      enddo
      meddy = seddy(1)/rhopbl
      fconv = 1./(1. + ((vk/(-hoverl))**0.3333)/(0.72*vk))
      mbar = meddy*fconv
      do k = 1,kpbl-1
        mbarks(k) = mbar
        rhosum = 0.
        do kk = k,kpbl
          rhosum = rhosum + dz(kk)
        enddo
        mdwn(k) = mbar*rhosum*dzinv(k)
      enddo
      mbarks(kpbl) = mbar
      mdwn(kpbl) = mbarks(kpbl)
c
c-----Modify timestep for ACM2
c
      rz = rhopbl*dzinv(1)
      dtacm = 1./(mbar*rz)
      dt = min(0.75*dtacm,deltat)
c
c-----Determine number of steps to take
c
      nstep = int(deltat/dt + 0.99)
      dt   = deltat/real(nstep)
      dfacp = theta*dt
      dfacq = thbar*dt

      do k = 1,nz
        dfsp(k) = dfacp*dzinv(k)
        dfsq(k) = dfacq*dzinv(k)
      enddo
c
c-----Initialize matrix elements
c
      do k = 1,nz
        aa(k) = 0.
        bb(k) = 0.
        ee(k) = 0.
      enddo
      bb(1) = 1. + rhopbl*dfsp(1)*mbarks(1)
      lfac1 = dfsq(1)*rhopbl*mbarks(1)
      lfac2 = dfsq(1)*mdwn(2)*dz(2)
      do k = 2,kpbl
        aa (k) = -dfacp*mbarks(k)
        bb(k) = 1. + dfacp*mdwn(k)
        ee(k) = -dfsp(k-1)*dz(k)*mdwn(k)
        mfac(k) = dz(k+1)*dzinv(k)*mdwn(k+1)
      enddo
c
c-----Loop over species
c
      do 301 l = 1,ns
        do m = 1,1+nddm
          do k = 1,nz
            cc(k,m) = con((m-1)*nz+k,l)/rho(k)
          enddo
        enddo
c
c-----Loop over sub-steps
c
        do 300 n = 1,nstep
          do m = 1,1+nddm
            do k = 1,nz
              dd(k) = 0.
              uu(k) = 0.
            enddo
c
c-----Call MATRIX1 solver for ACM convection
c
            dd(1) = cc(1,m) - lfac1*cc(1,m) + lfac2*cc(2,m)
            do k = 2,kpbl
              delc = mbarks(k)*cc(1,m)
     &               - mdwn(k)*cc(k,m)
     &               + mfac(k)*cc(k+1,m)
              dd(k) = cc(k,m) + dfacq*delc
            enddo
            call matrix1(nz,kpbl,aa,bb,ee,dd,uu)
            do k = 1,kpbl
              cc(k,m) = uu(k)
            enddo
          enddo
300     continue

        do m = 1,1+nddm
          do k = 1,nz
            con((m-1)*nz+k,l) = cc(k,m)*rho(k)
          enddo
        enddo

301   continue

      return
      end
c
c-------------------------------------------------------------------------------
c
      SUBROUTINE MATRIX1 ( nz, KL, A, B, E, D, X )
c
c This routine taken from:
c Models-3/CMAQ matrix.F,v 1.5 2011/10/21 16:11:45
c
C-------------------------------------------------------------------------------
C Rather than solving the ACM2 banded tridiagonal matrix using LU decomposition,
C it is much faster to split the solution into the ACM1 convective solver
C followed by the tridiagonal solver.
C MATRIX1 is the ACM1 solver. When the PBL is convective, this solver is called
C followed by TRI. If not convective, only TRI is called.
C
C-- ACM1 Matrix is in this form (there is no subdiagonal:
C   B1 E2                     <- note E2 (flux from layer above), not E1
C   A2 B2 E3
C   A3    B3 E4
C   A4       B4 E5
C   A5          B5 E6
C   A6             B6
C-------------------------------------------------------------------------------

      IMPLICIT NONE

C Arguments:

      integer nz
      INTEGER KL     ! CBL layer height
      REAL A( nz )   ! matrix column one
      REAL B( nz )   ! diagnonal
      REAL E( nz )   ! superdiagonal
      REAL D( nz )   ! R.H.S.
      REAL X( nz )   ! returned solution

C Locals:

      INTEGER NLAYS
      INTEGER L
      REAL BETA
      REAL ALPHA, GAMA

      NLAYS = nz

C-- ACM1 matrix solver

      BETA = D( 1 )
      GAMA = B( 1 )
      ALPHA = 1.0

      DO L = 2, KL
         ALPHA = -ALPHA * E( L ) / B( L )
         BETA  = ALPHA * D( L ) + BETA
         GAMA = GAMA + ALPHA * A( L )
      END DO

      X( 1 )  = BETA / GAMA
      X( KL ) = ( D( KL ) - A( KL ) * X( 1 ) ) / B( KL )

C-- Back sub for Ux=y

      DO L = KL-1, 2, -1
         X( L ) = ( D( L ) - A( L ) * X( 1 ) - E( L+1 ) * X( L+1 ) )
     &            / B( L )
      END DO

      RETURN
      END
