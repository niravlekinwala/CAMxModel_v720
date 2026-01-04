      subroutine vadvppm(nn,igrd,dtin,dz,con,vel,dil,sens,fluxlay,
     &                   fc1,fc2,fc3,ldoipts)
      use tracer
c  
c----CAMx v7.20 220430
c
c     VADVPPM performs concentration adjustments due to vertical advection
c     and changes in layer structure that result from the time- and
c     space-varying vertical coordinate system.  This version performs works
c     on both concentrations and DDM sensitivities.
c
c     It uses the one-dimensional implementation of the piecewise parabolic
c     method of Colella and Woodward (1984).
c     A piecewise continuous parabola is used as the intepolation polynomial.
c     The slope of the parabola at cell edges is computed from a cumulative
c     function of the advected quantity.  These slopes are further modified
c     so that the interpolation function is monotone.
c
c     Upper boundary conditions:
c        This routine accomodates non-zero vertical velocities at the top of
c        the column; the top boundary condition is held in the conc vector
c        at location (N+1).
c     Lower boundary conditions:
c        Zero flux is specified at the ground
c
c     This version adapted from CAMx/HADVPPM.F for variable vertical grid
c     spacing following CMAQ/VPPM.F
c
c     The following definitions are used:  
c  
c              |-----------> Positive direction (upward)
c  
c     |<--------------------------Domain----------------->|Boundary
c  
c     | CON(1) | CON(2) |  ...  | CON(I) |  ...  | CON(N) | CON(N+1) 
c  
c     VEL(1)-->|     VEL(I-1)-->|        |-->VEL(I)       |-->VEL(N)  
c  
c      FP(1)-->|      FP(I-1)-->|        |-->FP(I)        |-->FP(N)  
c  
c      FM(2)<--|        FM(I)<--|        |<--FM(I+1)      |<--FM(N+1)  
c             
c                            -->|  DZ(I) |<-- 
c 
c      Copyright 1996 - 2022
c     Ramboll
c  
c     Modifications:  
c        4/02/20   Original development
c        4/28/21   Added DDM sensitivities
c  
c     Input arguments:  
c        nn                  Number of cells  
c        igrd                index of current grid
c        dtin                Time step (s)  
c        dz                  Layer depth vector (m)
c        con                 Concentration vector (umol/m3) 
c        vel                 Wind speed+entrainment vector (m/s)
c                                              -- positive = upward 
c        dil                 Layer depth dilution rate (m/s)
c                                              -- positive = increasing/diluting
c        sens                sensitivity coefficients(conc/parameter unit)
c        ldoipts             flag to calculate Process Analysis data
c  
c     Output arguments:  
c        con                 Concentration vector (umol/m3)  
c        sens                sensitivity coefficients(conc/parameter unit)
c        fluxlay             Flux across top of each layer (umol/m2/s)
c        fc1                 Conc entrained through bottom of layer (umol/m3)
c        fc2                 Conc entrained through top of layer (umol/m3)
c        fc3                 Conc diluted by layer expansion/contraction
c                            (umol/m3)
c  
c     Routines called:  
c        PPM
c  
c     Called by:  
c        ZADVEC  
c 
      implicit none
      include "camx.prm"
c
      integer nn
      integer igrd
      real dtin
      real dz(nn),con(nn+1),vel(nn),dil(nn),fluxlay(nn)
c
      integer i,n,nstep,isp
      integer nsp
      real x,y,cfl,dtmx,tmpz,vwind,dt
      real con1(MXLAYER+1)
      real con2(MXLAYER)
      real flxarr(MXLAYER)
      real fm(MXLAYER+1)
      real fp(MXLAYER)
cae   real cm(MXLAYER)
      real cl(MXLAYER)
      real cr(MXLAYER)
      real dc(MXLAYER)
      real c6(MXLAYER)
c
c-----Local parameters
c
      real, parameter :: TWO3RDS=2./3.
c
c===========================DDM Begin=================================
c
      real sens(nn+1,nddmsp)
c
c     Note:  nddmsp is the total number of DDM parameters for which
c            sensitivities are tracked for the given species.
c
c===========================DDM End===================================
c
c
c================== Process Analysis Begin ===========================
c
      logical ldoipts
      real fc1(*)
      real fc2(*)
      real fc3(*)
c
c=================== Process Analysis End ============================
c
c
c-----Entry point 
c
      nsp = 1
      if( (lddm .OR. lhddm) .AND. lddmcalc(igrd) ) nsp = nddmsp + 1
      do i = 1,nn
        fluxlay(i) = 0.
        fc1(i) = 0.
        fc2(i) = 0.
        fc3(i) = 0.
      enddo
c
c-----Time split based on CFL
c
      cfl = 0.9
      dtmx = dtin
      tmpz = cfl*dz(1)/(abs(vel(1)) + 1.e-10)
      dtmx = min(dtmx,tmpz)
      do i = 2,nn
        if (vel(i).gt.0. .and. vel(i-1).lt.0.) then
          vwind = abs(vel(i)) + abs(vel(i-1))
        else
          vwind = max(abs(vel(i)),abs(vel(i-1)))
        endif
        tmpz = cfl*dz(i)/(vwind + 1.e-10)
        dtmx = min(dtmx,tmpz)
      enddo
      nstep = INT( 0.999*dtin/dtmx ) + 1
      dt = dtin/FLOAT(nstep)
c
c-----Loop over timesteps
c
      do n = 1,nstep
c
c-----Loop over this species and all its DDM sensitivity paraemters
c
      do isp = 1,nsp

      if (isp.eq.1) then
        do i = 1,nn+1
          con1(i) = con(i)
        enddo
      else
c
c============================DDM Begin================================
c
c-----If sensitivities are requested, load additional arrays into con,
c     one array for each parameter sensitivity coefficient.
c     Both the concentration and sensitivities are loaded into the same
c     array and the equations are solved together.
c
c     Load all layers 1 through nn, and the upper boundary condition.
c
        if( (lddm .OR. lhddm) .AND. lddmcalc(igrd) ) then
          do i = 1,nn+1
            con1(i) = sens(i,isp-1)
          enddo
        endif
c
c=============================DDM End=================================
c
      endif
c
c-----Set all fluxes to zero. Either positive or negative flux will
c     remain zero depending on the sign of the velocity
c
      do i = 1,nn
        fm(i+1) = 0.
        fp(i) = 0.
        flxarr(i) = 0.
      enddo
c
      call ppm(nn,dz,con1,cr,cl,dc,c6)
c
c-----Zero order polynomial at the boundary cells
c
cae   cm(2)  = con(2)
cae   cm(nn) = con(nn-1)
c
c-----First order polynomial at the next cells, no monotonicity constraint
c     needed
c
cae   cm(3)    = (con(3) + con(2))/2.
cae   cm(nn-1) = (con(nn-1) + con(nn-2))/2.
c
c-----Second order polynomial inside the domain
c
cae   do i = 3,nn-2
c
c-----Compute average slope in the i'th cell
c
cae     dc(i) = 0.5*(con(i+1) - con(i-1))
c      
c-----Guarantee that CM lies between CON(I) and CON(I+1)
c     monotonicity constraint
c     
cae     if ((con(i+1) - con(i))*(con(i) - con(i-1)).gt.0.) then
cae       dc(i) = sign(1.,dc(i))*min(
cae  &                               abs(dc(i)),
cae  &                               2.*abs(con(i+1) - con(i)),
cae  &                               2.*abs(con(i) - con(i-1)))
cae     else
cae       dc(i) = 0.
cae     endif
cae   enddo
c
cae   do i = 3,nn-3
cae     cm(i+1) = con(i) + 
cae  &            0.5*(con(i+1) - con(i)) + (dc(i) - dc(i+1))/6.
cae   enddo
c
cae   do i = 2,nn-1
cae     cr(i) = cm(i+1)
cae     cl(i) = cm(i)
cae   enddo
c
c-----Generate piecewise parabolic distributions
c
cae   do i = 2,nn-1
c
c-----Monotonicity
c 
cae     if ((cr(i) - con(i))*(con(i) - cl(i)).gt.0.) then
cae       dc(i) = cr(i) - cl(i)
cae       c6(i) = 6.*(con(i) - 0.5*(cl(i) + cr(i)))
c
c-----Overshoot cases
c
cae       if (dc(i)*c6(i) .gt. dc(i)*dc(i)) then
cae         cl(i) = 3.*con(i) - 2.*cr(i)
cae       elseif (-dc(i)*dc(i) .gt. dc(i)*c6(i)) then
cae         cr(i) = 3.*con(i) - 2.*cl(i)
cae       endif
cae     else
cae       cl(i) = con(i)
cae       cr(i) = con(i)
cae     endif
cae     dc(i) = cr(i) - cl(i)
cae     c6(i) = 6.*(con(i) - 0.5*(cl(i) + cr(i)))
cae   enddo
c
c-----Compute fluxes from the parabolic distribution
c
      do i = 2,nn
        y = max(0., -vel(i-1)*dt)
        x = y/dz(i)
        fm(i) = y*(cl(i) + 0.5*x*(dc(i) + c6(i)*(1. - TWO3RDS*x)))
        y = max(0., vel(i)*dt)
        x = y/dz(i)
        fp(i) = y*(cr(i) - 0.5*x*(dc(i) - c6(i)*(1. - TWO3RDS*x)))
      enddo
c
c-----Compute fluxes from boundary cells assuming uniform distribution
c
      if (vel(1).gt.0.) then
        y = vel(1)*dt
        fp(1) = y*con1(1)
      endif
c
      if (vel(nn).lt.0.) then
        y = -vel(nn)*dt
        fm(nn+1) = y*con1(nn+1)
      endif
c
c-----Update concentrations due to z-wind & cross-layer entrainment
c
      flxarr(1) = (fp(1) - fm(2))
      con2(1) = con1(1) - flxarr(1)/dz(1)
      do i = 2,nn
        flxarr(i) = (fp(i) - fm(i+1))
        con2(i) = con1(i) - (flxarr(i) - flxarr(i-1))/dz(i)
      enddo
c
c-----Account for layer dilution rate over this step, accumulate fluxes
c
      do i = 1,nn
        con2(i) = con2(i)*dz(i)/(dz(i) + dil(i)*dt)
      enddo
c
      if (isp.eq.1) then
        do i = 1,nn
          con(i) = con2(i)
          fluxlay(i) = fluxlay(i) + flxarr(i)/dtin
        enddo
c
c===================== Process Analysis Begin =========================
c
        if (ldoipts) then
          do i = 1,nn
            fc3(i) = fc3(i) - dt/dz(i)*dil(i)*con1(i)
          enddo
          fc1(1) = 0.
          do i = 2,nn
            fc1(i) = fc1(i) + flxarr(i-1)/dz(i)
          enddo
          do i = 1,nn
            fc2(i) = fc2(i) - flxarr(i)/dz(i)
          enddo
        endif
c
c====================== Process Analysis End ==========================
c
      else
c
c==============================DDM Begin===============================
c
        do i = 1,nn
          sens(i,isp-1) = con2(i)
        enddo
c
c===============================DDM End================================
c
      endif
c
      enddo ! End of species + DDM sensitivity loop

      enddo  !End timestep loop
c
      return
      end
c
c-------------------------------------------------------------------------------
c
      subroutine ppm(ni,ds,cn,cr,cl,dc,c6)

c-----Get PPM coefficients CR, CL, DC, and C6

      implicit none
      include "camx.prm"
c
c-----Arguments
c
      integer :: ni       ! number of layers
      real    :: ds(ni)   ! layer thickness
      real    :: cn(ni+1) ! column concentrations
      real    :: cr(ni)   ! layer r.h. (top) intercept
      real    :: cl(ni)   ! layer l.h. (bottom) intercept
      real    :: dc(ni)   ! cr - cl
      real    :: c6(ni)   ! coefficient of second-order term
c
c-----Local variables:
c
      integer i
      real a, b, c             ! temp lattice vars.
      real :: alpha (MXLAYER)  ! temp lattice var.
      real :: beta
      real :: chi   (MXLAYER)  ! lattice var. for dc
      real :: psi   (MXLAYER)  ! lattice var. for dc
      real :: mu    (MXLAYER)  ! lattice var. for cm
      real :: nu    (MXLAYER)  ! lattice var. for cm
      real :: lambda(MXLAYER)  ! lattice var. for cm
      real :: cm    (MXLAYER+1)! zone r.h. trial intercept
c
c-----Entry point
c
      do i = 2,ni-1
cae   do i = 3,ni-2
        alpha(i) = ds(i) + ds(i+1)
        beta     = ds(i-1) + ds(i)
        c        = ds(i)/(beta + ds(i+1))
        chi(i)   = c*(ds(i-1) + beta)/alpha(i)
        psi(i)   = c*(alpha(i) + ds(i+1))/beta
      enddo
      do i = 2,ni-2
cae   do i = 3,ni-3
        a     = ds(i)/alpha(i)
        b     = 2.0*ds(i+1)/alpha(i)
        c     = 1.0/(ds(i-1) + alpha(i) + ds(i+2))
        mu(i) = c*ds(i)*(ds(i-1) + ds(i))/(ds(i) + alpha(i))
        nu(i) = c*ds(i+1)*(ds(i+1) + ds(i+2))/(ds(i+1) + alpha(i))
        lambda(i) = a + mu(i)*b - 2.0*nu(i)*a
      enddo
c
c-----Zeroth order polynomial at the boundary cells
c
      cm(1)    = cn(1)
cae   cm(ni+1) = cn(ni)
      cm(ni+1) = cn(ni+1)
cae   cm(2)    = cn(2)
cae   cm(ni)   = cn(ni-1)
c
c-----First order polynomial at the next cells, no monotonicity constraint needed
c
      cm(2)  = (ds(1)*cn(2) + ds(2)*cn(1))/(ds(1) + ds(2))
      cm(ni) = (ds(ni-1)*cn(ni) + ds(ni)*cn(ni-1))/(ds(ni-1) + ds(ni))
cae   cm(3)  = (ds(2)*cn(3) + ds(3)*cn(2))/(ds(2) + ds(3))
cae   cm(ni-1) = (ds(ni-2)*cn(ni-1) + ds(ni-1)*cn(ni-2))/(ds(ni-2) + ds(ni-1))
c
c-----Second order polynomial inside the domain
c
      do i = 2,ni-1
cae   do i = 3,ni-2

c-----Compute average slope in i-th layer

        dc(i) = chi(i)*(cn(i+1) - cn(i)) + psi(i)*(cn(i) - cn(i-1)) ! equation (1.7)
c
c-----Guarantee that CM lies between CON(i) and CON(i+1) - monotonicity constraint
c
        if ((cn(i+1) - cn(i))*(cn(i) - cn(i-1)) .gt. 0.0) then
          dc(i) = sign(1.0,dc(i))* min(
     &                                 abs(dc(i)),
     &                                 2.0*abs(cn(i+1) - cn(i)),
     &                                 2.0*abs(cn(i) - cn(i-1)))
        else
          dc( i ) = 0.0
        endif           ! equation (1.8)
      enddo
c
      do i = 2,ni-2     ! equation (1.6)
cae   do i = 3,ni-3     ! equation (1.6)
        cm(i+1) = cn(i) + lambda(i)*(cn(i+1) - cn(i))
     &                  - mu(i)*dc(i+1) + nu(i)*dc(i)
      enddo
c
c-----Generate piecewise parabolic distributions
c
      do i = 1,ni
cae   do i = 2,ni-1
        cr(i) = cm(i+1) ! equation (1.15)
        cl(i) = cm(i)
c
c-----Monotonicity
c
        if ((cr(i) - cn(i))*(cn(i) - cl(i)) .gt. 0.0) then
          dc(i) = cr(i) - cl(i)  ! temporary computation of dc and c6
          c6(i) = 6.0*(cn(i) - 0.5*(cl(i) + cr(i)))
c
c-----Overshoot cases
c
          if (dc(i)*c6(i) .gt. dc(i)*dc(i)) then
            cl(i) = 3.0*cn(i) - 2.0*cr(i)
          elseif (-dc(i)*dc(i) .gt. dc(i)*c6(i)) then
            cr(i) = 3.0*cn(i) - 2.0*cl(i)
          endif
        else    ! local extremum: interpolation function is set to be a constant
          cl(i) = cn(i)
          cr(i) = cn(i)
        endif
        dc(i) = cr(i) - cl(i)    ! equation (1.5)
        c6(i) = 6.0*(cn(i) - 0.5*(cl(i) + cr(i)))
      enddo

      return
      end
