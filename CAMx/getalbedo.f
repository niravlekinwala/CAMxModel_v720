      subroutine getalbedo(igrd,ncol,nrow,snow,snowage,tsurf,fsurf,
     &                     albedo,snowfrc,snowalb)
c
c----CAMx v7.20 220430
c 
c     GETALB calcualtes surface albedo for the specified grid, according
c     to landuse and snow cover.
c 
c      Copyright 1996 - 2022
c     Ramboll
c           
c     Modifications:
c        11/4/16        Reductions to critical SWE values for some LU
c        04/23/21       Updated water and bareen UV albedo
c
c     Input arguments:
c        igrd                grid index
c        ncol                number of columns
c        nrow                number of rows
c        snow                snow cover water equivalent (m)
c        snowage             snow cover age (hr)
c        tsurf               surface temperature (K)
c        fsurf               fractional landuse cover field (fraction)
c             
c     Output arguments: 
c        albedo              surface net albedo (fraction)
c        snowfrc             net snow cover (fraction)
c        snowalb             surface snow albedo (fraction)
c             
c     Routines called: 
c        None
c             
c     Called by: 
c        EMISTRNS
c 
      implicit none
      include 'camx.prm'
      include 'camx.inc'
      include 'deposit.inc'
c
      integer :: igrd,ncol,nrow
      real, dimension(ncol,nrow) :: snow,snowage
      real, dimension(ncol,nrow) :: tsurf
      real, dimension(ncol,nrow) :: albedo
      real, dimension(ncol,nrow) :: snowfrc
      real, dimension(ncol,nrow) :: snowalb
      real, dimension(ncol,nrow,nlu) :: fsurf
c
      integer i,j,m
      real fsnow,age,c1,c2,albsno,snowc,delalb,old
      real albmax,albc1(2),albc2(2),albtot
      real snowcw(NLUW89),snowcz(NLUZ03)
      real albcw(NLUW89),albcz(NLUZ03)
c
      data albmax /0.85/
      data albc1 /0.94,0.82/
      data albc2 /0.58,0.46/
      data snowcw /.05,.01,.01,.20,.20,.20,.01,.01,.01,.01,.01/
      data snowcz /.01,.01,.01,.20,.20,.20,.20,.20,.20,.02,
     &             .02,.02,.01,.02,.01,.01,.01,.01,.01,.01,
     &             .05,.01,.01,.01,.20,.20/
      data albcw /.08,.05,.05,.05,.05,.05,.07,.10,.05,.05,.05/
      data albcz /.07,.50,.07,.05,.05,.05,.05,.05,.05,.05,
     &            .05,.05,.05,.05,.05,.05,.05,.05,.05,.05,
     &            .08,.05,.05,.10,.05,.05/
c
c-----Entry point
c
c-----Loop over rows and columns
c
      do 30 j = 2,nrow-1 
        do 20 i = 2,ncol-1
          snowfrc(i,j) = 0.
          albedo(i,j) = 0.
c
c-----Determine snow albedo based on snow age
c
          age = snowage(i,j)/24.  ! Convert from hr to days
          old = max(0.,age - dtinp/60./24.)
          c1 = albc1(1)
          c2 = albc2(1)
          if (tsurf(i,j).gt.272.) then
            c1 = albc1(2)
            c2 = albc2(2)
          endif
          if (snowalb(i,j).eq.0.) then
            albsno = albmax*(c1**(age**c2))
          else
            delalb = albmax*(c1**(age**c2) - c1**(old**c2))
            albsno = snowalb(i,j) + delalb
          endif
          albsno = max(0.4,albsno)
          snowalb(i,j) = albsno
c
c-----Determine snow cover fraction and snow albedo for each landuse
c
          do m = 1,nlu
            if (nlu .EQ. NLUW89) then
              snowc = snowcw(m)
            else
              snowc = snowcz(m)
            endif
            fsnow = 1. - exp(-2.6*snow(i,j)/snowc) + 
     &                   (snow(i,j)/snowc)*exp(-2.6)
            fsnow = max(0.,min(1.,fsnow))
c
c----Blend with default terrestrial (non-snow) albedo
c
            if (nlu .EQ. NLUW89) then
              if (m.eq.7) fsnow = 0.
              albtot = fsnow*albsno + (1.-fsnow)*albcw(m)
            else
              if (m.eq.1 .or. m.eq.3) fsnow = 0.
              albtot = fsnow*albsno + (1.-fsnow)*albcz(m)
            endif
            snowfrc(i,j) = snowfrc(i,j) + fsnow*fsurf(i,j,m)
            albedo(i,j) = albedo(i,j) + albtot*fsurf(i,j,m)
          enddo
c
 20     continue
 30   continue
c
      return
      end
