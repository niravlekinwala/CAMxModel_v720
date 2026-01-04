c----CAMx v7.20 220430
c  
c     DEPOSIT.COM contains arrays for dry deposition 
c                            
c      Copyright 1996 - 2022
c     Ramboll
c            
c     Modifications:  
c        7/21/03     removed istress as it can now be optionally specified
c                    via input file
c        5/15/13     Added surface model parameters (from RTRACSRF.INC)
c        8/22/14     Extended Zhang LU to surface model
c 
c-----------------------------------------------------------------------
c     Arrays defining the dry deposition model:
c
c     nlu     -- number of landuse categories for choses dry dep model
c     z0lu    -- surface roughness length by landuse and season (meter)
c     lai_ref -- reference leaf area index (LAI) by landuse/season
c     iseason -- season index map by latitude band and month
c     rj      -- baseline minimum stomatal resistance (s/m) 
c     rlu     -- baseline upper canopy (cuticle) resistance (s/m) 
c     rac     -- baseline canopy height/density resistance (s/m) 
c     rlcs    -- baseline SO2 lower canopy resistance (s/m) 
c     rlco    -- baseline O3 lower canopy resistance (s/m) 
c     rgss    -- baseline SO2 ground surface resistance (s/m) 
c     rgso    -- baseline O3 ground surface resistance (s/m)
c     dstress -- drought stress adjustment factors to stomatal resistance
c     lai_ref -- reference LAI for Zhang LU by month (and max/min)
c     z01     -- minimum surface roughness for Zhang LU
c     z02     -- maximum surface roughness for Zhang LU
c     fsoil   -- surface model soil fraction as f(Wesely LU)
c     fsoilz  -- surface model soil fraction as f(Zhang LU)
c     fsoiloc -- fraction soil organic matter f(Wesely LU)
c                (RTRAC surface model only)
c-----------------------------------------------------------------------
c
      integer nlu
      integer iseason(5,12)
      real    z0lu(NLUW89,5)
      real    rj(NLUW89,5)
      real    rlu(NLUW89,5)
      real    rac(NLUW89,5)
      real    rlcs(NLUW89,5)
      real    rlco(NLUW89,5)
      real    rgss(NLUW89,5)
      real    rgso(NLUW89,5)
      real    dstress(0:5)
      real    lai_ref(NLUZ03,15)
      real    z01(NLUZ03)
      real    z02(NLUZ03)
      real    fsoil(NLUW89)
      real    fsoilz(NLUZ03)
      real    fsoiloc(NLUW89)
c
      common /deposit/ nlu, z0lu, rj, rlu, rac, rlcs, rlco, rgss, rgso,
     &                 iseason, dstress, lai_ref, z01, z02, fsoil, 
     &                 fsoilz, fsoiloc
