c
c----CAMx v7.20 220430
c     
c      Copyright 1996 - 2022
c     Ramboll
c
c     Modifications:
c        09/16/16   --pk--       Added variables for SOA-forming species in aqueous-phase (GLY, MGLY & GLYD)
c                                Increased NWTMOL to include in-cloud SOA MW
c        11/22/17   --bkoo--     Increased NWTMOL to have Ca & Mg separate from CO3
c
      real nacl
      real co2
      real foa
      real mhp
      real paa
      real caco3
      real mgco3
      real a3fe
      real b2mn
      real potcl
      real gly
      real mgly
      real glyd
c
      common /radmx/ nacl, co2, foa, mhp, paa, caco3, mgco3, 
     &               a3fe, b2mn, potcl, gly, mgly, glyd

      real eps
      real e0
      real lv
      real rv
c
      common /rhcalx/ eps, e0, lv, rv
c
      integer NWTMOL
      parameter (NWTMOL = 12)
      real wtmol(NWTMOL)
      common /wtmolcom/ wtmol
c
c-----Pointers to CAMx species requird by PM routines
c
      integer kso2_c
      integer kh2o2_c
      integer kform_c
      integer khono_c
      integer ko3_c
      integer koh_c
      integer kho2_c
      integer kno3_c
      integer kno_c
      integer kno2_c
      integer kpan_c
      integer kcg1_c
      integer kcg2_c
      integer kcg3_c
      integer kcg4_c
      integer khno3_c
      integer knh3_c
      integer kh2so4_c
      integer khcl_c
      integer ksoa1_c
      integer ksoa2_c
      integer ksoa3_c
      integer ksoa4_c
      integer ksopa_c
      integer ksopb_c
      integer kcrst_c
      integer kpoa_c
      integer kpec_c
      integer kph2o_c
      integer kpcl_c
      integer kna_c
      integer kpnh4_c
      integer kpno3_c
      integer kpso4_c
      integer kn2o5_c
      integer khpo_c
      integer kfoa_c
      integer kmhp_c
      integer kpaa_c
      integer kohp_c
      integer kopa_c
      integer kgly_c
      integer kmgly_c
      integer kglyd_c
      integer ksoac_c
      integer ngas_c
      integer nspec_c
c
      common /pm_ptrs/ kso2_c,  kh2o2_c, kform_c, khono_c, ko3_c,
     &                 koh_c,   kho2_c,  kno3_c,  kno_c,   kno2_c,
     &                 kpan_c,  kcg1_c,  kcg2_c,  kcg3_c,  kcg4_c,
     &                 khno3_c, knh3_c,  kh2so4_c,khcl_c,  ksoa1_c,
     &                 ksoa2_c, ksoa3_c, ksoa4_c, ksopa_c, ksopb_c,
     &                 kcrst_c, kpoa_c,  kpec_c,  kph2o_c, kpcl_c,
     &                 kna_c,   kpnh4_c, kpno3_c, kpso4_c, kn2o5_c,
     &                 khpo_c,  kfoa_c,  kmhp_c,  kpaa_c,  kohp_c,
     &                 kopa_c,  kgly_c,  kmgly_c, kglyd_c, ksoac_c,
     &                 ngas_c,  nspec_c
