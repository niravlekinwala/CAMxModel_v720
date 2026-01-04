c----CAMx v7.20 220430
c  
c     CHMDAT.INC contains all chemistry variables 
c                            
c      Copyright 1996 - 2022
c     Ramboll 
c            
c     Modifications:  
c       4/4/00     Added aerosol deposition variables to /aerochm/
c       1/9/02     Aerosol size cut points and density now a function of
c                  species
c       8/20/02    Added minimum CWC to define presence of clouds
c      12/12/02    Expanded species list for Mechanism 4
c       1/10/03    Added array for deposition output species names
c      10/14/04    Modifications for mechanism 10
c      11/30/06    Created by taking Static variables out of chmstry.com
c      12/29/06 -bkoo-     Expanded species list for the updated SOA scheme
c      01/08/07 -bkoo-     Expanded species list for Mechanism 6 (CB05)
c                          Now ipigrxn is a scalar pointer to NO self-reaction
c      07/04/07 -bkoo-     Added pointer to hydrolysis of N2O5
c      04/28/08 -gyarwood- Added EBI solver
c      04/12/12 -cemery-   Added T and P adjustments to photolysis rates
c      07/20/12 -ou-       Added Mechanism 1 (CB06 with Iodine chemistry)
c      10/08/12 -jjung-    Added Mechanism 2 (CB6r1)
c      04/30/13 -cemery-   Added surface model
c      06/28/13 -bkoo-     Added indexes for new GREASD PiG chemistry
c      10/14/13 -bkoo-     Added pointer to hydrolysis of organic nitrate
c      10/24/13 -bkoo-     Revised Mechanism 2 (CB6r1 -> CB6r2)
c      03/18/14 -bkoo-     Added species for benzene SOA
c      04/07/14 -cemery-   Added top con file
c      06/27/14 -bkoo-     Added Mechanism 3 (CB6r2h)
c      11/11/14 -bkoo-     Added Mechanism 4 (CB6r3)
c      12/07/14 -bkoo-     Added species for VBS
c      01/08/16 -bkoo-     Updated for SAPRC99->SAPRC07TC
c                                      Revised VBS
c                                      Removal of IEH solver
c      02/29/16 -bkoo-     Added pointer to heterogeneous rxn of N2O5 + HCL
c      07/20/16 -bkoo-     Added pointer to heterogeneous hydrolysis of INTR
c      09/02/16 -bkoo-     Revised heterogeneous rxn pointers
c      11/28/16 -bkoo-     Revised VBS precursors
c      11/22/17 -bkoo-     Added PFE/PMN/PK/PCA/PMG
c      01/12/18 -bkoo-     Removed BNZA/TOLA/XYLA/ISP/TRP
c      01/06/19 -cemery-   Added PAL/PSI/PTI
c      01/09/19 -cemery-   Added DMS to CB6r4 (remains Mech 4)
c      11/17/20 -cemery-   Added definition of TUV photolysis file labels
c      01/20/21 -gwilson-  Added "blank" pointers for traditional tracer species
c 
c-----------------------------------------------------------------------
c     Parameters for some of the switches:
c
c     CDCMC  -- string for requesting the CMC chemistry solver
c     CDEBI  -- string for requesting the EBI chemistry solver
c     CDIEH  -- string for requesting the IEH chemistry solver
c     CDLSOD -- string for requesting the LSODE chemistry solver
c     IDCMC  -- code for using the CMC chemistry solver
c     IDEBI  -- code for using the EBI chemistry solver
c     IDIEH  -- code for using the IEH chemistry solver
c     IDLSOD -- code for using the LSODE chemistry solver
c
c-----------------------------------------------------------------------
c
      character*10 CDCMC
      character*10 CDEBI
      character*10 CDIEH
      character*10 CDLSOD
      integer      IDCMC
      integer      IDEBI
      integer      IDIEH
      integer      IDLSOD
c
      parameter( CDCMC  = "CMC       " )
      parameter( CDEBI  = "EBI       " )
      parameter( CDIEH  = "IEH       " )
      parameter( CDLSOD = "LSODE     " )
      parameter( IDCMC  = 1 )
      parameter( IDEBI  = 2 )
      parameter( IDIEH  = 3 )
      parameter( IDLSOD = 4 )
c 
c-----------------------------------------------------------------------
c    Variables for the number of species in input files:
c
c    ngas    --  number of gas species being modeled
c    naero   --  number of aersol species being modeled
c    nspec   --  total number of modeled species
c    nrad    --  number of radical species being modeled
c    nreact  --  number of chemical reactions
c    idmech  --  the code which determines which chemical mechanism is used
c    idsolv  --  the code which determines which chemstry solver to use
c    navspc  --  number of species to write to output average file
c    ndepspc --  number of species to write to deposition file
c    nicspc  --  number of species in the initial conditions file
c    nbcspc  --  number of species in the boundary conditions file
c    ntcspc  --  number of species in the top concentrations file
c    nemspc  --  number of species in the ALL emissions files
c    nsmspc  --  number of species in surface model
c    nsmrxn  --  number of reactions in surface model
c
c-----------------------------------------------------------------------
c
       integer      ngas
       integer      naero
       integer      nspec
       integer      nrad
       integer      nreact
       integer      idmech
       integer      idsolv
       integer      navspc
       integer      ndepspc
       integer      nicspc
       integer      nbcspc
       integer      ntcspc
       integer      nemspc
       integer      nsmspc
       integer      nsmrxn
c
      common /chm1/ ngas, naero, nspec, nrad, nreact,
     &                         idmech, idsolv, navspc, ndepspc, nicspc,
     &                           nbcspc, ntcspc, nsmspc, nsmrxn
c
c-----------------------------------------------------------------------
c
c    aeropt --  keyword for aerosol scheme ('NONE', 'INERT', 'CF' or 'CMU')
c    leqsam --  logical flag to use EQSAM chemistry
c
c-----------------------------------------------------------------------
c
      character*10 aeropt
      logical leqsam
c
      common /chm2/ aeropt
      common /chm2dat/ leqsam
c
c-----------------------------------------------------------------------
c     Variables for photolysis data:
c
c     nphot1   -- number of primary photolysis reactions
c     nphot2   -- number of secondary (scaled) photolysis reactions
c
c-----------------------------------------------------------------------
c
      integer nphot1
      integer nphot2
      character*20 tuvlabl(10)
c
      common /photmap/ nphot1, nphot2, tuvlabl
c 
c-----------------------------------------------------------------------
c     Mercury chemistry variables 
c     Variables for vertical profiles of Cl2, Br and BrO for Hg chemistry
c
c     NHTHAL   --  number of heights for halogen concentrations 
c     hthal    --  heights (km) for the top of each concentration range
c     cl2day  --  profile of Cl2 (ppm) with height during day over ocean
c     cl2nite --  profile of Cl2 (ppm) with height at night over ocean
c     brlprof  --  profile of Br (ppm) with height over land
c     brwprof  --  profile of Br (ppm) with height over ocean
c     brolprof --  profile of BrO (ppm) with height over land
c     browprof --  profile of BrO (ppm) with height over ocean
c
c     MWHG     --  Hg molecular weight 
c-----------------------------------------------------------------------
c
      integer, parameter :: NHTHAL = 6
      real, parameter :: MWHG = 200.59
c
      real, dimension(NHTHAL) :: hthal,brlprof,brwprof,brolprof,
     &                           browprof,cl2day,cl2nite
c
      common /halprofs/ hthal, brlprof, brwprof, brolprof, browprof,
     &                  cl2day, cl2nite
c 
c-----------------------------------------------------------------------
c     Parameters for mapping input species to internal model order:
c
c      NSPNAM  -- number of recognizable state species
c
      integer   NSPNAM
c
      parameter ( NSPNAM  =   538 )
c 
c-----------------------------------------------------------------------
c     Variables for mapping input species to internal model order:
c
c     kmap     -- mapping of species on chemistry parameters file to
c                 internal order
c
c-----------------------------------------------------------------------
c
      integer kmap(NSPNAM)
c
      common /kname/ kmap
c
      integer      kapo2,   kbzo2,   kc2o3
      integer       kcro,   kcxo3,   kepx2
      integer       kho2,      ki,     kio
      integer      kiso2,   kmeo2,      ko
      integer       ko1d,     koh,    koio
      integer      kopo3,    kro2,    kror
      integer       kto2,   ktpo2,   kxlo2
      integer       kxo2,   kxo2h,   kxo2n
      integer      kxpar,   kxprp,     kbr
      integer       kbro,     kcl,    kclo
      integer      khco3,   karpx,   khact
      integer      ktprd,   kaacd,   kacet
      integer      kald2,   kaldx,   kbenz
      integer      kcat1,     kco,   kcres
      integer      kcron,    kdms,   kech4
      integer      kepox,    keth,   ketha
      integer      kethy,   ketoh,   kfacd
      integer      kform,    kgly,   kglyd
      integer      kh2o2,    khcl,     khi
      integer      khio3,   khno3,    khoi
      integer      khono,   khpld,     ki2
      integer      ki2o2,   kino2,   kino3
      integer      kintr,   kiole,   kisop
      integer      kispd,   kispx,   kixoy
      integer       kket,   kmeoh,   kmepx
      integer      kmgly,   kn2o5,    knh3
      integer        kno,    kno2,    kno3
      integer      kntr1,   kntr2,     ko3
      integer       kole,   kopan,   kopen
      integer      kpacd,    kpan,   kpanx
      integer       kpar,    kpna,   kprpa
      integer      krooh,    kso2,    ksqt
      integer      ksulf,   kterp,    ktol
      integer      kxopn,    kxyl,   kccrs
      integer       kcg1,    kcg2,    kcg3
      integer       kcg4,   kcprm,   kfcrs
      integer      kfprm,   kivoa,   kivob
      integer        kna,    kpal,    kpca
      integer       kpcl,    kpec,    kpfe
      integer      kph2o,     kpk,    kpmg
      integer       kpmn,   kpnh4,   kpno3
      integer       kpoa,    kpsi,   kpso4
      integer       kpti,   ksoa1,   ksoa2
      integer      ksoa3,   ksoa4,   ksopa
      integer      ksopb,    kntr,    kbr2
      integer      kbrn2,   kbrn3,   kch3i
      integer       kcl2,   kcln2,   kcln3
      integer      kfmbr,   kfmcl,    khbr
      integer      khobr,   khocl,   kialk
      integer       kmb2,   kmb2c,    kmb3
      integer       kmbc,   kmbc2,    kmi2
      integer       kmib,    kmic,   kssbr
      integer      ksscl,   kssn3,    khg0
      integer       khg2,  khgiip, khgiipc
      integer       khgp,   kivod,   kivog
      integer      kpap0,   kpap1,   kpap2
      integer      kpap3,   kpap4,   kpas0
      integer      kpas1,   kpas2,   kpas3
      integer      kpas4,   kpbs0,   kpbs1
      integer      kpbs2,   kpbs3,   kpbs4
      integer      kpcp0,   kpcp1,   kpcp2
      integer      kpcp3,   kpcp4,   kpfp0
      integer      kpfp1,   kpfp2,   kpfp3
      integer      kpfp4,   kvap1,   kvap2
      integer      kvap3,   kvap4,   kvas1
      integer      kvas2,   kvas3,   kvas4
      integer      kvbs1,   kvbs2,   kvbs3
      integer      kvbs4,   kvcp1,   kvcp2
      integer      kvcp3,   kvcp4,   kvfp1
      integer      kvfp2,   kvfp3,   kvfp4
      integer      kaco3,   kactp,   kadcn
      integer      kaddc,   kapip,   kbal1
      integer      kbal2,   kbalp,   kbenp
      integer       kcho,   ketep,   kethp
      integer      khc3p,   khc5p,   khc8p
      integer      kketp,   klimp,   kmacp
      integer       kmcp,   kmcto,   kmctp
      integer      kmekp,    kmo2,   kmvkp
      integer      kolip,   kolnd,   kolnn
      integer      koltp,   korap,   kper1
      integer      kper2,   ktlp1,   ktolp
      integer       ktr2,   kualp,    kxy2
      integer      kxyl1,   kxylp,   kxyo2
      integer      kxyop,    kacd,    kace
      integer       kact,    kald,    kapi
      integer       kben,    kcsl,   kdcb1
      integer      kdcb2,   kdcb3,   kdien
      integer       keoh,    kepx,    kete
      integer      keteg,    khc3,    khc5
      integer       khc8,   khket,   khno4
      integer      kishp,    kiso,   kison
      integer       klim,   kmahp,    kmct
      integer       kmoh,   knald,    koli
      integer       kolt,   konit,    kop1
      integer       kop2,   kora1,   kora2
      integer       kpaa,   kphen,    kppn
      integer       kroh,   kuald,    kxym
      integer       kxyo,    kxyp,   kbzc3
      integer       kbzo,   kmac3,   kmco3
      integer       ko3p,   krco3,   kro2c
      integer      kro2x,   ktbuo,   kxace
      integer      kxacr,   kxaf1,   kxaf2
      integer      kxaf3,   kxbac,   kxbal
      integer      kxcch,    kxco,   kxgld
      integer      kxgly,   kxhch,   kxho2
      integer      kxipr,   kxma3,   kxmac
      integer      kxmc3,   kxmek,   kxmeo
      integer      kxmgl,   kxmvk,   kxno2
      integer       kxoh,   kxpd2,   kxrc3
      integer      kxrch,   kxrn3,   kxtbu
      integer      ky6px,   kyapx,   kyrpx
      integer      kzrn3,   kacro,   kacye
      integer      kafg1,   kafg2,   kafg3
      integer      kalk1,   kalk2,   kalk3
      integer      kalk4,   kalk5,   kapin
      integer      karo1,   karo2,   kb124
      integer      kbacl,   kbald,   kbd13
      integer      kccho,   kco3h,   kcooh
      integer      kethe,   khcho,   kiprd
      integer      kmacr,    kmek,   kmpan
      integer       kmvk,   kmxyl,   knphe
      integer      kole1,   kole2,   koxyl
      integer      kpan2,   kpbzn,   kprd2
      integer      kprpe,   kpxyl,   kr6px
      integer      krapx,   krcho,   krno3
      integer      kro3h,   ksesq,   ktolu
      integer        kxn,   ktr01,   ktr02
      integer      ktr03,   ktr04,   ktr05
      integer      ktr06,   ktr07,   ktr08
      integer      ktr09, kcrst_1,kcrst_10
      integer    kcrst_2, kcrst_3, kcrst_4
      integer    kcrst_5, kcrst_6, kcrst_7
      integer    kcrst_8, kcrst_9,   kna_1
      integer     kna_10,   kna_2,   kna_3
      integer      kna_4,   kna_5,   kna_6
      integer      kna_7,   kna_8,   kna_9
      integer     kpcl_1, kpcl_10,  kpcl_2
      integer     kpcl_3,  kpcl_4,  kpcl_5
      integer     kpcl_6,  kpcl_7,  kpcl_8
      integer     kpcl_9,  kpec_1, kpec_10
      integer     kpec_2,  kpec_3,  kpec_4
      integer     kpec_5,  kpec_6,  kpec_7
      integer     kpec_8,  kpec_9, kph2o_1
      integer   kph2o_10, kph2o_2, kph2o_3
      integer    kph2o_4, kph2o_5, kph2o_6
      integer    kph2o_7, kph2o_8, kph2o_9
      integer    kpnh4_1,kpnh4_10, kpnh4_2
      integer    kpnh4_3, kpnh4_4, kpnh4_5
      integer    kpnh4_6, kpnh4_7, kpnh4_8
      integer    kpnh4_9, kpno3_1,kpno3_10
      integer    kpno3_2, kpno3_3, kpno3_4
      integer    kpno3_5, kpno3_6, kpno3_7
      integer    kpno3_8, kpno3_9,  kpoa_1
      integer    kpoa_10,  kpoa_2,  kpoa_3
      integer     kpoa_4,  kpoa_5,  kpoa_6
      integer     kpoa_7,  kpoa_8,  kpoa_9
      integer    kpso4_1,kpso4_10, kpso4_2
      integer    kpso4_3, kpso4_4, kpso4_5
      integer    kpso4_6, kpso4_7, kpso4_8
      integer    kpso4_9, ksoa1_1,ksoa1_10
      integer    ksoa1_2, ksoa1_3, ksoa1_4
      integer    ksoa1_5, ksoa1_6, ksoa1_7
      integer    ksoa1_8, ksoa1_9, ksoa2_1
      integer   ksoa2_10, ksoa2_2, ksoa2_3
      integer    ksoa2_4, ksoa2_5, ksoa2_6
      integer    ksoa2_7, ksoa2_8, ksoa2_9
      integer    ksoa3_1,ksoa3_10, ksoa3_2
      integer    ksoa3_3, ksoa3_4, ksoa3_5
      integer    ksoa3_6, ksoa3_7, ksoa3_8
      integer    ksoa3_9, ksoa4_1,ksoa4_10
      integer    ksoa4_2, ksoa4_3, ksoa4_4
      integer    ksoa4_5, ksoa4_6, ksoa4_7
      integer    ksoa4_8, ksoa4_9, ksopa_1
      integer   ksopa_10, ksopa_2, ksopa_3
      integer    ksopa_4, ksopa_5, ksopa_6
      integer    ksopa_7, ksopa_8, ksopa_9
      integer    ksopb_1,ksopb_10, ksopb_2
      integer    ksopb_3, ksopb_4, ksopb_5
      integer    ksopb_6, ksopb_7, ksopb_8
      integer    ksopb_9
c
      equivalence
     &   (kmap(  1),   kapo2),(kmap(  2),   kbzo2),(kmap(  3),   kc2o3),
     &   (kmap(  4),    kcro),(kmap(  5),   kcxo3),(kmap(  6),   kepx2),
     &   (kmap(  7),    kho2),(kmap(  8),      ki),(kmap(  9),     kio),
     &   (kmap( 10),   kiso2),(kmap( 11),   kmeo2),(kmap( 12),      ko),
     &   (kmap( 13),    ko1d),(kmap( 14),     koh),(kmap( 15),    koio),
     &   (kmap( 16),   kopo3),(kmap( 17),    kro2),(kmap( 18),    kror),
     &   (kmap( 19),    kto2),(kmap( 20),   ktpo2),(kmap( 21),   kxlo2),
     &   (kmap( 22),    kxo2),(kmap( 23),   kxo2h),(kmap( 24),   kxo2n),
     &   (kmap( 25),   kxpar),(kmap( 26),   kxprp),(kmap( 27),     kbr),
     &   (kmap( 28),    kbro),(kmap( 29),     kcl),(kmap( 30),    kclo),
     &   (kmap( 31),   khco3),(kmap( 32),   karpx),(kmap( 33),   khact),
     &   (kmap( 34),   ktprd),(kmap( 35),   kaacd),(kmap( 36),   kacet),
     &   (kmap( 37),   kald2),(kmap( 38),   kaldx),(kmap( 39),   kbenz),
     &   (kmap( 40),   kcat1),(kmap( 41),     kco),(kmap( 42),   kcres),
     &   (kmap( 43),   kcron),(kmap( 44),    kdms),(kmap( 45),   kech4),
     &   (kmap( 46),   kepox),(kmap( 47),    keth),(kmap( 48),   ketha),
     &   (kmap( 49),   kethy),(kmap( 50),   ketoh),(kmap( 51),   kfacd),
     &   (kmap( 52),   kform),(kmap( 53),    kgly),(kmap( 54),   kglyd),
     &   (kmap( 55),   kh2o2),(kmap( 56),    khcl),(kmap( 57),     khi),
     &   (kmap( 58),   khio3),(kmap( 59),   khno3),(kmap( 60),    khoi),
     &   (kmap( 61),   khono),(kmap( 62),   khpld),(kmap( 63),     ki2),
     &   (kmap( 64),   ki2o2),(kmap( 65),   kino2),(kmap( 66),   kino3),
     &   (kmap( 67),   kintr),(kmap( 68),   kiole),(kmap( 69),   kisop),
     &   (kmap( 70),   kispd),(kmap( 71),   kispx),(kmap( 72),   kixoy),
     &   (kmap( 73),    kket),(kmap( 74),   kmeoh),(kmap( 75),   kmepx),
     &   (kmap( 76),   kmgly),(kmap( 77),   kn2o5),(kmap( 78),    knh3),
     &   (kmap( 79),     kno),(kmap( 80),    kno2),(kmap( 81),    kno3),
     &   (kmap( 82),   kntr1),(kmap( 83),   kntr2),(kmap( 84),     ko3),
     &   (kmap( 85),    kole),(kmap( 86),   kopan),(kmap( 87),   kopen),
     &   (kmap( 88),   kpacd),(kmap( 89),    kpan),(kmap( 90),   kpanx),
     &   (kmap( 91),    kpar),(kmap( 92),    kpna),(kmap( 93),   kprpa),
     &   (kmap( 94),   krooh),(kmap( 95),    kso2),(kmap( 96),    ksqt),
     &   (kmap( 97),   ksulf),(kmap( 98),   kterp),(kmap( 99),    ktol),
     &   (kmap(100),   kxopn),(kmap(101),    kxyl),(kmap(102),   kccrs),
     &   (kmap(103),    kcg1),(kmap(104),    kcg2),(kmap(105),    kcg3),
     &   (kmap(106),    kcg4),(kmap(107),   kcprm),(kmap(108),   kfcrs),
     &   (kmap(109),   kfprm),(kmap(110),   kivoa),(kmap(111),   kivob),
     &   (kmap(112),     kna),(kmap(113),    kpal),(kmap(114),    kpca),
     &   (kmap(115),    kpcl),(kmap(116),    kpec),(kmap(117),    kpfe),
     &   (kmap(118),   kph2o),(kmap(119),     kpk),(kmap(120),    kpmg),
     &   (kmap(121),    kpmn),(kmap(122),   kpnh4),(kmap(123),   kpno3),
     &   (kmap(124),    kpoa),(kmap(125),    kpsi),(kmap(126),   kpso4),
     &   (kmap(127),    kpti),(kmap(128),   ksoa1),(kmap(129),   ksoa2),
     &   (kmap(130),   ksoa3),(kmap(131),   ksoa4),(kmap(132),   ksopa),
     &   (kmap(133),   ksopb),(kmap(134),    kntr),(kmap(135),    kbr2),
     &   (kmap(136),   kbrn2),(kmap(137),   kbrn3),(kmap(138),   kch3i),
     &   (kmap(139),    kcl2),(kmap(140),   kcln2),(kmap(141),   kcln3),
     &   (kmap(142),   kfmbr),(kmap(143),   kfmcl),(kmap(144),    khbr),
     &   (kmap(145),   khobr),(kmap(146),   khocl),(kmap(147),   kialk),
     &   (kmap(148),    kmb2),(kmap(149),   kmb2c),(kmap(150),    kmb3),
     &   (kmap(151),    kmbc),(kmap(152),   kmbc2),(kmap(153),    kmi2),
     &   (kmap(154),    kmib),(kmap(155),    kmic),(kmap(156),   kssbr),
     &   (kmap(157),   ksscl),(kmap(158),   kssn3),(kmap(159),    khg0),
     &   (kmap(160),    khg2),(kmap(161),  khgiip),(kmap(162), khgiipc),
     &   (kmap(163),    khgp),(kmap(164),   kivod),(kmap(165),   kivog),
     &   (kmap(166),   kpap0),(kmap(167),   kpap1),(kmap(168),   kpap2),
     &   (kmap(169),   kpap3),(kmap(170),   kpap4),(kmap(171),   kpas0),
     &   (kmap(172),   kpas1),(kmap(173),   kpas2),(kmap(174),   kpas3),
     &   (kmap(175),   kpas4),(kmap(176),   kpbs0),(kmap(177),   kpbs1),
     &   (kmap(178),   kpbs2),(kmap(179),   kpbs3),(kmap(180),   kpbs4),
     &   (kmap(181),   kpcp0),(kmap(182),   kpcp1),(kmap(183),   kpcp2),
     &   (kmap(184),   kpcp3),(kmap(185),   kpcp4),(kmap(186),   kpfp0),
     &   (kmap(187),   kpfp1),(kmap(188),   kpfp2),(kmap(189),   kpfp3),
     &   (kmap(190),   kpfp4),(kmap(191),   kvap1),(kmap(192),   kvap2),
     &   (kmap(193),   kvap3),(kmap(194),   kvap4),(kmap(195),   kvas1),
     &   (kmap(196),   kvas2),(kmap(197),   kvas3),(kmap(198),   kvas4),
     &   (kmap(199),   kvbs1),(kmap(200),   kvbs2),(kmap(201),   kvbs3),
     &   (kmap(202),   kvbs4),(kmap(203),   kvcp1),(kmap(204),   kvcp2),
     &   (kmap(205),   kvcp3),(kmap(206),   kvcp4),(kmap(207),   kvfp1),
     &   (kmap(208),   kvfp2),(kmap(209),   kvfp3),(kmap(210),   kvfp4),
     &   (kmap(211),   kaco3),(kmap(212),   kactp),(kmap(213),   kadcn),
     &   (kmap(214),   kaddc),(kmap(215),   kapip),(kmap(216),   kbal1),
     &   (kmap(217),   kbal2),(kmap(218),   kbalp),(kmap(219),   kbenp),
     &   (kmap(220),    kcho),(kmap(221),   ketep),(kmap(222),   kethp),
     &   (kmap(223),   khc3p),(kmap(224),   khc5p),(kmap(225),   khc8p),
     &   (kmap(226),   kketp),(kmap(227),   klimp),(kmap(228),   kmacp),
     &   (kmap(229),    kmcp),(kmap(230),   kmcto),(kmap(231),   kmctp),
     &   (kmap(232),   kmekp),(kmap(233),    kmo2),(kmap(234),   kmvkp),
     &   (kmap(235),   kolip),(kmap(236),   kolnd),(kmap(237),   kolnn),
     &   (kmap(238),   koltp),(kmap(239),   korap),(kmap(240),   kper1),
     &   (kmap(241),   kper2),(kmap(242),   ktlp1),(kmap(243),   ktolp),
     &   (kmap(244),    ktr2),(kmap(245),   kualp),(kmap(246),    kxy2),
     &   (kmap(247),   kxyl1),(kmap(248),   kxylp),(kmap(249),   kxyo2),
     &   (kmap(250),   kxyop),(kmap(251),    kacd),(kmap(252),    kace),
     &   (kmap(253),    kact),(kmap(254),    kald),(kmap(255),    kapi),
     &   (kmap(256),    kben),(kmap(257),    kcsl),(kmap(258),   kdcb1),
     &   (kmap(259),   kdcb2),(kmap(260),   kdcb3),(kmap(261),   kdien),
     &   (kmap(262),    keoh),(kmap(263),    kepx),(kmap(264),    kete),
     &   (kmap(265),   keteg),(kmap(266),    khc3),(kmap(267),    khc5),
     &   (kmap(268),    khc8),(kmap(269),   khket),(kmap(270),   khno4)
c
      equivalence
     &   (kmap(271),   kishp),(kmap(272),    kiso),(kmap(273),   kison),
     &   (kmap(274),    klim),(kmap(275),   kmahp),(kmap(276),    kmct),
     &   (kmap(277),    kmoh),(kmap(278),   knald),(kmap(279),    koli),
     &   (kmap(280),    kolt),(kmap(281),   konit),(kmap(282),    kop1),
     &   (kmap(283),    kop2),(kmap(284),   kora1),(kmap(285),   kora2),
     &   (kmap(286),    kpaa),(kmap(287),   kphen),(kmap(288),    kppn),
     &   (kmap(289),    kroh),(kmap(290),   kuald),(kmap(291),    kxym),
     &   (kmap(292),    kxyo),(kmap(293),    kxyp),(kmap(294),   kbzc3),
     &   (kmap(295),    kbzo),(kmap(296),   kmac3),(kmap(297),   kmco3),
     &   (kmap(298),    ko3p),(kmap(299),   krco3),(kmap(300),   kro2c),
     &   (kmap(301),   kro2x),(kmap(302),   ktbuo),(kmap(303),   kxace),
     &   (kmap(304),   kxacr),(kmap(305),   kxaf1),(kmap(306),   kxaf2),
     &   (kmap(307),   kxaf3),(kmap(308),   kxbac),(kmap(309),   kxbal),
     &   (kmap(310),   kxcch),(kmap(311),    kxco),(kmap(312),   kxgld),
     &   (kmap(313),   kxgly),(kmap(314),   kxhch),(kmap(315),   kxho2),
     &   (kmap(316),   kxipr),(kmap(317),   kxma3),(kmap(318),   kxmac),
     &   (kmap(319),   kxmc3),(kmap(320),   kxmek),(kmap(321),   kxmeo),
     &   (kmap(322),   kxmgl),(kmap(323),   kxmvk),(kmap(324),   kxno2),
     &   (kmap(325),    kxoh),(kmap(326),   kxpd2),(kmap(327),   kxrc3),
     &   (kmap(328),   kxrch),(kmap(329),   kxrn3),(kmap(330),   kxtbu),
     &   (kmap(331),   ky6px),(kmap(332),   kyapx),(kmap(333),   kyrpx),
     &   (kmap(334),   kzrn3),(kmap(335),   kacro),(kmap(336),   kacye),
     &   (kmap(337),   kafg1),(kmap(338),   kafg2),(kmap(339),   kafg3),
     &   (kmap(340),   kalk1),(kmap(341),   kalk2),(kmap(342),   kalk3),
     &   (kmap(343),   kalk4),(kmap(344),   kalk5),(kmap(345),   kapin),
     &   (kmap(346),   karo1),(kmap(347),   karo2),(kmap(348),   kb124),
     &   (kmap(349),   kbacl),(kmap(350),   kbald),(kmap(351),   kbd13),
     &   (kmap(352),   kccho),(kmap(353),   kco3h),(kmap(354),   kcooh),
     &   (kmap(355),   kethe),(kmap(356),   khcho),(kmap(357),   kiprd),
     &   (kmap(358),   kmacr),(kmap(359),    kmek),(kmap(360),   kmpan),
     &   (kmap(361),    kmvk),(kmap(362),   kmxyl),(kmap(363),   knphe),
     &   (kmap(364),   kole1),(kmap(365),   kole2),(kmap(366),   koxyl),
     &   (kmap(367),   kpan2),(kmap(368),   kpbzn),(kmap(369),   kprd2),
     &   (kmap(370),   kprpe),(kmap(371),   kpxyl),(kmap(372),   kr6px),
     &   (kmap(373),   krapx),(kmap(374),   krcho),(kmap(375),   krno3),
     &   (kmap(376),   kro3h),(kmap(377),   ksesq),(kmap(378),   ktolu),
     &   (kmap(379),     kxn),(kmap(380),   ktr01),(kmap(381),   ktr02),
     &   (kmap(382),   ktr03),(kmap(383),   ktr04),(kmap(384),   ktr05),
     &   (kmap(385),   ktr06),(kmap(386),   ktr07),(kmap(387),   ktr08),
     &   (kmap(388),   ktr09),(kmap(389), kcrst_1),(kmap(390),kcrst_10),
     &   (kmap(391), kcrst_2),(kmap(392), kcrst_3),(kmap(393), kcrst_4),
     &   (kmap(394), kcrst_5),(kmap(395), kcrst_6),(kmap(396), kcrst_7),
     &   (kmap(397), kcrst_8),(kmap(398), kcrst_9),(kmap(399),   kna_1),
     &   (kmap(400),  kna_10),(kmap(401),   kna_2),(kmap(402),   kna_3),
     &   (kmap(403),   kna_4),(kmap(404),   kna_5),(kmap(405),   kna_6),
     &   (kmap(406),   kna_7),(kmap(407),   kna_8),(kmap(408),   kna_9),
     &   (kmap(409),  kpcl_1),(kmap(410), kpcl_10),(kmap(411),  kpcl_2),
     &   (kmap(412),  kpcl_3),(kmap(413),  kpcl_4),(kmap(414),  kpcl_5),
     &   (kmap(415),  kpcl_6),(kmap(416),  kpcl_7),(kmap(417),  kpcl_8),
     &   (kmap(418),  kpcl_9),(kmap(419),  kpec_1),(kmap(420), kpec_10),
     &   (kmap(421),  kpec_2),(kmap(422),  kpec_3),(kmap(423),  kpec_4),
     &   (kmap(424),  kpec_5),(kmap(425),  kpec_6),(kmap(426),  kpec_7),
     &   (kmap(427),  kpec_8),(kmap(428),  kpec_9),(kmap(429), kph2o_1),
     &   (kmap(430),kph2o_10),(kmap(431), kph2o_2),(kmap(432), kph2o_3),
     &   (kmap(433), kph2o_4),(kmap(434), kph2o_5),(kmap(435), kph2o_6),
     &   (kmap(436), kph2o_7),(kmap(437), kph2o_8),(kmap(438), kph2o_9),
     &   (kmap(439), kpnh4_1),(kmap(440),kpnh4_10),(kmap(441), kpnh4_2),
     &   (kmap(442), kpnh4_3),(kmap(443), kpnh4_4),(kmap(444), kpnh4_5),
     &   (kmap(445), kpnh4_6),(kmap(446), kpnh4_7),(kmap(447), kpnh4_8),
     &   (kmap(448), kpnh4_9),(kmap(449), kpno3_1),(kmap(450),kpno3_10),
     &   (kmap(451), kpno3_2),(kmap(452), kpno3_3),(kmap(453), kpno3_4),
     &   (kmap(454), kpno3_5),(kmap(455), kpno3_6),(kmap(456), kpno3_7),
     &   (kmap(457), kpno3_8),(kmap(458), kpno3_9),(kmap(459),  kpoa_1),
     &   (kmap(460), kpoa_10),(kmap(461),  kpoa_2),(kmap(462),  kpoa_3),
     &   (kmap(463),  kpoa_4),(kmap(464),  kpoa_5),(kmap(465),  kpoa_6),
     &   (kmap(466),  kpoa_7),(kmap(467),  kpoa_8),(kmap(468),  kpoa_9),
     &   (kmap(469), kpso4_1),(kmap(470),kpso4_10),(kmap(471), kpso4_2),
     &   (kmap(472), kpso4_3),(kmap(473), kpso4_4),(kmap(474), kpso4_5),
     &   (kmap(475), kpso4_6),(kmap(476), kpso4_7),(kmap(477), kpso4_8),
     &   (kmap(478), kpso4_9),(kmap(479), ksoa1_1),(kmap(480),ksoa1_10),
     &   (kmap(481), ksoa1_2),(kmap(482), ksoa1_3),(kmap(483), ksoa1_4),
     &   (kmap(484), ksoa1_5),(kmap(485), ksoa1_6),(kmap(486), ksoa1_7),
     &   (kmap(487), ksoa1_8),(kmap(488), ksoa1_9),(kmap(489), ksoa2_1),
     &   (kmap(490),ksoa2_10),(kmap(491), ksoa2_2),(kmap(492), ksoa2_3),
     &   (kmap(493), ksoa2_4),(kmap(494), ksoa2_5),(kmap(495), ksoa2_6),
     &   (kmap(496), ksoa2_7),(kmap(497), ksoa2_8),(kmap(498), ksoa2_9),
     &   (kmap(499), ksoa3_1),(kmap(500),ksoa3_10),(kmap(501), ksoa3_2),
     &   (kmap(502), ksoa3_3),(kmap(503), ksoa3_4),(kmap(504), ksoa3_5),
     &   (kmap(505), ksoa3_6),(kmap(506), ksoa3_7),(kmap(507), ksoa3_8),
     &   (kmap(508), ksoa3_9),(kmap(509), ksoa4_1),(kmap(510),ksoa4_10),
     &   (kmap(511), ksoa4_2),(kmap(512), ksoa4_3),(kmap(513), ksoa4_4),
     &   (kmap(514), ksoa4_5),(kmap(515), ksoa4_6),(kmap(516), ksoa4_7),
     &   (kmap(517), ksoa4_8),(kmap(518), ksoa4_9),(kmap(519), ksopa_1),
     &   (kmap(520),ksopa_10),(kmap(521), ksopa_2),(kmap(522), ksopa_3),
     &   (kmap(523), ksopa_4),(kmap(524), ksopa_5),(kmap(525), ksopa_6),
     &   (kmap(526), ksopa_7),(kmap(527), ksopa_8),(kmap(528), ksopa_9),
     &   (kmap(529), ksopb_1),(kmap(530),ksopb_10),(kmap(531), ksopb_2),
     &   (kmap(532), ksopb_3),(kmap(533), ksopb_4),(kmap(534), ksopb_5),
     &   (kmap(535), ksopb_6),(kmap(536), ksopb_7),(kmap(537), ksopb_8),
     &   (kmap(538), ksopb_9)
c
c-----------------------------------------------------------------------
c     Variables to define parameters for each chemical species:
c
c     henso20  -- Henry's Law constant at STP for SO2 (molar/atm)
c     tfactso2 -- Temperature dependence of SO2 Henry's Law constant (1/K)
c     nbin     -- Number of aerosol size bins
c     cwmin    -- Minimum cloud water threshold (g/m3)
c     tamin    -- Cloud water freezing threshold (K)
c
c-----------------------------------------------------------------------
c
      real henso20
      real tfactso2
      real cwmin
      real tamin
c
      common /depchm/ henso20, tfactso2, cwmin, tamin
c
      integer nbin
c
      common /aerochm/ nbin
c
c-----------------------------------------------------------------------
c     Pointers to photolysis reactions that receive T,P adjustments
c
      integer jno2rxn,jo3rxn,jhcho1rxn,jhcho2rxn,jch3chorxn
      common /phtrxn/ jno2rxn,jo3rxn,jhcho1rxn,jhcho2rxn,jch3chorxn
c
c-----------------------------------------------------------------------
c     Pointers to heterogeneous reactions to be adjusted in KHETERO
c                    N2O5 + H2O -> 2 HNO3       (1)
c                    N2O5 + HCL -> CLN2 + HNO3  (2)
c                    NTR2 + H2O -> HNO3         (3)
c                    INTR + H2O -> HNO3         (4)
c                    SO2        -> SULF         (5)
c
      integer, parameter :: NHETRXN = 5
      integer, dimension(NHETRXN) :: ihetrxn
      common /hetrxn/ ihetrxn
c
c-----------------------------------------------------------------------
c     Pointer used to lookup pig chemistry rate constant
c
c     ipigrxn  -- pointer to NO self-reaction
c                    NO + NO -> 2 NO2
c
      integer ipigrxn
c
      common /pigrxn/ ipigrxn

      integer, parameter :: l_O    = 1,
     &                      l_O1D  = 2,
     &                      l_OH   = 3,
     &                      l_HO2  = 4,
     &                      l_NO3  = 5,
     &                      l_N2O5 = 6,
     &                      l_NO   = 7,
     &                      l_NO2  = 8,
     &                      l_O3   = 9,
     &                      l_HONO = 10,
     &                      l_HNO3 = 11,
     &                      l_CO   = 12,
     &                      l_FORM = 13,
     &                      l_SO2  = 14,
     &                      l_SULF = 15
      integer, parameter :: NRPIGCHEM = 23
      integer, parameter :: NSPIGCHEM = 15
      integer, dimension(NSPIGCHEM) :: lsmap
      integer, dimension(NRPIGCHEM) :: lrmap
      common /pigchemmap/ lsmap, lrmap
c
c----------------------------------------------------------------------
c    Variables for controlling calls to aerosol routines
c
c     dtaero     -- user input time interval (min) between calls to aerosol 
c                   routines
c     dt_aero    -- adjusted (for I/O frequency) time interval (min) between
c                   calls to aerosol routines
c
      real dtaero
      real dt_aero
c
      common /aero_t/ dtaero, dt_aero
