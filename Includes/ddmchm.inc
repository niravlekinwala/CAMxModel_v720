c----CAMx v7.20 220430
c
c     DDMCHM.INC sets species pointers for the DDM chemistry
c     These equivalences must be consistent with the internal
c     species lists defined in data statements in READCHM
c
c      Copyright 1996 - 2022
c     Ramboll
c
c     Modifications:
c      01/20/21 -gwilson-  Added "blank" pointers for traditional tracer species
c
      integer lmap(NSPNAM)
c
      common /lname/ lmap
c
      integer      lapo2,   lbzo2,   lc2o3
      integer       lcro,   lcxo3,   lepx2
      integer       lho2,      li,     lio
      integer      liso2,   lmeo2,      lo
      integer       lo1d,     loh,    loio
      integer      lopo3,    lro2,    lror
      integer       lto2,   ltpo2,   lxlo2
      integer       lxo2,   lxo2h,   lxo2n
      integer      lxpar,   lxprp,     lbr
      integer       lbro,     lcl,    lclo
      integer      lhco3,   larpx,   lhact
      integer      ltprd,   laacd,   lacet
      integer      lald2,   laldx,   lbenz
      integer      lcat1,     lco,   lcres
      integer      lcron,    ldms,   lech4
      integer      lepox,    leth,   letha
      integer      lethy,   letoh,   lfacd
      integer      lform,    lgly,   lglyd
      integer      lh2o2,    lhcl,     lhi
      integer      lhio3,   lhno3,    lhoi
      integer      lhono,   lhpld,     li2
      integer      li2o2,   lino2,   lino3
      integer      lintr,   liole,   lisop
      integer      lispd,   lispx,   lixoy
      integer       lket,   lmeoh,   lmepx
      integer      lmgly,   ln2o5,    lnh3
      integer        lno,    lno2,    lno3
      integer      lntr1,   lntr2,     lo3
      integer       lole,   lopan,   lopen
      integer      lpacd,    lpan,   lpanx
      integer       lpar,    lpna,   lprpa
      integer      lrooh,    lso2,    lsqt
      integer      lsulf,   lterp,    ltol
      integer      lxopn,    lxyl,   lccrs
      integer       lcg1,    lcg2,    lcg3
      integer       lcg4,   lcprm,   lfcrs
      integer      lfprm,   livoa,   livob
      integer        lna,    lpal,    lpca
      integer       lpcl,    lpec,    lpfe
      integer      lph2o,     lpk,    lpmg
      integer       lpmn,   lpnh4,   lpno3
      integer       lpoa,    lpsi,   lpso4
      integer       lpti,   lsoa1,   lsoa2
      integer      lsoa3,   lsoa4,   lsopa
      integer      lsopb,    lntr,    lbr2
      integer      lbrn2,   lbrn3,   lch3i
      integer       lcl2,   lcln2,   lcln3
      integer      lfmbr,   lfmcl,    lhbr
      integer      lhobr,   lhocl,   lialk
      integer       lmb2,   lmb2c,    lmb3
      integer       lmbc,   lmbc2,    lmi2
      integer       lmib,    lmic,   lssbr
      integer      lsscl,   lssn3,    lhg0
      integer       lhg2,  lhgiip, lhgiipc
      integer       lhgp,   livod,   livog
      integer      lpap0,   lpap1,   lpap2
      integer      lpap3,   lpap4,   lpas0
      integer      lpas1,   lpas2,   lpas3
      integer      lpas4,   lpbs0,   lpbs1
      integer      lpbs2,   lpbs3,   lpbs4
      integer      lpcp0,   lpcp1,   lpcp2
      integer      lpcp3,   lpcp4,   lpfp0
      integer      lpfp1,   lpfp2,   lpfp3
      integer      lpfp4,   lvap1,   lvap2
      integer      lvap3,   lvap4,   lvas1
      integer      lvas2,   lvas3,   lvas4
      integer      lvbs1,   lvbs2,   lvbs3
      integer      lvbs4,   lvcp1,   lvcp2
      integer      lvcp3,   lvcp4,   lvfp1
      integer      lvfp2,   lvfp3,   lvfp4
      integer      laco3,   lactp,   ladcn
      integer      laddc,   lapip,   lbal1
      integer      lbal2,   lbalp,   lbenp
      integer       lcho,   letep,   lethp
      integer      lhc3p,   lhc5p,   lhc8p
      integer      lketp,   llimp,   lmacp
      integer       lmcp,   lmcto,   lmctp
      integer      lmekp,    lmo2,   lmvkp
      integer      lolip,   lolnd,   lolnn
      integer      loltp,   lorap,   lper1
      integer      lper2,   ltlp1,   ltolp
      integer       ltr2,   lualp,    lxy2
      integer      lxyl1,   lxylp,   lxyo2
      integer      lxyop,    lacd,    lace
      integer       lact,    lald,    lapi
      integer       lben,    lcsl,   ldcb1
      integer      ldcb2,   ldcb3,   ldien
      integer       leoh,    lepx,    lete
      integer      leteg,    lhc3,    lhc5
      integer       lhc8,   lhket,   lhno4
      integer      lishp,    liso,   lison
      integer       llim,   lmahp,    lmct
      integer       lmoh,   lnald,    loli
      integer       lolt,   lonit,    lop1
      integer       lop2,   lora1,   lora2
      integer       lpaa,   lphen,    lppn
      integer       lroh,   luald,    lxym
      integer       lxyo,    lxyp,   lbzc3
      integer       lbzo,   lmac3,   lmco3
      integer       lo3p,   lrco3,   lro2c
      integer      lro2x,   ltbuo,   lxace
      integer      lxacr,   lxaf1,   lxaf2
      integer      lxaf3,   lxbac,   lxbal
      integer      lxcch,    lxco,   lxgld
      integer      lxgly,   lxhch,   lxho2
      integer      lxipr,   lxma3,   lxmac
      integer      lxmc3,   lxmek,   lxmeo
      integer      lxmgl,   lxmvk,   lxno2
      integer       lxoh,   lxpd2,   lxrc3
      integer      lxrch,   lxrn3,   lxtbu
      integer      ly6px,   lyapx,   lyrpx
      integer      lzrn3,   lacro,   lacye
      integer      lafg1,   lafg2,   lafg3
      integer      lalk1,   lalk2,   lalk3
      integer      lalk4,   lalk5,   lapin
      integer      laro1,   laro2,   lb124
      integer      lbacl,   lbald,   lbd13
      integer      lccho,   lco3h,   lcooh
      integer      lethe,   lhcho,   liprd
      integer      lmacr,    lmek,   lmpan
      integer       lmvk,   lmxyl,   lnphe
      integer      lole1,   lole2,   loxyl
      integer      lpan2,   lpbzn,   lprd2
      integer      lprpe,   lpxyl,   lr6px
      integer      lrapx,   lrcho,   lrno3
      integer      lro3h,   lsesq,   ltolu
      integer        lxn,   ltr01,   ltr02
      integer      ltr03,   ltr04,   ltr05
      integer      ltr06,   ltr07,   ltr08
      integer      ltr09, lcrst_1,lcrst_10
      integer    lcrst_2, lcrst_3, lcrst_4
      integer    lcrst_5, lcrst_6, lcrst_7
      integer    lcrst_8, lcrst_9,   lna_1
      integer     lna_10,   lna_2,   lna_3
      integer      lna_4,   lna_5,   lna_6
      integer      lna_7,   lna_8,   lna_9
      integer     lpcl_1, lpcl_10,  lpcl_2
      integer     lpcl_3,  lpcl_4,  lpcl_5
      integer     lpcl_6,  lpcl_7,  lpcl_8
      integer     lpcl_9,  lpec_1, lpec_10
      integer     lpec_2,  lpec_3,  lpec_4
      integer     lpec_5,  lpec_6,  lpec_7
      integer     lpec_8,  lpec_9, lph2o_1
      integer   lph2o_10, lph2o_2, lph2o_3
      integer    lph2o_4, lph2o_5, lph2o_6
      integer    lph2o_7, lph2o_8, lph2o_9
      integer    lpnh4_1,lpnh4_10, lpnh4_2
      integer    lpnh4_3, lpnh4_4, lpnh4_5
      integer    lpnh4_6, lpnh4_7, lpnh4_8
      integer    lpnh4_9, lpno3_1,lpno3_10
      integer    lpno3_2, lpno3_3, lpno3_4
      integer    lpno3_5, lpno3_6, lpno3_7
      integer    lpno3_8, lpno3_9,  lpoa_1
      integer    lpoa_10,  lpoa_2,  lpoa_3
      integer     lpoa_4,  lpoa_5,  lpoa_6
      integer     lpoa_7,  lpoa_8,  lpoa_9
      integer    lpso4_1,lpso4_10, lpso4_2
      integer    lpso4_3, lpso4_4, lpso4_5
      integer    lpso4_6, lpso4_7, lpso4_8
      integer    lpso4_9, lsoa1_1,lsoa1_10
      integer    lsoa1_2, lsoa1_3, lsoa1_4
      integer    lsoa1_5, lsoa1_6, lsoa1_7
      integer    lsoa1_8, lsoa1_9, lsoa2_1
      integer   lsoa2_10, lsoa2_2, lsoa2_3
      integer    lsoa2_4, lsoa2_5, lsoa2_6
      integer    lsoa2_7, lsoa2_8, lsoa2_9
      integer    lsoa3_1,lsoa3_10, lsoa3_2
      integer    lsoa3_3, lsoa3_4, lsoa3_5
      integer    lsoa3_6, lsoa3_7, lsoa3_8
      integer    lsoa3_9, lsoa4_1,lsoa4_10
      integer    lsoa4_2, lsoa4_3, lsoa4_4
      integer    lsoa4_5, lsoa4_6, lsoa4_7
      integer    lsoa4_8, lsoa4_9, lsopa_1
      integer   lsopa_10, lsopa_2, lsopa_3
      integer    lsopa_4, lsopa_5, lsopa_6
      integer    lsopa_7, lsopa_8, lsopa_9
      integer    lsopb_1,lsopb_10, lsopb_2
      integer    lsopb_3, lsopb_4, lsopb_5
      integer    lsopb_6, lsopb_7, lsopb_8
      integer    lsopb_9
c
      equivalence
     &   (lmap(  1),   lapo2),(lmap(  2),   lbzo2),(lmap(  3),   lc2o3),
     &   (lmap(  4),    lcro),(lmap(  5),   lcxo3),(lmap(  6),   lepx2),
     &   (lmap(  7),    lho2),(lmap(  8),      li),(lmap(  9),     lio),
     &   (lmap( 10),   liso2),(lmap( 11),   lmeo2),(lmap( 12),      lo),
     &   (lmap( 13),    lo1d),(lmap( 14),     loh),(lmap( 15),    loio),
     &   (lmap( 16),   lopo3),(lmap( 17),    lro2),(lmap( 18),    lror),
     &   (lmap( 19),    lto2),(lmap( 20),   ltpo2),(lmap( 21),   lxlo2),
     &   (lmap( 22),    lxo2),(lmap( 23),   lxo2h),(lmap( 24),   lxo2n),
     &   (lmap( 25),   lxpar),(lmap( 26),   lxprp),(lmap( 27),     lbr),
     &   (lmap( 28),    lbro),(lmap( 29),     lcl),(lmap( 30),    lclo),
     &   (lmap( 31),   lhco3),(lmap( 32),   larpx),(lmap( 33),   lhact),
     &   (lmap( 34),   ltprd),(lmap( 35),   laacd),(lmap( 36),   lacet),
     &   (lmap( 37),   lald2),(lmap( 38),   laldx),(lmap( 39),   lbenz),
     &   (lmap( 40),   lcat1),(lmap( 41),     lco),(lmap( 42),   lcres),
     &   (lmap( 43),   lcron),(lmap( 44),    ldms),(lmap( 45),   lech4),
     &   (lmap( 46),   lepox),(lmap( 47),    leth),(lmap( 48),   letha),
     &   (lmap( 49),   lethy),(lmap( 50),   letoh),(lmap( 51),   lfacd),
     &   (lmap( 52),   lform),(lmap( 53),    lgly),(lmap( 54),   lglyd),
     &   (lmap( 55),   lh2o2),(lmap( 56),    lhcl),(lmap( 57),     lhi),
     &   (lmap( 58),   lhio3),(lmap( 59),   lhno3),(lmap( 60),    lhoi),
     &   (lmap( 61),   lhono),(lmap( 62),   lhpld),(lmap( 63),     li2),
     &   (lmap( 64),   li2o2),(lmap( 65),   lino2),(lmap( 66),   lino3),
     &   (lmap( 67),   lintr),(lmap( 68),   liole),(lmap( 69),   lisop),
     &   (lmap( 70),   lispd),(lmap( 71),   lispx),(lmap( 72),   lixoy),
     &   (lmap( 73),    lket),(lmap( 74),   lmeoh),(lmap( 75),   lmepx),
     &   (lmap( 76),   lmgly),(lmap( 77),   ln2o5),(lmap( 78),    lnh3),
     &   (lmap( 79),     lno),(lmap( 80),    lno2),(lmap( 81),    lno3),
     &   (lmap( 82),   lntr1),(lmap( 83),   lntr2),(lmap( 84),     lo3),
     &   (lmap( 85),    lole),(lmap( 86),   lopan),(lmap( 87),   lopen),
     &   (lmap( 88),   lpacd),(lmap( 89),    lpan),(lmap( 90),   lpanx),
     &   (lmap( 91),    lpar),(lmap( 92),    lpna),(lmap( 93),   lprpa),
     &   (lmap( 94),   lrooh),(lmap( 95),    lso2),(lmap( 96),    lsqt),
     &   (lmap( 97),   lsulf),(lmap( 98),   lterp),(lmap( 99),    ltol),
     &   (lmap(100),   lxopn),(lmap(101),    lxyl),(lmap(102),   lccrs),
     &   (lmap(103),    lcg1),(lmap(104),    lcg2),(lmap(105),    lcg3),
     &   (lmap(106),    lcg4),(lmap(107),   lcprm),(lmap(108),   lfcrs),
     &   (lmap(109),   lfprm),(lmap(110),   livoa),(lmap(111),   livob),
     &   (lmap(112),     lna),(lmap(113),    lpal),(lmap(114),    lpca),
     &   (lmap(115),    lpcl),(lmap(116),    lpec),(lmap(117),    lpfe),
     &   (lmap(118),   lph2o),(lmap(119),     lpk),(lmap(120),    lpmg),
     &   (lmap(121),    lpmn),(lmap(122),   lpnh4),(lmap(123),   lpno3),
     &   (lmap(124),    lpoa),(lmap(125),    lpsi),(lmap(126),   lpso4),
     &   (lmap(127),    lpti),(lmap(128),   lsoa1),(lmap(129),   lsoa2),
     &   (lmap(130),   lsoa3),(lmap(131),   lsoa4),(lmap(132),   lsopa),
     &   (lmap(133),   lsopb),(lmap(134),    lntr),(lmap(135),    lbr2),
     &   (lmap(136),   lbrn2),(lmap(137),   lbrn3),(lmap(138),   lch3i),
     &   (lmap(139),    lcl2),(lmap(140),   lcln2),(lmap(141),   lcln3),
     &   (lmap(142),   lfmbr),(lmap(143),   lfmcl),(lmap(144),    lhbr),
     &   (lmap(145),   lhobr),(lmap(146),   lhocl),(lmap(147),   lialk),
     &   (lmap(148),    lmb2),(lmap(149),   lmb2c),(lmap(150),    lmb3),
     &   (lmap(151),    lmbc),(lmap(152),   lmbc2),(lmap(153),    lmi2),
     &   (lmap(154),    lmib),(lmap(155),    lmic),(lmap(156),   lssbr),
     &   (lmap(157),   lsscl),(lmap(158),   lssn3),(lmap(159),    lhg0),
     &   (lmap(160),    lhg2),(lmap(161),  lhgiip),(lmap(162), lhgiipc),
     &   (lmap(163),    lhgp),(lmap(164),   livod),(lmap(165),   livog),
     &   (lmap(166),   lpap0),(lmap(167),   lpap1),(lmap(168),   lpap2),
     &   (lmap(169),   lpap3),(lmap(170),   lpap4),(lmap(171),   lpas0),
     &   (lmap(172),   lpas1),(lmap(173),   lpas2),(lmap(174),   lpas3),
     &   (lmap(175),   lpas4),(lmap(176),   lpbs0),(lmap(177),   lpbs1),
     &   (lmap(178),   lpbs2),(lmap(179),   lpbs3),(lmap(180),   lpbs4),
     &   (lmap(181),   lpcp0),(lmap(182),   lpcp1),(lmap(183),   lpcp2),
     &   (lmap(184),   lpcp3),(lmap(185),   lpcp4),(lmap(186),   lpfp0),
     &   (lmap(187),   lpfp1),(lmap(188),   lpfp2),(lmap(189),   lpfp3),
     &   (lmap(190),   lpfp4),(lmap(191),   lvap1),(lmap(192),   lvap2),
     &   (lmap(193),   lvap3),(lmap(194),   lvap4),(lmap(195),   lvas1),
     &   (lmap(196),   lvas2),(lmap(197),   lvas3),(lmap(198),   lvas4),
     &   (lmap(199),   lvbs1),(lmap(200),   lvbs2),(lmap(201),   lvbs3),
     &   (lmap(202),   lvbs4),(lmap(203),   lvcp1),(lmap(204),   lvcp2),
     &   (lmap(205),   lvcp3),(lmap(206),   lvcp4),(lmap(207),   lvfp1),
     &   (lmap(208),   lvfp2),(lmap(209),   lvfp3),(lmap(210),   lvfp4),
     &   (lmap(211),   laco3),(lmap(212),   lactp),(lmap(213),   ladcn),
     &   (lmap(214),   laddc),(lmap(215),   lapip),(lmap(216),   lbal1),
     &   (lmap(217),   lbal2),(lmap(218),   lbalp),(lmap(219),   lbenp),
     &   (lmap(220),    lcho),(lmap(221),   letep),(lmap(222),   lethp),
     &   (lmap(223),   lhc3p),(lmap(224),   lhc5p),(lmap(225),   lhc8p),
     &   (lmap(226),   lketp),(lmap(227),   llimp),(lmap(228),   lmacp),
     &   (lmap(229),    lmcp),(lmap(230),   lmcto),(lmap(231),   lmctp),
     &   (lmap(232),   lmekp),(lmap(233),    lmo2),(lmap(234),   lmvkp),
     &   (lmap(235),   lolip),(lmap(236),   lolnd),(lmap(237),   lolnn),
     &   (lmap(238),   loltp),(lmap(239),   lorap),(lmap(240),   lper1),
     &   (lmap(241),   lper2),(lmap(242),   ltlp1),(lmap(243),   ltolp),
     &   (lmap(244),    ltr2),(lmap(245),   lualp),(lmap(246),    lxy2),
     &   (lmap(247),   lxyl1),(lmap(248),   lxylp),(lmap(249),   lxyo2),
     &   (lmap(250),   lxyop),(lmap(251),    lacd),(lmap(252),    lace),
     &   (lmap(253),    lact),(lmap(254),    lald),(lmap(255),    lapi),
     &   (lmap(256),    lben),(lmap(257),    lcsl),(lmap(258),   ldcb1),
     &   (lmap(259),   ldcb2),(lmap(260),   ldcb3),(lmap(261),   ldien),
     &   (lmap(262),    leoh),(lmap(263),    lepx),(lmap(264),    lete),
     &   (lmap(265),   leteg),(lmap(266),    lhc3),(lmap(267),    lhc5),
     &   (lmap(268),    lhc8),(lmap(269),   lhket),(lmap(270),   lhno4)
c
      equivalence
     &   (lmap(271),   lishp),(lmap(272),    liso),(lmap(273),   lison),
     &   (lmap(274),    llim),(lmap(275),   lmahp),(lmap(276),    lmct),
     &   (lmap(277),    lmoh),(lmap(278),   lnald),(lmap(279),    loli),
     &   (lmap(280),    lolt),(lmap(281),   lonit),(lmap(282),    lop1),
     &   (lmap(283),    lop2),(lmap(284),   lora1),(lmap(285),   lora2),
     &   (lmap(286),    lpaa),(lmap(287),   lphen),(lmap(288),    lppn),
     &   (lmap(289),    lroh),(lmap(290),   luald),(lmap(291),    lxym),
     &   (lmap(292),    lxyo),(lmap(293),    lxyp),(lmap(294),   lbzc3),
     &   (lmap(295),    lbzo),(lmap(296),   lmac3),(lmap(297),   lmco3),
     &   (lmap(298),    lo3p),(lmap(299),   lrco3),(lmap(300),   lro2c),
     &   (lmap(301),   lro2x),(lmap(302),   ltbuo),(lmap(303),   lxace),
     &   (lmap(304),   lxacr),(lmap(305),   lxaf1),(lmap(306),   lxaf2),
     &   (lmap(307),   lxaf3),(lmap(308),   lxbac),(lmap(309),   lxbal),
     &   (lmap(310),   lxcch),(lmap(311),    lxco),(lmap(312),   lxgld),
     &   (lmap(313),   lxgly),(lmap(314),   lxhch),(lmap(315),   lxho2),
     &   (lmap(316),   lxipr),(lmap(317),   lxma3),(lmap(318),   lxmac),
     &   (lmap(319),   lxmc3),(lmap(320),   lxmek),(lmap(321),   lxmeo),
     &   (lmap(322),   lxmgl),(lmap(323),   lxmvk),(lmap(324),   lxno2),
     &   (lmap(325),    lxoh),(lmap(326),   lxpd2),(lmap(327),   lxrc3),
     &   (lmap(328),   lxrch),(lmap(329),   lxrn3),(lmap(330),   lxtbu),
     &   (lmap(331),   ly6px),(lmap(332),   lyapx),(lmap(333),   lyrpx),
     &   (lmap(334),   lzrn3),(lmap(335),   lacro),(lmap(336),   lacye),
     &   (lmap(337),   lafg1),(lmap(338),   lafg2),(lmap(339),   lafg3),
     &   (lmap(340),   lalk1),(lmap(341),   lalk2),(lmap(342),   lalk3),
     &   (lmap(343),   lalk4),(lmap(344),   lalk5),(lmap(345),   lapin),
     &   (lmap(346),   laro1),(lmap(347),   laro2),(lmap(348),   lb124),
     &   (lmap(349),   lbacl),(lmap(350),   lbald),(lmap(351),   lbd13),
     &   (lmap(352),   lccho),(lmap(353),   lco3h),(lmap(354),   lcooh),
     &   (lmap(355),   lethe),(lmap(356),   lhcho),(lmap(357),   liprd),
     &   (lmap(358),   lmacr),(lmap(359),    lmek),(lmap(360),   lmpan),
     &   (lmap(361),    lmvk),(lmap(362),   lmxyl),(lmap(363),   lnphe),
     &   (lmap(364),   lole1),(lmap(365),   lole2),(lmap(366),   loxyl),
     &   (lmap(367),   lpan2),(lmap(368),   lpbzn),(lmap(369),   lprd2),
     &   (lmap(370),   lprpe),(lmap(371),   lpxyl),(lmap(372),   lr6px),
     &   (lmap(373),   lrapx),(lmap(374),   lrcho),(lmap(375),   lrno3),
     &   (lmap(376),   lro3h),(lmap(377),   lsesq),(lmap(378),   ltolu),
     &   (lmap(379),     lxn),(lmap(380),   ltr01),(lmap(381),   ltr02),
     &   (lmap(382),   ltr03),(lmap(383),   ltr04),(lmap(384),   ltr05),
     &   (lmap(385),   ltr06),(lmap(386),   ltr07),(lmap(387),   ltr08),
     &   (lmap(388),   ltr09),(lmap(389), lcrst_1),(lmap(390),lcrst_10),
     &   (lmap(391), lcrst_2),(lmap(392), lcrst_3),(lmap(393), lcrst_4),
     &   (lmap(394), lcrst_5),(lmap(395), lcrst_6),(lmap(396), lcrst_7),
     &   (lmap(397), lcrst_8),(lmap(398), lcrst_9),(lmap(399),   lna_1),
     &   (lmap(400),  lna_10),(lmap(401),   lna_2),(lmap(402),   lna_3),
     &   (lmap(403),   lna_4),(lmap(404),   lna_5),(lmap(405),   lna_6),
     &   (lmap(406),   lna_7),(lmap(407),   lna_8),(lmap(408),   lna_9),
     &   (lmap(409),  lpcl_1),(lmap(410), lpcl_10),(lmap(411),  lpcl_2),
     &   (lmap(412),  lpcl_3),(lmap(413),  lpcl_4),(lmap(414),  lpcl_5),
     &   (lmap(415),  lpcl_6),(lmap(416),  lpcl_7),(lmap(417),  lpcl_8),
     &   (lmap(418),  lpcl_9),(lmap(419),  lpec_1),(lmap(420), lpec_10),
     &   (lmap(421),  lpec_2),(lmap(422),  lpec_3),(lmap(423),  lpec_4),
     &   (lmap(424),  lpec_5),(lmap(425),  lpec_6),(lmap(426),  lpec_7),
     &   (lmap(427),  lpec_8),(lmap(428),  lpec_9),(lmap(429), lph2o_1),
     &   (lmap(430),lph2o_10),(lmap(431), lph2o_2),(lmap(432), lph2o_3),
     &   (lmap(433), lph2o_4),(lmap(434), lph2o_5),(lmap(435), lph2o_6),
     &   (lmap(436), lph2o_7),(lmap(437), lph2o_8),(lmap(438), lph2o_9),
     &   (lmap(439), lpnh4_1),(lmap(440),lpnh4_10),(lmap(441), lpnh4_2),
     &   (lmap(442), lpnh4_3),(lmap(443), lpnh4_4),(lmap(444), lpnh4_5),
     &   (lmap(445), lpnh4_6),(lmap(446), lpnh4_7),(lmap(447), lpnh4_8),
     &   (lmap(448), lpnh4_9),(lmap(449), lpno3_1),(lmap(450),lpno3_10),
     &   (lmap(451), lpno3_2),(lmap(452), lpno3_3),(lmap(453), lpno3_4),
     &   (lmap(454), lpno3_5),(lmap(455), lpno3_6),(lmap(456), lpno3_7),
     &   (lmap(457), lpno3_8),(lmap(458), lpno3_9),(lmap(459),  lpoa_1),
     &   (lmap(460), lpoa_10),(lmap(461),  lpoa_2),(lmap(462),  lpoa_3),
     &   (lmap(463),  lpoa_4),(lmap(464),  lpoa_5),(lmap(465),  lpoa_6),
     &   (lmap(466),  lpoa_7),(lmap(467),  lpoa_8),(lmap(468),  lpoa_9),
     &   (lmap(469), lpso4_1),(lmap(470),lpso4_10),(lmap(471), lpso4_2),
     &   (lmap(472), lpso4_3),(lmap(473), lpso4_4),(lmap(474), lpso4_5),
     &   (lmap(475), lpso4_6),(lmap(476), lpso4_7),(lmap(477), lpso4_8),
     &   (lmap(478), lpso4_9),(lmap(479), lsoa1_1),(lmap(480),lsoa1_10),
     &   (lmap(481), lsoa1_2),(lmap(482), lsoa1_3),(lmap(483), lsoa1_4),
     &   (lmap(484), lsoa1_5),(lmap(485), lsoa1_6),(lmap(486), lsoa1_7),
     &   (lmap(487), lsoa1_8),(lmap(488), lsoa1_9),(lmap(489), lsoa2_1),
     &   (lmap(490),lsoa2_10),(lmap(491), lsoa2_2),(lmap(492), lsoa2_3),
     &   (lmap(493), lsoa2_4),(lmap(494), lsoa2_5),(lmap(495), lsoa2_6),
     &   (lmap(496), lsoa2_7),(lmap(497), lsoa2_8),(lmap(498), lsoa2_9),
     &   (lmap(499), lsoa3_1),(lmap(500),lsoa3_10),(lmap(501), lsoa3_2),
     &   (lmap(502), lsoa3_3),(lmap(503), lsoa3_4),(lmap(504), lsoa3_5),
     &   (lmap(505), lsoa3_6),(lmap(506), lsoa3_7),(lmap(507), lsoa3_8),
     &   (lmap(508), lsoa3_9),(lmap(509), lsoa4_1),(lmap(510),lsoa4_10),
     &   (lmap(511), lsoa4_2),(lmap(512), lsoa4_3),(lmap(513), lsoa4_4),
     &   (lmap(514), lsoa4_5),(lmap(515), lsoa4_6),(lmap(516), lsoa4_7),
     &   (lmap(517), lsoa4_8),(lmap(518), lsoa4_9),(lmap(519), lsopa_1),
     &   (lmap(520),lsopa_10),(lmap(521), lsopa_2),(lmap(522), lsopa_3),
     &   (lmap(523), lsopa_4),(lmap(524), lsopa_5),(lmap(525), lsopa_6),
     &   (lmap(526), lsopa_7),(lmap(527), lsopa_8),(lmap(528), lsopa_9),
     &   (lmap(529), lsopb_1),(lmap(530),lsopb_10),(lmap(531), lsopb_2),
     &   (lmap(532), lsopb_3),(lmap(533), lsopb_4),(lmap(534), lsopb_5),
     &   (lmap(535), lsopb_6),(lmap(536), lsopb_7),(lmap(537), lsopb_8),
     &   (lmap(538), lsopb_9)
