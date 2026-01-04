c**** SPCCB7SA
c
      subroutine spccb7sa(namecls,numcls,coefcon,coeflux,
     &                   nameyld,numyld,yieldH,yieldL,molwt_mod)
      use chmstry
      use tracer
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routine initializes the species variables that will be used
c   to determine the species that will be included in each tracer
c   species class
c    Argument descriptions:
c     Input:  
c       molwt_mod R  molecular weight of model species
c     Output:  
c       namecls   C  array of regular model species in each tracer class
c       numcls    I  the number of species contributing to this class
c       coefcon   R  2-D array of coefficients for making linear combo
c                   this one for concentrations and emissions
c       coeflux   R  2-D array of coefficients for making linear combo
c                   this one for fluxes
c       nameyld   C  array of regular model species for the yields of
c                   each tracer class
c       numyld   I  number of regular model species for the yields of
c                   each tracer class
c       yieldH    R  high-NOx yield rates for each species in each class
c       yieldL    R  low-NOx yield rates for each species in each class
c
c      Copyright 1996 - 2022
c     Ramboll
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     03/24/21   --gy--         Created based on spccb6r4sa.f
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'soap.inc'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      real         molwt_mod(MXSPEC)
      character*10 namecls(MXALCLS,MXSPEC)
      integer      numcls(MXALCLS)
      real         coefcon(MXALCLS,MXSPEC)
      real         coeflux(MXALCLS,MXSPEC)
      character*10 nameyld(MXALCLS,MXSPEC)
      integer      numyld(MXALCLS)
      real         yieldH(MXALCLS,MXSPEC),yieldL(MXALCLS,MXSPEC)
c
c-----------------------------------------------------------------------
c    Local parameters:
c-----------------------------------------------------------------------
c
c   NUMCB7        I  number of CB7 species used by PSAT treatment
c
      integer NUMCB7
c
      parameter( NUMCB7 = 62 )
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      character*10 cb7nam(NUMCB7)
      integer*4    i, j
      real         wtcb7(NUMCB7)
c
c-----------------------------------------------------------------------
c    Data statements:
c-----------------------------------------------------------------------
c
      data cb7nam /
     &         'NO        ', 'NO2       ', 'PAR       ', 'ETHA      ',
     &         'MEOH      ', 'ETOH      ', 'ETH       ', 'OLE       ',
     &         'IOLE      ', 'ISOP      ', 'TERP      ', 'FORM      ',
     &         'ALD2      ', 'ALDX      ', 'TOL       ', 'XYL       ',
     &         'ACET      ', 'BENZ      ', 'ETHY      ', 'KET       ',
     &         'PRPA      ', 'SO2       ', 'HONO      ', 'NO3       ',
     &         'N2O5      ', 'PAN       ', 'PANX      ', 'PNA       ',
     &         'NTR1      ', 'NTR2      ', 'HNO3      ', 'CRON      ',
     &         'INTR      ', 'OPAN      ', 'INO3      ', 'NH3       ',
     &         'IVOA      ', 'SQT       ', 'CG1       ', 'CG2       ',
     &         'CG3       ', 'CG4       ', 'HG0       ', 'HG2       ',
     &         'PSO4      ', 'PNO3      ', 'PNH4      ', 'PEC       ',
     &         'POA       ', 'FCRS      ', 'FPRM      ', 'CCRS      ',
     &         'CPRM      ', 'SOA1      ', 'SOA2      ', 'SOA3      ',
     &         'SOA4      ', 'SOPA      ', 'SOPB      ', 'HGP       ',
     &         'DMS       ', 'APIN      '/

      data wtcb7 /
     &          46.0       ,  46.0       ,  14.43      ,  30.07,
     &          32.04      ,  46.07      ,  28.05      ,  27.65,
     &          56.11      ,  68.12      , 136.23      ,  30.03,
     &          44.05      ,  43.65      ,  92.14      , 106.17,
     &          58.08      ,  78.11      ,  26.04      ,  28.82,
     &          44.10      ,  64.0       ,  46.0       ,  62.0,
     &         108.0       , 121.0       , 135.0       ,  79.0,
     &         119.1       , 135.1       ,  63.0       , 153.1,
     &         147.1       , 161.0       , 188.9       ,  17.0,
     &         212.0       , 204.0       , 150.0       , 150.0,
     &         180.0       , 180.0       , 200.6       , 200.6,
     &           1.0       ,   1.0       ,   1.0       ,   1.0,
     &           1.0       ,   1.0       ,   1.0       ,   1.0,
     &           1.0       ,   1.0       ,   1.0       ,   1.0,
     &           1.0       ,   1.0       ,   1.0       ,   1.0,
     &          62.1       , 136.23 /
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c  ---- initialize all arrays ---
c
      do i=1,MXALCLS
        numcls(i) = 0
        numyld(i) = 0
        do j=1,MXSPEC
          namecls(i,j) = ' '
          nameyld(i,j) = ' '
          coefcon(i,j) = 0.0
          coeflux(i,j) = 0.0
          yieldH(i,j) = 0.0
          yieldL(i,j) = 0.0
        enddo
      enddo
c
c  --- set the molecular weight for each modeled species ---
c
      do i=1,nspec
        molwt_mod(i) = 1.0
        do j=1,NUMCB7
           if( spname(i) .EQ. cb7nam(j) ) molwt_mod(i) = wtcb7(j)
        enddo
      enddo
c
c   --- if doing the OZONE or NITRATE species ---
c
      if( lozone .OR. lnitrate ) then
c
c  ---- NIT species ---
c
         numcls(ITRNIT) = 2
         namecls(ITRNIT,1) = 'NO'
         coefcon(ITRNIT,1) = 1.0
         coeflux(ITRNIT,1) = 1.0
         namecls(ITRNIT,2) = 'HONO'
         coefcon(ITRNIT,2) = 1.0
         coeflux(ITRNIT,2) = 1.0
c
c  ---- RGN species ---
c
         numcls(ITRRGN) = 4
         namecls(ITRRGN,1) = 'NO2'
         coefcon(ITRRGN,1) = 1.0
         coeflux(ITRRGN,1) = 1.0
         namecls(ITRRGN,2) = 'NO3'
         coefcon(ITRRGN,2) = 1.0
         coeflux(ITRRGN,2) = 1.0
         namecls(ITRRGN,3) = 'N2O5'
         coefcon(ITRRGN,3) = 2.0
         coeflux(ITRRGN,3) = 2.0
         namecls(ITRRGN,4) = 'INO3'
         coefcon(ITRRGN,4) = 1.0
         coeflux(ITRRGN,4) = 1.0
c
c  ---- TPN species ---
c
         numcls(ITRTPN) = 5
         namecls(ITRTPN,1) = 'PAN'
         coefcon(ITRTPN,1) = 1.0
         coeflux(ITRTPN,1) = 1.0
         namecls(ITRTPN,2) = 'PANX'
         coefcon(ITRTPN,2) = 1.0
         coeflux(ITRTPN,2) = 1.0
         namecls(ITRTPN,3) = 'PNA'
         coefcon(ITRTPN,3) = 1.0
         coeflux(ITRTPN,3) = 1.0
         namecls(ITRTPN,4) = 'OPAN'
         coefcon(ITRTPN,4) = 1.0
         coeflux(ITRTPN,4) = 1.0
         namecls(ITRTPN,5) = 'INTR'
         coefcon(ITRTPN,5) = 1.0
         coeflux(ITRTPN,5) = 1.0
c
c  ---- NTR species ---
c
         numcls(ITRNTR) = 3
         namecls(ITRNTR,1) = 'NTR1'
         coefcon(ITRNTR,1) = 1.0
         coeflux(ITRNTR,1) = 1.0
         namecls(ITRNTR,2) = 'NTR2'
         coefcon(ITRNTR,2) = 1.0
         coeflux(ITRNTR,2) = 1.0
         namecls(ITRNTR,3) = 'CRON'
         coefcon(ITRNTR,3) = 1.0
         coeflux(ITRNTR,3) = 1.0
c
c  ---- HN3 species ---
c
         numcls(ITRHN3) = 1
         namecls(ITRHN3,1) = 'HNO3'
         coefcon(ITRHN3,1) = 1.0
         coeflux(ITRHN3,1) = 1.0
      endif
c
c   --- if doing the OZONE species ---
c
      if( lozone ) then
c
c  ---- VOC species ---
c
         numcls(ITRVOC) = 21
         namecls(ITRVOC,1) = 'PAR'
         coefcon(ITRVOC,1) = 1.0
         coeflux(ITRVOC,1) = 1.0
         namecls(ITRVOC,2) = 'ETHA'
         coefcon(ITRVOC,2) = 2.0
         coeflux(ITRVOC,2) = 2.0
         namecls(ITRVOC,3) = 'MEOH'
         coefcon(ITRVOC,3) = 1.0
         coeflux(ITRVOC,3) = 1.0
         namecls(ITRVOC,4) = 'ETOH'
         coefcon(ITRVOC,4) = 2.0
         coeflux(ITRVOC,4) = 2.0
         namecls(ITRVOC,5) = 'ETH'
         coefcon(ITRVOC,5) = 2.0
         coeflux(ITRVOC,5) = 2.0
         namecls(ITRVOC,6) = 'OLE'
         coefcon(ITRVOC,6) = 2.0
         coeflux(ITRVOC,6) = 2.0
         namecls(ITRVOC,7) = 'IOLE'
         coefcon(ITRVOC,7) = 4.0
         coeflux(ITRVOC,7) = 4.0
         namecls(ITRVOC,8) = 'ISOP'
         coefcon(ITRVOC,8) = 5.0
         coeflux(ITRVOC,8) = 5.0
         namecls(ITRVOC,9) = 'TERP'
         coefcon(ITRVOC,9) = 10.0
         coeflux(ITRVOC,9) = 10.0
         namecls(ITRVOC,10) = 'FORM'
         coefcon(ITRVOC,10) = 1.0
         coeflux(ITRVOC,10) = 1.0
         namecls(ITRVOC,11) = 'ALD2'
         coefcon(ITRVOC,11) = 2.0
         coeflux(ITRVOC,11) = 2.0
         namecls(ITRVOC,12) = 'ALDX'
         coefcon(ITRVOC,12) = 2.0
         coeflux(ITRVOC,12) = 2.0
         namecls(ITRVOC,13) = 'TOL'
         coefcon(ITRVOC,13) = 7.0
         coeflux(ITRVOC,13) = 7.0
         namecls(ITRVOC,14) = 'XYL'
         coefcon(ITRVOC,14) = 8.0
         coeflux(ITRVOC,14) = 8.0
         namecls(ITRVOC,15) = 'PRPA'
         coefcon(ITRVOC,15) = 3.0
         coeflux(ITRVOC,15) = 3.0
         namecls(ITRVOC,16) = 'BENZ'
         coefcon(ITRVOC,16) = 6.0
         coeflux(ITRVOC,16) = 6.0
         namecls(ITRVOC,17) = 'ETHY'
         coefcon(ITRVOC,17) = 2.0
         coeflux(ITRVOC,17) = 2.0
         namecls(ITRVOC,18) = 'ACET'
         coefcon(ITRVOC,18) = 3.0
         coeflux(ITRVOC,18) = 3.0
         namecls(ITRVOC,19) = 'KET'
         coefcon(ITRVOC,19) = 1.0
         coeflux(ITRVOC,19) = 1.0
         namecls(ITRVOC,20) = 'APIN'
         coefcon(ITRVOC,20) = 10.0
         coeflux(ITRVOC,20) = 10.0
         namecls(ITRVOC,21) = 'SQT'
         coefcon(ITRVOC,21) = 15.0
         coeflux(ITRVOC,21) = 15.0
c
c  ---- O3-NOx species ---
c
         numcls(ITRO3N) = 1
         namecls(ITRO3N,1) = 'O3'
         coefcon(ITRO3N,1) = 0.5
         coeflux(ITRO3N,1) = 1.0
c
c  ---- O3-VOC species ---
c
         numcls(ITRO3V) = 1
         namecls(ITRO3V,1) = 'O3'
         coefcon(ITRO3V,1) = 0.5
         coeflux(ITRO3V,1) = 1.0
c
c  ---- Odd-oxygen in RGN from O3N ---
c
         numcls(ITROON) = 4
         namecls(ITROON,1) = 'NO2'
         coefcon(ITROON,1) = 0.5
         coeflux(ITROON,1) = 1.0
         namecls(ITROON,2) = 'NO3'
         coefcon(ITROON,2) = 0.5
         coeflux(ITROON,2) = 1.0
         namecls(ITROON,3) = 'N2O5'
         coefcon(ITROON,3) = 1.0
         coeflux(ITROON,3) = 2.0
         namecls(ITROON,4) = 'INO3'
         coefcon(ITROON,4) = 0.5
         coeflux(ITROON,4) = 1.0
c
c  ---- Odd-oxygen in RGN from O3V ---
c
         numcls(ITROOV) = 4
         namecls(ITROOV,1) = 'NO2'
         coefcon(ITROOV,1) = 0.5
         coeflux(ITROOV,1) = 1.0
         namecls(ITROOV,2) = 'NO3'
         coefcon(ITROOV,2) = 0.5
         coeflux(ITROOV,2) = 1.0
         namecls(ITROOV,3) = 'N2O5'
         coefcon(ITROOV,3) = 1.0
         coeflux(ITROOV,3) = 2.0
         namecls(ITROOV,4) = 'INO3'
         coefcon(ITROOV,4) = 0.5
         coeflux(ITROOV,4) = 1.0
      endif
c
c   --- if doing the SULFATE species ---
c
      if( lsulfate ) then
c
c  ---- SO2 species ---
c
        numcls(ITRSO2) = 1
        namecls(ITRSO2,1) = 'SO2'
        coefcon(ITRSO2,1) = 1.0
        coeflux(ITRSO2,1) = 1.0
c
c  ---- PS4 species ---
c
        numcls(ITRPS4) = 1
        namecls(ITRPS4,1) = 'PSO4'
        coefcon(ITRPS4,1) = 1.0
        coeflux(ITRPS4,1) = 1.0
c
c  ---- DMS species ---
c
        if( ldmschm ) then
          numcls(ITRDMS) = 1
          namecls(ITRDMS,1) = 'DMS'
          coefcon(ITRDMS,1) = 1.0
          coeflux(ITRDMS,1) = 1.0
        endif
      endif
c
c   --- if doing the NITRATE species ---
c
      if( lnitrate ) then
c
c  ---- NIT species ---
c
c         defined above
c
c  ---- RGN species ---
c
c         defined above
c
c  ---- TPN species ---
c
c         defined above
c
c  ---- NTR species ---
c
c         defined above
c
c  ---- HN3 species ---
c
c         defined above
c
c  ---- PN3 species ---
c
          numcls(ITRPN3) = 1
          namecls(ITRPN3,1) = 'PNO3'
          coefcon(ITRPN3,1) = 1.0
          coeflux(ITRPN3,1) = 1.0
c
c  ---- NH3 species ---
c
          numcls(ITRNH3) = 1
          namecls(ITRNH3,1) = 'NH3'
          coefcon(ITRNH3,1) = 1.0
          coeflux(ITRNH3,1) = 1.0
c
c  ---- PN4 species ---
c
          numcls(ITRPN4) = 1
          namecls(ITRPN4,1) = 'PNH4'
          coefcon(ITRPN4,1) = 1.0
          coeflux(ITRPN4,1) = 1.0
      endif
c
c   --- if doing the SOA species ---
c
      if( lsoa ) then
c
c  ---- ARO species ---
c
          numcls(ITRARO) = 4
          namecls(ITRARO,1) = 'BENZ'
          coefcon(ITRARO,1) = 1.0
          coeflux(ITRARO,1) = 1.0
          namecls(ITRARO,2) = 'TOL'
          coefcon(ITRARO,2) = 1.0
          coeflux(ITRARO,2) = 1.0
          namecls(ITRARO,3) = 'XYL'
          coefcon(ITRARO,3) = 1.0
          coeflux(ITRARO,3) = 1.0
          namecls(ITRARO,4) = 'IVOA'
          coefcon(ITRARO,4) = 1.0
          coeflux(ITRARO,4) = 1.0
c
c  ---- ISP species ---
c
          numcls(ITRISP) = 1
          namecls(ITRISP,1) = 'ISOP'
          coefcon(ITRISP,1) = 1.0
          coeflux(ITRISP,1) = 1.0
c
c  ---- TRP species ---
c
          numcls(ITRTRP) = 2
          namecls(ITRTRP,1) = 'TERP'
          coefcon(ITRTRP,1) = 1.0
          coeflux(ITRTRP,1) = 1.0
          namecls(ITRTRP,2) = 'APIN'
          coefcon(ITRTRP,2) = 1.0
          coeflux(ITRTRP,2) = 1.0
c
c  ---- SQT species ---
c
          numcls(ITRSQT) = 1
          namecls(ITRSQT,1) = 'SQT'
          coefcon(ITRSQT,1) = 1.0
          coeflux(ITRSQT,1) = 1.0
c
c  ---- CG1 species ---
c
          numcls(ITRCG1) = 1
          namecls(ITRCG1,1) = 'CG1'
          coefcon(ITRCG1,1) = 1.0
          coeflux(ITRCG1,1) = 1.0
          numyld(ITRCG1) = 4
          nameyld(ITRCG1,1) = 'BENZ'
          yieldH(ITRCG1,1) = MAX(EPSYLD, y_h(1,1))
          yieldL(ITRCG1,1) = MAX(EPSYLD, y_l(1,1))
          nameyld(ITRCG1,2) = 'TOL'
          yieldH(ITRCG1,2) = MAX(EPSYLD, y_h(1,2))
          yieldL(ITRCG1,2) = MAX(EPSYLD, y_l(1,2))
          nameyld(ITRCG1,3) = 'XYL'
          yieldH(ITRCG1,3) = MAX(EPSYLD, y_h(1,3))
          yieldL(ITRCG1,3) = MAX(EPSYLD, y_l(1,3))
          nameyld(ITRCG1,4) = 'IVOA'
          yieldH(ITRCG1,4) = MAX(EPSYLD, y_h(1,4))
          yieldL(ITRCG1,4) = MAX(EPSYLD, y_l(1,4))
c
c  ---- CG2 species ---
c
          numcls(ITRCG2) = 1
          namecls(ITRCG2,1) = 'CG2'
          coefcon(ITRCG2,1) = 1.0
          coeflux(ITRCG2,1) = 1.0
          numyld(ITRCG2) = 4
          nameyld(ITRCG2,1) = 'BENZ'
          yieldH(ITRCG2,1) = MAX(EPSYLD, y_h(2,1))
          yieldL(ITRCG2,1) = MAX(EPSYLD, y_l(2,1))
          nameyld(ITRCG2,2) = 'TOL'
          yieldH(ITRCG2,2) = MAX(EPSYLD, y_h(2,2))
          yieldL(ITRCG2,2) = MAX(EPSYLD, y_l(2,2))
          nameyld(ITRCG2,3) = 'XYL'
          yieldH(ITRCG2,3) = MAX(EPSYLD, y_h(2,3))
          yieldL(ITRCG2,3) = MAX(EPSYLD, y_l(2,3))
          nameyld(ITRCG2,4) = 'IVOA'
          yieldH(ITRCG2,4) = MAX(EPSYLD, y_h(2,4))
          yieldL(ITRCG2,4) = MAX(EPSYLD, y_l(2,4))
c
c  ---- CG3 species ---
c
          numcls(ITRCG3) = 1
          namecls(ITRCG3,1) = 'CG3'
          coefcon(ITRCG3,1) = 1.0
          coeflux(ITRCG3,1) = 1.0
c
c  ---- CG4 species ---
c
          numcls(ITRCG4) = 1
          namecls(ITRCG4,1) = 'CG4'
          coefcon(ITRCG4,1) = 1.0
          coeflux(ITRCG4,1) = 1.0
c
c  ---- PO1 species ---
c
          numcls(ITRPO1) = 1
          namecls(ITRPO1,1) = 'SOA1'
          coefcon(ITRPO1,1) = 1.0
          coeflux(ITRPO1,1) = 1.0
c
c  ---- PO2 species ---
c
          numcls(ITRPO2) = 1
          namecls(ITRPO2,1) = 'SOA2'
          coefcon(ITRPO2,1) = 1.0
          coeflux(ITRPO2,1) = 1.0
c
c  ---- PO3 species ---
c
          numcls(ITRPO3) = 1
          namecls(ITRPO3,1) = 'SOA3'
          coefcon(ITRPO3,1) = 1.0
          coeflux(ITRPO3,1) = 1.0
c
c  ---- PO4 species ---
c
          numcls(ITRPO4) = 1
          namecls(ITRPO4,1) = 'SOA4'
          coefcon(ITRPO4,1) = 1.0
          coeflux(ITRPO4,1) = 1.0
c
c  ---- PPA species ---
c
          numcls(ITRPPA) = 1
          namecls(ITRPPA,1) = 'SOPA'
          coefcon(ITRPPA,1) = 1.0
          coeflux(ITRPPA,1) = 1.0
          ! SPECIAL CASE: non-volatile CG; ARO -> SOPA directly (skip CG)
          numyld(ITRPPA) = 4
          nameyld(ITRPPA,1) = 'BENZ'
          yieldH(ITRPPA,1) = MAX(EPSYLD, y_h(3,1))
          yieldL(ITRPPA,1) = MAX(EPSYLD, y_l(3,1))
          nameyld(ITRPPA,2) = 'TOL'
          yieldH(ITRPPA,2) = MAX(EPSYLD, y_h(3,2))
          yieldL(ITRPPA,2) = MAX(EPSYLD, y_l(3,2))
          nameyld(ITRPPA,3) = 'XYL'
          yieldH(ITRPPA,3) = MAX(EPSYLD, y_h(3,3))
          yieldL(ITRPPA,3) = MAX(EPSYLD, y_l(3,3))
          nameyld(ITRPPA,4) = 'IVOA'
          yieldH(ITRPPA,4) = MAX(EPSYLD, y_h(3,4))
          yieldL(ITRPPA,4) = MAX(EPSYLD, y_l(3,4))
c
c  ---- PPB species ---
c
          numcls(ITRPPB) = 1
          namecls(ITRPPB,1) = 'SOPB'
          coefcon(ITRPPB,1) = 1.0
          coeflux(ITRPPB,1) = 1.0
      endif
c
c   --- if doing the PRIMARY species ---
c
      if( lprimary ) then
c
c  ---- PEC species ---
c
          numcls(ITRPEC) = 1
          namecls(ITRPEC,1) = 'PEC'
          coefcon(ITRPEC,1) = 1.0
          coeflux(ITRPEC,1) = 1.0
c
c  ---- POA species ---
c
          numcls(ITRPOA) = 1
          namecls(ITRPOA,1) = 'POA'
          coefcon(ITRPOA,1) = 1.0
          coeflux(ITRPOA,1) = 1.0
c
c  ---- PFC species ---
c
          numcls(ITRPFC) = 1
          namecls(ITRPFC,1) = 'FCRS'
          coefcon(ITRPFC,1) = 1.0
          coeflux(ITRPFC,1) = 1.0
c
c  ---- PFN species ---
c
          numcls(ITRPFN) = 1
          namecls(ITRPFN,1) = 'FPRM'
          coeflux(ITRPFN,1) = 1.0
          coefcon(ITRPFN,1) = 1.0
c
c  ---- PCC species ---
c
          numcls(ITRPCC) = 1
          namecls(ITRPCC,1) = 'CCRS'
          coefcon(ITRPCC,1) = 1.0
          coeflux(ITRPCC,1) = 1.0
c
c  ---- PCS species ---
c
          numcls(ITRPCS) = 1
          namecls(ITRPCS,1) = 'CPRM'
          coefcon(ITRPCS,1) = 1.0
          coeflux(ITRPCS,1) = 1.0

          if( lmineral ) then
c
c  ---- PFE species ---
c
              numcls(ITRPFE) = 1
              namecls(ITRPFE,1) = 'PFE'
              coefcon(ITRPFE,1) = 1.0
              coeflux(ITRPFE,1) = 1.0
c
c  ---- PMN species ---
c
              numcls(ITRPMN) = 1
              namecls(ITRPMN,1) = 'PMN'
              coefcon(ITRPMN,1) = 1.0
              coeflux(ITRPMN,1) = 1.0
c
c  ---- PMG species ---
c
              numcls(ITRPMG) = 1
              namecls(ITRPMG,1) = 'PMG'
              coefcon(ITRPMG,1) = 1.0
              coeflux(ITRPMG,1) = 1.0
c
c  ---- PK species ---
c
              numcls(ITRPK) = 1
              namecls(ITRPK,1) = 'PK'
              coefcon(ITRPK,1) = 1.0
              coeflux(ITRPK,1) = 1.0
c
c  ---- PCA species ---
c
              numcls(ITRPCA) = 1
              namecls(ITRPCA,1) = 'PCA'
              coefcon(ITRPCA,1) = 1.0
              coeflux(ITRPCA,1) = 1.0
c
c  ---- PAL species ---
c
              numcls(ITRPAL) = 1
              namecls(ITRPAL,1) = 'PAL'
              coefcon(ITRPAL,1) = 1.0
              coeflux(ITRPAL,1) = 1.0
c
c  ---- PSI species ---
c
              numcls(ITRPSI) = 1
              namecls(ITRPSI,1) = 'PSI'
              coefcon(ITRPSI,1) = 1.0
              coeflux(ITRPSI,1) = 1.0
c
c  ---- PTI species ---
c
              numcls(ITRPTI) = 1
              namecls(ITRPTI,1) = 'PTI'
              coefcon(ITRPTI,1) = 1.0
              coeflux(ITRPTI,1) = 1.0
          endif
      endif
c
c   --- if doing the MERCURY species ---
c
      if( lmercury ) then
c       
c  ---- HG0 species ---
c
          numcls(ITRHG0) = 1
          namecls(ITRHG0,1) = 'HG0'
          coefcon(ITRHG0,1) = 1.0
          coeflux(ITRHG0,1) = 1.0
c       
c  ---- HG2 species ---
c
          numcls(ITRHG2) = 1
          namecls(ITRHG2,1) = 'HG2'
          coefcon(ITRHG2,1) = 1.0
          coeflux(ITRHG2,1) = 1.0
c       
c  ---- PHG species ---
c
          numcls(ITRPHG) = 1
          namecls(ITRPHG,1) = 'HGP'
          coefcon(ITRPHG,1) = 1.0
          coeflux(ITRPHG,1) = 1.0
      endif
c       
      goto 9999
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
