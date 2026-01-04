c**** SPCSPRCSA
c
      subroutine spcsprcsa(namecls,numcls,coefcon,coeflux,
     &                     nameyld,numyld,yieldH,yieldL,molwt_in)
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
c     Output:  
c       namecls  C  array of regular model species in each tracer class
c       numcls   I  the number of species contributing to this class
c       coefcon  R  2-D array of coefficients for making linear combo
c                   this one for concentrations and emissions
c       coeflux  R  2-D array of coefficients for making linear combo
c                   this one for fluxes
c       nameyld  C  array of regular model species for the yields of
c                   each tracer class
c       numyld   I  number of regular model species for the yields of
c                   each tracer class
c       yieldH   R  high-NOx yield rates for each species in each class
c       yieldL   R  low-NOx yield rates for each species in each class
c       molwt_in R  molecular weight of model species
c
c      Copyright 1996 - 2022
c     Ramboll
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     09/20/03   --gwilson--    Original development
c     08/20/06   --bkoo--       Added ETOH, MTBE & MBUT for updated SAPRC99
c                               Fixed wts99
c     12/29/06   --bkoo--       Revised for the updated SOA scheme
c     03/18/14   --bkoo--       Revised for benzene SOA
c     08/25/16   --bkoo--       Updated for new SOAP
c     01/12/18   --bkoo--       Removed SOA class (to be updated)
c     01/16/18   --gyarwood-    Change from SAPRC99 to SAPRC07T
c     01/06/19   --cemery--     Added PFE/PMN/PMG/PK/PCA/PAL/PSI/PTI
c     06/18/19   --bkoo--       Updated SOA class for SAPRC07T
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
      character*10 namecls(MXALCLS,MXSPEC)
      integer      numcls(MXALCLS)
      real         coefcon(MXALCLS,MXSPEC)
      real         coeflux(MXALCLS,MXSPEC)
      character*10 nameyld(MXALCLS,MXSPEC)
      integer      numyld(MXALCLS)
      real         yieldH(MXALCLS,MXSPEC),yieldL(MXALCLS,MXSPEC)
      real         molwt_in(MXSPEC)
c
c-----------------------------------------------------------------------
c    Local parameters:
c-----------------------------------------------------------------------
c
c   NUMS07    I   number of SAPRC07T species used by PSAT treatment
c
      integer NUMS07
c
      parameter( NUMS07 = 69 )
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      character*10 s07nam(NUMS07)
      integer*4    i, j
      real         wts07(NUMS07)
c
c-----------------------------------------------------------------------
c    Data statements:
c-----------------------------------------------------------------------
c

      data s07nam /
     &        'NO        ', 'NO2       ', 'HCHO      ', 'CCHO      ',
     &        'RCHO      ', 'ACET      ', 'MEK       ', 'MEOH      ',
     &        'ETHE      ', 'PRPE      ', 'BD13      ', 'ISOP      ',
     &        'APIN      ', 'ACYE      ', 'BENZ      ', 'TOLU      ',
     &        'MXYL      ', 'OXYL      ', 'PXYL      ', 'B124      ',
     &        'ETOH      ', 'ALK1      ', 'ALK2      ', 'ALK3      ',
     &        'ALK4      ', 'ALK5      ', 'OLE1      ', 'OLE2      ',
     &        'ARO1      ', 'ARO2      ', 'TERP      ', 'SESQ      ',
     &        'SO2       ', 'HONO      ', 'NO3       ', 'N2O5      ',
     &        'PAN       ', 'PAN2      ', 'MPAN      ', 'PBZN      ',
     &        'PNA       ', 'RNO3      ', 'XN        ', 'NPHE      ',
     &        'HNO3      ', 'NH3       ', 'IVOA      ', 'CG1       ',
     &        'CG2       ', 'CG3       ', 'CG4       ', 'HG0       ',
     &        'HG2       ', 'PSO4      ', 'PNO3      ', 'PNH4      ',
     &        'PEC       ', 'POA       ', 'FCRS      ', 'FPRM      ',
     &        'CCRS      ', 'CPRM      ', 'SOA1      ', 'SOA2      ',
     &        'SOA3      ', 'SOA4      ', 'SOPA      ', 'SOPB      ',
     &        'HGP       '/

      data wts07 /
     &         46.01      ,  46.01      ,  30.03      ,  44.05      ,
     &         58.08      ,  58.08      ,  72.11      ,  32.04      ,
     &         28.05      ,  42.08      ,  54.09      ,  68.12      ,
     &        136.23      ,  26.04      ,  78.11      ,  92.14      ,
     &        106.17      , 106.17      , 106.17      , 120.19      ,
     &         46.07      ,  30.07      ,  44.10      ,  58.12      ,
     &         72.15      , 113.56      ,  70.13      ,  70.13      ,
     &        113.18      , 120.19      , 136.23      , 204.35      ,
     &         64.06      ,  46.01      ,  62.00      , 108.01      ,
     &        121.05      , 135.08      , 147.09      , 183.12      ,
     &         79.01      , 119.1       ,  14.01      , 139.11      ,
     &         63.01      ,  17.03      , 212.0       , 150.0       ,
     &        150.0       , 180.0       , 180.0       , 200.59      ,
     &        200.59      ,   1.0       ,   1.0       ,   1.0       ,
     &          1.0       ,   1.0       ,   1.0       ,   1.0       ,
     &          1.0       ,   1.0       ,   1.0       ,   1.0       ,
     &          1.0       ,   1.0       ,   1.0       ,   1.0       ,
     &          1.0/
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
        molwt_in(i) = 1.0
        do j=1,NUMS07
           if( spname(i) .EQ. s07nam(j) ) molwt_in(i) = wts07(j)
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
         numcls(ITRRGN) = 3
         namecls(ITRRGN,1) = 'NO2'
         coefcon(ITRRGN,1) = 1.0
         coeflux(ITRRGN,1) = 1.0
         namecls(ITRRGN,2) = 'NO3'
         coefcon(ITRRGN,2) = 1.0
         coeflux(ITRRGN,2) = 1.0
         namecls(ITRRGN,3) = 'N2O5'
         coefcon(ITRRGN,3) = 2.0
         coeflux(ITRRGN,3) = 2.0
c
c  ---- TPN species ---
c
         numcls(ITRTPN) = 5
         namecls(ITRTPN,1) = 'PAN'
         coefcon(ITRTPN,1) = 1.0
         coeflux(ITRTPN,1) = 1.0
         namecls(ITRTPN,2) = 'PAN2'
         coefcon(ITRTPN,2) = 1.0
         coeflux(ITRTPN,2) = 1.0
         namecls(ITRTPN,3) = 'MPAN'
         coefcon(ITRTPN,3) = 1.0
         coeflux(ITRTPN,3) = 1.0
         namecls(ITRTPN,4) = 'PBZN'
         coefcon(ITRTPN,4) = 1.0
         coeflux(ITRTPN,4) = 1.0
         namecls(ITRTPN,5) = 'PNA'
         coefcon(ITRTPN,5) = 1.0
         coeflux(ITRTPN,5) = 1.0
c
c  ---- NTR species ---
c
         numcls(ITRNTR) = 3
         namecls(ITRNTR,1) = 'RNO3'
         coefcon(ITRNTR,1) = 1.0
         coeflux(ITRNTR,1) = 1.0
         namecls(ITRNTR,2) = 'XN'
         coefcon(ITRNTR,2) = 1.0
         coeflux(ITRNTR,2) = 1.0
         namecls(ITRNTR,3) = 'NPHE'
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
         numcls(ITRVOC) = 30
         namecls(ITRVOC,1) = 'HCHO'
         coefcon(ITRVOC,1) = 1.0
         coeflux(ITRVOC,1) = 1.0
         namecls(ITRVOC,2) = 'CCHO'
         coefcon(ITRVOC,2) = 2.0
         coeflux(ITRVOC,2) = 2.0
         namecls(ITRVOC,3) = 'RCHO'
         coefcon(ITRVOC,3) = 3.0
         coeflux(ITRVOC,3) = 3.0
         namecls(ITRVOC,4) = 'ACET'
         coefcon(ITRVOC,4) = 3.0
         coeflux(ITRVOC,4) = 3.0
         namecls(ITRVOC,5) = 'MEK'
         coefcon(ITRVOC,5) = 4.0
         coeflux(ITRVOC,5) = 4.0
         namecls(ITRVOC,6) = 'MEOH'
         coefcon(ITRVOC,6) = 1.0
         coeflux(ITRVOC,6) = 1.0
         namecls(ITRVOC,7) = 'ETHE'
         coefcon(ITRVOC,7) = 2.0
         coeflux(ITRVOC,7) = 2.0
         namecls(ITRVOC,8) = 'PRPE'
         coefcon(ITRVOC,8) = 3.0
         coeflux(ITRVOC,8) = 3.0
         namecls(ITRVOC,9) = 'BD13'
         coefcon(ITRVOC,9) = 4.0
         coeflux(ITRVOC,9) = 4.0
         namecls(ITRVOC,10) = 'ISOP'
         coefcon(ITRVOC,10) = 5.0
         coeflux(ITRVOC,10) = 5.0
         namecls(ITRVOC,11) = 'APIN'
         coefcon(ITRVOC,11) = 10.0
         coeflux(ITRVOC,11) = 10.0
         namecls(ITRVOC,12) = 'ACYE'
         coefcon(ITRVOC,12) = 2.0
         coeflux(ITRVOC,12) = 2.0
         namecls(ITRVOC,13) = 'BENZ'
         coefcon(ITRVOC,13) = 6.0
         coeflux(ITRVOC,13) = 6.0
         namecls(ITRVOC,14) = 'TOLU'
         coefcon(ITRVOC,14) = 7.0
         coeflux(ITRVOC,14) = 7.0
         namecls(ITRVOC,15) = 'MXYL'
         coefcon(ITRVOC,15) = 8.0
         coeflux(ITRVOC,15) = 8.0
         namecls(ITRVOC,16) = 'OXYL'
         coefcon(ITRVOC,16) = 8.0
         coeflux(ITRVOC,16) = 8.0
         namecls(ITRVOC,17) = 'PXYL'
         coefcon(ITRVOC,17) = 8.0
         coeflux(ITRVOC,17) = 8.0
         namecls(ITRVOC,18) = 'B124'
         coefcon(ITRVOC,18) = 9.0
         coeflux(ITRVOC,18) = 9.0
         namecls(ITRVOC,19) = 'ETOH'
         coefcon(ITRVOC,19) = 2.0
         coeflux(ITRVOC,19) = 2.0
         namecls(ITRVOC,20) = 'ALK1'
         coefcon(ITRVOC,20) = 2.0
         coeflux(ITRVOC,20) = 2.0
         namecls(ITRVOC,21) = 'ALK2'
         coefcon(ITRVOC,21) = 2.5
         coeflux(ITRVOC,21) = 2.5
         namecls(ITRVOC,22) = 'ALK3'
         coefcon(ITRVOC,22) = 4.0
         coeflux(ITRVOC,22) = 4.0
         namecls(ITRVOC,23) = 'ALK4'
         coefcon(ITRVOC,23) = 5.4
         coeflux(ITRVOC,23) = 5.4
         namecls(ITRVOC,24) = 'ALK5'
         coefcon(ITRVOC,24) = 8.3
         coeflux(ITRVOC,24) = 8.3
         namecls(ITRVOC,25) = 'OLE1'
         coefcon(ITRVOC,25) = 5.2
         coeflux(ITRVOC,25) = 5.2
         namecls(ITRVOC,26) = 'OLE2'
         coefcon(ITRVOC,26) = 5.4
         coeflux(ITRVOC,26) = 5.4
         namecls(ITRVOC,27) = 'ARO1'
         coefcon(ITRVOC,27) = 7.2
         coeflux(ITRVOC,27) = 7.2
         namecls(ITRVOC,28) = 'ARO2'
         coefcon(ITRVOC,28) = 8.8
         coeflux(ITRVOC,28) = 8.8
         namecls(ITRVOC,29) = 'TERP'
         coefcon(ITRVOC,29) = 10.0
         coeflux(ITRVOC,29) = 10.0
         namecls(ITRVOC,30) = 'SESQ'
         coefcon(ITRVOC,30) = 15.0
         coeflux(ITRVOC,30) = 15.0
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
         numcls(ITROON) = 3
         namecls(ITROON,1) = 'NO2'
         coefcon(ITROON,1) = 0.5
         coeflux(ITROON,1) = 1.0
         namecls(ITROON,2) = 'NO3'
         coefcon(ITROON,2) = 0.5
         coeflux(ITROON,2) = 1.0
         namecls(ITROON,3) = 'N2O5'
         coefcon(ITROON,3) = 1.0
         coeflux(ITROON,3) = 2.0
c
c  ---- Odd-oxygen in RGN from O3V ---
c
         numcls(ITROOV) = 3
         namecls(ITROOV,1) = 'NO2'
         coefcon(ITROOV,1) = 0.5
         coeflux(ITROOV,1) = 1.0
         namecls(ITROOV,2) = 'NO3'
         coefcon(ITROOV,2) = 0.5
         coeflux(ITROOV,2) = 1.0
         namecls(ITROOV,3) = 'N2O5'
         coefcon(ITROOV,3) = 1.0
         coeflux(ITROOV,3) = 2.0
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
          numcls(ITRARO) = 9
          namecls(ITRARO,1) = 'BENZ'
          coefcon(ITRARO,1) = 1.0
          coeflux(ITRARO,1) = 1.0
          namecls(ITRARO,2) = 'TOLU'
          coefcon(ITRARO,2) = 1.0
          coeflux(ITRARO,2) = 1.0
          namecls(ITRARO,3) = 'ARO1'
          coefcon(ITRARO,3) = 1.0
          coeflux(ITRARO,3) = 1.0
          namecls(ITRARO,4) = 'OXYL'
          coefcon(ITRARO,4) = 1.0
          coeflux(ITRARO,4) = 1.0
          namecls(ITRARO,5) = 'MXYL'
          coefcon(ITRARO,5) = 1.0
          coeflux(ITRARO,5) = 1.0
          namecls(ITRARO,6) = 'PXYL'
          coefcon(ITRARO,6) = 1.0
          coeflux(ITRARO,6) = 1.0
          namecls(ITRARO,7) = 'B124'
          coefcon(ITRARO,7) = 1.0
          coeflux(ITRARO,7) = 1.0
          namecls(ITRARO,8) = 'ARO2'
          coefcon(ITRARO,8) = 1.0
          coeflux(ITRARO,8) = 1.0
          namecls(ITRARO,9) = 'IVOA'
          coefcon(ITRARO,9) = 1.0
          coeflux(ITRARO,9) = 1.0
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
          namecls(ITRSQT,1) = 'SESQ'
          coefcon(ITRSQT,1) = 1.0
          coeflux(ITRSQT,1) = 1.0
c
c  ---- CG1 species ---
c
          numcls(ITRCG1) = 1
          namecls(ITRCG1,1) = 'CG1'
          coefcon(ITRCG1,1) = 1.0
          coeflux(ITRCG1,1) = 1.0
          numyld(ITRCG1) = 9
          nameyld(ITRCG1,1) = 'BENZ'
          yieldH(ITRCG1,1) = MAX(EPSYLD, y_h(1,1))
          yieldL(ITRCG1,1) = MAX(EPSYLD, y_l(1,1))
          nameyld(ITRCG1,2) = 'TOLU'
          yieldH(ITRCG1,2) = MAX(EPSYLD, y_h(1,2))
          yieldL(ITRCG1,2) = MAX(EPSYLD, y_l(1,2))
          nameyld(ITRCG1,3) = 'ARO1'
          yieldH(ITRCG1,3) = MAX(EPSYLD, y_h(1,2))
          yieldL(ITRCG1,3) = MAX(EPSYLD, y_l(1,2))
          nameyld(ITRCG1,4) = 'OXYL'
          yieldH(ITRCG1,4) = MAX(EPSYLD, y_h(1,3))
          yieldL(ITRCG1,4) = MAX(EPSYLD, y_l(1,3))
          nameyld(ITRCG1,5) = 'MXYL'
          yieldH(ITRCG1,5) = MAX(EPSYLD, y_h(1,3))
          yieldL(ITRCG1,5) = MAX(EPSYLD, y_l(1,3))
          nameyld(ITRCG1,6) = 'PXYL'
          yieldH(ITRCG1,6) = MAX(EPSYLD, y_h(1,3))
          yieldL(ITRCG1,6) = MAX(EPSYLD, y_l(1,3))
          nameyld(ITRCG1,7) = 'B124'
          yieldH(ITRCG1,7) = MAX(EPSYLD, y_h(1,3))
          yieldL(ITRCG1,7) = MAX(EPSYLD, y_l(1,3))
          nameyld(ITRCG1,8) = 'ARO2'
          yieldH(ITRCG1,8) = MAX(EPSYLD, y_h(1,3))
          yieldL(ITRCG1,8) = MAX(EPSYLD, y_l(1,3))
          nameyld(ITRCG1,9) = 'IVOA'
          yieldH(ITRCG1,9) = MAX(EPSYLD, y_h(1,4))
          yieldL(ITRCG1,9) = MAX(EPSYLD, y_l(1,4))
c
c  ---- CG2 species ---
c
          numcls(ITRCG2) = 1
          namecls(ITRCG2,1) = 'CG2'
          coefcon(ITRCG2,1) = 1.0
          coeflux(ITRCG2,1) = 1.0
          numyld(ITRCG2) = 9
          nameyld(ITRCG2,1) = 'BENZ'
          yieldH(ITRCG2,1) = MAX(EPSYLD, y_h(2,1))
          yieldL(ITRCG2,1) = MAX(EPSYLD, y_l(2,1))
          nameyld(ITRCG2,2) = 'TOLU'
          yieldH(ITRCG2,2) = MAX(EPSYLD, y_h(2,2))
          yieldL(ITRCG2,2) = MAX(EPSYLD, y_l(2,2))
          nameyld(ITRCG2,3) = 'ARO1'
          yieldH(ITRCG2,3) = MAX(EPSYLD, y_h(2,2))
          yieldL(ITRCG2,3) = MAX(EPSYLD, y_l(2,2))
          nameyld(ITRCG2,4) = 'OXYL'
          yieldH(ITRCG2,4) = MAX(EPSYLD, y_h(2,3))
          yieldL(ITRCG2,4) = MAX(EPSYLD, y_l(2,3))
          nameyld(ITRCG2,5) = 'MXYL'
          yieldH(ITRCG2,5) = MAX(EPSYLD, y_h(2,3))
          yieldL(ITRCG2,5) = MAX(EPSYLD, y_l(2,3))
          nameyld(ITRCG2,6) = 'PXYL'
          yieldH(ITRCG2,6) = MAX(EPSYLD, y_h(2,3))
          yieldL(ITRCG2,6) = MAX(EPSYLD, y_l(2,3))
          nameyld(ITRCG2,7) = 'B124'
          yieldH(ITRCG2,7) = MAX(EPSYLD, y_h(2,3))
          yieldL(ITRCG2,7) = MAX(EPSYLD, y_l(2,3))
          nameyld(ITRCG2,8) = 'ARO2'
          yieldH(ITRCG2,8) = MAX(EPSYLD, y_h(2,3))
          yieldL(ITRCG2,8) = MAX(EPSYLD, y_l(2,3))
          nameyld(ITRCG2,9) = 'IVOA'
          yieldH(ITRCG2,9) = MAX(EPSYLD, y_h(2,4))
          yieldL(ITRCG2,9) = MAX(EPSYLD, y_l(2,4))
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
          ! SPECIAL CASE: non-volatile CG; anthro -> SOPA directly (skip CG)
          numyld(ITRPPA) = 9
          nameyld(ITRPPA,1) = 'BENZ'
          yieldH(ITRPPA,1) = MAX(EPSYLD, y_h(3,1))
          yieldL(ITRPPA,1) = MAX(EPSYLD, y_l(3,1))
          nameyld(ITRPPA,2) = 'TOLU'
          yieldH(ITRPPA,2) = MAX(EPSYLD, y_h(3,2))
          yieldL(ITRPPA,2) = MAX(EPSYLD, y_l(3,2))
          nameyld(ITRPPA,3) = 'ARO1'
          yieldH(ITRPPA,3) = MAX(EPSYLD, y_h(3,2))
          yieldL(ITRPPA,3) = MAX(EPSYLD, y_l(3,2))
          nameyld(ITRPPA,4) = 'OXYL'
          yieldH(ITRPPA,4) = MAX(EPSYLD, y_h(3,3))
          yieldL(ITRPPA,4) = MAX(EPSYLD, y_l(3,3))
          nameyld(ITRPPA,5) = 'MXYL'
          yieldH(ITRPPA,5) = MAX(EPSYLD, y_h(3,3))
          yieldL(ITRPPA,5) = MAX(EPSYLD, y_l(3,3))
          nameyld(ITRPPA,6) = 'PXYL'
          yieldH(ITRPPA,6) = MAX(EPSYLD, y_h(3,3))
          yieldL(ITRPPA,6) = MAX(EPSYLD, y_l(3,3))
          nameyld(ITRPPA,7) = 'B124'
          yieldH(ITRPPA,7) = MAX(EPSYLD, y_h(3,3))
          yieldL(ITRPPA,7) = MAX(EPSYLD, y_l(3,3))
          nameyld(ITRPPA,8) = 'ARO2'
          yieldH(ITRPPA,8) = MAX(EPSYLD, y_h(3,3))
          yieldL(ITRPPA,8) = MAX(EPSYLD, y_l(3,3))
          nameyld(ITRPPA,9) = 'IVOA'
          yieldH(ITRPPA,9) = MAX(EPSYLD, y_h(3,4))
          yieldL(ITRPPA,9) = MAX(EPSYLD, y_l(3,4))
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
