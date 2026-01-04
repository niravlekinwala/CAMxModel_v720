      Module chmstry
      include 'chmdat.inc'
      include 'chmstry.inc'
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c        This allocates the dynamic memory arrays in the CHMSTRY.COM
c        include file.
c-----------------------------------------------------------------------
c
c    Argument descriptions:
c     Input:  
c     Output:  
c
c      Copyright 1996 - 2022
c     Ramboll
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c        12/15/08    --gwilson--  Added code to handle averaging of
c                                 radicals
c        03/29/11    --cemery--   Support in-line TUV with aerosol optical depth
c        04/30/13    --cemery--   Added surface model
c
c-----------------------------------------------------------------------
c
      Contains
c
c-----------------------------------------------------------------------
c   BEGIN SUBROUTINE ALLOC_CHMSTRY
c-----------------------------------------------------------------------
c
         subroutine alloc_chmstry(numgrds,numspcs,
     &                numemissfiles,numpntfiles,numrxns,numpht1,numpht2)
c
c-----------------------------------------------------------------------
c    Modules used:
c-----------------------------------------------------------------------
c
         use camx_includes
         implicit none
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c    Argument descriptions:
c     Input:  
c        numgrds        I  number of grids
c        numspcs        I  number of species
c        numemissfiles  I  number of surface emissions files
c        numpntfiles    I  number of point emissions files
c        numrxns        I  number of reactions
c        numpht1        I  number of primary photalysis reactions
c        numpht2        I  number of secondary photalysis reactions
c     Output:  
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
         integer :: numgrds
         integer :: numspcs
         integer :: numemissfiles(*)
         integer :: numpntfiles
         integer :: numrxns
         integer :: numpht1
         integer :: numpht2
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c         
         integer :: max_emissfiles
         integer :: i
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c         
         max_emissfiles = maxval( numemissfiles(1:numgrds) )
         allocate( narspc        (numgrds,max_emissfiles) )
         allocate( idxems_files  (numgrds,max_emissfiles,MXSPEC) )
         allocate( idxpnt_files  (numpntfiles,MXSPEC) )
         idxpnt_files = 0
         allocate( lgas          (numspcs) )
         allocate( emspcname     (numspcs) )
         allocate( nptspc        (numpntfiles) )
         allocate( nptsrc_files  (numpntfiles) )
         allocate( idx_start_pts (numpntfiles) )
c
         allocate( ltdep   (numrxns) )
         allocate( lpdep   (numrxns) )
         allocate( bdnl    (numspcs+1) )
c
         allocate( idphot1 (numpht1) )
         allocate( idphot2 (numpht2) )
         allocate( idphot3 (numpht2) )
         allocate( phtscl  (numpht2) )
c
         allocate( spname  (numspcs+1) )
         allocate( depsp   (4*numspcs+2) )
c
         allocate( lbcmap  (numspcs) )
         allocate( ltcmap  (numspcs) )
         allocate( lavmap  (numspcs) )
         allocate( licmap  (numspcs) )
         allocate( lemmap  (numspcs) )
         allocate( ldepmap (numspcs) )
         lbcmap = 0
         ltcmap = 0
         lavmap = 0
         licmap = 0
         lemmap = 0
         ldepmap = 0
c
         allocate( rktbl(numrxns,NTEMPR,NPRESR) )
         allocate( prkn(NZEN,numpht1,NHGHT,NTRN,NALB,NOZN) )
c
         allocate( tempr  (NTEMPR) )
         allocate( presr  (NPRESR) )
         allocate( htint  (NHGHT)  )
         allocate( zenint (NZEN)   )
c
         allocate( henry0      (numspcs) )
         allocate( tfact       (numspcs) )
         allocate( diffrat     (numspcs) )
         allocate( f0          (numspcs) )
         allocate( rscale      (numspcs) )
         allocate( mole_weight (numspcs) )
c
         allocate( roprt    (numspcs)   )
         allocate( dcut     (numspcs,2) )
         allocate( bext     (numspcs)   )
         allocate( ssa      (numspcs)   )
         allocate( rhadj    (numspcs)   )
c
         allocate( time_aero (numgrds) )
         allocate( aero_dt   (numgrds) )
         allocate( date_aero (numgrds) )

         zenint(1) = 0.
         zenint(2) = 10.
         zenint(3) = 20.
         zenint(4) = 30.
         zenint(5) = 40.
         zenint(6) = 50.
         zenint(7) = 60.
         zenint(8) = 70.
         zenint(9) = 78.
         zenint(10) = 86.
c
         narspc = 0
c
         return
         end subroutine
c
c-----------------------------------------------------------------------
c   END SUBROUTINE ALLOC_CHMSTRY
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c   BEGIN SUBROUTINE ALLOC_CHMSTRY_AVG
c-----------------------------------------------------------------------
c
         subroutine alloc_chmstry_avg(numspcs)
c
c-----------------------------------------------------------------------
c    Modules used:
c-----------------------------------------------------------------------
c
         implicit none
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c    Argument descriptions:
c     Input:  
c        numspcs    I  number of species
c     Output:  
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
         integer :: numspcs
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
         allocate( spavg   (numspcs)   )
c
         return
         end subroutine
c
c-----------------------------------------------------------------------
c   END SUBROUTINE ALLOC_CHMSTRY_AVG
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c   BEGIN SUBROUTINE ALLOC_SRFMOD
c-----------------------------------------------------------------------
c
         subroutine alloc_srfmod(nsmspc,nsmrxn)
c
c-----------------------------------------------------------------------
c    Modules used:
c-----------------------------------------------------------------------
c
         implicit none
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c    Argument descriptions:
c     Input:
c        nsmspc     I  number of surface model species
c        nsmrxn     I  number of surface model reactions
c     Output:
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
         integer :: nsmspc
         integer :: nsmrxn
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
         allocate( smspc   (nsmspc)   )
         allocate( idsmsp  (nsmspc)   )
         allocate( smpre   (nsmrxn)   )
         allocate( smprd   (nsmrxn)   )
         allocate( smskrat  (nsmrxn)   )
         allocate( smsjrat  (nsmrxn)   )
         allocate( smvkrat  (nsmrxn)   )
         allocate( smvjrat  (nsmrxn)   )
         allocate( smikrat  (nsmrxn)   )
         allocate( smijrat  (nsmrxn)   )
         allocate( idsmpre (nsmrxn)   )
         allocate( idsmprd (nsmrxn)   )
         allocate( smssrb  (nsmspc)   )
         allocate( smvsrb  (nsmspc)   )
         allocate( smisrb  (nsmspc)   )
         allocate( smlch   (nsmspc)   )
         allocate( smpen   (nsmspc)   )
         allocate( smmlt   (nsmspc)   )
c
         return
         end subroutine
c
c-----------------------------------------------------------------------
c   END SUBROUTINE ALLOC_SRFMOD
c-----------------------------------------------------------------------
c
      end Module
