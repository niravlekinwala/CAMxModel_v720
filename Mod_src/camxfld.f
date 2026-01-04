C*** CAMXFLD
c
      Module camxfld
      include 'camxfld.inc'
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c        This allocates the dynamic memory arrays in the CAMXFLD.COM
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
c       08/08/14    Added new snow cover input fields
c       10/13/17    Added aerosol pH field (APH)
c       10/16/17    Added time-weighted NO2 photolysis rate field (AJNO2)
c       07/23/18    Added emiss flux of Bi-Di NH3 drydep (eflxnh3)
c-----------------------------------------------------------------------
c
      Contains
c
c-----------------------------------------------------------------------
c    BEGIN SUBROUTINE ALLOC_CAMXFLD
c-----------------------------------------------------------------------
c      
         subroutine alloc_camxfld(numgrds,numcols,numrows,numlays,
     &                       numspcs,numavg,numdeps,numsmspc,avg_3d)                            
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
c    This routine allocates the met and concentration fields that depend
c    on grid size and number of species. This version is called by the
c    master node, which needs space for the entire domain.    

c    Argument descriptions:
c     Input:  
c        numgrds    I  number of grids
c        numcols    I  number of cells in the X direction
c        numrows    I  number of cells in the Y direction
c        numlays    I  number of cells in the Z direction
c        numspcs    I  number of modeled species
c        numavg     I  number of species requested for output
c        numdeps    I  number of species in deposition array
c        numsmspc   I  number of surface model species
c        avg_3d     L  flag that determines if doing 3-D averages
c     Output:  
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
         integer :: numgrds
         integer :: numcols(numgrds)
         integer :: numrows(numgrds)
         integer :: numlays(numgrds)
         integer :: numspcs
         integer :: numavg
         integer :: numdeps
         integer :: numsmspc
         logical :: avg_3d(numgrds)
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
         integer :: mvec2d
         integer :: mvec3d
         integer :: mvec1lay
         integer :: mvec4d
         integer :: mvec2a
         integer :: mvec3a
         integer :: mveclu
         integer :: mvecdp
         integer :: mvecsm
         integer :: mvecavg
         integer :: mvecedge
         integer :: mveccig
         integer :: i, j
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c  -- calculate the size of the arrays ---
c
         mvec2d = 0
         mvec3d = 0
         mvecavg = 0
         mveccig = 0
         do i=1,numgrds
            mvec2d = mvec2d + numrows(i) * numcols(i)
            mvec3d = mvec3d + numrows(i) * numcols(i) * numlays(i)
            mveccig = mveccig + 2*numlays(i)*numlays(i)*numrows(i)*numcols(i)
            if( avg_3d(i) ) then
               mvecavg = mvecavg + numrows(i) * numcols(i) *
     &                                          numlays(i) * numavg
            else
               mvecavg = mvecavg + numrows(i) * numcols(i) * numavg
            endif
         enddo
         mvec1lay = mvec2d * numspcs
         mvec4d = mvec3d * numspcs
         mvec2a = MAXVAL(numrows(1:numgrds))*MAXVAL(numcols(1:numgrds))
         mvec3a = MAXVAL(numrows(1:numgrds))*
     &              MAXVAL(numcols(1:numgrds))*MAXVAL(numlays(1:numgrds))
         mveclu = MAX(mvec2d * nlu,1)
         mvecdp = mvec2d * (numdeps*3 + 2)
         mvecsm = mvec2d * numsmspc
         mvecscr=mvec3a*numspcs+ex_scratch
         mvecscr_dp=mvec2a*(3*numdeps+2) + ex_scratch
         mvecedge = MAX(numrows(1),numcols(1))
c
c ---- allocate the arrays that are 2-D fields ---
c
         allocate( cellon (mvec2d) )
         allocate( cellat (mvec2d) )
         allocate( mapscl (mvec2d) )
         allocate( tsurf  (mvec2d) )
         allocate( topo   (mvec2d) )
         allocate( lai    (mvec2d) )
         allocate( pspt   (mvec2d) )
         allocate( sfcz0  (mvec2d) )
         allocate( snow   (mvec2d) )
         allocate( snowage(mvec2d) )
         allocate( snowrat(mvec2d) )
         allocate( snowfrc(mvec2d) )
         allocate( snowalb(mvec2d) )
         allocate( albedo (mvec2d) )
c
c ---- allocate the logical arrays related with LAI
c
         allocate( lrdlai (numgrds) )
c
c --- allocate the 3-D fields that do not depend on species ---
c
         allocate( windu  (mvec3d) )
         allocate( windv  (mvec3d) )
         allocate( pupt   (mvec3d) )
         allocate( pvpt   (mvec3d) )
         allocate( tempk  (mvec3d) )
         allocate( ptpt   (mvec3d) )
         allocate( press  (mvec3d) )
         allocate( pppt   (mvec3d) )
         allocate( height (mvec3d) )
         allocate( phpt   (mvec3d) )
         allocate( rkv    (mvec3d) )
         allocate( pkpt   (mvec3d) )
         allocate( water  (mvec3d) )
         allocate( pwpt   (mvec3d) )
cae         allocate( fcloud (mvec3d) )
         allocate( depth  (mvec3d) )
         allocate( rkx    (mvec3d) )
         allocate( rky    (mvec3d) )
         allocate( cwc    (mvec3d) )
         allocate( pwr    (mvec3d) )
         allocate( pws    (mvec3d) )
         allocate( pwg    (mvec3d) )
         allocate( cod    (mvec3d) )
         allocate( cldtrns(mvec3d) )
         allocate( cph    (mvec3d) )
         allocate( aph    (mvec3d) )
         allocate( ajno2  (mvec3d) )
         allocate( cigfrc (mvec2d) )
         allocate( cigtim (mvec2d) )
         allocate( cigwtr (mvec3d) )
         allocate( cigph  (mvec3d) )
         allocate( cigpcp (mvec3d) )
         allocate( cigent (mvec3d) )
         allocate( cigdet (mvec3d) )
         allocate( cigmas (mveccig) )
c
c  --- allocate arrays for edge cells concentratons ---
c
         allocate( bndconc(4,mvecedge,numlays(1),numspcs) )
c
c  --- allocate arrays that depend on some kind if species number ---
c
         allocate( ctop   (mvec1lay) )
c
         allocate( conc   (mvec4d) )
c
         allocate( fluxtmp( MAXVAL(numrows(1:numgrds)),
     &                      MAXVAL(numcols(1:numgrds)),
     &                      MAXVAL(numlays(1:numgrds)), numspcs ) )
c
         allocate( avcnc  (mvecavg) )
c
         allocate( entrn  (mvec3d) )
         allocate( dilut  (mvec3d) )
c
         allocate( vdep   (mvec1lay) )
         allocate( fsurf  (mveclu) )
         allocate( depfld (mvecdp) )
         allocate( eflxnh3(mvec2d) )
c
         if (lsrfmod) then
           allocate( solmas (mvecsm) )
           allocate( vegmas (mvecsm) )
           allocate( reemis (mvecsm) )
         else
           allocate( solmas (1) )
           allocate( vegmas (1) )
           allocate( reemis (1) )
         endif
c
c  --- allocate the arrays used for scratch ---
c
         allocate( scr1(mvecscr) )
         allocate( scr1_dp(mvecscr_dp) )
c
c   --- allocate arrays used for mass summary ---
c
         allocate( xmass   (numspcs,numgrds)    )
         allocate( xmass0  (numspcs,numgrds)    )
         allocate( armass  (numspcs,numgrds)    )
         allocate( ptmass  (numspcs,numgrds)    )
         allocate( fluxes  (numspcs*11,numgrds) )
         allocate( xmschem (numspcs,numgrds)    )
         allocate( xmsold  (numspcs,numgrds)    )
         allocate( resid   (numspcs,numgrds)    )
         allocate( xmsfin  (numspcs,numgrds)    )
         allocate( xmstmp  (numspcs,numgrds)    )
         allocate( pigdump (numspcs,numgrds)    )
         allocate( pgmserr (numspcs,numgrds)    )
c
         allocate( tarmass (numspcs,numgrds)    )
         allocate( tptmass (numspcs,numgrds)    )
         allocate( tfluxes (numspcs*12,numgrds) )
         allocate( tresid  (numspcs,numgrds)    )
         allocate( txmschem(numspcs,numgrds)    )
         allocate( txmsfin (numspcs,numgrds)    )
c
c  --- initialize to zero ---
c
         do j=1,numgrds
           do i=1,numspcs
             xmass(i,j) = 0.
             xmass0(i,j) = 0.
             armass(i,j) = 0.
             ptmass(i,j) = 0.
             xmschem(i,j) = 0.
             xmsold(i,j) = 0.
             resid(i,j) = 0.
             xmsfin(i,j) = 0.
             xmstmp(i,j) = 0.
             pigdump(i,j) = 0.
             pgmserr(i,j) = 0.
             tarmass(i,j) = 0.
             tptmass(i,j) = 0.
             tresid(i,j) = 0.
             txmschem(i,j) = 0.
           enddo
           do i=1,numspcs*11
             fluxes(i,j) = 0.
           enddo
           do i=1,numspcs*12
             tfluxes(i,j) = 0.
           enddo
         enddo

         aph = 0.0
         ajno2 = 0.0
         eflxnh3 = 0.0
c
c  --- allocate the fields used to walk the met ---
c
         allocate( hnxt    (mvec3d) )
         allocate( pnxt    (mvec3d) )
         allocate( unxt    (mvec3d) )
         allocate( vnxt    (mvec3d) )
         allocate( tnxt    (mvec3d) )
         allocate( tsnxt   (mvec3d) )
         allocate( knxt    (mvec3d) )
         allocate( wnxt    (mvec3d) )
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
         return
         end subroutine
c
c-----------------------------------------------------------------------
c    END SUBROUTINE ALLOC_CAMXFLD
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c    BEGIN SUBROUTINE ALLOC_CAMXFLD_NODE
c-----------------------------------------------------------------------
c
         subroutine alloc_camxfld_node(numspcs, numgrds, numavg,
     &            numdeps, numsmspc, numcols, numrows, numlays, avg_3d)
c
c-----------------------------------------------------------------------
c    Modules used:
c-----------------------------------------------------------------------
c
         use camx_includes
         use node_mod
c
         implicit none
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c    This routine allocates the met and concentration fields that depend
c    on grid size and number of species. This version is called by the
c    compute nodes, which need space for the just the slice.
c
c    Argument descriptions:
c     Input:  
c        numspcs    I  number of modeled species
c        numgrds    I  number of grids
c        numavg     I  number of species requested for output
c        numdeps    I  number of depositions species
c        numsmspc   I  number of surface model species
c        numcols    I  number of cells in the X direction
c        numrows    I  number of cells in the Y direction
c        numlays    I  number of cells in the Z direction
c        avg_3d     L  flag that determines if doing 3-D averages
c     Output:  
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c                    
         integer :: numspcs
         integer :: numgrds
         integer :: numavg
         integer :: numdeps
         integer :: numsmspc
         integer :: numcols(numgrds)
         integer :: numrows(numgrds)
         integer :: numlays(numgrds)
         logical :: avg_3d(numgrds)
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
         integer :: mvec2d
         integer :: mvec3d
         integer :: mvec1lay
         integer :: mvec4d
         integer :: mvec2a
         integer :: mvec3a
         integer :: mveclu
         integer :: mvecdp
         integer :: mvecsm
         integer :: mvecavg
         integer :: mvecedge
         integer :: mveccig
         integer :: i, j
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c  -- calculate the size of the arrays ---
c
         mvec2d = 0
         mvec3d = 0
         mvecavg = 0
         mveccig = 0
         do i=1,numgrds
            mvec2d = mvec2d + mmyp(i) * mmxp(i)
            mvec3d = mvec3d + mmyp(i) * mmxp(i) * mmzp(i)
            mveccig = mveccig + 2*mmzp(i)*mmzp(i)*mmyp(i)*mmxp(i)
            if( avg_3d(i) ) then
               mvecavg = mvecavg + mmyp(i) * mmxp(i) * mmzp(i) * numavg
            else
               mvecavg = mvecavg + mmyp(i) * mmxp(i) * numavg
            endif
         enddo
         mvec1lay = mvec2d * numspcs
         mvec4d = mvec3d * numspcs
         mvec2a = MAXVAL(mmxp(1:numgrds))*MAXVAL(mmyp(1:numgrds))
         mvec3a = MAXVAL(mmxp(1:numgrds))*
     &                MAXVAL(mmyp(1:numgrds))*MAXVAL(mmzp(1:numgrds))
         mveclu = MAX(mvec2d * nlu,1)
         mvecdp = mvec2d * (numdeps*3 + 2)
         mvecsm = mvec2d * numsmspc
         mvecscr=mvec3a*numspcs+ex_scratch
         mvecscr_dp=mvec2a*(3*numdeps+2) + ex_scratch
         mvecedge = MAX(numrows(1),numcols(1))
c
c ---- allocate the arrays that are 2-D fields ---
c
         allocate( cellon (mvec2d) )
         allocate( cellat (mvec2d) )
         allocate( mapscl (mvec2d) )
         allocate( tsurf  (mvec2d) )
         allocate( topo   (mvec2d) )
         allocate( lai    (mvec2d) )
         allocate( pspt   (mvec2d) )
         allocate( sfcz0  (mvec2d) )
         allocate( snow   (mvec2d) )
         allocate( snowage(mvec2d) )
         allocate( snowrat(mvec2d) )
         allocate( snowfrc(mvec2d) )
         allocate( snowalb(mvec2d) )
         allocate( albedo (mvec2d) )
c
c ---- allocate the logical arrays related with LAI
c
         allocate( lrdlai (numgrds) )
c
c --- allocate the 3-D fields that do not depend on species ---
c
         allocate( windu  (mvec3d) )
         allocate( windv  (mvec3d) )
         allocate( pupt   (mvec3d) )
         allocate( pvpt   (mvec3d) )
         allocate( tempk  (mvec3d) )
         allocate( ptpt   (mvec3d) )
         allocate( press  (mvec3d) )
         allocate( pppt   (mvec3d) )
         allocate( height (mvec3d) )
         allocate( phpt   (mvec3d) )
         allocate( rkv    (mvec3d) )
         allocate( pkpt   (mvec3d) )
         allocate( water  (mvec3d) )
         allocate( pwpt   (mvec3d) )
cae         allocate( fcloud (mvec3d) )
         allocate( depth  (mvec3d) )
         allocate( rkx    (mvec3d) )
         allocate( rky    (mvec3d) )
         allocate( cwc    (mvec3d) )
         allocate( pwr    (mvec3d) )
         allocate( pws    (mvec3d) )
         allocate( pwg    (mvec3d) )
         allocate( cod    (mvec3d) )
         allocate( cldtrns(mvec3d) )
         allocate( cph    (mvec3d) )
         allocate( aph    (mvec3d) )
         allocate( ajno2  (mvec3d) )
         allocate( cigfrc (mvec2d) )
         allocate( cigtim (mvec2d) )
         allocate( cigwtr (mvec3d) )
         allocate( cigph  (mvec3d) )
         allocate( cigpcp (mvec3d) )
         allocate( cigent (mvec3d) )
         allocate( cigdet (mvec3d) )
         allocate( cigmas (mveccig) )
c
c  --- allocate arrays for edge cells concentratons ---
c
         allocate( bndconc(4,mvecedge,numlays(1),numspcs) )
c
c  --- allocate arrays that depend on some kine if species number ---
c
         allocate( ctop (mvec1lay) )
c
         allocate( conc   (mvec4d) )
c
         allocate( fluxtmp (MAXVAL(mmxp(1:numgrds)),
     &                      MAXVAL(mmyp(1:numgrds)),
     &                      MAXVAL(mmzp(1:numgrds)), numspcs ) )
c
         allocate( avcnc  (mvecavg) )
c
         allocate( entrn  (mvec3d) )
         allocate( dilut  (mvec3d) )
c
         allocate( vdep   (mvec1lay) )
         allocate( fsurf  (mveclu) )
         allocate( depfld (mvecdp) )
         allocate( eflxnh3(mvec2d) )
c
         if (lsrfmod) then
           allocate( solmas (mvecsm) )
           allocate( vegmas (mvecsm) )
           allocate( reemis (mvecsm) )
         else
           allocate( solmas (1) )
           allocate( vegmas (1) )
           allocate( reemis (1) )
         endif
c
c  --- allocate the arrays used for scratch ---
c
         allocate( scr1(mvecscr) )
         allocate( scr1_dp(mvecscr_dp) )
c
c   --- allocate arrays used for mass summary ---
c
         allocate( xmass   (numspcs,numgrds)    )
         allocate( xmass0  (numspcs,numgrds)    )
         allocate( armass  (numspcs,numgrds)    )
         allocate( ptmass  (numspcs,numgrds)    )
         allocate( fluxes  (numspcs*11,numgrds) )
         allocate( xmschem (numspcs,numgrds)    )
         allocate( xmsold  (numspcs,numgrds)    )
         allocate( resid   (numspcs,numgrds)    )
         allocate( xmsfin  (numspcs,numgrds)    )
         allocate( xmstmp  (numspcs,numgrds)    )
         allocate( pigdump (numspcs,numgrds)    )
         allocate( pgmserr (numspcs,numgrds)    )
c
         allocate( tarmass (numspcs,numgrds)    )
         allocate( tptmass (numspcs,numgrds)    )
         allocate( tfluxes (numspcs*12,numgrds) )
         allocate( tresid  (numspcs,numgrds)    )
         allocate( txmschem(numspcs,numgrds)    )
         allocate( txmsfin (numspcs,numgrds)    )
c
c  --- initialize to zero ---
c
         do j=1,numgrds
           do i=1,numspcs
             xmass(i,j) = 0.
             xmass0(i,j) = 0.
             armass(i,j) = 0.
             ptmass(i,j) = 0.
             xmschem(i,j) = 0.
             xmsold(i,j) = 0.
             resid(i,j) = 0.
             xmsfin(i,j) = 0.
             xmstmp(i,j) = 0.
             pigdump(i,j) = 0.
             pgmserr(i,j) = 0.
             tarmass(i,j) = 0.
             tptmass(i,j) = 0.
             tresid(i,j) = 0.
             txmschem(i,j) = 0.
           enddo
           do i=1,numspcs*11
             fluxes(i,j) = 0.
           enddo
           do i=1,numspcs*12
             tfluxes(i,j) = 0.
           enddo
         enddo

         aph = 0.0
         ajno2 = 0.0
         eflxnh3 = 0.0
c
c  --- allocate the fields used to walk the met ---
c
         allocate( hnxt    (mvec3d) )
         allocate( pnxt    (mvec3d) )
         allocate( unxt    (mvec3d) )
         allocate( vnxt    (mvec3d) )
         allocate( tnxt    (mvec3d) )
         allocate( tsnxt   (mvec3d) )
         allocate( knxt    (mvec3d) )
         allocate( wnxt    (mvec3d) )
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
         return
         end subroutine
c
c-----------------------------------------------------------------------
c    END SUBROUTINE ALLOC_CAMXFLD_NODE
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c    BEGIN SUBROUTINE ALLOC_EMISS_ARRAY
c-----------------------------------------------------------------------
c      
         subroutine alloc_emiss_array(numgrds,numcols,numrows,
     &                                             numlays_ems,numemiss)
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
c    This routine allocates the emissions array

c    Argument descriptions:
c     Input:  
c        numgrds    I  number of grids
c        numcols    I  number of cells in the X direction
c        numrows    I  number of cells in the Y direction
c        numlays    I  number of cells in the Z direction
c        numemiss   I  number of emissions species
c     Output:  
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
         integer :: numgrds
         integer :: numcols(numgrds)
         integer :: numrows(numgrds)
         integer :: numlays_ems
         integer :: numemiss
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
         integer :: mvec3d
         integer :: mvecem
         integer :: i
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c  -- calculate the size of the arrays ---
c
         mvec3d = 0
         do i=1,numgrds
            mvec3d = mvec3d + numrows(i) * numcols(i) * numlays_ems
         enddo
         mvecem = mvec3d * numemiss
c
c ---- allocate the arrays that are 2-D fields ---
c
         allocate( aremis (mvecem) )
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
         return
         end subroutine
c
c-----------------------------------------------------------------------
c    END SUBROUTINE ALLOC_EMISS_ARRAY
c-----------------------------------------------------------------------
c
      end Module
