c**** SPECDDM
c
      subroutine specddm( )
      use filunit
      use chmstry
      use grid
      use tracer
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routine sets up the species names and pointers into the species
c   for all of the DDM species.  Pointers will be set up for both the
c   concentration array and the emissions array.
c
c      Copyright 1996 - 2022
c     Ramboll
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     03/23/99   --gwilson--    Original development
c     10/03/03   --gwilson--    Added flag for gaseous species
c     07/16/07   --bkoo--       Revised for HDDM
c     06/11/08   --bkoo--       Added rate constant sensitivity
c     11/12/09   --gwilson--    Added initialization of factor for
c                               applying new type of top boundary
c     08/23/13   --bkoo--       Added PM species for DDM
c     03/18/14   --bkoo--       Added species for benzene SOA
c     08/25/16   --bkoo--       Updated DDM-PM species mapping for new SOAP
c     11/28/16   --bkoo--       Corrected top fac argument in filspddm call
c                               Added pttop initialization
c     08/09/18   --bkoo--       Added rate term sensitivity
c     10/26/18   --bkoo--       Updated DDM-PM species mapping for in-cloud
c                               SOA formation & revised CA scaling
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'flags.inc'
      include 'camx_aero.inc'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      character*10 caffec, cinflu, srcnam
      character*3  edgnam(5)
      integer      ispc, iaffec, iddm, mvecedge, nedge, i, l
      integer*8    idxsa, mvec4d
      logical      lout, lsns
c
c-----------------------------------------------------------------------
c    Data statements:
c-----------------------------------------------------------------------
c
      data nedge /5/
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c  --- load the names of the boundary edges ---
c
      edgnam(IDXBWS) = 'WST'
      edgnam(IDXBES) = 'EST'
      edgnam(IDXBST) = 'STH'
      edgnam(IDXBNT) = 'NTH'
      edgnam(IDXBTP) = 'TOP'
c
c  --- echo the table of species to the output file ---
c
      write(idiag,*) ' '
      write(idiag,*) ' Affected   Influencing   Source',
     &                   '                       Long            Short' 
      write(idiag,*) ' Species      Species      Type ',
     &                   '      Group   Region   Name            Name' 
      write(idiag,'(1X,79A)') ('-',i=1,79)
c
c  --- calculate the number of DDM families per model species ---
c
      if( lbndry ) then
          nbdic = 5 * nbcddm + nicddm
      else
          nbdic = nbcddm + nicddm
      endif

      nddmsp = nbdic + nemddm * nregin * ngroup + nrateddm
     &                                          + ntermddm + nhddm
c
c  --- calculate the number of DDM species needed and allocate
c      some arrays ---
c
      ntotsp = nspec * nddmsp

      if ( nddmsp .GT. MXFDDM ) goto 7000
      if ( ntotsp .GT. MXTRSP ) goto 7001
      
      lsns = .NOT. lmpi
      mvecedge = MAX(ncol(1),nrow(1))
      call alloc_ddm(lsns,ngrid,ncol,nrow,nlay,nlayers_ems,nspec,
     &                                 nddmsp,nhddm,mvecedge,iout)
c
c  --- loop over the modeled species and setup the names ---
c
      do ispc=1,nspec
c
c  --- set the pointer into the array for this family ---
c
          iptddm(ispc) = (ispc-1)*nddmsp + 1
          iaffec = ispc
          caffec = spname(ispc)
c
c  --- set the flag for determining if this species should
c      be output to average file ---
c
          lout = .FALSE.
          do i=1,navspc
            if( lavmap(i) .EQ. ispc ) lout = .TRUE.
          enddo
c
c  --- intialize the index into the array ---
c
          idxddm = iptddm(ispc) - 1
c
c  --- loop over all of the DDM initial condition species ---
c
          do iddm = 1,nicddm
             cinflu = icddmsp(iddm)
             srcnam = 'IC'
c
c  --- call routine to fill the family of species names ---
c
             call filspddm(iaffec,caffec,cinflu,
     &                                    srcnam,idxddm,0,0,lout,0.0)
          enddo
c
c  --- loop over all of the DDM boundary condition species ---
c
          do iddm = 1,nbcddm
             cinflu = bcddmsp(iddm)
c
c  --- if the stratify boundary is off, just fill the one species name ---
c
             if( .NOT. lbndry ) then
                srcnam = 'BCALL'
c
c  --- call routine to fill the family of species names ---
c
                call filspddm(iaffec,caffec,cinflu,
     &                                     srcnam,idxddm,0,0,lout,1.0)
c
c  --- otherwise, loop over all of the boundary edges ---
c
             else
                do i=1,nedge
                   srcnam = 'BC'//edgnam(i)(1:3)
c
c  --- call routine to fill the family of species names ---
c
                   if ( i.eq.IDXBTP ) then
                     call filspddm(iaffec,caffec,cinflu,
     &                                       srcnam,idxddm,0,0,lout,1.0)
                   else
                     call filspddm(iaffec,caffec,cinflu,
     &                                       srcnam,idxddm,0,0,lout,0.0)
                   endif
                enddo
             endif
          enddo
c
c  --- loop over all of the DDM emissions condition species ---
c
          do iddm = 1,nemddm
             cinflu = emddmsp(iddm)
             srcnam = 'EM'
c
c  --- call routine to fill the family of species names ---
c
              call filspddm(iaffec,caffec,cinflu,srcnam,idxddm,
     &                                    MAX(1,ngroup),nregin,lout,0.0)
          enddo
c
c  --- loop over all of the rate constant sensitivity groups ---
c
          do iddm = 1,nrateddm
             cinflu = rateddm(iddm)
             srcnam = 'RATE'
c
c  --- call routine to fill the family of species names ---
c
             call filspddm(iaffec,caffec,cinflu,
     &                                    srcnam,idxddm,iddm,0,lout,0.0)
          enddo
c
c  --- loop over all of the rate term sensitivity groups ---
c
          do iddm = 1,ntermddm
             cinflu = termddm(iddm)
             srcnam = 'TERM'
c
c  --- call routine to fill the family of species names ---
c
             call filspddm(iaffec,caffec,cinflu,
     &                                    srcnam,idxddm,iddm,0,lout,0.0)
          enddo
c
c  --- loop over all of the HDDM sensitivity groups ---
c
          do iddm = 1,nhddm
             cinflu = ' '
             srcnam = 'HDDM'
c
c  --- call routine to fill the family of species names ---
c
             call filspddm(iaffec,caffec,cinflu,
     &                                    srcnam,idxddm,iddm,0,lout,0.0)
          enddo
c
c  --- get then next affect species ---
c
      enddo
      ntotsp = idxddm
      ipttim = ntotsp + 1
      nsaspc = ntotsp
      write(idiag,'(1X,79A)') ('-',i=1,79)
      write(idiag,*) 
c
c  --- set the flag for gaseous species ---
c
      lsagas = .FALSE.
      do i=iptddm(nrad+1),iptddm(ngas)+nddmsp-1
         lsagas(i) = .TRUE.
      enddo
c
c  --- initialize all of the tracers concs to zero to start off ---
c
      mvec4d = 0
      do i=1,ngrid
         mvec4d = mvec4d + DBLE(ncol(i)) * DBLE(nrow(i)) * DBLE(nlay(i))
      enddo
      mvec4d = mvec4d * DBLE(ntotsp)
      do idxsa=1,mvec4d
         ptconc(idxsa) = 0.
      enddo
      pttop = 0.
      do l=1,ntotsp
         do i=1,MXRECP
            conrcp(l,i) = 0.
         enddo
      enddo
c
c  --- set pointers for DDM-PM species
c
      if (lchem .AND. (aeropt.EQ.'CMU' .OR. aeropt.EQ.'CF')) then
         jdpmap(jdpSO2 ) = kso2
         jdpmap(jdpH2O2) = khpo_c
         jdpmap(jdpO3  ) = ko3
         jdpmap(jdpFOA ) = kfoa_c
         jdpmap(jdpMHP ) = kmhp_c
         jdpmap(jdpOHP ) = kohp_c
         jdpmap(jdpPAA ) = kpaa_c
         jdpmap(jdpOPA ) = kopa_c
         jdpmap(jdpN2O5) = kn2o5
         jdpmap(jdpHNO3) = khno3
         jdpmap(jdpNH3 ) = knh3
         jdpmap(jdpSULF) = ksulf
         jdpmap(jdpHCL ) = khcl
         jdpmap(jdpCG1 ) = kcg1
         jdpmap(jdpCG2 ) = kcg2
         jdpmap(jdpCG3 ) = kcg3
         jdpmap(jdpCG4 ) = kcg4
         jdpmap(jdpGLY ) = kgly_c
         jdpmap(jdpMGLY) = kmgly_c
         jdpmap(jdpGLYD) = kglyd_c
         jdpmap(jdpHO  ) = koh
         jdpmap(jdpNA  ) = kna
         jdpmap(jdpPSO4) = kpso4
         jdpmap(jdpPNH4) = kpnh4
         jdpmap(jdpPNO3) = kpno3
         jdpmap(jdpPCL ) = kpcl
         jdpmap(jdpSOA1) = ksoa1
         jdpmap(jdpSOA2) = ksoa2
         jdpmap(jdpSOA3) = ksoa3
         jdpmap(jdpSOA4) = ksoa4
         jdpmap(jdpSOPA) = ksopa
         jdpmap(jdpSOPB) = ksopb
         jdpmap(jdpPOA ) = kpoa
         jdpmap(jdpPH2O) = kph2o
         jdpmap(jdpFPRM) = kfprm
         jdpmap(jdpFCRS) = kfcrs
         do i = 1, NSDDMPM
            if ( jdpmap(i).eq.nspec+1 ) then
               jdpmap(i) = 0
            else
               jdpmap(i) = iptddm( jdpmap(i) )
            endif
         enddo
      endif
c
c  --- return to calling routine ---
c
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,A)') 'ERROR in SPECDDM:'
      write(iout,'(/,1X,A,I5)')
     &            'ERROR: Number of DDM families exceeds max: ',MXFDDM
      write(iout,'(1X,A,I5)')
     &        'You need room for at least this many families: ',nddmsp
      write(iout,'(1X,A)') 'Increase parameter MXFDDM in camx.prm'
      write(iout,'(1X,2A,I5)') 'Also increase parameter MXTRSP',
     &                     ' in camx.prm if it is less than ',ntotsp
      call camxerr()
c
 7001 continue
      write(iout,'(//,a)') 'ERROR in SPECDDM:'
      write(iout,'(/,1X,A,I5)')
     &                 'Number of tracer species exceeds max: ',MXTRSP
      write(iout,'(1X,A,I5)')
     &         'You need room for at least this many tracers: ',ntotsp
      write(iout,'(1X,A)') 'Increase parameter MXTRSP in camx.prm'
      call camxerr()
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
