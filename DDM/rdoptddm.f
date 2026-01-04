c**** RDOPTDDM
c
      subroutine rdoptddm()
      use filunit
      use grid
      use chmstry
      use procan
      use tracer
      implicit none
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c     This routine loads all of the user options and flags for the
c     source apportionment algorithm.  This version is for the
c     DDM technology type, which has different options.
c
c      Copyright 1996 - 2022
c     Ramboll
c
c     Argument description:
c           none
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     03/23/99   --gwilson--    Original development
c     07/19/01   --gyarwood-    Initialize ngroup, nregin, lbndry
c     10/06/04   --cemery  -    Restructured for namelist input
c     07/16/07   --bkoo--       Revised for HDDM
c     06/11/08   --bkoo--       Added rate constant sensitivity
c     07/16/08   --bkoo--       Added DDM turn-off flag
c     07/23/18   --bkoo--       Removed memory allocation for lddmcalc
c     08/09/18   --bkoo--       Added rate term sensitivity
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'namelist.inc'
c
c-----------------------------------------------------------------------
c    Local Variables:
c-----------------------------------------------------------------------
c
      character*200 strtmp
      character*10  tmpnam
      integer i, j, k1, k2, k3, n
      integer numrxn,numreac,numprod,idxrxn,idxspc
      integer, allocatable :: reactmp(:,:), prodtmp(:,:)
      real,    allocatable :: coeftmp(:,:)
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c   --- if option is turned off, just return ---
c
      if( .NOT. lddm .AND. .NOT. lhddm ) goto 9999
c
c  --- check that switches are compatible ---
c
      if( naero .GT. 0 .AND. lhddm ) then
         write(iout,'(//,a)') 'ERROR in RDOPTDDM:'
         write(iout,*) 'Cannot use HDDM with Aerosol chemistry.'
         write(iout,*)
     &   '  Turn off HDDM or use a different chemical mechanism.'
         call camxerr()
      endif
!bk      if( naero .GT. 0 ) then
!bk         write(iout,'(//,a)') 'ERROR in RDOPTDDM:'
!bk         write(iout,*) 'Cannot use DDM with Aerosol chemistry'
!bk         write(iout,*) 'in this version of CAMx.'
!bk         write(iout,*)
!bk     &   '  Turn off DDM or run a gas-only chemical mechanism.'
!bk         call camxerr()
!bk      endif
c
c   --- get the number of Initial conditions species ---
c
      nicddm =  Number_of_IC_Species_Groups
      if( nicddm .GT. 0 ) then
         do i = 1,nicddm
           icddmsp(i) = IC_Species_Groups(i)
           call toupper( icddmsp(i) )
           call jstlft( icddmsp(i) )
         enddo
      endif
c
c   --- get the number of boundary conditions species ---
c
      lbndry = .FALSE.
      nbcddm = Number_of_BC_Species_Groups
      if( nbcddm .GT. 0 ) then
         do i = 1,nbcddm
           bcddmsp(i) = BC_Species_Groups(i)
           call toupper( bcddmsp(i) )
           call jstlft( bcddmsp(i) )
         enddo
c
c   --- flag for stratifying the boundary by edge ---
c
         lbndry = DDM_Stratify_Boundary
      endif
c
c   --- get the number of emissions species ---
c
      ngroup = 0
      nregin = 0
      nemddm = Number_of_EM_Species_Groups
      if( nemddm .GT. 0 ) then
         do i=1,nemddm
           emddmsp(i) = Emis_Species_Groups(i)
           call toupper( emddmsp(i) )
           call jstlft( emddmsp(i) )
         enddo
c
c   --- number of source regions ---
c
        nregin = DDM_Number_of_Source_Regions
        if( nregin .EQ. 0 ) then
           write(iout,'(//,A)') 'ERROR in RDOPTDDM:'
           write(iout,'(/,1X,2A,/,A)') 'When requesting ',
     &                'emissions species for DDM you must ',
     &                'provide at least 1 emissions region.'
           call camxerr()
        endif
c
c   --- number of source emissions groupings ---
c
        ngroup = DDM_Number_of_Source_Groups
        if( ngroup .EQ. 0 ) then
           write(iout,'(//,A)') 'ERROR in RDOPTDDM:'
           write(iout,'(/,1X,2A,/,A)') 'When requesting ',
     &                'emissions species for DDM you must ',
     &                'provide at least 1 emissions group.'
           call camxerr()
        endif
      endif
c
c   --- leftover and 3D average is always false with DDM ---
c
      leftovr_area = .FALSE.
      lsa_3davrg = .FALSE.
c
c   --- check for array overflow ---
c
      if( ngroup .GT. MXTEMF-1 ) then
         write(iout,'(//,A)') 'ERROR in RDOPTDDM:'
         write(iout,'(/,1X,A,I4,A)') 
     &                  'Number of source groupings ',ngroup,
     &                  ' exceeds maximum.  Increase MXTEMF.'
         call camxerr()
      endif
      if( ngroup .LT. 0 ) then
         write(iout,'(//,A)') 'ERROR in RDOPTDDM:'
         write(iout,'(1X,A,I4,A)') 
     &                  'Number of emissions groups ',ngroup,
     &                                         ' is invalid.'
         call camxerr()
      endif
c
c   --- number of timing realeases is always zero for DDM ---
c
      ntrtim = 0
c
c   --- number of rate constant sensitivity groups ---
c
      nrateddm = Number_of_Rate_Const_Groups
      allocate( rateddm(nrateddm) )
      allocate( iprate(0:nreact, nrateddm ) )
      iprate = 0
      do i = 1, nrateddm
        strtmp = ADJUSTL(Rate_Const_Groups(i))
        j = INDEX(strtmp,':')
        if (j.eq.0) then
          write(iout,'(//,A)') 'ERROR in RDOPTDDM:'
          write(iout,'(1X,A,I3,A)')
     &        'Delimiter (:) is not found in Rate_Const_Groups(',i,').'
          call camxerr()
        endif
        rateddm(i) = strtmp(:j-1)
        do
          strtmp = strtmp(j+1:)
          j = INDEX(strtmp,',')
          if (j.eq.0) then
            j = LEN_TRIM(strtmp) + 1
            if (j.eq.1) EXIT
          endif
          iprate(0,i) = iprate(0,i) + 1
          if (iprate(0,i).gt.nreact) then
            write(iout,'(//,A)') 'ERROR in RDOPTDDM:'
            write(iout,'(1X,A,I3,A)')
     &            'Too many rxns (>nreact) in Rate_Const_Groups(',i,').'
            call camxerr()
          endif
          if ( VERIFY( TRIM( ADJUSTL( strtmp(:j-1) ) ), '1234567890' )
     &                                                    .ne. 0 ) then
            write(iout,'(//,A)') 'ERROR in RDOPTDDM:'
            write(iout,'(1X,A,I3,A)')
     &                 'Invalid rxn index in Rate_Const_Groups(',i,').'
            call camxerr()
          endif
          read(strtmp(:j-1),'(I5)') iprate(iprate(0,i),i)
          if ( iprate(iprate(0,i),i).lt.1 .or.
     &         iprate(iprate(0,i),i).gt.nreact ) then
            write(iout,'(//,A)') 'ERROR in RDOPTDDM:'
            write(iout,'(1X,A,I3,A)')
     &            'Out-of-range rxn index in Rate_Const_Groups(',i,').'
            call camxerr()
          endif
        enddo
        if (iprate(0,i).eq.0) then
          write(iout,'(//,A)') 'ERROR in RDOPTDDM:'
          write(iout,'(1X,A,I3,A)')
     &                   'No rxn is found in Rate_Const_Groups(',i,').'
          call camxerr()
        endif
      enddo
c
c   --- number of rate term sensitivity groups ---
c
      ntermddm = Number_of_Rate_Term_Groups
      allocate( termddm(ntermddm) )
      termddm = ' '
      allocate( ipterm(0:ngas, nreact, ntermddm ) )
      ipterm = 0
      allocate( wfterm(  ngas, nreact, ntermddm ) )
      wfterm = 1.0
      allocate( reactmp(ngas,nreact) )
      allocate( prodtmp(ngas,nreact) )
      allocate( coeftmp(ngas,nreact) )
      do i = 1, ntermddm
        numrxn = 0
        reactmp = 0
        prodtmp = 0
        coeftmp = 1.0
        do n = 1, MXNAM
          strtmp = ADJUSTL(Rate_Term_Groups(i,n))
          if ( LEN_TRIM(strtmp).eq.0 ) CYCLE
          j = INDEX(strtmp,':')
          if (j.eq.0) then
            write(iout,'(//,A)') 'ERROR in RDOPTDDM:'
            write(iout,'(1X,A,I3,A)')
     &          'Delimiter (:) is not found in Rate_Term_Groups(',i,').'
            call camxerr()
          endif
          if ( LEN_TRIM(termddm(i)).eq.0 ) then
            termddm(i) = strtmp(:j-1)
          else
            if ( termddm(i).ne.strtmp(:j-1) ) then
              write(iout,'(//,A)') 'ERROR in RDOPTDDM:'
              write(iout,'(1X,A,I3,A)')
     &                   'Inconsistent name of Rate_Term_Groups(',i,').'
              call camxerr()
            endif
          endif
          do
            strtmp = strtmp(j+1:)
            j = INDEX(strtmp,',')
            if (j.eq.0) then
              j = LEN_TRIM(strtmp) + 1
              if (j.eq.1) EXIT
            endif
            k1 = INDEX(strtmp(:j-1),'[R')
            k2 = INDEX(strtmp(:j-1),'[P')
            k3 = INDEX(strtmp(:j-1),']')
            if (k1+k2.eq.0 .or. k3.eq.0 .or. k3.lt.k1+k2+2) then
              write(iout,'(//,A)') 'ERROR in RDOPTDDM:'
              write(iout,'(1X,A,I3,A)')
     &                     'Invalid keyword in Rate_Term_Groups(',i,').'
              call camxerr()
            endif
            if ( VERIFY( TRIM( ADJUSTL( strtmp(:k1+k2-1) ) ),
     &                                     '1234567890' ) .ne. 0 ) then
              write(iout,'(//,A)') 'ERROR in RDOPTDDM:'
              write(iout,'(1X,A,I3,A)')
     &                   'Invalid rxn index in Rate_Term_Groups(',i,').'
              call camxerr()
            endif
            read(strtmp(:k1+k2-1),*) idxrxn
            if ( idxrxn.lt.1 .or. idxrxn.gt.nreact ) then
              write(iout,'(//,A)') 'ERROR in RDOPTDDM:'
              write(iout,'(1X,A,I3,A)')
     &              'Out-of-range rxn index in Rate_Term_Groups(',i,').'
              call camxerr()
            endif
            numrxn = numrxn + 1
            ipterm(0,idxrxn,i) = 1
            tmpnam = strtmp(k3+1:j-1)
            do idxspc = 1, ngas
              if (tmpnam.eq.spname(idxspc)) EXIT
            enddo
            if (idxspc.gt.ngas) then
              write(iout,'(//,A)') 'ERROR in RDOPTDDM:'
              write(iout,'(1X,A,I3,A,I3,A)')
     &          'Invalid spc name in rxn(',idxrxn,
     &                                   ') of Rate_Term_Groups(',i,').'
              call camxerr()
            endif
            if (k1.gt.0) then ! Reactant
              reactmp(idxspc,idxrxn) = 1
            else              ! Product
              prodtmp(idxspc,idxrxn) = 1
            endif
            if (strtmp(k1+k2+2:k1+k2+2).eq.'#') then
              if ( k2.ne.0 ) then ! Check if a weighting factor is assigned to a product; may be allowed later
                write(iout,'(//,A)') 'ERROR in RDOPTDDM:'
                write(iout,'(1X,A)')
     &                        'Products must not have weighting factors'
                call camxerr()
              endif
              if ( VERIFY( TRIM( ADJUSTL( strtmp(k1+k2+3:k3-1) ) ),
     &                                '1234567890.+-Ee' ) .ne. 0 ) then
                write(iout,'(//,A)') 'ERROR in RDOPTDDM:'
                write(iout,'(1X,A,I3,A)')
     &                 'Invalid coefficient in Rate_Term_Groups(',i,').'
                call camxerr()
              endif
              read(strtmp(k1+k2+3:k3-1),*) coeftmp(idxspc,idxrxn)
            endif
          enddo
        enddo ! loop n
        if (numrxn.eq.0) then
          write(iout,'(//,A)') 'ERROR in RDOPTDDM:'
          write(iout,'(1X,A,I3,A)')
     &                     'No rxn is found in Rate_Term_Groups(',i,').'
          call camxerr()
        endif
        do idxrxn = 1, nreact
          numreac = 0
          do idxspc = 1, ngas
            if ( reactmp(idxspc,idxrxn).eq.0 ) CYCLE
            numreac = numreac + 1
            ipterm(1+numreac,idxrxn,i) = idxspc
            wfterm(1+numreac,idxrxn,i) = coeftmp(idxspc,idxrxn)
          enddo
          ipterm(1,idxrxn,i) = numreac
          numprod = 0
          do idxspc = 1, ngas
            if ( prodtmp(idxspc,idxrxn).eq.0 ) CYCLE
            numprod = numprod + 1
            if ( 2+numreac+numprod.gt.ngas ) then
              write(iout,'(//,A)') 'ERROR in RDOPTDDM:'
              write(iout,'(1X,A,I3,A,I3,A)')
     &          'Too many spc selected in rxn(',idxrxn,
     &                                  ') of Rate_Const_Groups(',i,').'
              call camxerr()
            endif
            ipterm(2+numreac+numprod,idxrxn,i) = idxspc
          enddo
          ipterm(2+numreac,idxrxn,i) = numprod
        enddo ! loop idxrxn
      enddo ! loop i
      deallocate( reactmp )
      deallocate( prodtmp )
      deallocate( coeftmp )
c
c   --- number of HDDM sensitivity groups ---
c
      nhddm = Number_of_HDDM_Sens_Groups
      if ( .NOT. lhddm .AND. nhddm .GT. 0 ) then
        write(iout,'(//,A)') 'ERROR in RDOPTDDM:'
        write(iout,'(1X,A)')
     &              'Number_of_HDDM_Sens_Groups must be 0 if not HDDM.'
        call camxerr()
      endif
      allocate( hddmsp(2,nhddm) )
      do i = 1,nhddm
        do j = 1, 2
          hddmsp(j,i) = HDDM_parameters(i,j)
          call toupper( hddmsp(j,i) )
          call jstlft( hddmsp(j,i) )
          if ( hddmsp(j,i)(:4).eq.'TERM' ) then !bk_dbg=>
            write(iout,'(//,A)') 'ERROR in RDOPTDDM:'
            write(iout,'(1X,A)')
     &        'Rate term sensitivity is currently limited to 1st-order.'
            call camxerr()
          endif                                 !bk_dbg<=
        enddo
      enddo
c
c   --- set DDM turn-off flags ---
c
!bk - now allocated in ALLOC_LDDMCALC
!bk      allocate( lddmcalc(ngrid) )
      do i = 1, ngrid
        lddmcalc(i) = DDM_Calc_Grid(i)
        if ( lddmcalc(i) .OR. nbcddm.EQ.0 ) CYCLE
        write(iout,'(//,A)') 'ERROR in RDOPTDDM:'
        write(iout,'(1X,A)')
     &               'Flexi-DDM is not allowed with sensitivity to BC.'
        call camxerr()
      enddo
c
c   --- dummy array for tracer wetdep field ---
c
      allocate( ptwetfld(1) )
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
