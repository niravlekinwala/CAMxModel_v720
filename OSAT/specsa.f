c**** SPECSA
c
      subroutine specsa(idate,begtim,jdate,endtim)
      use filunit
      use tracer
      use grid
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routine sets up the species names and pointers into the species
c   for all of the tracer species.  Pointers will be set up for both the
c   concentration array and the emissions array.
c
c      Copyright 1996 - 2022
c     Ramboll
c
c     Argument descriptions:
c      Inputs:
c        idate  I   date of the beginning of the simulation (YYJJJ)
c        begtim R   hour of the begining of simulation
c        jdate  I   date of the ending of the simulation (YYJJJ)
c        endtim R   hour of the endng of simulation
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     05/26/96   --gwilson--    Original development
c     12/12/97   --gwilson--    Fixed bug in initializing the timing
c                               tracers
c     11/06/01   --cemery--     Input dates are now Julian
c     12/29/06   --bkoo--       Revised for the updated SOA scheme
c     11/4/09    -cemery-       Removed reference to input top conc array
c     11/12/09   --gwilson--    Added initialization of factor for
c                               applying new type of top boundary 
c     03/18/14   --bkoo--       Added tracer for SOAH
c     08/25/16   --bkoo--       Updated for new SOAP
c     01/06/19   --cemery--     Added PFE/PMN/PMG/PK/PCA/PAL/PSI/PTI
c     01/09/19   --cemery-      Added DMS
c     
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      integer   idate
      integer   jdate
      real      begtim
      real      endtim     
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      character*10 name
      integer*8    mvsa3d, idxsa
      integer      ibegdt, ienddt
      integer      ncount, ioff, idtnow, nhours, numic_spcs, numbc_spcs
      integer      numrt_spcs, num_icbc_specs, i, j, k, l
      real         molwt_cls(MXALCLS), timnow, btim, etim
      logical      lgasflg(MXALCLS), lwrtcls(MXALCLS), lfound_class
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
      ibegdt = idate
      btim = begtim/100.
      ienddt = jdate
      etim = endtim/100.
c
c   --- set flags for whether this class should be 
c       written to average file ---
c
      do i=1,MXALCLS
        lwrtcls(i) = .FALSE.
        molwt_cls(i) = 1.0
      enddo
      lwrtcls(ITRO3N) = .TRUE.
      lwrtcls(ITRO3V) = .TRUE.
      lwrtcls(ITRPS4) = .TRUE.
      lwrtcls(ITRPN3) = .TRUE.
      lwrtcls(ITRPN4) = .TRUE.
      lwrtcls(ITRPO1) = .TRUE.
      lwrtcls(ITRPO2) = .TRUE.
      lwrtcls(ITRPO3) = .TRUE.
      lwrtcls(ITRPO4) = .TRUE.
      lwrtcls(ITRPPA) = .TRUE.
      lwrtcls(ITRPPB) = .TRUE.
      lwrtcls(ITRPEC) = .TRUE.
      lwrtcls(ITRPOA) = .TRUE.
      lwrtcls(ITRPFC) = .TRUE.
      lwrtcls(ITRPFN) = .TRUE.
      lwrtcls(ITRPCC) = .TRUE.
      lwrtcls(ITRPCS) = .TRUE.
      lwrtcls(ITRPFE) = .TRUE.
      lwrtcls(ITRPMN) = .TRUE.
      lwrtcls(ITRPMG) = .TRUE.
      lwrtcls(ITRPK)  = .TRUE.
      lwrtcls(ITRPCA) = .TRUE.
      lwrtcls(ITRPAL) = .TRUE.
      lwrtcls(ITRPSI) = .TRUE.
      lwrtcls(ITRPTI) = .TRUE.
      lwrtcls(ITRHG0) = .TRUE.
      lwrtcls(ITRHG2) = .TRUE.
      lwrtcls(ITRPHG) = .TRUE.
c
c  --- define the names of the tracer clasess ---
c
      clsnam(ITRVOC) = 'VOC'
      lgasflg(ITRVOC) = .TRUE.
      molwt_cls(ITRVOC) = 14.5
c
      clsnam(ITRO3N) = 'O3N'
      lgasflg(ITRO3N) = .TRUE.
      molwt_cls(ITRO3N) = 48.0
c
      clsnam(ITRO3V) = 'O3V'
      lgasflg(ITRO3V) = .TRUE.
      molwt_cls(ITRO3V) = 48.0
c
      clsnam(ITROON) = 'OON'
      lgasflg(ITROON) = .TRUE.
      molwt_cls(ITROON) = 48.0
c
      clsnam(ITROOV) = 'OOV'
      lgasflg(ITROOV) = .TRUE.
      molwt_cls(ITROOV) = 48.0
c
      clsnam(ITRSO2) = 'SO2'
      lgasflg(ITRSO2) = .TRUE.
      molwt_cls(ITRSO2) = 64.0
c
      clsnam(ITRPS4) = 'PS4'
      lgasflg(ITRPS4) = .FALSE.
c
      clsnam(ITRDMS) = 'DMS'
      lgasflg(ITRDMS) = .TRUE.
      molwt_cls(ITRDMS) = 62.1
c
      clsnam(ITRNIT) = 'NIT'
      lgasflg(ITRNIT) = .TRUE.
      molwt_cls(ITRNIT) = 30.0
c
      clsnam(ITRRGN) = 'RGN'
      lgasflg(ITRRGN) = .TRUE.
      molwt_cls(ITRRGN) = 46.0
c
      clsnam(ITRTPN) = 'TPN'
      lgasflg(ITRTPN) = .TRUE.
      molwt_cls(ITRTPN) = 121.0
c
      clsnam(ITRNTR) = 'NTR'
      lgasflg(ITRNTR) = .TRUE.
      molwt_cls(ITRNTR) = 147.0
c
      clsnam(ITRHN3) = 'HN3'
      lgasflg(ITRHN3) = .TRUE.
      molwt_cls(ITRHN3) = 63.0
c
      clsnam(ITRPN3) = 'PN3'
      lgasflg(ITRPN3) = .FALSE.
c
      clsnam(ITRNH3) = 'NH3'
      lgasflg(ITRNH3) = .TRUE.
      molwt_cls(ITRNH3) = 17.0
c
      clsnam(ITRPN4) = 'PN4'
      lgasflg(ITRPN4) = .FALSE.
c
      clsnam(ITRARO) = 'ARO'
      lgasflg(ITRARO) = .TRUE.
      molwt_cls(ITRARO) = 99.2
c
      clsnam(ITRISP) = 'ISP'
      lgasflg(ITRISP) = .TRUE.
c
      clsnam(ITRTRP) = 'TRP'
      lgasflg(ITRTRP) = .TRUE.
      molwt_cls(ITRTRP) = 136.2
c
      clsnam(ITRSQT) = 'SQT'
      lgasflg(ITRSQT) = .TRUE.
      molwt_cls(ITRSQT) = 204.4
c
      clsnam(ITRCG1) = 'CG1'
      lgasflg(ITRCG1) = .TRUE.
      molwt_cls(ITRCG1) = 150.0
c
      clsnam(ITRCG2) = 'CG2'
      lgasflg(ITRCG2) = .TRUE.
      molwt_cls(ITRCG2) = 150.0
c
      clsnam(ITRCG3) = 'CG3'
      lgasflg(ITRCG3) = .TRUE.
      molwt_cls(ITRCG3) = 180.0
c
      clsnam(ITRCG4) = 'CG4'
      lgasflg(ITRCG4) = .TRUE.
      molwt_cls(ITRCG4) = 180.0
c
      clsnam(ITRPO1) = 'PO1'
      lgasflg(ITRPO1) = .FALSE.
c
      clsnam(ITRPO2) = 'PO2'
      lgasflg(ITRPO2) = .FALSE.
c
      clsnam(ITRPO3) = 'PO3'
      lgasflg(ITRPO3) = .FALSE.
c
      clsnam(ITRPO4) = 'PO4'
      lgasflg(ITRPO4) = .FALSE.
c
      clsnam(ITRPPA) = 'PPA'
      lgasflg(ITRPPA) = .FALSE.
c
      clsnam(ITRPPB) = 'PPB'
      lgasflg(ITRPPB) = .FALSE.
c
      clsnam(ITRPEC) = 'PEC'
      lgasflg(ITRPEC) = .FALSE.
c
      clsnam(ITRPOA) = 'POA'
      lgasflg(ITRPOA) = .FALSE.
c
      clsnam(ITRPFC) = 'PFC'
      lgasflg(ITRPFC) = .FALSE.
c
      clsnam(ITRPFN) = 'PFN'
      lgasflg(ITRPFN) = .FALSE.
c
      clsnam(ITRPCC) = 'PCC'
      lgasflg(ITRPCC) = .FALSE.
c
      clsnam(ITRPCS) = 'PCS'
      lgasflg(ITRPCS) = .FALSE.
c
      clsnam(ITRPFE) = 'PFE'
      lgasflg(ITRPFE) = .FALSE.
c
      clsnam(ITRPMN) = 'PMN'
      lgasflg(ITRPMN) = .FALSE.
c
      clsnam(ITRPMG) = 'PMG'
      lgasflg(ITRPMG) = .FALSE.
c
      clsnam(ITRPK) = 'PK_'
      lgasflg(ITRPK) = .FALSE.
c
      clsnam(ITRPCA) = 'PCA'
      lgasflg(ITRPCA) = .FALSE.
c
      clsnam(ITRPAL) = 'PAL'
      lgasflg(ITRPAL) = .FALSE.
c
      clsnam(ITRPSI) = 'PSI'
      lgasflg(ITRPSI) = .FALSE.
c
      clsnam(ITRPTI) = 'PTI'
      lgasflg(ITRPTI) = .FALSE.
c
      clsnam(ITRHG0) = 'HG0'
      lgasflg(ITRHG0) = .TRUE.
      molwt_cls(ITRHG0) = 200.6
c
      clsnam(ITRHG2) = 'HG2'
      lgasflg(ITRHG2) = .TRUE.
      molwt_cls(ITRHG2) = 253.1
c
      clsnam(ITRPHG) = 'PHG'
      lgasflg(ITRPHG) = .FALSE.
c
c  --- calculate the beginning of the various tracer types ---
c      there will be (ngroup+1) if there is an extra group for the 
c      "leftover" group  ----
c
      if( lsa_ioric ) then
          do icls=1,ntrcls
            lfound_class = .FALSE.
            do ispc=1,num_ioric
               if( spc_ioric(ispc)(1:3) .EQ. clsnam(idxcls(icls)) ) 
     &                                        lfound_class = .TRUE.
            enddo
            if( .NOT. lfound_class ) goto 7000
          enddo
          numic_spcs = ncls_ioric
      else
          numic_spcs = 1
          ncls_ioric = 1
      endif
      if( lsa_iorbc ) then
          do icls=1,ntrcls
            lfound_class = .FALSE.
            do ispc=1,num_iorbc
               if( spc_iorbc(ispc)(1:3) .EQ. clsnam(idxcls(icls)) ) 
     &                                        lfound_class = .TRUE.
            enddo
            if( .NOT. lfound_class ) goto 7001
          enddo
          numbc_spcs = ncls_iorbc
          numtc_spcs = 0
          if( .NOT. lsa_iortc ) numtc_spcs = 1
      else
          numtc_spcs = 0
          if( lbndry ) then
              numbc_spcs = 5 
          else
              numbc_spcs = 1 
          endif
      endif
      num_icbc_specs = numic_spcs + numbc_spcs + numtc_spcs
      if( ngroup .EQ. 0 ) then
          ncount = nregin
      else
          if( leftovr_area ) then
             ncount = (ngroup + 1) * nregin
          else
             ncount = ngroup * nregin
          endif
      endif
c
c  --- set the flag for gaseous species ---
c
      do i=1,ntrcls
        do j=iptcls(i),nptcls(i)
           lsagas(j) = lgasflg(idxcls(i))
           sa_mole_weight(j) = molwt_cls(idxcls(i))
        enddo
      enddo
c
c   --- load the species from the initial conditions file ---
c
      if( lsa_ioric ) then
         do icls=1,ntrcls
           ioff = 0
           do ispc=1,num_ioric
             if( spc_ioric(ispc)(1:3) .EQ. clsnam(idxcls(icls)) ) then
                ptname(iptcls(icls)+ioff) = spc_ioric(ispc)
                ptop_fac(iptcls(icls)+ioff) = 0.0
                lsagas(iptcls(icls)+ioff) = lgasflg(idxcls(icls))
                sa_mole_weight(iptcls(icls)+ioff) = molwt_cls(idxcls(icls))
                ioff = ioff + 1
             endif
           enddo
         enddo
c
c   --- set the names for the initial condition tracers ---
c
      else
         do i=1,ntrcls
            ptname(iptcls(i)) = clsnam(idxcls(i))//'000IC  '
            ptop_fac(iptcls(i)) = 0.0
         enddo
         num_ioric = 1
      endif
c
c   --- load the species from the boundary conditions file ---
c
      if( lsa_iorbc ) then
         do icls=1,ntrcls
           ioff = numic_spcs
           do ispc=1,num_iorbc
             if( spc_iorbc(ispc)(1:3) .EQ. clsnam(idxcls(icls)) ) then
                ptname(iptcls(icls)+ioff) = spc_iorbc(ispc)
                lsagas(iptcls(icls)+ioff) = lgasflg(idxcls(icls))
                sa_mole_weight(iptcls(icls)+ioff) = molwt_cls(idxcls(icls))
                if( lsa_iortc ) then 
                    ptop_fac(iptcls(icls)+ioff) = 1.0
                else
                    ptop_fac(iptcls(icls)+ioff) = 0.0
                endif
                ioff = ioff + 1
             endif
           enddo
           if( .NOT. lsa_iortc ) then
             ptname(iptcls(icls)+ioff) = clsnam(idxcls(icls))//'TOPBC  '
             ptop_fac(iptcls(icls)+ioff) = 1.0
           endif
         enddo
c
c   --- if stratifying by boundary there will be 5 boundary condition
c       tracers, otherwise there will be only one ---
c
      else
          if( lbndry ) then
              do i=1,ntrcls
                 ptname(iptcls(i) + numic_spcs-1 + IDXBNT) = clsnam(idxcls(i))//'NTHBC  '
                 ptname(iptcls(i) + numic_spcs-1 + IDXBES) = clsnam(idxcls(i))//'ESTBC  '
                 ptname(iptcls(i) + numic_spcs-1 + IDXBST) = clsnam(idxcls(i))//'STHBC  '
                 ptname(iptcls(i) + numic_spcs-1 + IDXBWS) = clsnam(idxcls(i))//'WSTBC  '
                 ptname(iptcls(i) + numic_spcs-1 + IDXBTP) = clsnam(idxcls(i))//'TOPBC  '
                 ptop_fac(iptcls(i) + numic_spcs-1 + IDXBNT) = 0.0
                 ptop_fac(iptcls(i) + numic_spcs-1 + IDXBES) = 0.0
                 ptop_fac(iptcls(i) + numic_spcs-1 + IDXBST) = 0.0
                 ptop_fac(iptcls(i) + numic_spcs-1 + IDXBWS) = 0.0
                 ptop_fac(iptcls(i) + numic_spcs-1 + IDXBTP) = 1.0
              enddo
          else
              do i=1,ntrcls
                 ptname(iptcls(i)+numic_spcs) = clsnam(idxcls(i))//'000BC  '
                 ptop_fac(iptcls(i)+numic_spcs) = 1.0
              enddo
          endif
      endif
c
c  --- construct the tracer names and put into names array ---
c
      if( ngroup .EQ. 0 ) then
          ioff = num_icbc_specs
          do i=1,nregin
             do k=1,ntrcls
                 write(name,'(A,I3.3,I3.3)') clsnam(idxcls(k)),1,i
                 ptname(iptcls(k)+ioff) = name
                 ptop_fac(iptcls(k)+ioff) = 0.0
             enddo
             ioff = ioff + 1
          enddo
      else
          ioff = num_icbc_specs
          if( leftovr_area ) then
             ncount = ngroup + 1
          else
             ncount = ngroup 
          endif
          do j=1,ncount
             do i=1,nregin
                do k=1,ntrcls
                   write(name,'(A,I3.3,I3.3)') clsnam(idxcls(k)),j,i
                   ptname(iptcls(k)+ioff) = name
                   ptop_fac(iptcls(k)+ioff) = 0.0
                enddo
                ioff = ioff + 1
             enddo
          enddo
      endif
c
c  --- calculate the number of timing tracers there will be and put
c      the names into the names array ---
c
      ntotsp = ipttim - 1 
      if( ntrtim .GT. 0 ) then
        if( etim .EQ. 0. ) then
            etim = 24.
            ienddt = ienddt - 1
        endif 
        timnow = btim
        idtnow = ibegdt
        nhours = (ienddt-ibegdt)*24 + INT( etim - btim ) 
        npttim = 1
        do i=1,nhours
           if( MOD( INT(timnow), 24/ntrtim ) .EQ. 0 .OR. i .EQ. 1) then
              do j=1,nregin
                  write(name,'(A,I3.3,I2.2,I3.3)') 'I',MOD(idtnow,1000),
     &                                                  INT(timnow),j
                  ptname(ntotsp+1) = name
                  ptop_fac(ntotsp+1) = 0.0
                  write(name,'(A,I3.3,I2.2,I3.3)') 'D',MOD(idtnow,1000),
     &                                                  INT(timnow),j
                  ptname(ntotsp+2) = name
                  ptop_fac(ntotsp+2) = 0.0
                  npttim = npttim + 2
                  ntotsp = ntotsp + 2
              enddo
           endif
           timnow = timnow + 1.0
           if( timnow .EQ. 24.0 ) then
               timnow = 0.
               idtnow = idtnow + 1
           endif
        enddo
      endif
c
c  --- initialize all of the tracers concs to zero to start off ---
c
      mvsa3d = 0
      do i=1,ngrid
         mvsa3d = mvsa3d + DBLE(ncol(i)) * DBLE(nrow(i)) * DBLE(nlay(i))
      enddo
      mvsa3d = DBLE(mvsa3d) * DBLE(ntotsp)
      do idxsa=1,mvsa3d
        ptconc(idxsa) = 0.
      enddo
      do l=1,ntotsp
         lsamap(l) = l
         do i=1,MXRECP
            conrcp(l,i) = 0.        
         enddo
      enddo
c
c  --- set the flag for outputting the species to average file
c      to true automatically for tracer species ---
c
      if( lallout ) then
         do i=1,ntotsp
            loutsa(i) = .TRUE.
         enddo
      else
         do i=1,ntotsp
            loutsa(i) = .FALSE.
         enddo
         do j=1,ntrcls
           if( lwrtcls(idxcls(j)) ) then
              do i=iptcls(j),nptcls(j)
                 loutsa(i) = .TRUE.
              enddo
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
      write(iout,'(//,a)') 'ERROR in SPECSA:'
      write(iout,'(2A)',ERR=9999) 'You requested a tracer species ',
     &       'that is not included in your SA IC file.'
      write(iout,'(2A)',ERR=9999) 'Either disable the class or ',
     &       'supply the data in the SA IC file.'
      write(iout,'(10X,2A)',ERR=9999) 'SA species: ',clsnam(idxcls(icls))
      call camxerr()
c
 7001 continue
      write(iout,'(//,a)') 'ERROR in SPECSA:'
      write(iout,'(2A)',ERR=9999) 'You requested a tracer species ',
     &       'that is not included in your SA BC file.'
      write(iout,'(2A)',ERR=9999) 'Either disable the class or ',
     &       'supply the data in the SA BC file.'
      write(iout,'(10X,2A)',ERR=9999) 'SA species: ',clsnam(idxcls(icls))
      call camxerr()
c
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
