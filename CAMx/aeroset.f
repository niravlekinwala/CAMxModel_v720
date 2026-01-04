      subroutine aeroset(dsec_i,ierr)
      use camxcom
      use grid
      use chmstry
      use filunit
c
c----CAMx v7.20 220430
c
c     AEROSET prepares for the AERO chemistry routines in CAMx
c
c      Copyright 1996 - 2022
c     Ramboll
c  
c     Modifications:
c        12/29/06   --bkoo--     Added species mapping for the updated SOA scheme
c        01/10/07   --bkoo--     Included camx_aero.inc
c                                Set pointers for aqueous PM
c        03/18/14   --bkoo--     Added pointers for benzene SOA
c        12/07/14   --bkoo--     Set pointers for VBS species
c        01/08/16   --bkoo--     Updated for SAPRC99->SAPRC07TC
c                                            Revised VBS
c        08/25/16   --bkoo--     Updated for new SOAP
c        09/16/16   --pk--       Added pointers for GLY, MGLY & GLYD
c        09/16/16   --bkoo--     Set WTMOL here (moved from CHMDAT.F)
c        10/16/17   --bkoo--     Added sfacJSOA; removed LSOAP2
c        11/22/17   --bkoo--     Updated WTMOL and separated Ca & Mg from CO3
c        02/21/20   --rlb--      Added CB6r5 (idmech=1)
c        03/24/21   --gy--       Added CB7 (idmech=7)
c
c     Input arguments: 
c        dsec_i              User supplied section cutpoints 
c 
c     Output arguments: 
c        ierr                Error flag
c            
c     Routines called: 
c        UPTIME
c            
c     Called by: 
c        READCHM
c 
      include 'camx.prm'
      include 'soap.inc'
      include 'vbs.inc'
      include 'section.inc'
      include 'section_aq.inc'
      include 'camx_aero.inc'
      include 'camx_aero_cmu.inc'
c
      real*8  dsec_i(nsecp1)
      integer ierr
c
      real    dtio
      integer nsteps, i, n
c     
c-----Entry point 
c
      ierr = 0
      lfrst = .true.
c
c-----Reset dtaero to hit output interval exactly
c
      dtio = amin1(60.,dtinp,dtems,dtout)
      nsteps = INT( 0.999*dtio/dtaero ) + 1
      dt_aero = dtio/nsteps
      do n = 1,ngrid
        time_aero(n) = time
        date_aero(n) = date
        call uptime(time_aero(n),date_aero(n),60.*dt_aero)
        aero_dt(n)  = 0.0
      enddo
c
c-----Set pointers for aqueous PM
c
      if (idmech.EQ.5) then ! SAPRC07TC
        khpo_c = kh2o2
        kfoa_c = kfacd
        kmhp_c = kcooh
        kpaa_c = kco3h
        kohp_c = krooh
        kopa_c = kro3h
        kgly_c = kgly
        kmgly_c = kmgly
        kglyd_c = kglyd
      elseif (idmech.EQ.6) then ! CB05
        khpo_c = kh2o2
        kfoa_c = kfacd
        kmhp_c = kmepx
        kpaa_c = kpacd
        kohp_c = krooh
        kopa_c = nspec+1
        kgly_c = nspec+1
        kmgly_c = kmgly
        kglyd_c = nspec+1
      elseif (idmech.EQ.1 .or. idmech.EQ.3 
     &         .or. idmech.EQ.4 .or. idmech.EQ.7 ) then ! CB6 or CB7
!bkoo - to be updated: Add ISPX to hydroperoxides for CB6?
        khpo_c = kh2o2
        kfoa_c = kfacd
        kmhp_c = kmepx
        kpaa_c = kpacd
        kohp_c = krooh
        kopa_c = nspec+1
        kgly_c = kgly
        kmgly_c = kmgly
        kglyd_c = kglyd
      endif
c
      if (lvbs .and. LVBSPREPROC) then
        kpap_c(0) = kpap0
        kpap_c(1) = kpap1
        kpap_c(2) = kpap2
        kpap_c(3) = kpap3
        kpap_c(4) = kpap4
        kpcp_c(0) = kpcp0
        kpcp_c(1) = kpcp1
        kpcp_c(2) = kpcp2
        kpcp_c(3) = kpcp3
        kpcp_c(4) = kpcp4
        kpfp_c(0) = kpfp0
        kpfp_c(1) = kpfp1
        kpfp_c(2) = kpfp2
        kpfp_c(3) = kpfp3
        kpfp_c(4) = kpfp4
      endif
c
c-----SOAP parameters based on old chamber data
c
!bk      y_h      = yh_1
!bk      y_l      = yl_1
!bk      csatS    = csatS_1
!bk      cstempS  = cstempS_1
!bk      deltahS  = deltahS_1
!bk      kpolya   = kpolya_1
!bk      kpolyb   = kpolyb_1
!bk      sfacJSOA = sfacJSOA_1
c
c-----SOAP parameters based on new chamber data (Hodzic et al., 2016, ACP)
c
      y_h      = yh_2
      y_l      = yl_2
      csatS    = csatS_2
      cstempS  = cstempS_2
      deltahS  = deltahS_2
      kpolya   = kpolya_2
      kpolyb   = kpolyb_2
      sfacJSOA = sfacJSOA_2
c
c-----Molecular weights for PM chemistry modules;
c     Used in AEROCHEM_CF, AEROCHEM_CMU, and DDM routines
c
      wtmol( 1) =  96.1 ! SO4
      wtmol( 2) =  18.  ! NH4
      wtmol( 3) =  62.  ! NO3
      wtmol( 4) =  35.5 ! Cl
      wtmol( 5) =  23.  ! Na
      wtmol( 6) =  39.1 ! K
      wtmol( 7) =  60.  ! CO3
      wtmol( 8) =  40.1 ! Ca
      wtmol( 9) =  24.3 ! Mg
      wtmol(10) =  55.8 ! Fe
      wtmol(11) =  54.9 ! Mn
      if (lvbs) then
        ksoac_c = kPBS0
        wtmol(12) = mwPBS0
      else
        ksoac_c = ksopb
        wtmol(12) = mwsopb
      endif
c
      if (aeropt .ne. 'CMU') return
c
c-----Check nbin against nsec parameters used in AERO modules
c
      if ( nbin .ne. nsec ) then
        write(iout,*) ' ERROR:  nbin and nsec must be equal!'
        write(iout,*) ' NBIN =',nbin, '; NSEC =',nsec
        write(iout,*) ' Set NSEC in SECTION.INC to ',nbin
        ierr = 1
        return
      endif
      if ( nsec .ne. nsect ) then
        write(iout,*) ' ERROR:  nsec and nsect must be equal!'
        write(iout,*) ' NSEC =',nsec, '; NSECT =',nsect
        write(iout,*) ' Set NSECT in SECTION_AQ.INC to ',nsec
        ierr = 1
        return
      endif
c
c-----The following specifies which AERO Modules to use 
c
      chaero = 'EQUI'
c      chaero = 'MADM'
c      chaero = 'HYBR'
      chaq = 'RADM'
c      chaq = 'VSRM' ! DO NOT USE THIS !!!
c      chaq = 'OVSR' ! DO NOT USE THIS !!!
c
c-----Set cwmin & tamin
c
      aqcwmin = cwmin
      aqtamin = tamin
c
c-----Calculate the initial diameters
c
c-----Check that user supplied section cutpoints are monotonically increasing
c
      do i = 2,nsecp1
        if ( dsec_i(i) .le. dsec_i(i-1) ) then
          write(iout,'(//,a)') 'ERROR in AEROSET:'
          write(iout,'(1X,2a)') 'Invalid section cut-points:',
     &                         ' they must be monotonically increasing!'
          write(iout,'(1X,a,/)') 'Input cut-points are :'
          do n = 1,nsecp1
            write(iout,'(4X,i2,1X,D9.3)') n,dsec_i(n)
          enddo
          call camxerr()
        endif
      enddo
c
c-----For the MOVING sectional approach dsec is the SECTIONAL DIAMETER
c
      dmin = dsec_i(1)
      dmax = dsec_i(nsecp1)
      do i = 1,nsecp1
        dsecf_c(i) = dsec_i(i)
      enddo

      write(idiag,*) ' '
      write(idiag,*) 'Particle section cut-points:'
      do i = 1,nsecp1
        write(idiag,'(1x,i2,1x,d9.3)') i,dsecf_c(i)
      enddo
      write(idiag,*) ' '
c
c-----Set moving diameters to logarithmic mean of fixed section diameters 
c     (tmg,01/25/02)
c
      do i = 1,nsec
        dsec_c(i) = sqrt(dsecf_c(i)*dsecf_c(i+1))
      enddo
c
c-----Additional pointers for RADM
c
      kso2_c   = kso2
      ko3_c    = ko3
      kn2o5_c  = kn2o5
      khno3_c  = khno3
      knh3_c   = knh3
      kh2so4_c = ksulf
      kpso4_c  = kpso4_1
      kpnh4_c  = kpnh4_1
      kpno3_c  = kpno3_1
      kna_c    = kna_1
      kpcl_c   = kpcl_1
c
c-----Additional pointers for SOAP
c
      kcg1_c   = kcg1
      kcg2_c   = kcg2
      kcg3_c   = kcg3
      kcg4_c   = kcg4
      ksoa1_c  = ksoa1_1
      ksoa2_c  = ksoa2_1
      ksoa3_c  = ksoa3_1
      ksoa4_c  = ksoa4_1
      ksopa_c  = ksopa_1
      ksopb_c  = ksopb_1
      ksoac_c  = ksopb_1
      kpoa_c   = kpoa_1
c
c-----Additional pointers for ISORROPIA
c
      khcl_c   = khcl
      kpec_c   = kpec_1
      kph2o_c  = kph2o_1
c
c-----Additional pointers for total PM
c
      kcrst_c  = kcrst_1
c
c-----Pointers for VSRM only
c
      kh2o2_c  = kh2o2
      kform_c  = kform
      khono_c  = khono
      koh_c    = koh
      kho2_c   = kho2
      kno3_c   = kno3
      kno_c    = kno
      kno2_c   = kno2
      kpan_c   = kpan
c
c-----Number of modeling species
c
      ngas_c   = ngas
      nspec_c  = nspec
c
      return
      end
