      subroutine nodes_init(numprocs,iproc_id)
      use grid
      use chmstry
      use filunit
      use o3colmap
      use bndary
      use camxfld
      use camxcom
      use pigsty
      use ptemiss
      use procan
      use rtracchm
      use rtcmcchm
      use tracer
      use node_mod
c
      implicit none
c
c----CAMx v7.20 220430
c
c     This routine initializes all of the data for the compute
c     nodes when running in MPI mode.  All of the I/O is
c     handled by the master process (ID = 0). All other processes
c     will do computation but need to have data structures 
c     allocated and filled through message passing from the
c     master node.
c     
c
c      Copyright 1996 - 2022
c     Ramboll 
c
c
c     Modifications:
c       12/15/08    Added code to handle averaging of radicals
c       03/15/09    Added code for deposition output for tracers
c       10/29/09    Added code for RTRAC surface model
c       11/4/09     Removed input top concentrations
c       01/04/11    Revised for new met input format
c       11/06/12    Fixed Wall of Cells receptors for MPI
c       11/06/12    Pass receptor information for RTRAC
c       04/30/13    Added surface model
c       05/24/13    Pass RTCMC photolysis/dep pointers
c       09/02/14    Added subgrid convective model
c       12/07/14    Added VBS parameters
c       01/08/16    Updated for Revised VBS & Removal of IEH solver
c       03/01/15    Added partial source area map
c       08/25/16    Updated for new SOAP
c       09/02/16    Renamed pointers to het rxns (-> ihetrxn)
c       09/16/16    Added WTMOL & pointers for in-cloud SOA species
c       10/16/17    Added sfacJSOA; removed LSOAP2
c       07/23/18    Added Bi-Di NH3 drydep flag
c                   Added npartial
c       08/09/18    Added variables for rate term sensitivity
c       07/23/19    Staggered wind flag is now grid-specific
c       11/17/20    RTCMC dry/wet deposition is now consistent with
c                   RTRAC
c
c     Input arguments:
c        numprocs            the number of processes
c        iproc_id            process number for current process
c
c     Output arguments:
c        none
c
c     Routines called:
c     
c     Called by:
c        CAMx
c
      include 'camx.prm'
      include 'flags.inc'
      include 'deposit.inc'
      include 'camx_aero.inc'
      include 'soap.inc'
      include 'vbs.inc'
      include 'ddmchm.inc'
      include 'rtracsrf.inc'
      include 'mpif.h'
c
c----Argument declarations
c
      integer numprocs
      integer iproc_id
c     
c-----Local variables
c
      integer ncola, nrowa, nlaya, igrd, ispc, ilay, ierr, i
      integer nfilcnt, mvec2d, mvec3d, max_emissfiles, max_ptfiles
c     
c-----Entry point
c
      if( .NOT. lmpi .AND. iproc_id .GT. 0 ) return
c
c   --- calculate the number of bytes to send ---
c
      nrowa = maxval(nrow(1:ngrid) )
      ncola = maxval(ncol(1:ngrid) )
      nlaya = maxval(nlay(1:ngrid) )
c
c  --- find the number of output files --
c
      nfilcnt = 1
      if( ngrid .GT. 1 ) nfilcnt = 2
c
c  --- send the other variables in the filunit include file ---
c
      call nodes_pass(icur_unit,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(imass,    1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(iconc,    nfilcnt,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(ifconc,   nfilcnt,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(ipig,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(ichem,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(iphot,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(iic,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(is_netcdf_iic,1,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(ibc,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(is_netcdf_ibc,1,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(itc,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(is_netcdf_itc,1,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(io3col,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(iptem,npoint_files,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(is_netcdf_iptem,npoint_files,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(idxpnt_files,npoint_files*MXSPEC,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(irstc,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(irstf,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(irstp,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(ipr_unit,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(irr_unit,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(iavg,ngrid,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(idep,ngrid,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(ismin,ngrid,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(ismout,ngrid,MPI_INTEGER,itag,numprocs,iproc_id)
      max_emissfiles = MAXVAL( nemiss_files(1:ngrid) )
      call nodes_pass(iarem,ngrid*max_emissfiles,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(buffer_offset_iarem,ngrid*max_emissfiles,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(is_netcdf_iarem,ngrid*max_emissfiles,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(idxems_files,ngrid*max_emissfiles*MXSPEC,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(isurf,ngrid,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(is_netcdf_isurf,ngrid,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(i3dmet,ngrid,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(is_netcdf_i3dmet,ngrid,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(i2dmet,ngrid,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(is_netcdf_i2dmet,ngrid,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(ikv,ngrid,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(is_netcdf_ikv,ngrid,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(icld,ngrid,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(lcig,ngrid,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(n3dmet,ngrid,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(n2dmet,ngrid,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(nkvmet,ngrid,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(ncldmet,ngrid,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(nsrfvar,ngrid,MPI_INTEGER,itag,numprocs,iproc_id)
c
c  --- send the variables for the o3colmap include file ---
c
      call nodes_pass(lrdocn,ngrid,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(lrdlai,ngrid,MPI_LOGICAL,itag,numprocs,iproc_id)
c
c  --- send the variables in the flags include file ----
c
      call nodes_pass(iadvct,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(lrstrt,1,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(lflexi,1,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(idrydep,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(ldry,1,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(lwet,1,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(lutm,1,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(llatlon,1,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(lrpolar,1,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(lambrt,1,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(lpolar,1,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(lmerc,1,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(lstagw,MXGRID,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(larsrc,1,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(lptsrc,1,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(ltopcon,1,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(ltrace,1,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(lsa_ioric,1,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(num_ioric,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(ncls_ioric,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(is_netcdf_ioric,1,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(lsa_iorbc,1,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(num_iorbc,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(ncls_iorbc,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(is_netcdf_iorbc,1,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(lrd_sa_bc_hdr,1,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(lsa_iortc,1,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(ncls_iortc,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(is_netcdf_iortc,1,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(lrd_sa_tc_hdr,1,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(tectyp,  10,MPI_CHARACTER,itag,numprocs,iproc_id)
      call nodes_pass(lapca,1,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(lapcapt,1,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(lpartial,1,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(lptdepout,1,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(lproca,1,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(lipr,1,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(lirr,1,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(lcpacum,1,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(lddm,1,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(lhddm,1,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(lddmcalc,ngrid,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(lsrfmod,1,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(lsrfmodrt,1,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(lixemis,1,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(lixbypass,1,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(lbidinh3,1,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(lstrato3,1,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(lgas_out_ppm,1,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(lzppm,1,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(lparttn,1,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(le1day,1,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(lairqul,1,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(ldiag,1,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(lacm2,1,MPI_LOGICAL,itag,numprocs,iproc_id)
c
c  --- send the other variables in the chmstry include file ----
c
      call nodes_pass(idmech,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(idsolv,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(nicspc,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(nbcspc,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(ntcspc,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(aeropt,10,MPI_CHARACTER,itag,numprocs,iproc_id)
      call nodes_pass(leqsam,1,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(jno2rxn,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(jo3rxn,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(jhcho1rxn,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(jhcho2rxn,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(jch3chorxn,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(hthal,NHTHAL,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(brlprof,NHTHAL,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(brwprof,NHTHAL,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(brolprof,NHTHAL,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(browprof,NHTHAL,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(cl2day,NHTHAL,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(cl2nite,NHTHAL,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(kmap,NSPNAM,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(henso20,1,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(tfactso2,1,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(cwmin,1,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(tamin,1,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(nbin,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(ipigrxn,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(lsmap,NSPIGCHEM,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(lrmap,NRPIGCHEM,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(ihetrxn,NHETRXN,MPI_INTEGER,itag,numprocs,iproc_id)
c
c  --- send the variables in the soap include file ----
c
      call nodes_pass(kpolya,1,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(kpolyb,1,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(sfacJSOA,1,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(csatS,NSOAP,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(cstempS,NSOAP,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(deltahS,NSOAP,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(y_h,3*NPREC_S,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(y_l,3*NPREC_S,MPI_REAL,itag,numprocs,iproc_id)
c
c  --- send the variables in the vbs include file ----
c
      call nodes_pass(lvbs,1,MPI_LOGICAL,itag,numprocs,iproc_id)
      if ( lvbs .and. LVBSPREPROC ) then
        call nodes_pass(kpap_c,NVOLBIN+1,MPI_INTEGER,itag,numprocs,iproc_id)
        call nodes_pass(kpcp_c,NVOLBIN+1,MPI_INTEGER,itag,numprocs,iproc_id)
        call nodes_pass(kpfp_c,NVOLBIN+1,MPI_INTEGER,itag,numprocs,iproc_id)
      endif
c
c  --- send WTMOL for the PM chemistry ---
c
      call nodes_pass(wtmol,NWTMOL,MPI_REAL,itag,numprocs,iproc_id)
c
c  --- send the pointers for the PM chemistry ---
c
      call nodes_pass(kso2_c,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(kh2o2_c,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(kform_c,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(khono_c,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(ko3_c,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(koh_c,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(kho2_c,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(kno3_c,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(kno_c,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(kno2_c,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(kpan_c,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(kcg1_c,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(kcg2_c,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(kcg3_c,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(kcg4_c,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(khno3_c,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(knh3_c,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(kh2so4_c,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(khcl_c,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(ksoa1_c,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(ksoa2_c,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(ksoa3_c,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(ksoa4_c,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(ksopa_c,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(ksopb_c,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(kcrst_c,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(kpoa_c,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(kpec_c,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(kph2o_c,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(kpcl_c,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(kna_c,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(kpnh4_c,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(kpno3_c,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(kpso4_c,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(kn2o5_c,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(khpo_c,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(kfoa_c,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(kmhp_c,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(kpaa_c,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(kohp_c,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(kopa_c,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(kgly_c,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(kmgly_c,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(kglyd_c,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(ksoac_c,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(ngas_c,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(nspec_c,1,MPI_INTEGER,itag,numprocs,iproc_id)
c
c  --- send the other variables in the chmstry include file ----
c
      max_emissfiles = MAXVAL( nemiss_files(1:ngrid) )
      call nodes_pass(narspc,ngrid*max_emissfiles,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(lgas,nspec,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(emspcname,10*nspec,MPI_CHARACTER,itag,numprocs,iproc_id)
      call nodes_pass(ltdep,nreact,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(lpdep,nreact,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(bdnl,nspec+1,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(idphot1,nphot1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(idphot2,nphot2,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(idphot3,nphot2,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(phtscl,nphot2,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(spname,10*(nspec+1),MPI_CHARACTER,itag,
     &                                               numprocs,iproc_id)
      call nodes_pass(depsp,10*(4*nspec+2),MPI_CHARACTER,itag,
     &                                               numprocs,iproc_id)
      call nodes_pass(lbcmap,nspec,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(ltcmap,nspec,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(lavmap,nspec+nrad,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(licmap,nspec,MPI_INTEGER,itag,
     &                                              numprocs,iproc_id)
      call nodes_pass(ldepmap,nspec,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(lemmap,nspec,MPI_INTEGER,itag,
     &                                               numprocs,iproc_id)
      call nodes_pass(rktbl,nreact*NTEMPR*NPRESR,MPI_REAL,itag,
     &                                               numprocs,iproc_id)
      call nodes_pass(prkn,NZEN*nphot1*NHGHT*NTRN*NALB*NOZN,
     &                                MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(tempr,NTEMPR,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(presr,NPRESR,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(htint,NHGHT,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(zenint,NZEN,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(henry0,nspec,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(tfact,nspec,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(diffrat,nspec,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(f0,nspec,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(rscale,nspec,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(mole_weight,nspec,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(roprt,nspec,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(dcut,nspec*2,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(bext,nspec,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(ssa,nspec,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(rhadj,nspec,MPI_INTEGER,itag,numprocs,iproc_id)
c
      if( lsrfmod ) then
         if (iproc_id .GT. 0) call alloc_srfmod(nsmspc,nsmrxn)
         call nodes_pass(idsmsp,nsmspc,MPI_INTEGER,itag,numprocs,iproc_id)
         call nodes_pass(idsmpre,nsmrxn,MPI_INTEGER,itag,numprocs,iproc_id)
         call nodes_pass(idsmprd,nsmrxn,MPI_INTEGER,itag,numprocs,iproc_id)
         call nodes_pass(smskrat,nsmrxn,MPI_REAL,itag,numprocs,iproc_id)
         call nodes_pass(smsjrat,nsmrxn,MPI_REAL,itag,numprocs,iproc_id)
         call nodes_pass(smvkrat,nsmrxn,MPI_REAL,itag,numprocs,iproc_id)
         call nodes_pass(smvjrat,nsmrxn,MPI_REAL,itag,numprocs,iproc_id)
         call nodes_pass(smikrat,nsmrxn,MPI_REAL,itag,numprocs,iproc_id)
         call nodes_pass(smijrat,nsmrxn,MPI_REAL,itag,numprocs,iproc_id)
         call nodes_pass(smssrb,nsmspc,MPI_REAL,itag,numprocs,iproc_id)
         call nodes_pass(smvsrb,nsmspc,MPI_REAL,itag,numprocs,iproc_id)
         call nodes_pass(smisrb,nsmspc,MPI_REAL,itag,numprocs,iproc_id)
         call nodes_pass(smlch,nsmspc,MPI_REAL,itag,numprocs,iproc_id)
         call nodes_pass(smpen,nsmspc,MPI_REAL,itag,numprocs,iproc_id)
         call nodes_pass(smmlt,nsmspc,MPI_REAL,itag,numprocs,iproc_id)
      endif
c
c  --- send the other variables in the grid include file ----
c
      call nodes_pass(iuzon,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(itzon,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(nnest,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(xorg,1,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(yorg,1,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(delx,1,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(dely,1,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(polelon,1,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(polelat,1,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(tlat1,1,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(tlat2,1,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(deltay,ngrid,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(ntim,ngrid,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(ntimcrs,ngrid,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(i1,ngrid,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(j1,ngrid,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(i2,ngrid,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(j2,ngrid,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(nmesh,ngrid,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(nchdrn,ngrid,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(idchdrn,MXCHDRN*ngrid,MPI_INTEGER,itag,
     &                                               numprocs,iproc_id)
      call nodes_pass(nosrc,ngrid,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(inst1,ngrid,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(inst2,ngrid,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(jnst1,ngrid,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(jnst2,ngrid,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(meshold,ngrid,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(mapgrd,ngrid,MPI_INTEGER,itag,numprocs,iproc_id)
c
c  --- send the variables in the camx include file ----
c
      call nodes_pass(begdate,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(begtim,1,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(dtinp,1,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(dtems,1,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(dtout,1,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(dtmax,1,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(densfac,1,MPI_REAL,itag,numprocs,iproc_id)

c  --- send variables that depend on rows, cols or layers ---
c
      call nodes_pass(deltax,ngrid*nrowa,MPI_REAL,itag,
     &                                                numprocs,iproc_id)
      do igrd=1,ngrid
        itag = itag+1
        if (iproc_id .eq. 0) then
          CALL master_send_gridded_data(idfin(iptr2d(igrd)), igrd, 1, 1, itag)
        else
          CALL node_recv_gridded_data  (idfin(iptr2d(igrd)), igrd, 1, 1, itag)
        end if
        itag = itag+1
        if (iproc_id .eq. 0) then
          CALL master_send_gridded_data(ldark(iptr2d(igrd)), igrd, 1, 1, itag)
        else
          CALL node_recv_gridded_data(ldark(iptr2d(igrd)), igrd, 1, 1, itag)
        end if
        itag = itag+1
        if (iproc_id .eq. 0) then
          CALL master_send_gridded_data(cellon(iptr2d(igrd)), igrd, 1, 1, itag)
        else
          CALL node_recv_gridded_data(cellon(iptr2d(igrd)), igrd, 1, 1, itag)
        end if
        itag = itag+1
        if (iproc_id .eq. 0) then
          CALL master_send_gridded_data(cellat(iptr2d(igrd)), igrd, 1, 1, itag)
        else
          CALL node_recv_gridded_data(cellat(iptr2d(igrd)), igrd, 1, 1, itag)
        end if
        itag = itag+1
        if (iproc_id .eq. 0) then
          CALL master_send_gridded_data(lai(iptr2d(igrd)), igrd, 1, 1, itag)
        else
          CALL node_recv_gridded_data(lai(iptr2d(igrd)), igrd, 1, 1, itag)
        end if
        itag = itag+1
        if (iproc_id .eq. 0) then
          CALL master_send_gridded_data(mapscl(iptr2d(igrd)), igrd, 1, 1, itag)
        else
          CALL node_recv_gridded_data(mapscl(iptr2d(igrd)), igrd, 1, 1, itag)
        end if
        itag = itag+1
        if (iproc_id .eq. 0) then
          CALL master_send_gridded_data(fsurf(iptrlu(igrd)), igrd, nlu, 1, itag)
        else
          CALL node_recv_gridded_data(fsurf(iptrlu(igrd)), igrd, nlu, 1, itag)
        end if
        itag = itag+1
      enddo
c
c  --- send variables that depend on point sources ---
c
      call nodes_pass(idsrc,nptsrc*ngrid,MPI_INTEGER,itag,
     &                                               numprocs,iproc_id)
      call nodes_pass(isrc,nptsrc*ngrid,MPI_INTEGER,itag,
     &                                               numprocs,iproc_id)
      call nodes_pass(jsrc,nptsrc*ngrid,MPI_INTEGER,itag,
     &                                               numprocs,iproc_id)
c
c  --- pass the arrays for point source data ---
c
      call nodes_pass(hstk,nptsrc,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(dstk,nptsrc,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(tstk,nptsrc,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(lpiglet,nptsrc,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(xstk,nptsrc*ngrid,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(ystk,nptsrc*ngrid,MPI_REAL,itag,numprocs,iproc_id)
c
c  --- send the other variables in the o3colmap include file ----
c
      call nodes_pass(albcl,NALB,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(trncl,NTRN,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(ozcl,NOZN,MPI_REAL,itag,numprocs,iproc_id)
c
c  --- pass the variables in the deposit include file ---
c
      call nodes_pass(rgso,NLUW89*5,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(z0lu,NLUW89*5,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(rj,NLUW89*5,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(rlu,NLUW89*5,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(rac,NLUW89*5,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(rlcs,NLUW89*5,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(rlco,NLUW89*5,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(rgss,NLUW89*5,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(lai_ref,NLUZ03*15,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(iseason,5*12,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(dstress,6,MPI_REAL,itag,numprocs,iproc_id)
c
c  --- pass the variables in the camx_aero include file ---
c
      call nodes_pass(nacl,1,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(co2 ,1,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(foa,1,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(mhp,1,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(paa,1,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(caco3,1,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(mgco3,1,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(a3fe,1,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(b2mn,1,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(potcl,1,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(eps,1,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(e0,1,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(lv,1,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(rv,1,MPI_REAL,itag,numprocs,iproc_id)
c
c ----- Chemistry parameters ---
c
      call nodes_pass(lmap,NSPNAM,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(nirrrxn,1,MPI_INTEGER,itag,numprocs,iproc_id)
      if( lproca ) then
         mvec2d = 0
         mvec3d = 0
         do i=1,ngrid
            mvec2d = mvec2d + nrow(i) * ncol(i)
            mvec3d = mvec3d + nrow(i) * ncol(i) * nlay(i)
         enddo
         call nodes_pass(ptname,10*MXCPA,MPI_CHARACTER,itag,numprocs,iproc_id)
         call nodes_pass(cpadesc,60*MXCPA,MPI_CHARACTER,itag,numprocs,iproc_id)
         call nodes_pass(cpaunit,20*MXCPA,MPI_CHARACTER,itag,numprocs,iproc_id)
         call nodes_pass(ptop_fac,MXCPA,MPI_REAL,itag,numprocs,iproc_id)
         call nodes_pass(loutsa,MXCPA,MPI_LOGICAL,itag,numprocs,iproc_id)
         call nodes_pass(ipacl_2d,mvec2d,MPI_INTEGER,itag,numprocs,iproc_id)
         call nodes_pass(ipacl_3d,mvec3d,MPI_INTEGER,itag,numprocs,iproc_id)
         call nodes_pass(ipadom,npa_cels,MPI_INTEGER,itag,numprocs,iproc_id)
         call nodes_pass(ipax,npa_cels,MPI_INTEGER,itag,numprocs,iproc_id)
         call nodes_pass(ipay,npa_cels,MPI_INTEGER,itag,numprocs,iproc_id)
         call nodes_pass(ipaz,npa_cels,MPI_INTEGER,itag,numprocs,iproc_id)
         call nodes_pass(ipanst,npa_cels,MPI_INTEGER,itag,numprocs,iproc_id)
      endif
c
c ----- pass variables in the tracer include file ----
c       if not doing probing tools, just skip this ---
c
      if( .NOT. ltrace .and. .NOT. lddm .AND. .NOT. lhddm ) goto 111
      call nodes_pass(num_iortem,ngrid*(ngroup+2),MPI_INTEGER,
     &                                      itag,numprocs,iproc_id)
      call nodes_pass(num_iortpt,ngroup+2,MPI_INTEGER,
     &                                      itag,numprocs,iproc_id)
      if( iproc_id .GT. 0) call alloc_tracer_emfiles(ngroup,ngrid,
     &           nemiss_files,npoint_files,num_iortem,num_iortpt,nspec)
c
c========================= Probing Tools Begin =======================
c
c ----- call routine to allocate the RTRAC data structures ---
c
      if( iproc_id .GT. 0) then
         if( tectyp .EQ. RTRAC ) then
            call alloc_rtracchm(ngrid,mmxp,mmyp,nspec,ntotsp)
         elseif ( tectyp .EQ. RTCMC ) then
            call alloc_rtracchm(ngrid,mmxp,mmyp,nspec,ntotsp)
            call alloc_rtcmc(MXTRSP,MXSPEC,MXRX,
     &                       MXKPRM,MXPHOT,MXZEN,MXRCT,MXPRD,
     &                       MXJACTRM,MXEQM,MXSLO)
         endif
      endif
c
c ----- Variables for gridded tracer emissions data
c
      call nodes_pass(ipigmap,nptsrc,MPI_INTEGER,itag,
     &                                          numprocs,iproc_id)
      call nodes_pass(ipiggrp,nptsrc,MPI_INTEGER,itag,
     &                                          numprocs,iproc_id)
      do igrd=1,ngrid
        do ilay=1,nlayers_ems
          do ispc=1,ntotsp
            if (iproc_id .eq. 0) then
              CALL master_send_1species_data(saemis(ipsa3d_ems(igrd)), igrd, 
     &                               nlayers_ems, ilay, ntotsp, ispc, itag)
            else
              CALL node_recv_1species_data(saemis(ipsa3d_ems(igrd)), igrd, 
     &                               nlayers_ems, ilay, ntotsp, ispc, itag)
            endif
            itag = itag+1
          enddo
        enddo
      enddo
      if( lptsrc ) then
         call nodes_pass_sapnts(numprocs,iproc_id)
         call nodes_pass(xlocpt,nptsrc,MPI_REAL,itag,numprocs,iproc_id)
         call nodes_pass(ylocpt,nptsrc,MPI_REAL,itag,numprocs,iproc_id)
         call nodes_pass(lpigsa,nptsrc,MPI_LOGICAL,itag,numprocs,iproc_id)
      endif
      call nodes_pass(nreles,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(ntrtim,1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(lreles,1,MPI_LOGICAL,itag,numprocs,iproc_id)
c
c  ---- Variables for tracer names ----
c
      call nodes_pass(ptname,10*ntotsp,MPI_CHARACTER,itag,
     &                                               numprocs,iproc_id)
      call nodes_pass(ptop_fac,ntotsp,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(lsamap,ntotsp,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(lsagas,ntotsp,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(sa_mole_weight,ntotsp,MPI_REAL,itag,numprocs,iproc_id)
      if( tectyp .NE. RTRAC .AND. tectyp .NE. RTCMC ) then
          call nodes_pass(trspmap,nspec*ntrcls,MPI_REAL,itag,
     &                                              numprocs,iproc_id)
          call nodes_pass(fluxmap,nspec*ntrcls,MPI_REAL,itag,
     &                                               numprocs,iproc_id)
          call nodes_pass(yhratmap,nspec*ntrcls,MPI_REAL,itag,
     &                                               numprocs,iproc_id)
          call nodes_pass(ylratmap,nspec*ntrcls,MPI_REAL,itag,
     &                                               numprocs,iproc_id)
          call nodes_pass(lusespc,nspec,MPI_LOGICAL,itag,numprocs,iproc_id)
c
          call nodes_pass(clsnam,3*MXALCLS,MPI_CHARACTER,itag,
     &                                               numprocs,iproc_id)
          call nodes_pass(iptcls,MXALCLS,MPI_INTEGER,itag,
     &                                               numprocs,iproc_id)
          call nodes_pass(nptcls,MXALCLS,MPI_INTEGER,itag,
     &                                               numprocs,iproc_id)
          call nodes_pass(ipttrc,MXALCLS,MPI_INTEGER,itag,
     &                                               numprocs,iproc_id)
          call nodes_pass(npttrc,MXALCLS,MPI_INTEGER,itag,
     &                                               numprocs,iproc_id)
          call nodes_pass(idxcls,MXALCLS,MPI_INTEGER,itag,
     &                                               numprocs,iproc_id)
          call nodes_pass(idxipt,MXALCLS,MPI_INTEGER,itag,
     &                                               numprocs,iproc_id)
          call nodes_pass(iemcls,MXALCLS,MPI_INTEGER,itag,
     &                                               numprocs,iproc_id)
          call nodes_pass(nemcls,MXALCLS,MPI_INTEGER,itag,
     &                                               numprocs,iproc_id)
          call nodes_pass(nsaspc,MXALCLS,MPI_INTEGER,itag,
     &                                               numprocs,iproc_id)
          call nodes_pass(npttim,1,MPI_INTEGER,itag,numprocs,iproc_id)
          call nodes_pass(ipttim,1,MPI_INTEGER,itag,numprocs,iproc_id)
          call nodes_pass(iemtim,1,MPI_INTEGER,itag,numprocs,iproc_id)
          call nodes_pass(nbdic, 1,MPI_INTEGER,itag,numprocs,iproc_id)
      endif
c
c  ---- Variables for turning parts of SA ---
c
      call nodes_pass(lsulfate,1,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(lnitrate,1,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(lsoa,    1,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(lprimary,1,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(lmercury,1,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(lozone,  1,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(lmineral,1,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(ldmschm ,1,MPI_LOGICAL,itag,numprocs,iproc_id)
c
c   ---- Variables for user options and flags ---
c
      call nodes_pass(loutsa,ntotsp,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass (ngroup,  1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass (nchar,   1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(lrestrt, 1,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(leftovr_area, 1,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(lbndry,  1,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(lallout, 1,MPI_LOGICAL,itag,numprocs,iproc_id)
c
c  ---- Variables for region mapping ---
c
      call nodes_pass(nxcell,ngrid,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(nycell,ngrid,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(npartial,(ngroup+1)*ngrid,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(igrmap,(ngroup+1)*MXPARTIAL*ngrid*ncola*nrowa,
     &                             MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(frcmap,(ngroup+1)*MXPARTIAL*ngrid*ncola*nrowa,
     &                                MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(nregin,    1,MPI_INTEGER,itag,numprocs,iproc_id)
c
c  ---- Variables for receptor data ---
c
      call nodes_pass(conrcp,ntotsp*MXRECP,MPI_REAL,itag,
     &                                                numprocs,iproc_id)
      call nodes_pass(rcpnam,10*MXRECP,MPI_CHARACTER,itag,
     &                                                numprocs,iproc_id) 
      call nodes_pass(idrcp,  MXRECP,MPI_INTEGER,itag,
     &                                                numprocs,iproc_id)
      call nodes_pass(igrdrcp,MXRECP,MPI_INTEGER,itag,
     &                                                numprocs,iproc_id)
      call nodes_pass(irecep,MXRECP*MXCELR,MPI_INTEGER,itag,
     &                                                numprocs,iproc_id)
      call nodes_pass(jrecep,MXRECP*MXCELR,MPI_INTEGER,itag,
     &                                                numprocs,iproc_id)
      call nodes_pass(nclrcp,MXRECP,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(iwalbg,MXRECP,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(iwalnd,MXRECP,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(jwalbg,MXRECP,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(jwalnd,MXRECP,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(kwalbg,MXRECP,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(kwalnd,MXRECP,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(recepx,MXRECP,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(recepy,MXRECP,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(nrecep,     1,MPI_INTEGER,itag,numprocs,iproc_id)
      call nodes_pass(lwalls,     1,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(lrcpfil,    1,MPI_LOGICAL,itag,numprocs,iproc_id)
c
c  ---- Variables for species order and species flags ----
c
      max_ptfiles = MAXVAL(num_iortpt)
      call nodes_pass(idxpts,(ngroup+2)*max_ptfiles*nspec*100,MPI_INTEGER,itag,
     &                                               numprocs,iproc_id)
      call nodes_pass(idx_point_in_list,(ngroup+2)*max_ptfiles,MPI_INTEGER,itag,
     &                                               numprocs,iproc_id)
      call nodes_pass(nspcpt,(ngroup+2)*max_ptfiles,MPI_INTEGER,itag,
     &                                               numprocs,iproc_id)
      call nodes_pass(lvocsp,nspec,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(lvocsoa,nspec,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(lhrvoc,nspec,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(lnoxsp,nspec,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(lo3sp, nspec,MPI_LOGICAL,itag,numprocs,iproc_id)
      call nodes_pass(crbnum,nspec,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(mwspec,nspec,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(rkohrt,nspec,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(rmirrt,nspec,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(wtkoh,ntotsp,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(wtmir,ntotsp,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(yhrates,ntotsp,MPI_REAL,itag,numprocs,iproc_id)
      call nodes_pass(ylrates,ntotsp,MPI_REAL,itag,numprocs,iproc_id)
c
      if( tectyp .EQ. RTRAC ) then
          call nodes_pass(lsecnd,ntotsp,MPI_LOGICAL,itag,numprocs,iproc_id)
          call nodes_pass(lreg,ntotsp,MPI_LOGICAL,itag,numprocs,iproc_id)
          call nodes_pass(ksec,ntotsp,MPI_INTEGER,itag,numprocs,iproc_id)
          call nodes_pass(rtlbnd,ntotsp,MPI_REAL,itag,numprocs,iproc_id)
          call nodes_pass(rthlaw,ntotsp,MPI_REAL,itag,numprocs,iproc_id)
          call nodes_pass(rttfact,ntotsp,MPI_REAL,itag,numprocs,iproc_id)
          call nodes_pass(rtmolwt,ntotsp,MPI_REAL,itag,numprocs,iproc_id)
          call nodes_pass(rtdrate,ntotsp,MPI_REAL,itag,numprocs,iproc_id)
          call nodes_pass(rtreact,ntotsp,MPI_REAL,itag,numprocs,iproc_id)
          call nodes_pass(rtscale,ntotsp,MPI_REAL,itag,numprocs,iproc_id)
          call nodes_pass(rtdens,ntotsp,MPI_REAL,itag,numprocs,iproc_id)
          call nodes_pass(rtlcut,ntotsp,MPI_REAL,itag,numprocs,iproc_id)
          call nodes_pass(rtucut,ntotsp,MPI_REAL,itag,numprocs,iproc_id)
          call nodes_pass(jnum,ntotsp,MPI_INTEGER,itag,numprocs,iproc_id)
          call nodes_pass(rtjfact,ntotsp,MPI_REAL,itag,numprocs,iproc_id)
          call nodes_pass(aoh,ntotsp,MPI_REAL,itag,numprocs,iproc_id)
          call nodes_pass(eaoh,ntotsp,MPI_REAL,itag,numprocs,iproc_id)
          call nodes_pass(boh,ntotsp,MPI_REAL,itag,numprocs,iproc_id)
          call nodes_pass(troh,ntotsp,MPI_REAL,itag,numprocs,iproc_id)
          call nodes_pass(ano3,ntotsp,MPI_REAL,itag,numprocs,iproc_id)
          call nodes_pass(eano3,ntotsp,MPI_REAL,itag,numprocs,iproc_id)
          call nodes_pass(bno3,ntotsp,MPI_REAL,itag,numprocs,iproc_id)
          call nodes_pass(trno3,ntotsp,MPI_REAL,itag,numprocs,iproc_id)
          call nodes_pass(ao3,ntotsp,MPI_REAL,itag,numprocs,iproc_id)
          call nodes_pass(eao3,ntotsp,MPI_REAL,itag,numprocs,iproc_id)
          call nodes_pass(bo3,ntotsp,MPI_REAL,itag,numprocs,iproc_id)
          call nodes_pass(tro3,ntotsp,MPI_REAL,itag,numprocs,iproc_id)
          call nodes_pass(fsoil,NLUW89,MPI_REAL,itag,numprocs,iproc_id)
          call nodes_pass(fsoilz,NLUZ03,MPI_REAL,itag,numprocs,iproc_id)
          call nodes_pass(fsoiloc,NLUW89,MPI_REAL,itag,numprocs,iproc_id)
          call nodes_pass(eqkoa,ntotsp,MPI_REAL,itag,numprocs,iproc_id)
          call nodes_pass(khydro,ntotsp,MPI_REAL,itag,numprocs,iproc_id)
          call nodes_pass(kleach,ntotsp,MPI_REAL,itag,numprocs,iproc_id)
          call nodes_pass(kpen,ntotsp,MPI_REAL,itag,numprocs,iproc_id)
          call nodes_pass(kphot,ntotsp,MPI_REAL,itag,numprocs,iproc_id) 
          if( lrcpfil ) then
              mvec3d = 0
              do i=1,ngrid
                 mvec3d = mvec3d + nrow(i) * ncol(i) * nlay(i)
              enddo
              call nodes_pass(ipacl_3d,mvec3d,MPI_INTEGER,
     &                                         itag,numprocs,iproc_id)
          endif
      elseif( tectyp .EQ. RTCMC ) then
          call nodes_pass(ityprtc,MXRX,MPI_INTEGER,itag,numprocs,iproc_id)
          call nodes_pass(nrkprm,MXRX,MPI_INTEGER,itag,numprocs,iproc_id)
          call nodes_pass(rkprmrtc,MXRX*MXKPRM,MPI_REAL,itag,numprocs,iproc_id)
          call nodes_pass(srkrtc,MXRX,MPI_REAL,itag,numprocs,iproc_id)
          call nodes_pass(rkrtc,MXRX,MPI_DOUBLE_PRECISION,itag,numprocs,iproc_id)
          call nodes_pass(rthlaw,ntotsp,MPI_REAL,itag,numprocs,iproc_id)
          call nodes_pass(rttfact,ntotsp,MPI_REAL,itag,numprocs,iproc_id)
          call nodes_pass(rtdrate,ntotsp,MPI_REAL,itag,numprocs,iproc_id)
          call nodes_pass(rtreact,ntotsp,MPI_REAL,itag,numprocs,iproc_id)
          call nodes_pass(rtscale,ntotsp,MPI_REAL,itag,numprocs,iproc_id)
          call nodes_pass(rtdens,ntotsp,MPI_REAL,itag,numprocs,iproc_id)
          call nodes_pass(rtlcut,ntotsp,MPI_REAL,itag,numprocs,iproc_id)
          call nodes_pass(rtucut,ntotsp,MPI_REAL,itag,numprocs,iproc_id)
          call nodes_pass(ijschm,MXPHOT,MPI_INTEGER,itag,numprocs,iproc_id)
          call nodes_pass(ijpntr,MXPHOT,MPI_INTEGER,itag,numprocs,iproc_id)
          call nodes_pass(jscale,MXPHOT,MPI_REAL,itag,numprocs,iproc_id)
          call nodes_pass(zenschm,MXZEN,MPI_REAL,itag,numprocs,iproc_id)
          call nodes_pass(rjschm,MXZEN*MXPHOT,MPI_REAL,itag,numprocs,iproc_id)
          call nodes_pass(itypsp,MXTRSP+MXSPEC,MPI_INTEGER,itag,numprocs,iproc_id)
          call nodes_pass(itypschm,MXTRSP+MXSPEC,MPI_INTEGER,itag,numprocs,iproc_id)
          call nodes_pass(lblrxn,MXRX,MPI_INTEGER,itag,numprocs,iproc_id)
          call nodes_pass(nrct,MXRX,MPI_INTEGER,itag,numprocs,iproc_id)
          call nodes_pass(nprd,MXRX,MPI_INTEGER,itag,numprocs,iproc_id)
          call nodes_pass(idxrct,MXRX*MXRCT,MPI_INTEGER,itag,numprocs,iproc_id)
          call nodes_pass(idxprd,MXRX*MXPRD,MPI_INTEGER,itag,numprocs,iproc_id)
          call nodes_pass(conschm,MXTRSP+MXSPEC,MPI_REAL,itag,numprocs,iproc_id)
          call nodes_pass(spdcoef,MXRX*MXPRD,MPI_REAL,itag,numprocs,iproc_id)
          call nodes_pass(prdcoef,MXRX*MXPRD,MPI_DOUBLE_PRECISION,itag,numprocs,iproc_id)
          call nodes_pass(spnmrt,MXTRSP+MXSPEC*10,MPI_CHARACTER,itag,numprocs,iproc_id)
          call nodes_pass(spnmschm,MXTRSP+MXSPEC,MPI_CHARACTER,itag,numprocs,iproc_id)
          call nodes_pass(hlawschm,MXTRSP+MXSPEC,MPI_REAL,itag,numprocs,iproc_id)
          call nodes_pass(tfacschm,MXTRSP+MXSPEC,MPI_REAL,itag,numprocs,iproc_id)
          call nodes_pass(molwschm,MXTRSP+MXSPEC,MPI_REAL,itag,numprocs,iproc_id)
          call nodes_pass(reacschm,MXTRSP+MXSPEC,MPI_REAL,itag,numprocs,iproc_id)
          call nodes_pass(rsclschm,MXTRSP+MXSPEC,MPI_REAL,itag,numprocs,iproc_id)
          call nodes_pass(namrct,MXRX*MXRCT,MPI_CHARACTER,itag,numprocs,iproc_id)
          call nodes_pass(namprd,MXRX*MXPRD,MPI_CHARACTER,itag,numprocs,iproc_id)
          call nodes_pass(nrctfst,MXRX,MPI_INTEGER,itag,numprocs,iproc_id)
          call nodes_pass(nprdfst,MXRX,MPI_INTEGER,itag,numprocs,iproc_id)
          call nodes_pass(idxrctfst,MXRX*MXRCT,MPI_INTEGER,itag,numprocs,iproc_id)
          call nodes_pass(idxprdfst,MXRX*MXPRD,MPI_INTEGER,itag,numprocs,iproc_id)
          call nodes_pass(prdcofst,MXRX*MXPRD,MPI_DOUBLE_PRECISION,itag,numprocs,iproc_id)
          call nodes_pass(ipd,MXJACTRM,MPI_INTEGER,itag,numprocs,iproc_id)
          call nodes_pass(jpd,MXJACTRM,MPI_INTEGER,itag,numprocs,iproc_id)
          call nodes_pass(idrxjac,MXJACTRM,MPI_INTEGER,itag,numprocs,iproc_id)
          call nodes_pass(nspjac,MXJACTRM,MPI_INTEGER,itag,numprocs,iproc_id)
          call nodes_pass(idspjac,MXJACTRM*(MXRCT-1),MPI_INTEGER,itag,numprocs,iproc_id)
          call nodes_pass(scoefjac,MXJACTRM,MPI_REAL,itag,numprocs,iproc_id)
          call nodes_pass(coefjac,MXJACTRM,MPI_DOUBLE_PRECISION,itag,numprocs,iproc_id)
          call nodes_pass(idslo,MXSLO,MPI_INTEGER,itag,numprocs,iproc_id)
          call nodes_pass(nslgain,MXSLO,MPI_INTEGER,itag,numprocs,iproc_id)
          call nodes_pass(nslloss,MXSLO,MPI_INTEGER,itag,numprocs,iproc_id)
          call nodes_pass(islgain,MXSLO*MXRX,MPI_INTEGER,itag,numprocs,iproc_id)
          call nodes_pass(islloss,MXSLO*MXRX,MPI_INTEGER,itag,numprocs,iproc_id)
          call nodes_pass(spdcoslo,MXSLO*MXRX,MPI_DOUBLE_PRECISION,itag,numprocs,iproc_id)
          call nodes_pass(prdcoslo,MXSLO*MXRX,MPI_DOUBLE_PRECISION,itag,numprocs,iproc_id)
          call nodes_pass(ideqm,MXEQM,MPI_INTEGER,itag,numprocs,iproc_id)
          call nodes_pass(nrxgain,MXEQM,MPI_INTEGER,itag,numprocs,iproc_id)
          call nodes_pass(nrxloss,MXEQM,MPI_INTEGER,itag,numprocs,iproc_id)
          call nodes_pass(irxgain,MXEQM*MXRX,MPI_INTEGER,itag,numprocs,iproc_id)
          call nodes_pass(irxloss,MXEQM*MXRX,MPI_INTEGER,itag,numprocs,iproc_id)
          call nodes_pass(irxupdt,MXRX,MPI_INTEGER,itag,numprocs,iproc_id)
          call nodes_pass(spdcoeqm,MXEQM*MXRX,MPI_REAL,itag,numprocs,iproc_id)
          call nodes_pass(prdcoeqm,MXEQM*MXRX,MPI_DOUBLE_PRECISION,itag,numprocs,iproc_id)
          call nodes_pass(idxfix,MXTRSP+MXSPEC,MPI_INTEGER,itag,numprocs,iproc_id)
          if( lrcpfil ) then
              mvec3d = 0 
              do i=1,ngrid
                 mvec3d = mvec3d + nrow(i) * ncol(i) * nlay(i)
              enddo
              call nodes_pass(ipacl_3d,mvec3d,MPI_INTEGER,
     &                                         itag,numprocs,iproc_id)
          endif
      endif
c
c========================= Probing Tools End =======================
c
c======================== DDM Begin ================================
c
c  --- Tracer species map for DDM ----
c
      if( lddm .OR. lhddm ) then
c
         call nodes_pass(icddmsp,10*nspec,MPI_CHARACTER,itag,
     &                                                numprocs,iproc_id) 
c
         call nodes_pass(bcddmsp,10*nspec,MPI_CHARACTER,itag,
     &                                                numprocs,iproc_id) 
c
         call nodes_pass(emddmsp,10*nspec,MPI_CHARACTER,itag,
     &                                                numprocs,iproc_id) 
c
         call nodes_pass(nicddm,  1,MPI_INTEGER,itag,numprocs,iproc_id)
c
         call nodes_pass(nbcddm,  1,MPI_INTEGER,itag,numprocs,iproc_id)
c
         call nodes_pass(nemddm,  1,MPI_INTEGER,itag,numprocs,iproc_id)
c
c  ----- Species list and pointers into arrays for DDM species ---
c
         call nodes_pass(ptlong,14*ntotsp,MPI_CHARACTER,itag,
     &                                                numprocs,iproc_id)
         call nodes_pass(iptddm,  nspec,MPI_INTEGER,itag,
     &                                                numprocs,iproc_id)
         call nodes_pass(rateddm, 10*nrateddm, MPI_CHARACTER,
     &                                          itag,numprocs,iproc_id)
         call nodes_pass(termddm, 10*ntermddm, MPI_CHARACTER,
     &                                          itag,numprocs,iproc_id)
         call nodes_pass(hddmsp, 10*nhddm,MPI_CHARACTER,
     &                                          itag,numprocs,iproc_id)
         call nodes_pass(iprate,(nreact+1)*nrateddm, MPI_INTEGER,
     &                                          itag,numprocs,iproc_id)
         call nodes_pass(ipterm,(ngas+1)*nreact*ntermddm, MPI_INTEGER,
     &                                          itag,numprocs,iproc_id)
         call nodes_pass(wfterm,    ngas*nreact*ntermddm, MPI_REAL,
     &                                          itag,numprocs,iproc_id)
         call nodes_pass(iphddm,2*nhddm,MPI_INTEGER,
     &                                          itag,numprocs,iproc_id)
       endif
c
c======================== DDM End ==================================
c
c  --- return to the calling routine ---
c
  111 continue
      end
