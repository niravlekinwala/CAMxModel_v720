      subroutine pasetup
      use filunit
      use grid
      use chmstry
      use procan
      use tracer
c
c----CAMx v7.20 220430
c
c     This routine initializes some of the data strucutures for the
c     Process Analysis algorithm.  The "species" names array is filled 
c     and the arrays are checked to make sure they have been allocated 
c     properly.  The pointers into the gridded arrays are initialized.
c
c      Copyright 1996 - 2022
c     Ramboll
c
c     Modifications:
c        10/13/17 -bkoo-     Added aerosol pH to CPA variable list
c
c     Input arguments:
c        none
c
c     Subroutines Called:
c        CPAMECH
c
c     Called by:
c        STARTUP
c
      implicit none
      include 'camx.prm'
c
c-----Argument declarations
c
c-----Local variables
c
      real padum(MXCPA)
      real dtfact
      integer i,  npa_init
c
      real rdum(MXREACT)
      real rkdum(MXREACT)
c
c-----Entry point
c
c  --- set flag to determine whether CPA variables are cumulative ---
c
      lcpacum = .false.
c            
c  --- initialize the number of reactions ---
c
      nirrrxn =  nreact
c
c  --- set the number of species and initialize species names ---
c
c   --- dummy array set to zero for first call ---
c
      do i=1,nreact
       rdum(i) = 0.0
       rkdum(i) = 0.0
      enddo
      dtfact = 0.
c
c   --- dummy call sets number of CPA outputs and defines 
c       variable output names ---
c
      if(idmech .EQ. 1) then
         call cpamech1(rdum,rkdum,dtfact,nreact,padum,MXCPA,
     &                                           npa_init,.false.)
      elseif(idmech .EQ. 3) then
         call cpamech3(rdum,rkdum,dtfact,nreact,padum,MXCPA,
     &                                           npa_init,.false.)
      elseif(idmech .EQ. 4) then
         call cpamech4(rdum,rkdum,dtfact,nreact,padum,MXCPA,
     &                                           npa_init,.false.)
      elseif(idmech .EQ. 5) then
         call cpamech5(rdum,rkdum,dtfact,nreact,padum,MXCPA,
     &                                           npa_init,.false.)
      elseif(idmech .EQ. 6) then
         call cpamech6(rdum,rkdum,dtfact,nreact,padum,MXCPA,
     &                                           npa_init,.false.)
      elseif(idmech .EQ. 7) then
         call cpamech7(rdum,rkdum,dtfact,nreact,padum,MXCPA,
     &                                           npa_init,.false.)
      else
          goto 7001
      endif
c
c  --- add the cloud adjustment to photolysis rates ---
c
      npa_init = npa_init + 1
      ptname(npa_init)  = 'J_CLDADJ'
      cpadesc(npa_init) = 'Cloud/haze adjustment factor'
      cpaunit(npa_init) = 'hr^-1'
c
c  --- add the aerosol pH ---
c
      npa_init = npa_init + 1
      ptname(npa_init)  = 'AER_PH'
      cpadesc(npa_init) = 'aerosol pH'
      cpaunit(npa_init) = 'dimensionless'
c     
c  --- set names of radical concentrations to be saved by CPA 
c      don't save if CPA is set to accumulate values ---
c       
      if(.NOT.lcpacum)
     &     call cparad(rdum, nrad, padum, MXCPA, npa_init, 0.0)
c
      if( npa_init .GT. MXCPA ) then
         write(iout,'(//,a)') 'ERROR in PASETUP:'
         write(iout,*) 'Number of outputs requested exceeds limit.'
         write(iout,*)
     &     'Increase the parameter MXCPA in include file procan.inc'
         write(iout,'(1X,A,I5)')
     &         'You need room for at least this many species: ',npa_init
         call camxerr()
      endif
c
      ntotsp = npa_init
      nsaspc = npa_init
c
c  --- define the pointers for each grid ---
c
      ipsa3d(1) = 1
      do i=2,ngrid
         ipsa3d(i) = ipsa3d(i-1) + 
     &       DBLE(ncol(i-1))*DBLE(nrow(i-1))*DBLE(nlay(i-1))*DBLE(ntotsp)
      enddo
c
c   --- set the output flags for all species ----
c
      do i=1,MXCPA
        loutsa(i) = .TRUE.
      enddo
      lfirst = .TRUE.
c
c   --- report the CPA configuration to the *.diag file ---
c
      write(idiag,'(//,A)') ' Chemical Process Analysis (CPA) is active'
      write(idiag,'(A,I5)') ' Number of CPA parameters = ', npa_init
      write(idiag,'(A)')    ' The parameter names are:'
      write(idiag,'(10X,A)') (ptname(i), i=1,npa_init)
      if(lcpacum) write(idiag,'(A)')
     &  ' The flag to accumulate parameters is set true'
      write(idiag,'(//)')
c
c   --- done ---
c
      goto 9999
c
c-----------------------------------------------------------------------
c   Error messages:
c-----------------------------------------------------------------------
c
 7001 continue
      write(iout,'(//,a)') 'ERROR in PASETUP:'
      write(iout,'(2A,I2)') 'Chemical Process Analysis (CPA) does ',
     &                       'not support chemical mechanism: ',idmech
      write(iout,*) 'Please choose another mechanism and try again.'
      call camxerr()
      goto 9999
c
c-----------------------------------------------------------------------
c   Return point:
c-----------------------------------------------------------------------
c
 9999 continue
c
      return
      end
