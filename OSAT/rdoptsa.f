c**** RDOPTSA
c
      subroutine rdoptsa()
      use filunit
      use grid
      use procan
      use tracer
      use chmstry
      implicit none
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c     This routine loads all of the user options and flags for the
c     source apportionment algorithm.  
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
c     05/16/96   --gwilson--    Original development
c     02/10/97   --cemery--     Added read of fileroot for SA output files
c     04/14/97   --gwilson--    Changed the way number of source groups
c                               is specified
c     04/28/97   --gwilson--    Added flag for OPPAT
c     05/22/97   --gwilson--    Added flag for APCA
c     01/10/02   --cemery --    Eliminated the read for fine grid
c                               flag if no nests 
c     03/21/03   --gwilson--    Removed the OSAT technology type OPPAT
c     10/06/04   --cemery --    Restructured for namelist input
c     10/07/12   --gwilson--    Disabled timing tracers
c     09/11/15   --bkoo--       Revised for SA v3
c     03/01/16   --gwilson--    Added partial source area map
c     11/09/16   --cemery--     Added Baker APCA point source override option
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'namelist.inc'
      include 'flags.inc'
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      integer nemiss, ncount
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c
c   --- if option is turned off, just return ---
c
      if( .NOT. ltrace ) goto 9999
c
c   --- flag for stratifying the boundary by edge ---
c
      lbndry = SA_Stratify_Boundary
c
c   --- flag for writing 3D average ---
c
      lsa_3davrg = SA_3D_Average
c
c   --- number of source regions ---
c
      nregin = SA_Number_of_Source_Regions
c
c   --- number of source emissions groupings ---
c
      nemiss = SA_Number_of_Source_Groups
c
c   --- flag for determining if the leftover group should be used ---
c
      if( Use_Leftover_Group ) then
         write(iout,'(//,a)') 'ERROR in RDOPTSA:'
         write(iout,'(/,1X,A)') 'The Use_Leftover_Group flag is obsolete.'
         write(iout,'(1X,2A)') 'This version allows the leftover group for ',
     &                         'gridded sources only.'
         write(iout,'(1X,2A)') 'You must include all of the point sources ',
     &                         'in emissions groups.'
         write(iout,'(1X,2A,/,1X,A)') 'To turn on the leftover group for ',
     &                               'gridded sources only set the',
     &                               'flag Use_Gridded_Leftover_Group to true.'
         call camxerr()
      endif
      leftovr_area = Use_Gridded_Leftover_Group
c
c   --- flag for determining what type of output to produce ----
c
      lallout = (.NOT. SA_Summary_Output )
c
c   --- set the number of source groups from the emissions groups ---
c
      if( nemiss .EQ. 1 .AND. leftovr_area ) then
         write(iout,'(//,a)') 'ERROR in RDOPTSA:'
         write(iout,'(/,1X,2A)') 'Cannot have leftover group ',
     &                           'with only one source group.'
         write(iout,'(1X,2A)') 'Set number of groups to 2 or turn ',
     &                         'off leftover group.'
         call camxerr()
      endif
      if( leftovr_area ) then
          ngroup = nemiss - 1
      else
          ngroup = nemiss 
      endif
c
c   --- check for array overflow ---
c
      if( ngroup .GT. MXTEMF-1 ) then
         write(iout,'(//,a)') 'ERROR in RDOPTSA:'
         write(iout,'(/,1X,A,I4,A)')'Number of source groupings ',
     &                ngroup,' exceeds maximum.  Increase MXTEMF.'
         call camxerr()
      endif
      if( ngroup .LT. 0 ) then
         write(iout,'(//,a)') 'ERROR in RDOPTSA:'
         write(iout,'(1X,A,I4,A)')'Number of emissions groups ',nemiss,
     &                            ' is invalid.'
         call camxerr()
      endif
c
c   --- number of timing releases per day ----
c
      ntrtim = Number_of_Timing_Releases
c
c   --- disable timing tracers ---
c
      if( ntrtim .NE. 0 ) goto 7002
c
c   --- set flag for partial source are map ---
c
      if( SA_Use_Partial_SourceMap ) lpartial = .TRUE.
c
c   --- if doing SA, get the classes that should be treat ---
c	
      if( tectyp .EQ. SA ) then
         ncount = 0
         lsulfate = .FALSE.
         lnitrate = .FALSE.
         lsoa     = .FALSE.
         lprimary = .FALSE.
         lmercury = .FALSE.
         lozone   = .FALSE.
         lapca    = .FALSE.
         lapcapt  = .FALSE.
         if( SA_Treat_SULFATE_Class ) then
            lsulfate = .TRUE.
            ncount = ncount + 1
         endif
         if( SA_Treat_NITRATE_Class ) then
            lnitrate = .TRUE.
            if( SA_Use_APCA ) lapca = .TRUE.
            if( lbidinh3 .AND. .NOT. lapca ) goto 7004
            ncount = ncount + 1
         endif
         if( SA_Treat_SOA_Class ) then
            lsoa = .TRUE.
            ncount = ncount + 1
         endif
         if( SA_Treat_PRIMARY_Class ) then
            lprimary = .TRUE.
            ncount = ncount + 1
         endif
         if( SA_Treat_MERCURY_Class ) then
            lmercury = .TRUE.
            ncount = ncount + 1
         endif
         if( SA_Treat_OZONE_Class ) then
            lozone = .TRUE.
            ncount = ncount + 1
            if( SA_Use_APCA ) lapca = .TRUE.
            if( SA_Use_APCA_Ptoverride ) lapcapt = .TRUE.
         endif
         if( ncount .LE. 0 ) goto 7001
      endif
c
c-----Make sure that if you are treating PM classes with SA
c     that PM is used in this run
c
      if( aeropt(1:4) .EQ. 'NONE' ) then
         if( lsulfate ) goto 7003
         if( lnitrate ) goto 7003
         if( lsoa ) goto 7003
         if( lprimary ) goto 7003
      endif
c
c   --- Need the biogenics group if doing APCA ---
c
      if( lapca .AND. ngroup .EQ. 0 .AND. .NOT.lapcapt) then
         write(iout,'(//,a)') 'ERROR in RDOPTSA:'
         write(iout,'(/,1X,3A)')'Need biogenic sources as a separate ',
     &                          'emissions group when doing APCA.'
         call camxerr()
      endif
c
c   --- do not allow APCA with only one group ---
c
      if( ngroup .EQ. 1 .AND.  lapca .AND. .NOT. leftovr_area ) then
         write(iout,'(//,a)') 'ERROR in RDOPTSA:'
         write(iout,'(/,1X,A,A)')'You asked to use APCA for SA but only',
     &                                           ' provided 1 source group.'
         write(iout,'(1X,A)')'APCA treatment requires at least 2 source groups.'
         call camxerr()
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
 7001 continue
      write(iout,'(//,a)') 'ERROR in RDOPTSA:'
      write(iout,'(/,1X,2A)') 'The SA technology type was ',
     &     'specified but no classes are selected for treatment.'
      write(iout,'(1X,2A)') 'Please activate at least one SA class ',
     &                          'in the CAMx control file.'
      write(iout,'(1X,2A)') 'With this version the SA_Control namelist ',
     &                   'variables have changed significantly.'
      write(iout,'(1X,2A)') 'In particular, the OSAT/APCA options and ',
     &           'PSAT options have been consolidated.' 
      write(iout,'(1X,2A)') 'Include the following for all ',
     &           'source apportionment types: OSAT/APCA/PSAT'
      write(iout,'(1X,A)') '    SA_Treat_SULFATE_Class'
      write(iout,'(1X,A)') '    SA_Treat_NITRATE_Class'
      write(iout,'(1X,A)') '    SA_Treat_SOA_Class'
      write(iout,'(1X,A)') '    SA_Treat_PRIMARY_Class'
      write(iout,'(1X,A)') '    SA_Treat_MERCURY_Class'
      write(iout,'(1X,A)') '    SA_Treat_OZONE_Class'
      write(iout,'(1X,A)') '    SA_Use_APCA'
      write(iout,'(1X,2A)') "See the User's Guide or the namelist template ",
     &           'for details.'
      call camxerr()
c
 7002 continue
      write(iout,'(//,a)') 'ERROR in RDOPTSA:'
      write(iout,'(/,1X,2A)') 'The CAMx model no longer supports ',
     &                                             'timing tracers.'
      write(iout,'(1X,2A)') 'Please remove the Number_of_Timing_Releases ',
     &                          'variable from the SA_Control namelist.'
      call camxerr()
c
 7003 continue
      write(iout,'(//,a)') 'ERROR in RDOPTSA:'
      write(iout,'(2a)') 'At least one PM class is being treated by',
     &                        ' the Source Apportionment (SA).'
      write(iout,'(2a)') 'But the chemistry parameters file does not ',
     &                                        'include PM species.'
      write(iout,'(2a/,a)') 'Please run with a PM chemparam file or',
     &  ' turn off treatment of PM class in ','the SA_Control namelist.'
      call camxerr()
c
 7004 continue
      write(iout,'(//,a)') 'ERROR in RDOPTSA:'
      write(iout,'(2a)') 'In order to run BiDI with ',
     &                 'the Source Apportionment (SA) you must have'
      write(iout,'(a)') 'Biogenics sources as group 1.'
      write(iout,'(2a)') 'Make sure that the group 1 emissions ',
     &                           'files contain Biogenic sources '
      write(iout,'(a)') 'and turn on the SA_Use_APCA flag.'
      call camxerr()
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
