      subroutine rdpartial_map()
      use filunit
      use grid
      use tracer
      implicit none
c
c----CAMx v7.20 220430
c
c     This routione reads the partial source area map file and loads 
c     all of the data in the panels into the global arrays to be used
c     to aggregate the emissions into tracer regions.
c
c      Copyright 1996 - 2022
c     Ramboll
c          
c     Modifications:
c
c     Input arguments:
c
c     Output arguments:
c
c     Routines called:
c
c     Called by:
c
      include "camx.prm"
c
c-----External Functions
c
      integer istrln
c
c-----Local variables ----
c
      character*(KEYLEN) keywrd
      integer            irowvl(MXCELLS), npart
      integer            ierr, ipart, igrp, igrd, irow, irec, jerr, i
      real               frcrow(MXCELLS)
c
c-----Entry point
c
c  --- loop over grrids ---
c
      do igrd=1,ngrid
c
c  --- loop over groups, initializing the count ---
c
       do 10 igrp=1,ngroup
         npart = 0
c
c  --- skip if file is not supplied ---
c
         if( .NOT. lmapfl(igrp,igrd) ) goto 10
c  
c  ---- open the file ---
c        
         open(unit=iormap(igrp,igrd),file=mapfil(igrp,igrd),status='UNKNOWN')
c
c-----write some information to the diag file ---
c
         write(idiag,'(A,I3)') 'Parsing Partial Source Map for Grid:',igrd
c
c  --- loop back here until out of map panels ---
c
  111    continue   
         npart = npart + 1
         if( npart .GT. MXPARTIAL ) then
             write(iout,'(//,a)') 'ERROR in RDPARTIAL_MAP:'
             write(iout,*) 'The number of partial source maps exceeds the maximum: ',MXPARTIAL
             write(iout,*) 'Increase MXPARTIAL and recompile'
             call camxerr()
         endif
         ipart = npart
c
c  --- create the keyword for this group and map ---
c
         write(keywrd,'(A,I2.2,A,I2.2,A)') '/SRCMAP',igrp,'-',ipart,'/'
c
c  --- call routine to find the keyword for this map ---
c
         call fndkey(ierr,iormap(igrp,igrd),keywrd)
c
c  --- read error, exit ---
c
         if( ierr .EQ. IRDERR ) then
             write(iout,'(//,a)') 'ERROR in RDPARTIAL_MAP'
             write(iout,'(/,1X,2A,I2,A,I2)') 'Reading the partial source ',
     &                    'area map file for grid: ',igrd,' group: ',igrp
             write(iout,'(10X,A,/,A)') 'Source map filename: ',
     &                       mapfil(igrp,igrd)(:istrln(mapfil(igrp,igrd)))
            call camxerr()
c
c  --- end of file just means packet was not found, no more panels
c      for this group/grid ---
c
         else if( ierr .EQ. IEOF ) then
             npart = npart - 1
c
c  --- panel was found, read and store the data ---
c
         else if( ierr .EQ. ISUCES ) then
             write(idiag,'(2A)') 'Parsing panel: ',keywrd(:istrln(keywrd))
             irec = 0
             do irow=nycell(igrd),1,-1
                irec = irec + 1
                read(iormap(igrp,igrd),'(500(:,I3,1X,F5.0))',IOSTAT=jerr)
     &                                (irowvl(i),frcrow(i),i=1,nxcell(igrd))
                if( jerr .NE. 0 ) then
                   write(iout,'(//,a)') 'ERROR in RDPARTIAL_MAP:'
                   write(iout,'(/,1X,2A,I4)') 'Reading the partial source ',
     &                                     'area map file: '
                   write(iout,'(/,1X,3A,I5)') 'Parsing data in panel: ',
     &                         keywrd(:istrln(keywrd)),' at record: ',irec
                   write(iout,'(10X,A,/,A)') 'Source map filename: ',
     &                         mapfil(igrp,igrd)(:istrln(mapfil(igrp,igrd)))
                   call camxerr()
                endif
c
c  --- check that region values are valid ---
c
                do i=1,nxcell(igrd)
                   if( (irowvl(i) .LE. 0 .AND. frcrow(i) .GT. 0.) 
     &                                   .OR. irowvl(i) .GT. nregin ) then
                      write(iout,'(//,a)') 'ERROR in RDPARTIAL_MAP:'
                      write(iout,'(/,1X,A,A,I4)') 'Mapping value read',
     &                      ' from map file is out of range: ',irowvl(i)
                      write(iout,'(/,1X,3A,I5)') 'Parsing data in panel: ',
     &                         keywrd(:istrln(keywrd)),' at record: ',irec
                      write(iout,'(10X,A,/,A)') 'Source map filename: ',
     &                               mapfil(igrp,igrd)(:istrln(mapfil(igrp,igrd)))
                      call camxerr()
                   endif
c
c  --- load the data into global array ---
c
                   igrmap(igrp,ipart,igrd,i,irow) = irowvl(i)
                   frcmap(igrp,ipart,igrd,i,irow) = frcrow(i)/100.
                enddo
            enddo
c
c  --- get next panel ---
c
            goto 111
         endif
c
c  ---- set global array value ---
c
         if( npart .GT. 0 ) npartial(igrp,igrd) = npart
  10   continue
       write(idiag,'(/)')
c
c ---- next grid ---
c
      enddo
c
      return
      end
