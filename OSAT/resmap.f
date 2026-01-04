c**** RESMAP
c
      subroutine resmap()
      use filunit
      use grid
      use tracer
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routine reads the mapping file and sets the source region
c   for each grid cell.  The mapping file is an ASCII grid which
c   contains the source area number for each grid. 
c
c      Copyright 1996 - 2022
c     Ramboll
c
c     Argument description:
c      Inputs:
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     05/20/96   --gwilson--    Original development
c     07/19/02   --gwilson--    Added source area map for nests
c     07/16/07   --bkoo--       Added check for HDDM
c     03/01/16   --gwilson--    Added partial source area map
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
c
c-----------------------------------------------------------------------
c    Arguement declarations:
c-----------------------------------------------------------------------
c
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      integer   irec, irow, ierr, igrd, igroup, i
c
      integer irowvl(MXCELLS)
c
c-----------------------------------------------------------------------
c    Enxternal functions:
c-----------------------------------------------------------------------
c
      integer istrln
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c  --- initialize global variables --
c
      npartial = 1
      igrmap = 0
      frcmap = 1.0
c
c  --- loop over grids ---
c
      do 10 igrd=1,ngrid
c
c   --- if mapping file is not supplied, go to next grid ---
c
         if( .NOT. lmapfl(0,igrd) ) goto 10
c
c  ---- open the file ---
c
         open(unit=iormap(0,igrd),file=mapfil(0,igrd),status='UNKNOWN')
c
c  ---- loop over rows (Y direction), go backwards since the ASCII map
c       is a grid image ---
c
         irec = 0
         do 20 irow=nycell(igrd),1,-1
           irec = irec + 1
           read(iormap(0,igrd),'(500(:,I3))',IOSTAT=ierr) 
     &                                (irowvl(i),i=1,nxcell(igrd))
           if( ierr .NE. 0 ) then
               write(iout,'(//,a)') 'ERROR in RESMAP:'
               write(iout,'(/,1X,2A,I4)') 'Reading the source ',
     &                                     'map file at line: ',irec
               write(iout,'(10X,A,/,A)') 'Source map filename: ',
     &                              mapfil(0,igrd)(:istrln(mapfil(0,igrd)))
              call camxerr()
           endif
c
c  --- check each mapping value for validity and load data 
c      into the global array ---
c
           do 30 i=1,nxcell(igrd)
              if( (irowvl(i) .LE. 0 .AND. .NOT. (lddm.OR.lhddm)) .OR.
     &                                     irowvl(i) .GT. nregin ) then
                  write(iout,'(//,a)') 'ERROR in RESMAP:'
                  write(iout,'(/,1X,A,A,I4)') 'Mapping value read',
     &                      ' from map file is out of range: ',irowvl(i)
                  write(iout,'(10X,A,/,A)') 'Source map filename: ',
     &                               mapfil(0,igrd)(:istrln(mapfil(0,igrd)))
                  call camxerr()
              endif
              if( ltrace ) then
                 do igroup=0,ngroup
                    igrmap(igroup,1,igrd,i,irow) = irowvl(i)
                    frcmap(igroup,1,igrd,i,irow) = 1.0
                 enddo
              else
                 igrmap(0,1,igrd,i,irow) = irowvl(i)
                 frcmap(0,1,igrd,i,irow) = 1.0
              endif
   30      continue
c
c  --- read the next row ---
c
   20   continue
c
c  --- close file and return to calling routine ---
c
        close(iormap(0,igrd))
c
c   --- next nest ---
c
   10 continue
c
      goto 9999
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
c
      return
      end
