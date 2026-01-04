c***** RDARRT.F
c
      subroutine rdarrt(igrid,ndate,ttime,nox,noy,
     &                                  nlay_ems,nrtarsp,emisrt)
      use filunit
      use grid
      use chmstry
      use bndary
      use camxcom
      use rtracchm
      use tracer
c
c----CAMx v7.20 220430
c
c
c      Copyright 1996 - 2022
c     Ramboll
c
c-----------------------------------------------------------------------
c
c   This routine reads one hour of emissions for the RTRAC process
c   and fills the approproate arrays.  The emissions file for one grid
c   is read.
c    Argument descriptions:
c     Outputs:
c       emisrt    R    array of emissions for RTRAC speices
c     Inputs:
c       igrid     I    grid number
c       ndate     I    julian day of current hour
c       ttime     R    current hour
c       nox       I    number of columns in grid
c       noy       I    number of rows in grid
c       nlays_ems I    number of layers in emissions files
c       nrtarsp   I    number of RTRAC species
c       
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c       10/24/01  Removed BSWAP and converted integer strings to character*4
c
c-----------------------------------------------------------------------
c   Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'flags.inc'
c       
c-----------------------------------------------------------------------
c   Argument declarations:
c-----------------------------------------------------------------------
c
      integer igrid
      integer ndate
      real    ttime
      integer nox
      integer noy
      integer nlay_ems
      integer nrtarsp
      real    emisrt(nox,noy,nlay_ems,nrtarsp)
c
c-----------------------------------------------------------------------
c   Local variables:
c-----------------------------------------------------------------------
c
      character*200 fname
      character*4   iname(10)
      integer       ibgdat, iendat, iounit, idx, iseg
      integer       i, j, ispc, buffer_offset
      real          emstmp, bgtim, edtim
      logical       lfound, lpass
c
      real, allocatable, dimension(:,:) :: emsgrd
c
c-----------------------------------------------------------------------
c   Entry point:
c-----------------------------------------------------------------------
c
c   --- initialize emissions to zero ---
c
      do i=1,nox
         do j=1,noy
            do l=1,nrtarsp
               emisrt(i,j,1,l) = 0.
            enddo
         enddo
      enddo
c
c   --- skip if filename not supplied ---
c
      if( .NOT. ltemfl(igrid,1,1) .OR. .NOT. larsrc ) goto 9999
      if( is_netcdf_iortem(igrid,1,1) ) goto 9999
c
c   --- set the unit number for surface emissions file ---
c
      iounit = iortem(igrid,1,1)
      fname = temfil(igrid,1,1)
      buffer_offset = buffer_offset_iortem(igrid,1,1)
c
c   --- allocate the local array ---
c
      allocate( emsgrd(nox-2*buffer_offset,noy-2*buffer_offset) )
c
c   --- read the date and time, again ---
c
      lfound = .FALSE.
      lpass = .FALSE.
  111 continue
      read(iounit,END=222) ibgdat, bgtim, iendat, edtim
      bgtim = 100.*aint(bgtim) + 60.*amod(bgtim,1.)
      edtim = 100.*aint(edtim) + 60.*amod(edtim,1.)
c
c   --- read the emissions for this hour ---
c
      do 20 ispc=1,nspcem(igrid,1,1)
         read(iounit,ERR=7000) iseg, (iname(i),i=1,10),((emsgrd(i,j),
     &                   i=1,nox-2*buffer_offset),j=1,noy-2*buffer_offset)
         do i=1,nox
           emsgrd(i,1) = emsgrd(i,1+buffer_offset)
           emsgrd(i,noy) = emsgrd(i,noy-buffer_offset)
         enddo
         do j=1,noy
           emsgrd(1,j) = emsgrd(1+buffer_offset,j)
           emsgrd(nox,j) = emsgrd(nox-buffer_offset,j)
         enddo
c
c   --- if date and time does not match this hour, skip this record ---
c
         if( le1day ) then
           if( abs(bgtim-ttime).ge.0.01 ) goto 20
         else
           if( ndate .NE. ibgdat .OR. abs(bgtim-ttime).ge.0.01 ) goto 20
         endif
         lfound  = .TRUE.
c
c   --- if the species is a not modeled or not used for RTRAC ----
c
         idx = idxems(igrid,1,1,ispc)
         if( idx .LE. 0 ) goto 20
c
c  --- if this RTRAC species matches this modeled species, load it ---
c
         do j=1+buffer_offset,noy-buffer_offset
            jtmp = j-buffer_offset
            do i=1+buffer_offset,nox-buffer_offset
               itmp = i-buffer_offset
c
c  --- convert to emissions time and put into array ---
c
               emstmp = emsgrd(itmp,jtmp)/(60.*dtems)
               emisrt(i,j,1,idx) = emstmp
c
c  --- next cell ---
c
            enddo
         enddo
c
c   --- next species ---
c
  20  continue
c
c   --- if the correct hour has not been found, 
c       go back and read some more else read next file ---
c
      if( .NOT. lfound ) then
         goto 111
      else
         do i=1,nspcem(igrid,1,1) + 1
            backspace(iounit)
         enddo
         goto 9999
      endif
c
c   --- if using 1 day emissions, we need to rewind the file to
c       get the current hour ---
c
  222 continue
      if( le1day ) then
         rewind(iounit)
         read(iounit)
         read(iounit)
         read(iounit)
         read(iounit)
         if( lpass ) goto 7001
         if( .NOT. lpass ) lpass = .TRUE.
         goto 111
      else
         goto 7001
      endif
c
c  --- get the next file ---
c
       deallocate( emsgrd )
       goto 9999
c
c----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,A)') 'ERROR in RDARRT:'
      write(iout,'(/,1X,2A,I8.5,F8.1,2A)') 'ERROR: Reading emissions ',
     &   'after hour ',ibgdat, bgtim,' in file: ',fname(:istrln(fname))
      call camxerr()
c
 7001 continue
      write(iout,'(//,A)') 'ERROR in RDARRT:'
      write(iout,'(/,1X,2A,I8.5,F8.1,2A)') 
     &              'ERROR: Premature end-of-file',
     &              ' in emissions file after hour ',ibgdat, bgtim,
     &                              ' in file: ',fname(:istrln(fname))
      call camxerr()
c
c----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
c
      return
      end
