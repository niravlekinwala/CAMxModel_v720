      subroutine sumgrps(numcols,numrows,nspmod,nsptrac,ndlast,ttlast,
     &                      emstot,emslft,emsbas,emsoth,emssum,lemit)
      use filunit
      use grid
      use chmstry
      use bndary
      use camxcom
      use ptemiss
      use tracer
c
c      Copyright 1996 - 2022
c     Ramboll
c
c
c----CAMx v7.20 220430
c
c     
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c     
c     05/01/03  Time span of emissions must now match emiss update interval
c     10/28/09  Changed dimension of variables to accomodate the
c               dynamic memory allocation
c     12/20/13   -gwilson-  Added the compact point source file
c     03/01/16   -gwilson-  Added partial source area map
c     
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'flags.inc'
c
c-----------------------------------------------------------------------
c    Argument Declarations:
c-----------------------------------------------------------------------
c
      integer numcols
      integer numrows
      integer nspmod
      integer nsptrac
      integer ndlast
      real    ttlast
      real*8  emstot(numcols,numrows,nspmod)
      real*8  emslft(numcols,numrows,nspmod)
      real    emsbas(nspmod,nsptrac)
      real    emsoth(nspmod,nsptrac)
      real    emssum(nspmod,nsptrac)
      logical lemit(*)
c
c-----------------------------------------------------------------------
c    Variables:
c-----------------------------------------------------------------------
c
      character*200 fname
      character*10  cname
      character*4   iname(10)
      integer       ibgdat, iendat, iounit, igroup, idxfile, idx
      integer       npts_in_point, iseg, npoint, idum, ispc, i, j
      integer       igrid, num_emsfiles, num_ptsfiles, buffer_offset
      real          bgtim, edtim
      logical       lpass
c
      real,    allocatable, dimension(:,:,:) ::  emsgrd
      real,    allocatable, dimension(:)     ::  emspnt
      integer, allocatable, dimension(:)     ::  idcompact
      integer, allocatable, dimension(:)     ::  izcel_in
      integer, allocatable, dimension(:)     ::  izcel
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c   --- loop over all of the groups ----
c
      allocate( emspnt(nptsrc) )
      allocate( idcompact(nptsrc) )
      allocate( izcel(nptsrc) )
      allocate( izcel_in(nptsrc) )
      izcel = 0
c
      do 20 igroup=0,ngroup
c
c   --- only process if filename is supplied for group ---
c
         do 25 igrid=1,ngrid
            num_emsfiles = num_iortem(igrid,igroup)
            if( igroup .EQ. 0 ) num_emsfiles = nemiss_files(igrid)
            do 26 idxfile=1,num_emsfiles
              if( allocated(emsgrd) ) deallocate( emsgrd )
              if( .NOT. larsrc) goto 26
c
c   --- set the unit number for surface emissions file ---
c
              lpass = .FALSE.
              if( igroup .EQ. 0 ) then
                  if( is_netcdf_iarem(igrid,idxfile) ) goto 26
                  iounit = iarem(igrid,idxfile)
                  write(fname,'(A,I3)') 'EMISSIONS -- UNIT ',
     &                                                iarem(igrid,idxfile)
                  buffer_offset = buffer_offset_iarem(igrid,idxfile)
              else
                  if( is_netcdf_iortem(igrid,igroup,idxfile) ) goto 26
                  iounit = iortem(igrid,igroup,idxfile)
                  fname = temfil(igrid,igroup,idxfile)
                  buffer_offset = buffer_offset_iortem(igrid,igroup,idxfile)
              endif
              allocate( emsgrd(ncol(igrid),nrow(igrid),1) )
c
c   --- read the date and time, again ---
c
  111         continue
              read(iounit,END=333) ibgdat, bgtim, iendat, edtim
              if( INT(edtim) .EQ. 0 ) then
                edtim = 24.
                iendat = ibgdat
              endif
              if( le1day ) iendat = ndlast
              ichktm1 = NINT( 1000*(bgtim) )
              if( le1day ) then
                ichktm2 = NINT( 1000*(edtim) )
              else
              ichktm2 = NINT( 1000*(edtim)+24000*(iendat-ibgdat) )
              endif
              if( NINT(edtim) .EQ. 0 ) ichktm2 = 24000
              ichkems = NINT( 1000*(dtems/60.) )
              if((ichktm2 - ichktm1) .NE. ichkems ) then
                write(iout,'(//,a)')'ERROR in SUMGRPS:'
                write(iout,*) 'Time interval in surface emissions file does'
                write(iout,*) ' not match emissions update time interval.'
                write(iout,*) '  Beginning Date/Time (Hour): ',ibgdat,bgtim
                write(iout,*) '  Ending Date/Time    (Hour): ',iendat,edtim
                write(iout,*) '  Emiss Input interval(Hour): ',dtems/60.
                call camxerr()
              endif
              bgtim = 100.*aint(bgtim) + 60.*amod(bgtim,1.)
              edtim = 100.*aint(edtim) + 60.*amod(edtim,1.)
c
c   --- read the emissions for this hour ---
c
              do 30 ispc=1,nspcem(igrid,igroup,idxfile)
                  read(iounit) iseg, (iname(i),i=1,10), 
     &               ((emsgrd(i,j,1),i=1+buffer_offset,ncol(igrid)-buffer_offset),
     &                                 j=1+buffer_offset,nrow(igrid)-buffer_offset)
c
c   --- if the species is a not modelled or not a VOC species skip it ---
c
                  idx = idxems(igrid,igroup,idxfile,ispc)
                  if( idx .LE. 0 ) goto 30
                  call sum1grd(numcols,numrows,ncol(igrid),nrow(igrid),1,
     &                      nspmod,nsptrac,igroup,igrid,idx,emssum,emsgrd,
     &                                   emsbas,emsoth,emslft,emstot,lemit)
c
c   --- next species ---
c
  30          continue
c
c   --- check if end of simulation read ---
c
              if( iendat .LT. ndlast  .OR. (iendat .EQ. ndlast .AND. 
     &                          INT(edtim) .LT. INT(ttlast)) ) goto 111
              goto 26
c
c   --- if using 1 day emissions, we might need to keep going
c       to finish out the simulation ----
c
  333          continue
               if( le1day ) then
                    if( lpass ) goto 7005
                    lpass = .TRUE.
                    rewind(iounit)
                    read(iounit)
                    read(iounit)
                    read(iounit)
                    read(iounit)
                    goto 111
                 endif
   26       continue 
            if( allocated(emsgrd) ) deallocate( emsgrd )
   25    continue 
c
c   --- no longer a zero group for point sources ---
c
         if( igroup .EQ. 0 ) cycle
c
c   --- only process if filename is supplied for group ---
c
         do idxfile=1,num_iortpt(igroup)
            if( .NOT. lptsrc .OR. .NOT. ltptfl(igroup,idxfile) ) cycle
c
c   --- set the unit number for elevated points emissions file ---
c
            lpass = .FALSE.
            if( is_netcdf_iortpt(igroup,idxfile) ) cycle
            iounit = iortpt(igroup,idxfile)
            fname = tptfil(igroup,idxfile)
c
c   --- read the date and time, again ---
c
  222       continue
            read(iounit,END=444) ibgdat, bgtim, iendat, edtim
            if( INT(edtim) .EQ. 0 ) then
                edtim = 24.
                iendat = ibgdat
            endif
            if( le1day ) iendat = ndlast
            ichktm1 = NINT( 1000*(bgtim) )
            if( le1day ) then
              ichktm2 = NINT( 1000*(edtim) )
            else
              ichktm2 = NINT( 1000*(edtim)+24000*(iendat-ibgdat) )
            endif
            if( NINT(edtim) .EQ. 0 ) ichktm2 = 24000
            ichkems = NINT( 1000*(dtems/60.) )
            if( (ichktm2 - ichktm1) .NE. ichkems ) then
               write(iout,'(//,a)')'ERROR in SUMGRPS:'
               write(iout,*) 'Time interval in surface emissions file does'
               write(iout,*) ' not match emissions update time interval.'
               write(iout,*) '  Beginning Date/Time (Hour): ',ibgdat,bgtim
               write(iout,*) '  Ending Date/Time    (Hour): ',iendat,edtim
               write(iout,*) '  Emiss Input interval(Hour): ',dtems/60.
               call camxerr()
            endif
            bgtim = 100.*aint(bgtim) + 60.*amod(bgtim,1.)
            edtim = 100.*aint(edtim) + 60.*amod(edtim,1.)
c
c   --- read the emissions for this hour ---
c
            read(iounit,ERR=7000,END=7001) iseg, npoint
            if( npoint .GT. nptsrc ) goto 7002
            read(iounit,ERR=7000,END=7001)  (idcompact(i), idum,
     &                                izcel_in(i), rdum, rdum,i=1,npoint)
            if( .NOT. lcompactpt(igroup,idxfile) ) then
               do i=1,npoint
                  idcompact(i) = i
                  izcel(i) = izcel_in(i)
               enddo
            else
               do i=1,npoint
                  izcel(idcompact(i)) = izcel_in(i)
               enddo
            endif
            do 60 ispc=1,nspcpt(igroup,idxfile)
               read(iounit) iseg, (iname(i),i=1,10), 
     &                                        (emspnt(i),i=1,npoint)
c
c   --- if the species is a not modelled or not a VOC species skip it ---
c
               idx = idxpts(igroup,idxfile,ispc)
               if( idx .LE. 0 ) goto 60
c
c   --- sum up the emissions for each point ---
c
               call sum1pnt(numcols,numrows,nspmod,nsptrac,igroup,idxfile,
     &                       idx,nptsrc,emsbas,emsoth,emslft,
     &                          emstot,emspnt,emssum,izcel,idcompact,lemit)
c
c   --- next species ---
c
  60        continue
            if( iendat .LT. ndlast .OR. (iendat .EQ. ndlast .AND. 
     &                          INT(edtim) .LT. INT(ttlast)) ) goto 222
c
            cycle
c
c   --- if using 1 day emissions, we might need to keep going
c       to finish out the simulation ----
c
 444        continue
            if( le1day ) then
                if( lpass ) goto 7005
                lpass = .TRUE.
                rewind(iounit)
                read(iounit)
                read(iounit)
                read(iounit)
                read(iounit)
                goto 222
            endif
c
c   --- next species ---
c
         enddo
c
c  --- get the next file ---
c
         if( igroup .EQ. 0 ) npts_in_list = npts_in_list + npoint
  20  continue
      deallocate( emspnt )
      deallocate( idcompact )
      deallocate( izcel )
      deallocate( izcel_in )
c
      return
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,a)') 'ERROR in SUMGRPS:'
      write(iout,'(/,1X,A,I10.5,F10.1,2A)') 
     &      'Reading emissions after hour ',ibgdat, bgtim,
     &      ' in file: ',fname(:istrln(fname))
      call camxerr()
c
 7001 continue
      write(iout,'(//,a)') 'ERROR in SUMGRPS:'
      write(iout,'(/,1X,3A)') 'Premature end-of-file reading ',
     &                    'emissions from file: ',fname(:istrln(fname))
      call camxerr()
c
 7002 continue
      write(iout,'(//,a)') 'ERROR in SUMGRPS:'
      write(iout,'(/,1X,A,I10,2A)') 'Number of points: ',npoint,
     &                   ' exceeds max in file: ',fname(:istrln(fname))
      write(iout,'(1X,2A,/,1X,A)') 'Make sure the Probing Tools ',
     &              'point source files are consistent with regular',
     &                                                   ' model file.'
      call camxerr()
c
 7005 continue
      write(iout,'(//,a)') 'ERROR in SUMGRPS:'
      write(iout,'(/,1X,2A,/,10X,2A)') 'Emissions file does cover',
     &    ' entire simulation period.','File = ',fname(:istrln(fname))
      call camxerr()
c
      end
