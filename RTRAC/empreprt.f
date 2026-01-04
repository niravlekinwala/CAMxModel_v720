*** EMPREPRT
c
      subroutine empreprt(idate,btim,igrid)
      use filunit
      use grid
      use chmstry
      use ptemiss
      use rtracchm
      use tracer
      implicit none
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routine reads the emissions files for the RTRAC and moves the 
c   file pointer to the first hour to be used in this run.
c
c      Copyright 1996 - 2022
c     Ramboll
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c      01/16/02   --gwilson-- original development
c      11/12/03   --cemery -- Changed to handle optional area/pt files, 
c                             and to check Rtrac point source file against 
c                             regular model pt file
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'flags.inc'
c
c-----------------------------------------------------------------------
c    Arguement declarations:
c-----------------------------------------------------------------------
c
      integer   idate
      real      btim
      integer   igrid
c
c-----------------------------------------------------------------------
c    External functions:
c-----------------------------------------------------------------------
c
      integer istrln
c
c-----------------------------------------------------------------------
c    Local parameters:
c-----------------------------------------------------------------------
c
      integer MXOUT
c
      parameter (MXOUT = MAX(MXSPEC,MXTRSP))
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      character*200 fname
      character*10  namein
      character*4   ifname(10), ifnote(60)
      integer       nseg, nspect, ipt, numpts, numpts_lst, igrd, ii, jj
      integer       ispc, ibgdat, iendat, buffer_offset
      integer       idateb, idatee, nxcelt, nycelt, nzcelt
      integer       idum, izonet, iounit, isegm, npoint
      integer       i, j, ndate, seg4(4)
      real          tutmx, tutmy, tdxcel, tdycel
      real          tbegti, tendti, begtim, endtim, dum
      real          txorig, tyorig, ttime
      logical       luse
c
      character*4 iname(10,MXOUT)

      real,    allocatable, dimension(:) :: xlocin
      real,    allocatable, dimension(:) :: ylocin
      real,    allocatable, dimension(:) :: hstkin
      real,    allocatable, dimension(:) :: dstkin
      real,    allocatable, dimension(:) :: tstkin
      real,    allocatable, dimension(:) :: vstkin
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c   --- set the date and time ---
c
      ndate = idate
      ttime = btim/100.0
      if( ttime .EQ. 24.0 ) then
         ndate = ndate + 1
         ttime = 0.
         if( MOD(ndate,1000) .GT. 365 ) then
            if( MOD(INT(ndate/1000),4) .EQ. 0 ) then
               if( MOD(ndate,1000) .EQ. 367 )
     &                     ndate = (INT(ndate/1000)+1)*1000 + 1
            else
               ndate = (INT(ndate/1000)+1)*1000 + 1
            endif
         endif
      endif
c
c   --- only read if surface emissions file is supplied ---
c
      if( .NOT. is_netcdf_iortem(igrid,1,1) .AND. ltemfl(igrid,1,1) 
     &                                           .AND. larsrc ) then
        iounit = iortem(igrid,1,1)
        fname = temfil(igrid,1,1)
c
c   --- rewind the file ----
c
        rewind(iounit,ERR=7005)
c
c   --- read the header of the emissions file ---
c
        read(iounit,ERR=7000) ifname, ifnote, nseg, nspect, idateb, 
     &                                        tbegti, idatee, tendti
        if( nspect .GT. MXSPEC ) then
          write(iout,'(//,A)') 'ERROR in EMPREPRT:'
          write(iout,*) 'A parameter in the camx.prm is not ',
     &                                        'sufficiently large.'
          write(iout,*)
     &                    'Please change the value for parameter: MXTRSP'
          write(iout,*) 'It should be set to a value of at least: ',
     &                                                             nspect
          call flush(iout)
          call camxerr()
        endif
c
        read(iounit,ERR=7000) tutmx, tutmy, izonet, txorig,tyorig, 
     &                           tdxcel, tdycel, nxcelt, nycelt, nzcelt
        read(iounit,ERR=7000) seg4
        read(iounit,ERR=7000) ((iname(i,ispc),i=1,10),ispc=1,nspect)
        buffer_offset_iortem(igrid,1,1) = 0
        if( nxcelt .EQ. ncol(igrid)-2 .AND. nycelt .EQ. nrow(igrid)-2 )
     &                                buffer_offset_iortem(igrid,1,1) = 1
        buffer_offset = buffer_offset_iortem(igrid,1,1)
c
c   --- set up index for species names ---
c
        nspcem(igrid,1,1) = nspect
        do ispc=1,nspcem(igrid,1,1)
           write(namein,'(10A1)') (iname(i,ispc),i=1,10)
           luse = .FALSE.
           do j=1,ntotsp
             if(namein .EQ. ptname(j) ) then
                 idxems(igrid,1,1,ispc) = j
                 luse = .TRUE.
             endif
           enddo
c
c  --- if not found in model species list, write a message ----     
c   
           if( .NOT. luse ) then
                write(idiag,'(1X,4A)') 'Species in RTRAC ',
     &              'surface emissions file: ',namein(:istrln(namein)),
     &                      ' not found in species list ... Skipping.'
           endif
        enddo
c
c   --- read the data until the correct position in the file is reached ---
c
  111   continue
        read(iounit,ERR=7002,END=7003) ibgdat, begtim, iendat, endtim
        if( le1day ) then
           if( INT(begtim) .EQ. INT(ttime) ) goto 222
        else
           if( ibgdat .EQ. ndate .AND. 
     &                    INT(begtim) .EQ. INT(ttime) ) goto 222
        endif
        do ispc=1,nspcem(igrid,1,1)
            read(iounit,ERR=7002,END=7003) idum, 
     &          (dum,i=1,10),((dum,i=1+buffer_offset,ncol(igrid)-buffer_offset),
     &                                j=1+buffer_offset,nrow(igrid)-buffer_offset)
        enddo
        goto 111
  222   continue
        backspace(iounit)
c
      endif
c
c   --- point sources are only for the coarse grid -- grid #1 ---
c
      if( igrid .GT. 1 ) goto 9999
c
c   --- only read if elevated point source emissions file is supplied ---
c
      if( .NOT. ltptfl(1,1) ) goto 9999
      if( .NOT. lptsrc ) goto 9999
      if( is_netcdf_iortpt(1,1) ) goto 9999
      iounit = iortpt(1,1)
      fname = tptfil(1,1)
c
c   --- rewind the file ----
c
      rewind(iounit,ERR=7005)
c
c   --- read the header of the emissions file ---
c
      read(iounit,ERR=7000) ifname, ifnote, nseg, nspect, idateb, 
     &                                         tbegti, idatee, tendti
      if( nspect .GT. MXTRSP ) then
          write(iout,'(//,A)') 'ERROR in EMPREPRT:'
          write(iout,*) 'A parameter in the camx.prm is not ',
     &                                        'sufficiently large.'
          write(iout,*)
     &                    'Please change the value for parameter: MXTRSP'
          write(iout,*) 'It should be set to a value of at least: ',
     &                                                             nspect
          call flush(iout)
          call camxerr()
      endif
c
      read(iounit,ERR=7000) tutmx, tutmy, izonet, txorig, tyorig, 
     &                          tdxcel, tdycel, nxcelt, nycelt, nzcelt
      read(iounit,ERR=7000) seg4
      read(iounit,ERR=7000) ((iname(i,ispc),i=1,10),ispc=1,nspect)
c
c   --- set up index for species names ---
c
      nspcpt(1,1) = nspect
      do ispc=1,nspcpt(1,1)
         write(namein,'(10A1)') (iname(i,ispc),i=1,10)
         luse = .FALSE.
         do j=1,ntotsp
           if(namein .EQ. ptname(j) ) then
               idxpts(1,1,ispc) = j
               luse = .TRUE.
           endif
         enddo
c
c  --- if not found in model species list, write a message ----
c  
         if( .NOT. luse ) then
              write(idiag,'(1X,4A)') 'Species in RTRAC ',
     &            'point emissions file: ',namein(:istrln(namein)),
     &                    ' not found in species list ... Skipping.'
         endif
      enddo
      read(iounit,ERR=7000) isegm, npoint
c
c --- allocate local arrays ---
c
      allocate( xlocin(npoint) )
      allocate( ylocin(npoint) )
      allocate( hstkin(npoint) )
      allocate( dstkin(npoint) )
      allocate( tstkin(npoint) )
      allocate( vstkin(npoint) )
c
      read(iounit,ERR=7002) (xlocin(i),ylocin(i),hstkin(i),
     &                     dstkin(i),tstkin(i),vstkin(i),i=1,npoint)
c
c  --- make sure this source is in the base inventory and flag it ---
c
      numpts_lst = nptsrc - npoint
      do ipt=1,npoint
        vstkin(ipt) = vstkin(ipt)/3600.
        numpts_lst = numpts_lst + 1 
        if (llatlon) then
          xstk(numpts_lst,1) = xlocin(ipt) - xorg
          ystk(numpts_lst,1) = ylocin(ipt) - yorg
        else
          xstk(numpts_lst,1) = xlocin(ipt)/1000. - xorg
          ystk(numpts_lst,1) = ylocin(ipt)/1000. - yorg
        endif
        do igrd = 2,ngrid
           xstk(numpts_lst,igrd) = xstk(numpts_lst,1) - (inst1(igrd) - 1)*delx
           ystk(numpts_lst,igrd) = ystk(numpts_lst,1) - (jnst1(igrd) - 1)*dely
           ii = 2 + FLOOR(xstk(numpts_lst,igrd)/delx*FLOAT( meshold(igrd) ) )
           jj = 2 + FLOOR(ystk(numpts_lst,igrd)/dely*FLOAT( meshold(igrd) ) )
           if (ii .GT. 1 .AND. ii .LT. ncol(igrd) .AND.
     &                    jj .GT. 1 .AND. jj .LT. nrow(igrd)) then
              nosrc(igrd) = nosrc(igrd) + 1
              idsrc(nosrc(igrd),igrd) = numpts_lst
              isrc(nosrc(igrd),igrd) = ii
              jsrc(nosrc(igrd),igrd) = jj
           endif
        enddo
        xlocpt(numpts_lst) = xlocin(ipt)
        ylocpt(numpts_lst) = ylocin(ipt)
        xloc(numpts_lst) = xlocin(ipt)
        yloc(numpts_lst) = ylocin(ipt)
        hstk(numpts_lst) = hstkin(ipt)
        dstk(numpts_lst) = dstkin(ipt)
        tstk(numpts_lst) = tstkin(ipt)
        vstk(numpts_lst) = vstkin(ipt)
        lpiglet(numpts_lst) = .FALSE.
        if( dstk(numpts_lst) .LT. 0 ) lpiglet(numpts_lst) = .TRUE.
        idx_point_in_list(1,1,ipt) = numpts_lst
      enddo
c
c  --- deallocate the local arrays ---
c
      deallocate( xlocin )
      deallocate( ylocin )
      deallocate( hstkin )
      deallocate( dstkin )
      deallocate( tstkin )
      deallocate( vstkin )
c
c   --- read the data until the correct position in the file is reached ---
c
  333 continue
      read(iounit,ERR=7002,END=7003) ibgdat, begtim, iendat, endtim
      if( le1day ) then
         if( begtim .EQ. ttime ) goto 444
      else
         if( ibgdat .EQ. ndate .AND. begtim .EQ. ttime ) goto 444
      endif
      read(iounit,ERR=7002) isegm, npoint
      read(iounit,ERR=7002) (idum,idum,idum,dum,dum,i=1,npoint)
      do ispc=1,nspcpt(1,1)
          read(iounit,ERR=7002,END=7003) idum,(dum,i=1,10),
     &                                            (dum,i=1,npoint)
      enddo
      goto 333
  444 continue
      backspace(iounit)
c
c  ---- return to calling routine ---
c
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,A)') 'ERROR in EMPREPRT:'
      write(iout,'(/,1X,3A)') 'Reading header of ',
     &                     'emissions file: ',fname(:istrln(fname))
      call camxerr()
c
 7002 continue
      write(iout,'(//,A)') 'ERROR in EMPREPRT:'
      write(iout,'(/,1X,2A)') 'Reading emissions grid in file: ',
     &                                        fname(:istrln(fname))
      call camxerr()
c
 7003 continue
      write(iout,'(//,A)') 'ERROR in EMPREPRT:'
      write(iout,'(/,1X,3A)') 'Premature end-of-file reading ',
     &                     'emissions file: ',fname(:istrln(fname))
      call camxerr()
c
 7004 continue
      write(iout,'(//,A)') 'ERROR in EMPREPRT:'
      write(iout,'(1X,A,I10,2A)') 'Location for point: ',i,
     &     ' is not consistent with regular emissions in file: ',
     &                                        fname(:istrln(fname))
      call camxerr()
c
 7005 continue
      write(iout,'(//,A)') 'ERROR in EMPREPRT:'
      write(iout,'(/,1X,3A)') 'Cannot access emissions file',
     &     ' provided for OSAT processing: ',fname(:istrln(fname))
      call camxerr()
c
c-----------------------------------------------------------------------
c    Format statements:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
c
      return
      end
