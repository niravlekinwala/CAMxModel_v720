c**** RDTCDDM.F
c
      subroutine rdtcddm()
      use filunit
      use grid
      use chmstry
      use camxfld
      use camxcom
      use tracer
      implicit none
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routine fills one hour of top concentrations for the DDM
c   species. It then places these concentrations in the  appropriate 
c   place in the gridded array used for tracer concentrations.  
c
c      Copyright 1996 - 2022
c     Ramboll
c
c      Argument description:
c       Outputs:
c       Inputs:
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     11/28/16   --bkoo--       Created
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
c-----------------------------------------------------------------------
c  External functions:
c-----------------------------------------------------------------------
c
      integer ncf_chkfile
c
c-----------------------------------------------------------------------
c    Local variables:
c----------------------------------------------------------------------
c
      character*200 action
      character*10 aqfil, cfile, cspec
      character*4  ifile(10), note(60), ispec(10)
      integer      iseg, nspc, idat1, idat2, ierr
      integer      izone, nx, ny, nz, nedge, ioff
      integer      ispc, imod, iddm, iptr, n3d
      integer      ip, ic, ig
      integer      i, j, k, n
      real         tim1, tim2, orgx, orgy, utmx, utmy, dx, dy
      logical      lexist, luse
c
      real contop(MXCELLS,MXCELLS)
c
c-----------------------------------------------------------------------
c    Data statements:
c-----------------------------------------------------------------------
c
      data aqfil /'AIRQUALITY'/
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c  --- set the number of BC edges --
c
      if( lbndry ) then
        nedge = 5
        ioff = IDXBTP
      else
        nedge = 1
        ioff = 1
      endif
c
c  --- open the TC ---
c
      inquire(file=tcfil,exist=lexist)
      if( .NOT. lexist ) goto 7000
      action = 'Opening DDM Top Concentrations file.'
      ierr = ncf_chkfile(iortc,tcfil,action,aqfil)
      if( ierr .EQ. ISUCES ) then
         open(file=tcfil,unit=iortc,status='UNKNOWN',
     &                                   form='UNFORMATTED',ERR=7001)
      else
         return
      endif
c
c  --- read 1st TC header record and check inputs ---
c
      read(iortc,ERR=7002) ifile,note,iseg,nspc,idat1,tim1,idat2,tim2
      write(cfile,'(10A1)') (ifile(i),i=1,10)
      if( cfile .NE. aqfil ) goto 7003
      if( INT(tim2) .EQ. 24 ) then
        idat2 = idat2 + 1
        tim2 = 0.
        if( MOD(idat2,1000) .GT. 365 ) then
           if( MOD(INT(idat2/1000),4) .EQ. 0 ) then
              if( MOD(idat2,1000) .EQ. 367 )
     &                     idat2 = (INT(idat2/1000)+1)*1000 + 1
           else
              idat2 = (INT(idat2/1000)+1)*1000 + 1
           endif
        endif
      endif
      tim1 = 100.*tim1
      tim2 = 100.*tim2
      if( idat1 .GT. date .AND. .NOT. le1day ) goto 7004
      if( idat1 .EQ. date .AND. tim1 .GT. time
     &                              .AND. .NOT. le1day ) goto 7004
c
c  --- read 2nd BC header record and check inputs ---
c
      read(iortc,ERR=7002) orgx,orgy,izone,utmx,utmy,dx,dy,nx,ny,nz
      if( .NOT. llatlon ) then
        dx = dx/1000.
        dy = dy/1000.
      endif
      if( nx .NE. ncol(1) .OR. ny .NE. nrow(1) ) goto 7005
c
c  --- skip the next 2 records ---
c
      read(iortc,ERR=7002)
      read(iortc,ERR=7002)
c
c  --- read the date and time and make sure it is the correct hour ---
c
  111 continue
      read(iortc,ERR=7007,END=7006) idat1, tim1, idat2, tim2
      if( INT(tim2) .EQ. 24 ) then
        tim2 = 0.
        idat2 = idat2 + 1
        if( MOD(idat2,1000) .GT. 365 ) then
           if( MOD(INT(idat2/1000),4) .EQ. 0 ) then
              if( MOD(idat2,1000) .EQ. 367 )
     &                     idat2 = (INT(idat2/1000)+1)*1000 + 1
           else
              idat2 = (INT(idat2/1000)+1)*1000 + 1
           endif
        endif
      endif
      tim1 = tim1*100
      tim2 = tim2*100
      if( (idat1 .LT. date .OR. 
     &       (idat1 .EQ. date .AND. tim1 .LE. time)) .AND.
     &            (idat2 .GT. date .OR. 
     &                (idat2 .EQ. date .AND. tim2 .GT. time))) then
c
c  --- read the concentrations for each species ---
c
         do ispc=1,nspc
            read(iortc,ERR=7007,END=7006) iseg, (ispec(n),n=1,10),
     &                                     ((contop(i,j),i=1,nx),j=1,ny)
            write(cspec,'(10A1)') ispec
            if( cspec .EQ. 'HNO2      ') cspec = 'HONO      '
            if( cspec .EQ. 'HCHO      ' .and. kHCHO.eq.nspec+1 ) 
     &                                             cspec = 'FORM      '
c
c  --- find this species in the modeled species list ---
c
            do imod=1,nspec
               if( cspec .EQ. spname(imod) ) then
                   do iddm=1,nbcddm
                      luse = .FALSE.
                      if( bcddmsp(iddm) .EQ. cspec ) luse = .TRUE. 
                      if( bcddmsp(iddm) .EQ. NAMALL  ) luse = .TRUE. 
                      if( bcddmsp(iddm) .EQ. NAMVOC 
     &                               .AND. lvocsp(imod) ) luse = .TRUE. 
                      if( bcddmsp(iddm) .EQ. NAMNOX 
     &                               .AND. lnoxsp(imod) ) luse = .TRUE. 
                      if( bcddmsp(iddm) .EQ. NAMHRV
     &                               .AND. lhrvoc(imod) ) luse = .TRUE.
c
c  --- if this DDM species matches this modeled species, load it ---
c
                      if( luse ) then
                         iptr = iptddm(imod) + (iddm-1)*nedge + nicddm + ioff - 1
                         do j=1,ny
                         do i=1,nx
                            n3d = i + nx*(j - 1) + nx*ny*(iptr - 1)
c
c  --- load the top concentrations ---
c
                            if( contop(i,j) .GT. bdnl(imod) ) then
                               pttop(n3d) = contop(i,j)
                            else
                               pttop(n3d) = 0.0
                            endif
                         enddo
                         enddo
c
c  --- interpolate top con sens to all nested grids
c
                         if ( ngrid.gt.1 ) then
                            do ip = 1,ngrid
                            do ic = 1,nchdrn(ip)
                               ig = idchdrn(ic,ip)
                               call interp2d(ncol(ip),nrow(ip),1,
     &                                       i1(ig),j1(ig),nmesh(ig),
     &                                       ncol(ig),nrow(ig),
     &                                       pttop( ipsa2d(ip) + ncol(ip)*nrow(ip)*(iptr - 1) ),
     &                                       pttop( ipsa2d(ig) + ncol(ig)*nrow(ig)*(iptr - 1) ))
                            enddo
                            enddo
                         endif
                      endif
c
c  --- check the next BC DDM species ---
c
                   enddo
c
c  --- check the next modeled species ---
c
               endif
            enddo
c
c  --- read next species worth of data ---
c
         enddo
c
c  --- if not the right hour, read through it ---
c
      else
         do ispc=1,nspc
            read(iortc,ERR=7007) 
         enddo
         goto 111
      endif
c
c  --- finally got the right hour, close the file ---
c
      close(iortc)
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue 
      write(iout,'(//,A)') 'ERROR in RDTCDDM:' 
      write(iout,*) 'ERROR:  TC file for DDM does not exist: ',TRIM(tcfil)
      call camxerr()
c
 7001 continue 
      write(iout,'(//,A)') 'ERROR in RDTCDDM:' 
      write(iout,*) 'ERROR:  Opening TC file for DDM: ',TRIM(tcfil)
      call camxerr()
c
 7002 continue
      write(iout,'(//,A)') 'ERROR in RDTCDDM:' 
      write(iout,*) 'ERROR:  Reading header of TC file for DDM: ',TRIM(tcfil)
      call camxerr()
c
 7003 continue
      write(iout,'(//,A)') 'ERROR in RDTCDDM:' 
      write(iout,*) 'ERROR:  TC file for DDM is not labelled ',TRIM(aqfil)
      call camxerr()
c
 7004 continue
      write(iout,'(//,A)') 'ERROR in RDTCDDM:' 
      write(iout,*) 'ERROR:  TC file for DDM is not for correct time period.'
      write(iout,*) '   ***  Episode ***'
      write(iout,'(a,i10.5)') 'Date   : ',date
      write(iout,'(a,f10.1)') 'Time   : ',time
      write(iout,*) '   ***  TC File ***'
      write(iout,'(a,i10.5)') 'Date   : ',idat1
      write(iout,'(a,f10.1)') 'Time   : ',tim1
      call camxerr()
c
 7005 continue
      write(iout,'(//,A)') 'ERROR in RDTCDDM:' 
      write(iout,*) 'ERROR:  TC file for DDM is not for correct grid.'
      write(iout,*) '   ***  Coarse Grid ***'
      write(iout,*) 'No. of Cells : (',ncol(1),',',nrow(1),')'
      write(iout,*) '   ***  TC File ***'
      write(iout,*) 'No. of Cells : (',nx,',',ny,')'
      call camxerr()
c
 7006 continue
      write(iout,'(//,A)') 'ERROR in RDTCDDM:' 
      write(iout,*) 'ERROR:  Premature end-of-file reached in TC file for DDM.'
      write(iout,*) 'Make sure the file contains the correct date/time.'
      call camxerr()
c
 7007 continue
      write(iout,'(//,A)') 'ERROR in RDTCDDM:' 
      write(iout,*) 'ERROR:  Reading TC file for DDM: ',TRIM(tcfil)
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
