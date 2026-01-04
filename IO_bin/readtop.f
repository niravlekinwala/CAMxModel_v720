      subroutine readtop(toptim,topdate)
      use filunit
      use grid
      use chmstry
      use bndary
      use camxfld
      use camxcom
c 
c----CAMx v7.20 220430
c 
c     READTOP reads and cycles through the TOPCON file to the
c     current time/date, and loads master grid top concentrations.
c     It then interpolates to all nested grids.
c 
c      Copyright 1996 - 2022
c     Ramboll
c           
c     Modifications: 
c        none
c 
c     Input arguments: 
c        toptim                 model simulation time (HHMM)
c        topdate                model simulation date (YYJJJ)
c             
c     Output arguments: 
c        toptim                 next boundary update time (HHMM)
c        topdate                next boundary update date (YYJJJ)
c             
c     Routines Called: 
c        none
c             
c     Called by: 
c        CAMx
c
      include 'camx.prm'
c
      character*4 tcspec(10)
      integer topdate, ncells
c
      real, allocatable, dimension(:,:,:) :: tctmp
c
c-----Entry point
c
      allocate( tctmp(ncol(1),nrow(1),MXSPEC) )
c
c-----Read through records until current time/date
c
 100  read(itc,end=900) idat1,tim1,idat2,tim2
      if (INT(tim2) .EQ. 24) then
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
      tim1 = 100.*tim1
      tim2 = 100.*tim2
      do l = 1,ntcspc
         read(itc) idum,(tcspec(m),m=1,10),
     &             ((tctmp(i,j,l),i=1,ncol(1)),j=1,nrow(1)) 
      enddo
      write(iout,'(a40,2(f7.0,i8.5))')
     &  'Read top concentrations file at ',tim1,idat1,tim2,idat2
      call flush(iout)
      if ((idat1.lt.date .or. (idat1.eq.date .and. tim1.le.time)) .and.
     &    (idat2.gt.date .or. (idat2.eq.date .and. tim2.gt.time))) then
c
c-----Load top concentrations for master grid
c
        do 60 l = 1,nspec
          ltc = ltcmap(l)
          do 50 j = 1,nrow(1)
            do 40 i = 1,ncol(1)
              n3d = i + ncol(1)*(j - 1) + ncol(1)*nrow(1)*(l - 1)
              ctop(n3d) = bdnl(l)
              if (ltc.gt.0) then
                 if( tctmp(i,j,ltc).gt.bdnl(l) ) 
     &               ctop(n3d) = tctmp(i,j,ltc)
              endif
c
 40         continue
 50       continue
 60     continue
      else
        goto 100
      endif
c
c-----Interpolate topcon to all nested grids
c
      if (ngrid.gt.1) then
        do ip = 1,ngrid
          do ic = 1,nchdrn(ip)
            ig = idchdrn(ic,ip)
            call intrpcnc(nspec,ncol(ip),nrow(ip),1,i1(ig),
     &                    j1(ig),nmesh(ig),ncol(ig),nrow(ig),1,
     &                    ctop(iptr1lay(ip)),ctop(iptr1lay(ig)))
          enddo
        enddo
      endif
c
c-----Set next topcon update time
c
      toptim = tim2
      topdate = idat2
      if (toptim.ge.2400.) then
        toptim = toptim - 2400.
        topdate = topdate + 1
        if( MOD(topdate,1000) .GT. 365 ) then
           if( MOD(INT(topdate/1000),4) .EQ. 0 ) then
              if( MOD(topdate,1000) .EQ. 367 )
     &                     topdate = (INT(topdate/1000)+1)*1000 + 1
           else
              topdate = (INT(topdate/1000)+1)*1000 + 1
           endif
        endif
      endif
      goto 999
c
c-----End of TC file reached
c
 900  write(iout,'(//,a)') 'ERROR in READTOP:'
      write(iout,*)'Premature End of TC file reached.'
      write(iout,*)'Make sure topcon file contains simulation ',
     &                                               'time period.'
      call camxerr()
c
 999  continue
c
c  --- dealocate the local array ---
c
      deallocate( tctmp )
      return
      end
