      subroutine rd_sa_tcfile(bndtim,bnddate)
      use tracer
      use grid
      use camxfld
      use filunit
      implicit none
c 
c----CAMx v7.20 220430
c 
c     RD_SA_TCFILE reads and cycles through the source apportionment TC
c     file to the current time/date, and loads top concentrations into
c     the tracer array.
c 
c      Copyright 1996 - 2022
c     Ramboll
c           
c     Modifications: 
c        none
c 
c     Input arguments: 
c        bndtim                 model simulation time (HHMM)
c        bnddate                model simulation date (YYJJJ)
c             
c     Output arguments: 
c        bndtim                 next boundary update time (HHMM)
c        bnddate                next boundary update date (YYJJJ)
c             
c     Routines Called: 
c        none
c             
c     Called by: 
c        CAMx
c
      include 'camx.prm'
      include 'camx.inc'
c
c     Argument delcarations:
c
      integer bnddate
      real    bndtim
c
c     Local Variables
c
      character*10 spcname
      character*4  tcspec(10)
      integer      idat1, idat2, ncells, nz, i, j, k, ltc, itr, n
      integer      nc, idum, n3d
      real         tim1, tim2
c
      real,    allocatable, dimension(:,:) :: tctmp
c
c-----Entry point
c
      nz = nlay(1)
      allocate( tctmp(ncol(1),nrow(1)) )
c
c-----If this is the first time, read through the header records ---
c
      if( lrd_sa_tc_hdr ) then
        lrd_sa_tc_hdr = .FALSE.
        rewind(iortc)
        do i=1,8
          read(iortc)
        enddo
      endif
c
c-----Read through coarse grid concentration records until current time/date
c
 100  read(iortc,end=900) idat1,tim1,idat2,tim2
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
      write(iout,'(a40,2(f7.0,i8.5))')
     &  'Read SA top concentrations condition file at ',tim1,idat1,tim2,idat2
      call flush(iout)
      if ((idat1.lt.date .or. (idat1.eq.date .and. tim1.le.time)) .and.
     &    (idat2.gt.date .or. (idat2.eq.date .and. tim2.gt.time))) then
c
c-----read the data for each species and load into the proper place in SA \
c     grdded concnetration array ---
c
         do ltc = 1,num_iorbc
            read(iortc) idum,(tcspec(n),n=1,10),
     &                  ((tctmp(i,j),i=1,ncol(1)),j=1,nrow(1))
            write(spcname,'(10A1)') (tcspec(n),n=1,10)
            do itr = 1,ntotsp
               if( spcname .NE. ptname(itr) ) cycle
               do j=1,nrow(1)
                  do i = 1,ncol(1)
                    n3d = i + ncol(1)*(j-1) + ncol(1)*nrow(1)*(itr - 1)
                    pttop(n3d) = BNDLPT
                    pttop(n3d) = amax1(ptconc(n3d),tctmp(i,j))
                  enddo
                enddo
            enddo
         enddo
         goto 999
      else
        goto 100
      endif
c
c-----End of TC file reached
c
 900  continue
      write(iout,'(//,a)') 'ERROR in RD_TA_BCFILE :'
      write(iout,*)'Premature End of source apportionment TC file reached.'
      write(iout,*)'Make sure boundary file contains simulation ',
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
