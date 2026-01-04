      subroutine rd_sa_bcfile( )
      use tracer
      use grid
      use camxfld
      use filunit
      implicit none
c 
c----CAMx v7.20 220430
c 
c     RD_SA_BCFILE reads and cycles through the source apportionment BC
c     file to the current time/date, and loads coarse grid boundary 
c     concentrations. The first time it reads the through the header
c     records to get to the first hour.
c 
c      Copyright 1996 - 2022
c     Ramboll
c           
c     Modifications: 
c        none
c 
c     Input arguments: 
c             
c     Output arguments: 
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
c     Local Variables
c
      character*10 spcname
      integer      idat1, idat2, ncells, nz, i, j, k, lbc, itr, n
      integer      nc, idum, iedge, n3d, n4d
      real         tim1, tim2, convfac
c
      real,        allocatable, dimension(:,:,:,:) :: bctmp
      character*4, allocatable, dimension(:,:)     :: bcspec
c
c-----Entry point
c
      ncells = MAX(ncol(1),nrow(1))
      nz = nlay(1)
      allocate( bctmp(ncells,nz,4,num_iorbc) )
      allocate( bcspec(10,num_iorbc) )
c
c-----If this is the first time, read through the header records ---
c
      if( lrd_sa_bc_hdr ) then
        lrd_sa_bc_hdr = .FALSE.
        rewind(iorbc)
        do i=1,8
          read(iorbc)
        enddo
      endif
c
c-----Read through coarse grid concentration records until current time/date
c
 100  read(iorbc,end=900) idat1,tim1,idat2,tim2
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
      do lbc = 1,num_iorbc
        do n = 1,4
          nc = nrow(1)
          if (n.gt.2) nc = ncol(1)
          read(iorbc) idum,(bcspec(j,lbc),j=1,10),iedge,
     &              ((bctmp(i,k,n,lbc),k=1,nz),i=1,nc) 
        enddo
      enddo
      write(iout,'(a40,2(f7.0,i8.5))')
     &  'Read SA boundary condition file at ',tim1,idat1,tim2,idat2
      call flush(iout)
      if ((idat1.lt.date .or. (idat1.eq.date .and. tim1.le.time)) .and.
     &    (idat2.gt.date .or. (idat2.eq.date .and. tim2.gt.time))) then
c
c-----Load boundary concentrations; convert gasses from ppm to umol/m3,
c     PM stays at ug/m3
c
        do lbc = 1,num_iorbc
          write(spcname,'(10A1)') (bcspec(i,lbc),i=1,10)
          do itr = 1,ntotsp
             if( spcname .NE. ptname(itr) ) cycle
             do k = 1,nz
                do j = 2,nrow(1)-1
                  i = 1
                  n3d = i + ncol(1)*(j - 1) + ncol(1)*nrow(1)*(k - 1)
                  n4d = n3d + ncol(1)*nrow(1)*nlay(1)*(itr - 1)
                  if( lsagas(itr) ) then
                    convfac = densfac*273./tempk(n3d)*press(n3d)/1013.
                  else
                    convfac = 1.
                  endif
                  ptconc(n4d) = MAX(BNDLPT,bctmp(j,k,1,lbc))*convfac
                  i = ncol(1)
                  n3d = i + ncol(1)*(j - 1) + ncol(1)*nrow(1)*(k - 1)
                  n4d = n3d + ncol(1)*nrow(1)*nlay(1)*(itr - 1)
                  if( lsagas(itr) ) then
                    convfac = densfac*273./tempk(n3d)*press(n3d)/1013.
                  else
                    convfac = 1.
                  endif
                  ptconc(n4d) = MAX(BNDLPT,bctmp(j,k,2,lbc))*convfac
                enddo
c
                do i = 2,ncol(1)-1
                  j = 1
                  n3d = i + ncol(1)*(j - 1) + ncol(1)*nrow(1)*(k - 1) 
                  n4d = n3d + ncol(1)*nrow(1)*nlay(1)*(itr - 1)
                  if( lsagas(itr) ) then
                    convfac = densfac*273./tempk(n3d)*press(n3d)/1013.
                  else
                    convfac = 1.
                  endif
                  ptconc(n4d) = MAX(BNDLPT,bctmp(i,k,3,lbc))*convfac
c
                  j = nrow(1)
                  n3d = i + ncol(1)*(j - 1) + ncol(1)*nrow(1)*(k - 1) 
                  n4d = n3d + ncol(1)*nrow(1)*nlay(1)*(itr - 1)
                  if( lsagas(itr) ) then
                    convfac = densfac*273./tempk(n3d)*press(n3d)/1013.
                  else
                    convfac = 1.
                  endif
                  ptconc(n4d) = MAX(BNDLPT,bctmp(i,k,4,lbc))*convfac
                enddo
             enddo
           enddo
        enddo
        goto 999
      else
        goto 100
      endif
c
c-----Set next boundary update time
c
      sa_bndtim = tim2
      sa_bnddate = idat2
c
c-----End of BC file reached
c
 900  continue
      write(iout,'(//,a)') 'ERROR in RD_SA_BCFILE :'
      write(iout,*)'Premature End of source apportionment BC file reached.'
      write(iout,*)'Make sure boundary file contains simulation ',
     &                                               'time period.'
      call camxerr()
c
 999  continue
c
c  --- dealocate the local array ---
c
      deallocate( bctmp )
      deallocate( bcspec )
      return
      end
