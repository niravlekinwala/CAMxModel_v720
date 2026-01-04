      subroutine rd_sa_icfile()
      use tracer
      use grid
      use camxfld
      use chmstry
      use filunit
      implicit none
c     
c----CAMx v7.20 220430
c
c     Reads the source apportionment IC file and loads the values into
c     the tracer concentration array.
c 
c      Copyright 1996 - 2022
c     Ramboll
c           
c     Modifications: 
c  
c     Input arguments: 
c        none
c             
c     Output arguments: 
c        none
c             
c     Routines Called: 
c        INTRPCNC
c             
c     Called by: 
c        CAMx
c
      include 'camx.prm'
      include 'camx.inc'
c
c-----Local variables
c
      character*10 spcname
      character*4  icspec(10)
      integer      nx, ny, nz, idat1, idat2, i, j, k, n
      integer      n3d, n4d, ip, ic, ig, igrd, ispc, idum, itr
      real         tim1, tim2, convfac
c
      real cinit(MXCELLS,MXCELLS)
c
c-----Entry point
c
      nx = ncol(1)
      ny = nrow(1)
      nz = nlay(1)
c
c-----Read through the header records to get to the first timestep
c
      rewind(ioric)
      do i=1,4
        read(ioric) 
      enddo
c
c-----Read through coarse grid concentration records until current time/date
c
 100  continue
      read(ioric,end=900) idat1,tim1,idat2,tim2
      tim1 = 100.*tim1
      tim2 = 100.*tim2
      write(iout,'(a40,2(f7.0,i8.5))') 
     &      'Read SA initial condition file at ',tim1,idat1,tim2,idat2
      call flush(iout)
      if ((idat1.lt.date .or. (idat1.eq.date .and. tim1.le.time)) .and.
     &    (idat2.gt.date .or. (idat2.eq.date .and. tim2.ge.time))) then
        do ispc = 1,num_ioric
          do k = 1,nz
            read(ioric) idum,(icspec(n),n=1,10), 
     &                  ((cinit(i,j),i=1,nx),j=1,ny) 
            write(spcname,'(10A1)') (icspec(n),n=1,10)
            do itr = 1,ntotsp
              if( spcname .NE. ptname(itr) ) cycle
              do j = 2,ny-1
                do i = 2,nx-1
                  n3d = i + nx*(j - 1) + nx*ny*(k - 1)
                  n4d = n3d + nx*ny*nz*(itr - 1)
                  ptconc(n4d) = BNDLPT
                  ptconc(n4d) = amax1(ptconc(n4d),cinit(i,j))
                enddo
              enddo
            enddo
          enddo
        enddo
      else
        do ispc = 1,num_ioric
          do k = 1,nz
            read(ioric)
          enddo
        enddo
        goto 100
      endif
c
c-----If this is not a restart, interpolate coarse grid concentrations
c     to all fine grids
c
      if (ngrid.gt.1) then
          do ip = 1,ngrid
            do ic = 1,nchdrn(ip)
              ig = idchdrn(ic,ip)
              call intrpcnc(ntotsp,ncol(ip),nrow(ip),nlay(ip),i1(ig),
     &                      j1(ig),nmesh(ig),ncol(ig),nrow(ig),nlay(ig),
     &                      ptconc(ipsa3d(ip)),ptconc(ipsa3d(ig)) )
            enddo
          enddo
      endif
c
c-----Convert from ppm to umol/m3
c
      do igrd = 1,ngrid
        nx = ncol(igrd)
        ny = nrow(igrd)
        nz = nlay(igrd)
        do itr = 1,ntotsp
          if( .NOT. lsagas(itr) .OR. ptname(itr)(7:10) .NE. 'IC  ' ) cycle
          do k = 1,nz
            do j = 2,ny-1
              do i = 2,nx-1
                n3d = i + nx*(j - 1) + nx*ny*(k - 1)
                n4d = n3d + nx*ny*nz*(itr - 1)
                convfac = densfac*273./tempk(iptr3d(igrd)-1+n3d)*
     &                        press(iptr3d(igrd)-1+n3d)/1013.
                ptconc(ipsa3d(igrd)-1+n4d) = convfac*
     &                     AMAX1(BNDLPT, ptconc(ipsa3d(igrd)-1+n4d) )
              enddo
            enddo
          enddo
        enddo
      enddo
      goto 999

c-----End of IC file reached
c
 900  write(iout,'(//,a)') 'ERROR in RD_SA_ICFILE:'
      write(iout,*)'End of SA IC file'
      write(iout,*)'Make sure initial condition file contains the ',
     &                                   'simulation beginning hour.'
      call camxerr()
c
 999  continue
c
      return
      end
