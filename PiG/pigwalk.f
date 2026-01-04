      subroutine pigwalk(dt) 
      use grid
      use filunit
      use camxfld
      use pigsty
c
c----CAMx v7.20 220430
c
c     PIGWALK transports PiG puffs for the duration of one input timestep
c
c      Copyright 1996 - 2022
c     Ramboll
c
c     Modifications:
c        9/5/03              Revised to transport puff ends in "chained"
c                            puff approach
c        9/30/14             Added vertical advection of puff centerpoint
c
c     Input arguments:
c        dt                  time step (s)     
c
c     Output arguments:
c        none
c
c     Subroutines called:
c        PIGCOORD
c        WALK1PUF
c
c     Called by:
c        EMISTRNS
c
      implicit none
      include "camx.prm"
c
      real dt,dtcell,xpig,ypig,dzbot,dztop,puffb,pufft,ztop
      real ww(0:MXLAYER)
      integer n,kount,igrd,iipig,jjpig,ingrdf,ingrdb
      integer nflip
      logical lkick
c
c-----Entry point
c
c-----Loop over all PiGs; find active puffs to move
c
      do 20 n = 1,npig
        if (ingrd(n).eq.0) goto 20
c
        xpig = (xpigf(n) + xpigb(n))/2.
        ypig = (ypigf(n) + ypigb(n))/2.
        call pigcoord(xpig,ypig,iipig,jjpig,igrd)
c
c-----Walk the leading (front) edge
c     New born PiGs walk for the duration of their age;
c     Older puffs walk for the duration of the current timestep
c
        call pigcoord(xpigf(n),ypigf(n),iipig,jjpig,ingrdf)
c
        if (lnewt(n)) then
          dtcell = agepigf(n)
        else
          dtcell = dt
          agepigf(n) = agepigf(n) + dt
        endif
c
        kount = 0
        nflip = 0
        lkick = .false.
  10    kount = kount + 1
        igrd = ingrdf
        call walk1puf(dtcell,igrd,ncol(igrd),nrow(igrd),nlay(igrd),
     &                iipig,jjpig,ingrdf,lkick,xpigf(n),ypigf(n),
     &                pufftop(n),puffbot(n),windu(iptr3d(igrd)),
     &                windv(iptr3d(igrd)),height(iptr3d(igrd)),
     &                press(iptr3d(igrd)),tempk(iptr3d(igrd)))
        if (igrd.ne.ingrdf) nflip = nflip + 1
        if (kount.gt.50 .and. dtcell.gt.0.) then
          if (nflip.ge.10) then
            nflip = 0
            lkick = .true.
            goto 10
          endif
          write(iout,'(//,a)') 'ERROR in PIGWALK:'
          write(iout,*) 'Number of steps > 50'
          write(iout,*) 'Front edge: Puff#,grid#,timestep left:'
          write(iout,*) n,ingrdf,dtcell
          write(iout,*) 'Location (i,j), puff top & bottom'
          write(iout,*) iipig,jjpig,pufftop(n),puffbot(n)
          call camxerr()
        endif
        if (dtcell.gt.0.) goto 10
        if (ingrdf.eq.0) then
          ingrd(n) = 0
          goto 20
        endif
c
c-----Walk the trailing (back) edge
c
        call pigcoord(xpigb(n),ypigb(n),iipig,jjpig,ingrdb)
c
        if (lnewt(n)) then
          dtcell = agepigb(n)
        else
          dtcell = dt
          agepigb(n) = agepigb(n) + dt
        endif
        if (dtcell.eq.0.) goto 30
c  
        kount = 0
        nflip = 0
        lkick = .false.
  11    kount = kount + 1
        igrd = ingrdb
        call walk1puf(dtcell,igrd,ncol(igrd),nrow(igrd),nlay(igrd),
     &                iipig,jjpig,ingrdb,lkick,xpigb(n),ypigb(n),
     &                pufftop(n),puffbot(n),windu(iptr3d(igrd)),
     &                windv(iptr3d(igrd)),height(iptr3d(igrd)),
     &                press(iptr3d(igrd)),tempk(iptr3d(igrd)))
        if (igrd.ne.ingrdb) nflip = nflip + 1
        if (kount.gt.50 .and. dtcell.gt.0.) then
          if (nflip.ge.10) then
            nflip = 0
            lkick = .true.
            goto 11
          endif
          write(iout,'(//,a)') 'ERROR in PIGWALK:'
          write(iout,*) 'Number of steps > 50'
          write(iout,*) 'Back edge: Puff#,grid#,timestep left:'
          write(iout,*) n,ingrdb,dtcell
          write(iout,*) 'Location (i,j), puff top & bottom'
          write(iout,*) iipig,jjpig,pufftop(n),puffbot(n)
          write(iout,*) 
          call camxerr()
        endif
        if (dtcell.gt.0.) goto 11
        if (ingrdb.eq.0) then
          ingrd(n) = 0
          goto 20
        endif
c
c-----Determine new puff center point coords, host grid, and puff length
c
 30     xpig = (xpigf(n) + xpigb(n))/2.
        ypig = (ypigf(n) + ypigb(n))/2.
        call pigcoord(xpig,ypig,iipig,jjpig,ingrd(n))
        igrd = ingrd(n)
        if (lnewt(n)) then
          lnewt(n) = .false.
          goto 20
        endif
c
c-----Calculate vertical velocity profile and move puff for full timestep
c
        call getdwdz(ncol(igrd),nrow(igrd),nlay(igrd),iipig,jjpig,
     &               deltax(jjpig,igrd),deltay(igrd),
     &               height(iptr3d(igrd)),press(iptr3d(igrd)),
     &               tempk(iptr3d(igrd)),windu(iptr3d(igrd)),
     &               windv(iptr3d(igrd)),ztop,ww)
c
        dzbot = zpig(n) - puffbot(n)
        dztop = pufftop(n) - zpig(n)
c
        call walk1pufz(dt,ncol(igrd),nrow(igrd),nlay(igrd),iipig,jjpig,
     &                 zpig(n),height(iptr3d(igrd)),ww)
c
        puffb = amax1(0.,zpig(n) - dzbot)
        pufft = amin1(ztop,zpig(n) + dztop)
        if (puffb.eq.0.) pufft = puffb + axisz(n)
        if (pufft.eq.ztop) puffb = pufft - axisz(n)
        puffbot(n) = puffb
        pufftop(n) = pufft
        zpig(n) = (pufftop(n) + puffbot(n))/2.
        axisz(n) = pufftop(n) - puffbot(n)
c
  20  continue
c
      return
      end
c
c-------------------------------------------------------------------------------
c
      subroutine getdwdz(ncol,nrow,nlay,i,j,dx,dy,
     &                   height,press,tempk,windu,windv,ztop,ww)
c
c-----GETDWDZ calculates the vertical velocity profile
c     
      implicit none
      include "camx.prm"
      integer ncol,nrow,nlay,i,j
      real dx,dy,ztop,height(ncol,nrow,nlay),
     &     press(ncol,nrow,nlay),tempk(ncol,nrow,nlay),
     &     windu(ncol,nrow,nlay),windv(ncol,nrow,nlay)
      integer kk
      real deplyr,uup,uum,vvp,vvm,rho,rhoxp,rhoxm,rhoyp,rhoym
      real ww(0:MXLAYER)
c 
c-----Entry point
c
      ztop = height(i,j,nlay)
      ww(0) = 0.
      do kk = 1,nlay
        deplyr = height(i,j,kk)
        if (kk.gt.1) deplyr = height(i,j,kk) - height(i,j,kk-1)
        rho  = press(i,j,kk)/tempk(i,j,kk)
        rhoxp = press(i+1,j,kk)/tempk(i+1,j,kk)
        rhoxm = press(i-1,j,kk)/tempk(i-1,j,kk)
        rhoyp = press(i,j+1,kk)/tempk(i,j+1,kk)
        rhoym = press(i,j-1,kk)/tempk(i,j-1,kk)
        uup = windu(i,j,kk)*(rhoxp + rho)/2.
        uum = windu(i-1,j,kk)*(rhoxm + rho)/2.
        vvp = windv(i,j,kk)*(rhoyp + rho)/2.
        vvm = windv(i,j-1,kk)*(rhoym + rho)/2.
        ww(kk) = ww(kk-1) - (deplyr/rho)*
     &                      ((uup - uum)/dx + (vvp - vvm)/dy)
      enddo
c
      return
      end
