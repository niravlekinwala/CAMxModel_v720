      subroutine wrtdep(tim2,idat2,igrd,iunit,nox,noy,nsptmp,nspdry,
     &                  vdep,depfld)
      use grid
      use chmstry
      use camxcom
c
c----CAMx v7.20 220430
c 
c     WRTDEP writes deposition fields.
c 
c      Copyright 1996 - 2022
c     Ramboll
c           
c     Modifications: 
c        05/17/16       Added in-line Ix emissions to dep output
c 
c     Input arguments:
c        tim2                output time (HHMM)
c        idat2               output date (YYJJJ)
c        igrd                grid number
c        iunit               output unit
c        nox                 number of cells in x-direction
c        noy                 number of cells in y-direction
c        nsptmp              number of dep field species
c        nspdry              number of dry dep species
c        vdep                dry deposition velocities
c        depfld              deposition field to output
c             
c     Output arguments: 
c        none
c             
c     Routines Called: 
c        none
c             
c     Called by: 
c        CAMx
c 
      include 'camx.prm'
      include 'flags.inc'
c
      real depfld(nox,noy,nsptmp),vdep(nox,noy,nspdry)
c
      character*4 ispec(10,4*MXSPEC+2)
c
      data nseg /1/
c
c-----Entry point
c
      ibegcel = 1
      iendcel = nox
      jbegcel = 1
      jendcel = noy
      if( igrd .GT. 1 ) then
        ibegcel = 2
        iendcel = nox-1
        jbegcel = 2
        jendcel = noy-1
      endif
c
c-----Determine time/date range
c
      idat1 = idat2
      etim = AINT(ANINT(tim2)/100.) + amod(ANINT(tim2),100.)/60.
      btim = ANINT( 1000*(etim - ANINT(dtout)/60.) )/1000.
      if (btim.lt.0.) then
        btim = btim + 24.
        idat1 = idat1 - 1
        if (MOD(idat1,1000) .EQ. 0 ) then
          if ( MOD(INT(idat1/1000)-1,4) .EQ. 0 ) then
            idat1 = (INT(idat1/1000)-1)*1000 + 366
          else
            idat1 = (INT(idat1/1000)-1)*1000 + 365
          end if
        endif
      endif 
c
c-----Write gridded deposition field
c
      do l = 1,4*ndepspc
        read(depsp(l),'(10a1)') (ispec(n,l),n=1,10)
      enddo
      write(iunit) idat1,btim,idat2,etim
      do l = 1,ndepspc
        ll = ldepmap(l)
        write(iunit) nseg,(ispec(n,l),n=1,10),
     &               ((vdep(i,j,ll),i=ibegcel,iendcel),j=jbegcel,jendcel)
      enddo
      do l = 1,3*ndepspc
        ll = l + ndepspc
        write(iunit) nseg,(ispec(n,ll),n=1,10),
     &               ((depfld(i,j,l),i=ibegcel,iendcel),j=jbegcel,jendcel)
      enddo
      if (lixemis) then
        l = 4*ndepspc + 1
        ll = 3*ndepspc + 1
        read(depsp(l),'(10a1)') (ispec(n,l),n=1,10)
        write(iunit) nseg,(ispec(n,l),n=1,10),
     &               ((depfld(i,j,ll),i=ibegcel,iendcel),j=jbegcel,jendcel)
        l = 4*ndepspc + 2
        ll = 3*ndepspc + 2
        read(depsp(l),'(10a1)') (ispec(n,l),n=1,10)
        write(iunit) nseg,(ispec(n,l),n=1,10),
     &               ((depfld(i,j,ll),i=ibegcel,iendcel),j=jbegcel,jendcel)
      endif
c
      return
      end
