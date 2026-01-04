      subroutine wrtsrf(lsrfmod,ngrd,tim2,idat2,iunit,nox,noy,nspdep,spname,
     &                  solmas,vegmas)
c
c----CAMx v7.20 220430
c 
c     WRTSRF writes surface model mass fields
c 
c      Copyright 1996 - 2022
c     Ramboll
c           
c     Modifications: 
c     08/26/21 -cemery- RTRAC surface file doubles as (default) dep output
c                       or (flagged) surface model output
c     09/01/21 -cemery- Updated to not write nested grid buffer cells
c 
c     Input arguments:
c        lsrfmod             surface model flag
c        ngrd                grid index
c        tim2                output time (HHMM)
c        idat2               output date (YYJJJ)
c        iunit               output unit
c        nox                 number of cells in x-direction
c        noy                 number of cells in y-direction
c        nspdep              number of dep field species
c        spname              species name array
c        solmas              surface soil mass field to output
c        vegmas              surface veg mass field to output
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
      implicit none
c
      include 'camx.inc'
c
      logical lsrfmod
      integer ngrd
      integer idat2, iunit, nox, noy, nspdep
      real tim2
      real solmas(nox,noy,nspdep)
      real vegmas(nox,noy,nspdep)
      character*10 spname(nspdep)
c
      integer nseg, idat1, i, j, n, l
      integer ibegcel,iendcel,jbegcel,jendcel
      real tim1,btim,etim
      character*4 ispec(10)
c
      data nseg /1/
c
c-----Entry point
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
c-----Write the time stamp ---
c
      write(iunit) idat1,btim,idat2,etim
c
c-----Write gridded surface mass field
c
      ibegcel = 1
      iendcel = nox
      jbegcel = 1
      jendcel = noy
      if( ngrd .GT. 1 ) then
        ibegcel = 2
        iendcel = nox-1
        jbegcel = 2
        jendcel = noy-1
      endif

      do l = 1,nspdep
        if (lsrfmod) then
          ispec(1) = 'S'
          ispec(2) = '_'
        else
          ispec(1) = 'D'
          ispec(2) = '_'
        endif
        read(spname(l)(1:8),'(8a1)') (ispec(n),n=3,10)
        write(iunit) nseg,(ispec(n),n=1,10),
     &           ((solmas(i,j,l),i=ibegcel,iendcel),j=jbegcel,jendcel)
      enddo
      do l = 1,nspdep
        if (lsrfmod) then
          ispec(1) = 'V'
          ispec(2) = '_'
        else
          ispec(1) = 'W'
          ispec(2) = '_'
        endif
        read(spname(l)(1:8),'(8a1)') (ispec(n),n=3,10)
        write(iunit) nseg,(ispec(n),n=1,10),
     &           ((vegmas(i,j,l),i=ibegcel,iendcel),j=jbegcel,jendcel)
      enddo
c
      return
      end
