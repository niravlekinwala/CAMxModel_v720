      subroutine o3colprep
      use filunit
      use grid
      use o3colmap
      implicit none
c
c----CAMx v7.20 220430
c
c     O3COLPREP reads the header records of the ozone column file,
c     and any optional constant codes
c
c      Copyright 1996 - 2022
c     Ramboll
c
c     Modifications:
c        06/13/03   Added optional snow cover, land-ocean, drought, 
c                   and roughness maps; only snow cover is time-varying
c        02/11/11   Moved snow panel to 2D met file, removed roughness,
c                   removed albedo panel (now determined in SRFPREP)
c        03/30/12   Removed haze and drought stress codes
c        07/07/16   Removed ocean mask option
c
c     Input arguments:
c        none
c
c     Output arguments:
c        none
c
c     Routines called:
c        none
c
c     Called by:
c        STARTUP
c
      include 'camx.prm'
      include 'flags.inc'
c
      character*10  title, name
      character*80  action
      character*180 line
      real          rdozn(NOZN)
      integer       i, j, nhdro3col
c
      data name    /'OZONE COL '/
c
c-----Entry point
c
      read(io3col,*)
      nhdro3col = 1
c
c-----Read ozone class and check inputs
c
      action = 'looking for ozone column header'
      read(io3col,'(a10,5f10.0)',err=900) title,(rdozn(i),i=1,NOZN)
      nhdro3col = nhdro3col + 1
      if (title.ne.name) then
        write(iout,'(//,a)') 'ERROR in O3COLPREP:'
        write(iout,*) 'After reading line: ',nhdro3col
        write(iout,*) 'Expecting keyword: ',name
        write(iout,*) 'Read from ozone column file: ',title
        write(iout,'(/,2A)') 'This version uses a different format for ',
     &                             'what was the Albedo/Haze/Ozone file.'
        write(iout,'(/,2A,/,A)') 'For more information refer to the User Guide ',
     &                         'and Release Notes included with ',
     &                          'this CAMx distribution.'
        call camxerr()
      endif
      if (lchem) then
        do i = 1,NOZN
          if (rdozn(i).ne.ozcl(i)) then
            write(iout,'(//,a)') 'ERROR in O3COLPREP:'
            write(iout,*) 'After reading line: ',nhdro3col
            write(iout,*) 'Mismatch in ozone class'
            write(iout,*) 'Photolysis rates  file: ',(ozcl(j),j=1,NOZN)
            write(iout,*) 'Ozone column file: ',(rdozn(j),j=1,NOZN)
            write(iout,'(/,2A)') 'This version uses a different format for ',
     &                             'what was the Albedo/Haze/Ozone file.'
            write(iout,'(/,2A,/,A)') 'For more information refer to the User Guide ',
     &                         'and Release Notes included with ',
     &                          'this CAMx distribution.'
            call camxerr()
          endif
        enddo
      endif
c
c-----Read any optional records -- these should not exist
c
      action = 'looking for headers for any optional constant inputs'
      read(io3col,'(a)',err=900) line
      nhdro3col = nhdro3col + 1
      read(line,'(a10)') title
      if (title.eq.name) then
        backspace(io3col)
        nhdro3col = nhdro3col - 1
        return
      endif

      write(iout,'(//,a)') 'ERROR in O3COLPREP:'
      write(iout,*) 'After reading line: ',nhdro3col
      write(iout,*) 'Unrecognized keyword in O3MAP file!'
      write(iout,*) 'Acceptable keyword is: OZONE COL'
      write(iout,*) 'O3MAP file should not contain any optional fields'
      write(iout,*) 'Review your ozone column file.'
      write(iout,'(/,2A)') 'This version uses a different format for ',
     &                             'what was the Albedo/Haze/Ozone file.'
      write(iout,'(/,2A,/,A)') 'For more information refer to the User Guide ',
     &                         'and Release Notes included with ',
     &                          'this CAMx distribution.'
      call camxerr()
c
c --- Trap read errors
c
 900  write(iout,'(//,a,i4)') 'ERROR in O3COLPREP reading line ',
     &                                                     nhdro3col+1
      write(iout,'(2a)') 'While ', action
      write(iout,'(/,2A)') 'This version uses a different format for ',
     &                             'what was the Albedo/Haze/Ozone file.'
      write(iout,'(/,2A,/,A)') 'For more information refer to the User Guide ',
     &                         'and Release Notes included with ',
     &                          'this CAMx distribution.'
      call camxerr()
c
      end
