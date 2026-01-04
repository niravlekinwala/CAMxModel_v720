      subroutine banner( )
c
c----CAMx v7.20 220430
c
c     BANNER writes to the screen the license blurb for CAMx.
c
c      Copyright 1996 - 2022
c     Ramboll
c
c     Modifications:
c        none
c
c     Input arguments:
c        none
c
c     Output arguments:
c        iounit
c
c     Routines called:
c        none
c
c     Called by:
c         MAIN
c
c-----Entry point
c
      write(*,'(/,1X,2A)') 'By running CAMx you accept the terms ',
     &                                           'of the CAMx LICENSE.'
      write(*,'(1X,2A)') 'The CAMx LICENSE is included with the ',
     &                     'CAMx source code and also available here:'
      write(*,'(2X,A,//)') 'http://www.camx.com/Files/LICENSE_v720.txt'
      call system("sleep 5")
 
      return
      end
