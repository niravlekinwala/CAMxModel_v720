C**** FNDKEY.F 
c
c----CAMx v7.20 220430
c
c      Copyright 1996 - 2022
c     Ramboll
c
      subroutine fndkey( ierr,
     &                       iounit, keyin )
      implicit none
c
c-----------------------------------------------------------------------
c
c     This routine searches the file attached to logical unit number
c     "iounit" for a string matching "keyin".  The first 20 cahracters
c     of each line are read and converted to upper case.  It is assumed
c     that the "keyin" string is upper case.  An initial search is made,
c     and if not successful, the file is rewound and a second search is
c     made.  If a match is still not found, the error flag is set to
c     failure.
c
c   Arguments:
c
c     Outputs:
c       ierr     I   error flag (either ISUCES or IFAIL)
c     Inputs:
c       iounit   I   logical unit number of file to read
c       keyin    C   string to match (all upper case)
c
c-----------------------------------------------------------------------
c   LOG:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c   Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
c
c-----------------------------------------------------------------------
c   External functions:
c-----------------------------------------------------------------------
c
c       istrln   I   returns actual length of string (no trailing blanks)
c
      integer istrln
c
c-----------------------------------------------------------------------
c   Argument decleration:
c-----------------------------------------------------------------------
c
      integer            ierr
      integer            iounit
      character*(KEYLEN) keyin
c
c-----------------------------------------------------------------------
c   Local variables:
c-----------------------------------------------------------------------
c
c       keywrd   C   keyword to check against argument for match
c       ilen     I   actual length of string in argument list
c       lpass2   L   flag to determine if currently on second pass
c
      character*(KEYLEN) keywrd
      integer            ilen
      logical            lpass2
c
c-----------------------------------------------------------------------
c   Entry point:
c-----------------------------------------------------------------------
c
      ierr = IFAIL
      lpass2 = .FALSE.
c
c   ---- get the actual length of Argument ----
c
      ilen = istrln ( keyin )
c
c   ---- read next string from file ----
c
 111  continue
      read( iounit, 8000, ERR=7000, END=7001 ) keywrd
c
c   ---- convert ot upper case ---
c
      call toupper( keywrd )
c
c   ---- get rid of leading blanks ----
c
      call jstlft( keywrd )
c
c   ---- if match return success; otherwise read next record ---
c
      if( keywrd(1:ilen) .EQ. keyin(1:ilen) ) then
          ierr = ISUCES
          goto 9999
      else
          goto 111
      endif
c
c----------------------------------------------------------------------
c   Error statements:
c----------------------------------------------------------------------
c
 7000 continue
      ierr = IRDERR
      goto 9999
c
 7001 continue
      rewind( iounit, ERR=7000 )
      if( lpass2 ) then
         ierr = IEOF
         goto 9999
      else
         lpass2 = .TRUE.
         goto 111
      endif
c
c-----------------------------------------------------------------------
c   Format statements:
c-----------------------------------------------------------------------
c
 8000 format(A20)
c
c-----------------------------------------------------------------------
c   Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
