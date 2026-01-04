c**** NCF_CHKFILE
c
      function ncf_chkfile(iounit,fname,action,cname)
      use filunit
      implicit none
      integer ncf_chkfile
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routine tests the file to see if it is an old-style CAMx file
c   by reading the first 10 bytes and comparing against a string.
c
c      Copyright 1996 - 2022
c     Ramboll
c      Argument description:
c       Inputs:
c         iounit I NCF file ID
c         fname  C filename
c         action C string that describes file being read
c         cname  C type of file to compare against
c       Outputs:
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     02/20/17   --gwilson--    Original development
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      integer        iounit
      character*(*)  fname
      character*(*)  action
      character*(10) cname
c
c-----------------------------------------------------------------------
c    External functions:
c-----------------------------------------------------------------------
c
      integer istrln
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      character*10 this_name
      character*4  in_name(10)
      integer      istatus, rec_count
c
c-----------------------------------------------------------------------
c    Parameters:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c  --- assume success ---
c
      ncf_chkfile = ISUCES
c
c  --- open the file as type binary ---
c
      open(unit=iounit,file=fname,form='UNFORMATTED',access='STREAM',
     &                                              status='OLD',ERR=7000)
c 
c  --- read the first 10 bytes ---
c
      read(UNIT=iounit,IOSTAT=istatus) rec_count,in_name
      write(this_name,'(10A1)') in_name
c
c  --- check against input string ---
c
      if( cname(:istrln(cname)) .NE. this_name(:istrln(cname)) ) 
     &                                              ncf_chkfile = IFAIL 
c
c  --- close file ---
c
      close(iounit)
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,a)') 'ERROR in NCF_CHKFILE:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Could not open file: ',fname(:istrln(fname))
      call camxerr()
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
