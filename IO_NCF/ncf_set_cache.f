c**** NCF_SET_CACHE
c
      subroutine ncf_set_cache()
      use grid
      use filunit
      implicit none
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routine sets the cache parameters that will control how
c   the NetCDF output files will be chunked.
c
c      Copyright 1996 - 2022
c     Ramboll
c      Argument description:
c       Inputs:
c       Outputs:
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     04/20/17   --gwilson--    Original development
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'namelist.inc'
      include 'ncf_iodat.inc'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
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
c   primes I array of prime numbers 
c
      integer primes(168)
      logical lfound_prime
      integer ncolmax, nrowmax, i
c
c-----------------------------------------------------------------------
c    Data Statements:
c-----------------------------------------------------------------------
c
      data primes /2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,
     &             71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,
     &             149,151,157,163,167,173,179,181,191,193,197,199,211,
     &             223,227,229,233,239,241,251,257,263,269,271,277,281,
     &             283,293,307,311,313,317,331,337,347,349,353,359,367,
     &             373,379,383,389,397,401,409,419,421,431,433,439,443,
     &             449,457,461,463,467,479,487,491,499,503,509,521,523,
     &             541,547,557,563,569,571,577,587,593,599,601,607,613,
     &             617,619,631,641,643,647,653,659,661,673,677,683,691,
     &             701,709,719,727,733,739,743,751,757,761,769,773,787,
     &             797,809,811,821,823,827,829,839,853,857,859,863,877,
     &             881,883,887,907,911,919,929,937,941,947,953,967,971,
     &             977,983,991,997/
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
      if ( NCF_CACHESIZE .LE. 0 ) then
        write(iout,'(//,A)') 'ERROR in NCF_IODAT.INC:'
        write(iout,'(A)') 'Invalid value provided for parameter NCF_CACHESIZE.'
        write(iout,'(2A)') 'Value represents total size of the raw data ',
     &                                           'chunk cache in MegaBytes.'
        write(iout,'(A)') 'It must be a postive number.'
        call camxerr()
      endif
c
      lfound_prime = .FALSE.
      do i=1,166
         if( NCF_NELEMS .EQ. primes(i) ) lfound_prime = .TRUE.
      enddo
      if( NCF_NELEMS .LE. 0 .OR. .NOT. lfound_prime ) then
        write(iout,'(//,A)') 'ERROR in NCF_IODAT.INC:'
        write(iout,'(A)') 'Invalid value provided for parameter NCF_NELEMS.'
        write(iout,'(2A)') 'This value must be a prime number.'
        call camxerr()
      endif
c
      if( NCF_PREEMPTION .LT. 0 .OR. NCF_PREEMPTION .GT. 100 ) then
        write(iout,'(//,A)') 'ERROR in NCF_IODAT.INC:'
        write(iout,'(A)') 'Invalid value provided for parameter NCF_PREEMPTION.'
        write(iout,'(A)') 'This value must be between 0 and 100.'
        call camxerr()
      endif
c
      if( NCF_DEFLATE_LEVEL .LT. 0 .OR.  NCF_DEFLATE_LEVEL .GT. 9 ) then
        write(iout,'(//,A)') 'ERROR in NCF_IODAT.INC:'
        write(iout,'(A)') 'Invalid value provided for parameter NCF_DEFLATE_LEVEL.'
        write(iout,'(A)') 'This value must be between 0 and 9.'
        call camxerr()
      endif
c
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
 
