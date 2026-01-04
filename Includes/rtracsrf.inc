c**** RTRACSRF
c
c----CAMx v7.20 220430
c
c
c----------------------------------------------------------------------
c
c    Include file for variables used in the RTRAC surface model.
c
c      Copyright 1996 - 2022
c     Ramboll
c
c-----------------------------------------------------------------------
c   LOG:
c-----------------------------------------------------------------------
c
c     10/29/09   -cemery-   original development
c     05/15/13   -cemery-   moved soil and shade fractions to deposit.inc
c
c-----------------------------------------------------------------------
c   Variables for the RTRAC surface and partitioning models
c-----------------------------------------------------------------------
c
c   lsrfmodrt L   flag to invoke the RTRAC surface model 
c   lparttn   L   flag to invoke RTRAC gas/PM partitioning
c
      logical      lsrfmodrt
      logical      lparttn
c
      common /rtsrfpart/ lsrfmodrt, lparttn
c
