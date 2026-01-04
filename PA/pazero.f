      subroutine pazero
      use grid
      use chmstry
      use procan
      use tracer
      use node_mod
c
c----CAMx v7.20 220430
c
c      Copyright 1996 - 2022
c     Ramboll
c
c     Initializes or re-initializes the Process Analysis data structures
c     to zero.  This routine is called at the beginning of each time
c     step.
c
c     Subroutines Called:
c
c
c     Called by:
c        CAMX
c
      implicit none
      include 'camx.prm'
      include 'flags.inc'
c
      integer i, j, k, igrd
c
c   ---- initialize the integrated process rates array to zero ---
c
      if( lipr ) then
         do i=1,NPAPRC
           do j=1,npa_cels
              do k=1,nspec
                 cipr(i,j,k) = 0.
              enddo
           enddo
         enddo 
c
c  ---- initialize the the number of steps to zero ---
c
         do i=1,npa_cels
           do j=1,nspec
             npastep(i,j) = 0
           enddo
         enddo 
      endif
c
c  --- zero out data structures for the integrated reaction rate data ---
c
      if( lirr ) then
c
c  ---- initialize the integrated reaction rates array to zero ---
c
        do i=1,npa_cels
            do j=1,nreact
               cirr(i,j) = 0.
            enddo
        enddo
c
c  --- zero out the gridded chemical process analysis array ---
c
        ptconc = 0.
      endif
c
      return
      end
