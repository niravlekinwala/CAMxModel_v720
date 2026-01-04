      subroutine emassign(ncol,nrow,nlay_ems,io,jo,nmesh,ncolf,nrowf,
     &                                            nspcems,cgems,fgems)
c
c----CAMx v7.20 220430
c
c     EMASSIGN assigns emissions to fine grid values from coarse grid
c
c      Copyright 1996 - 2022
c     Ramboll
c          
c     Modifications:
c        none
c
c     Input arguments:
c        ncol              number of columns in the parent grid
c        nrow              number of rows in the parent grid
c        nlay_ems          number of layers in emissions filas
c        io                starting i index for the fine grid
c        jo                starting j index for the fine grid
c        nmesh             mesh number
c        ncolf             number of columns in fine grid
c        nrowf             number of rows in fine grid
c        nspcems           number of emissions species
c        cgems             emissions value on coarse grid
c
c     Output arguments:
c        fgems             emissions value on fine grid
c
c     Subroutine called:
c        none
c
c     Called by:
c
      real cgems(ncol,nrow,nlay_ems,nspcems)
      real fgems(ncolf,nrowf,nlay_ems,nspcems)
c
c-----Entry point
c
        do k=1,nspcems
           do 40 jfin = 2,nrowf-1
             j = (jfin - 2)/nmesh + jo
             do 30 ifin = 2,ncolf-1
               i = (ifin - 2)/nmesh + io
               do ilay=1,nlay_ems
                   fgems(ifin,jfin,ilay,k) = cgems(i,j,ilay,k)/(nmesh*nmesh)
               enddo
  30         continue
  40       continue
        enddo
c
      return
      end
