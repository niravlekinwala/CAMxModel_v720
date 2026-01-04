c**** NCF_WRT_GLOBAL_DDM
c
      subroutine ncf_wrt_global_ddm(action,iounit)
      use ncf_iomod
      use filunit
      use grid
      use tracer
      implicit none
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routine writes the Global attributes that are specific to
c   SA to the NetCDF file
c
c      Copyright 1996 - 2022
c     Ramboll
c      Argument description:
c       Inputs:
c           action    C  name of file to open
c           iounit    I  NetCDF file ID of file
c       Outputs:
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     02/20/17   --gwilson--    Original development
c     08/09/18   --bkoo--       Added ntermddm
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'netcdf.inc'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      character*(*) action
      integer       iounit
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
      integer ierr
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
      ierr = nf_put_att_text(iounit, NF_GLOBAL, 'PROBING_TOOL', 
     &                                              10, "DDM       ")
      if( ierr .NE. NF_NOERR ) goto 7000
c
      if( lbndry ) then
         ierr = nf_put_att_int(iounit, NF_GLOBAL, 'STRATIFY_BOUNDARY', 
     &                                                    NF_INT, 1, 1)
      else
         ierr = nf_put_att_int(iounit, NF_GLOBAL, 'STRATIFY_BOUNDARY', 
     &                                                    NF_INT, 1, 0)
      endif
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'NUMBER_SOURCE_REGIONS',
     &                                               NF_INT, 1, nregin)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'NUMBER_SOURCE_GROUPS',
     &                                               NF_INT, 1, ngroup)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'NUMBER_IC_GROUPS',
     &                                               NF_INT, 1, nicddm)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'NUMBER_BC_GROUPS',
     &                                               NF_INT, 1, nbcddm)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'NUMBER_EMIS_GROUPS',
     &                                               NF_INT, 1, nemddm)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'NUMBER_RATE_GROUPS',
     &                                               NF_INT, 1, nrateddm)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'NUMBER_TERM_GROUPS',
     &                                               NF_INT, 1, ntermddm)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      ierr = nf_put_att_int(iounit, NF_GLOBAL, 'NUMBER_HDDM_GROUPS',
     &                                               NF_INT, 1, nhddm)
      if( ierr .NE. NF_NOERR ) goto 7000
c
      if( lptoverride ) then
         ierr = nf_put_att_int(iounit, NF_GLOBAL, 'POINT_SOURCE_OVERRIDE',
     &                                                        NF_INT, 1, 0)
      else
         ierr = nf_put_att_int(iounit, NF_GLOBAL, 'POINT_SOURCE_OVERRIDE',
     &                                                        NF_INT, 1, 1)
      endif
      if( ierr .NE. NF_NOERR ) goto 7000
c
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,a)') 'ERROR in NCF_WRT_GLOBAL_SA:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(A)') 'Cannot write global atttributes to file.'
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
 
