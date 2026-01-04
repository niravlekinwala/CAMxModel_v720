      subroutine ncf_rdstacks( )
      use chmstry
      use filunit
      use ptemiss
      use tracer
      implicit none
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c     NCF_RDPTHDR reads the stack data from the time invariant data in 
c     the NetCDF point source file and loads the data into global varables.
c                           
c
c      Copyright 1996 - 2022
c     Ramboll
c
c      Argument description:
c       Inputs:
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
      include 'flags.inc'
      include 'netcdf.inc'
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
      character*200 action
      character*10  this_var
      integer       this_varid, idxfile, ierr, idx_start, ipt, n, i
      real          pi
c
      real*4,  allocatable, dimension(:) :: darray
      integer, allocatable, dimension(:) :: iarray
c
c-----------------------------------------------------------------------
c    Data statments:
c-----------------------------------------------------------------------
c
      data pi /3.1415927/
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c   --- loop over number of files ---
c
      do idxfile=1,npoint_files
         if( .NOT. is_netcdf_iptem(idxfile) ) cycle
         action = 'Reading stack data from NetCDF point source file.'
c
c  --- set the index into stacks array ---
c
         idx_start = idx_start_pts(idxfile)
c
c --- allocate the local array to read the stack data ---
c
         allocate( darray(nptsrc_files(idxfile)) ) 
         allocate( iarray(nptsrc_files(idxfile)) ) 
c
c  --- X coordinates ---
c
         this_var = "xcoord"
         ierr = nf_inq_varid(iptem(idxfile),this_var,this_varid)
         if( ierr .NE. NF_NOERR ) goto 7000
         ierr = nf_get_var_real(iptem(idxfile),this_varid,darray)
         if( ierr .NE. NF_NOERR ) goto 7001
         do ipt=1,nptsrc_files(idxfile)
             xloc(ipt+idx_start) = darray(ipt)
         enddo
c
c  --- Y coordinates ---
c
         this_var = "ycoord"
         ierr = nf_inq_varid(iptem(idxfile),this_var,this_varid)
         if( ierr .NE. NF_NOERR ) goto 7000
         ierr = nf_get_var_real(iptem(idxfile),this_varid,darray)
         if( ierr .NE. NF_NOERR ) goto 7001
         do ipt=1,nptsrc_files(idxfile)
             yloc(ipt+idx_start) = darray(ipt)
         enddo
c
c  --- stack height ---
c
         this_var = "stkheight"
         ierr = nf_inq_varid(iptem(idxfile),this_var,this_varid)
         if( ierr .NE. NF_NOERR ) goto 7000
         ierr = nf_get_var_real(iptem(idxfile),this_varid,darray)
         if( ierr .NE. NF_NOERR ) goto 7001
         do ipt=1,nptsrc_files(idxfile)
             hstk(ipt+idx_start) = darray(ipt)
         enddo
c
c  --- stack diameter ---
c
         this_var = "stkdiam"
         ierr = nf_inq_varid(iptem(idxfile),this_var,this_varid)
         if( ierr .NE. NF_NOERR ) goto 7000
         ierr = nf_get_var_real(iptem(idxfile),this_varid,darray)
         if( ierr .NE. NF_NOERR ) goto 7001
         do ipt=1,nptsrc_files(idxfile)
             dstk(ipt+idx_start) = ABS(darray(ipt))
         enddo
c
c  --- stack gas exit temperature ---
c
         this_var = "stktemp"
         ierr = nf_inq_varid(iptem(idxfile),this_var,this_varid)
         if( ierr .NE. NF_NOERR ) goto 7000
         ierr = nf_get_var_real(iptem(idxfile),this_varid,darray)
         if( ierr .NE. NF_NOERR ) goto 7001
         do ipt=1,nptsrc_files(idxfile)
             tstk(ipt+idx_start) = darray(ipt)
         enddo
c
c  --- stack gas exit velocity ---
c
         this_var = "stkspeed"
         ierr = nf_inq_varid(iptem(idxfile),this_var,this_varid)
         if( ierr .NE. NF_NOERR ) goto 7000
         ierr = nf_get_var_real(iptem(idxfile),this_varid,darray)
         if( ierr .NE. NF_NOERR ) goto 7001
         do ipt=1,nptsrc_files(idxfile)
             vstk(ipt+idx_start) = darray(ipt)/3600.
         enddo
c
c  --- PiG flag ---
c
         if( ipigflg .NE. 0 ) then
            this_var = "pigflag"
            ierr = nf_inq_varid(iptem(idxfile),this_var,this_varid)
            if( ierr .NE. NF_NOERR ) goto 7000
            ierr = nf_get_var_int(iptem(idxfile),this_varid,iarray)
            if( ierr .NE. NF_NOERR ) goto 7001
            do ipt=1,nptsrc_files(idxfile)
               lpiglet(idx_start+ipt) = .FALSE.
               if( iarray(ipt) .GT. 0 ) lpiglet(idx_start+ipt) = .TRUE.
            enddo         
         endif
c
c --- deallocate local arrays ---
c
         deallocate( darray )
         deallocate( iarray )
c
      enddo
c
c  --- load  coordinates into array for probing tools ---
c
      if( ltrace  .OR. lddm .OR. lhddm ) then
         call alloc_tracer_pts(MAX(nptsrc,1))
         if( lptsrc ) then
            do ipt = 1,nptsrc
               xlocpt(ipt) = xloc(ipt)
               ylocpt(ipt) = yloc(ipt)
            enddo
         endif
      endif
c
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,a)') 'ERROR in NCF_RDSTACKS:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot find variable id for: ',
     &                                      this_var(:istrln(this_var))
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7001 continue
      write(iout,'(//,a)') 'ERROR in NCF_RDSTACKS:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot read data for variable: ',
     &                                      this_var(:istrln(this_var))
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
