c**** NCF_RDBCRT.F
c
      subroutine ncf_rdbcrt(nox,noy,noz,nspsa,saconc,tpgrd,prgrd)
      use tracer
      use grid
      use chmstry
      use rtracchm
      use filunit
      implicit none
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routine fils one hour of boundary conditions for the RTRAC
c   species. It then places these concentrations in the  appropriate 
c   place in the gridded array used for tracer concentrations.  
c   The values are left in PPM so that they can be interpolated to
c   the nests in a later step. This version is for NetCDF files.
c
c      Copyright 1996 - 2022
c     Ramboll
c
c      Argument description:
c       Outputs:
c           saconc   R  tracer concentrations
c       Inputs:
c           noy      I  number of Y cells in the grid
c           noz      I  number of layers in the grid
c           nspsa    I  number of tracer species
c           tpgrd    I  3-D temperature field
c           prgrd    I  3-D pressure field
c       
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'camx.inc'
      include 'netcdf.inc'         
c
c-----------------------------------------------------------------------
c    Arguement declarations:
c-----------------------------------------------------------------------
c
      integer   nox
      integer   noy
      integer   noz
      integer   nspsa
      real      saconc(nox,noy,noz,nspsa)
      real      tpgrd(nox,noy,noz)
      real      prgrd(nox,noy,noz)
c
c-----------------------------------------------------------------------
c   External functions:
c-----------------------------------------------------------------------
c
      integer ncf_chkfile
      integer ncf_get_tstep
      integer istrln
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      character*200 action
      character*10  bcname, cfile, this_var, name_in
      integer       ierr, jerr, this_tstep, i, j, k, ispc
      integer       data_start(4), data_count(4), this_varid
      integer       ibcdate_time_tflag, ibcdate_time_etflag
      integer       nedge, ioff, itmp
      logical       lexist, luse
      real          convfac 
c
      real, allocatable, dimension(:,:,:) :: bcgrid
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c  --- if no BC file provided then we skip to the end ---
c
      if( .NOT. lbcfil ) goto 9999
c
c  --- set the number of BC edges --
c
      if( lbndry ) then
        nedge = 5
      else
        nedge = 1
      endif
c
c  --- open the BC ---
c
      inquire(file=bcfil,exist=lexist)
      if( .NOT. lexist ) goto 7000
      action = 'Opening RTRAC Boundary Conditions file.'
      bcname = 'BOUNDARY  '
      ierr = ncf_chkfile(iorbc,bcfil,action,bcname)
      if( ierr .NE. ISUCES ) then
         jerr = nf_open(bcfil,NF_NOWRITE,iorbc)
         if( jerr .NE. NF_NOERR ) goto 7001
      else
         return
      endif
c
c  --- initialize the boundary conditions to lower bound ---
c  
      do ispc=1,ntotsp
        do k=1,noz
          do j=1,noy
            saconc(1,j,k,ispc) = rtlbnd(ispc)
            saconc(nox,j,k,ispc) = rtlbnd(ispc)
          enddo
          do i=1,nox
            saconc(i,1,k,ispc) = rtlbnd(ispc)
            saconc(i,noy,k,ispc) = rtlbnd(ispc)
          enddo
        enddo
      enddo
c
c --- get the type of file to make sure it is a BC file ---
c
      this_var = 'CAMx_NAME'
      ierr = nf_get_att_text(iorbc, NF_GLOBAL, this_var, name_in)
      if( ierr .NE. NF_NOERR ) then
         this_var = 'NAME'
         ierr = nf_get_att_text(iorbc, NF_GLOBAL, this_var, name_in)
         if( ierr .NE. NF_NOERR ) goto 7003
      endif
      if( name_in(:8) .NE. bcname(:8) ) goto 7002
c
c --- call routine to make sure grid defintion is consistent ---
c
      call ncf_chk_griddef(iorbc,action,1,.TRUE.,.TRUE.,.FALSE.,.FALSE.,itmp)
c
c --- call routine to make sure file spans the episode ---
c
      call ncf_chk_tstep(iorbc,action,begdate,begtim,begdate,begtim,.FALSE.)
c
c  --- call routine to zero out the array ---
c
      allocate( bcgrid(ncol(1),nrow(1),nlay(1)) )
c
c  --- set the parameters for how to read the this data ---
c
      data_start(1) = 1
      data_count(1) = ncol(1)
      data_start(2) = 1
      data_count(2) = nrow(1)
      data_start(3) = 1
      data_count(3) = nlay(1)
c
c  --- get the timestep for the first hour ---
c
      this_tstep = ncf_get_tstep(iorbc,action,begdate,begtim,
     &          ibcdate_time_tflag,ibcdate_time_etflag,.FALSE.,.TRUE.)
      data_start(4) = this_tstep
      data_count(4) = 1
c
c  --- read the concentrations for each species ---
c
      do ispc=1,ntotsp
c
c  ---- load this species into local array ---
c
         this_var = ptname(ispc)
         ierr = nf_inq_varid(iorbc,this_var,this_varid)
         if( ierr .NE. NF_NOERR ) cycle
         ierr = nf_get_vara_real(iorbc,this_varid,data_start,data_count,bcgrid)
         if( ierr .NE. NF_NOERR) then
           write(iout,'(//,a)') 'ERROR in NCF_RDBCRT:'
           write(iout,*)'Cannot read boundary conditions data for speices: ',
     &                                   spname(ispc)(:istrln(spname(ispc)))
           write(iout,'(A,I5)') 'NetCDF error code: ',ierr
           call camxerr()
         endif
c
c  --- load into global array ---
c
         do k=1,noz
           do j=1,noy
c
c  --- get the gas conversion factor to umol/m3 ---
c
              if( ispc .LE. ngas ) then
                 convfac = densfac*273./
     &                               tpgrd(1,j,k)*prgrd(1,j,k)/1013.
              else
                 convfac = 1.
              endif

c  --- load the west boundary ---
c
              if( lbndry ) then
                 ioff = IDXBWS
              else
                 ioff = 1
              endif
              if( bcgrid(1,j,k) .GT. bdnl(ispc) ) then
                  saconc(1,j,k,ispc) = bcgrid(1,j,k) * convfac
              else
                  saconc(1,j,k,ispc) = 0.0
              endif
c
c  --- load the east boundary ---
c
c  --- get the gas conversion factor to umol/m3 ---
c
              if( ispc .LE. ngas ) then
                  convfac = densfac*273./
     &                          tpgrd(nox,j,k)*prgrd(nox,j,k)/1013.
              else
                  convfac = 1.
              endif
              if( lbndry ) then
                  ioff = IDXBES
              else
                  ioff = 1
              endif
              if( bcgrid(nox,j,k) .GT. bdnl(ispc) ) then
                  saconc(nox,j,k,ispc) = bcgrid(nox,j,k) * convfac
              else
                 saconc(nox,j,k,ispc) = 0.0
              endif
           enddo
c
c  --- load the south boundary ---
c
           do i=1,nox
c
c  --- get the gas conversion factor to umol/m3 ---
c
              if( ispc .LE. ngas ) then
                 convfac = densfac*273./
     &                           tpgrd(i,1,k)*prgrd(i,1,k)/1013.
              else
                 convfac = 1.
              endif
              if( lbndry ) then
                 ioff = IDXBST
              else
                 ioff = 1
              endif
              if( bcgrid(i,1,k) .GT. bdnl(ispc) ) then
                  saconc(i,1,k,ispc) = bcgrid(i,1,k) * convfac
              else
                  saconc(i,1,k,ispc) = 0.0
              endif
c
c  --- load the east boundary ---
c
c
c  --- get the gas conversion factor to umol/m3 ---
c
              if( ispc .LE. ngas ) then
                  convfac = densfac*273./
     &                                tpgrd(i,noy,k)*prgrd(i,noy,k)/1013.
              else
                  convfac = 1.
              endif
              if( lbndry ) then
                  ioff = IDXBNT
              else
                  ioff = 1
              endif
              if( bcgrid(i,noy,k) .GT. bdnl(ispc) ) then
                  saconc(i,noy,k,ispc) = bcgrid(i,noy,k) * convfac
              else
                  saconc(i,noy,k,ispc) = 0.0
              endif
           enddo
c
c  --- next layer --
c
         enddo
c
c  --- check the next modeled species ---
c
      enddo
c
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue 
      write(iout,'(//,A)') 'ERROR in NCF_RDBCRT:'
      write(iout,*) 'ERROR:  BC file for RTRAC does not exist: ',
     &                                             bcfil(:istrln(bcfil))
      call camxerr()
c
 7001 continue
      write(iout,'(//,A)') 'ERROR in NCF_RDBCRT:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Could not open file: ',
     &                                   bcfil(:istrln(bcfil))
      write(iout,'(3A)') 'If this is a NCF file, check that its',
     &                   ' format is consistent with the NCF library',
     &                   ' used to build CAMx.'
      write(iout,'(2A)') 'The CAMx makefile supports netCDF3-Classic,',
     &                   ' and netCDF4 compressed or uncompressed.'
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7002 continue
      write(iout,'(//,A)') 'ERROR in NCF_RDBCRT:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)',ERR=9999) 'Filetype is not correct.'
      write(iout,'(2A)',ERR=9999) 'Expecting    : ',bcname(:istrln(bcname))
      write(iout,'(2A)',ERR=9999) 'Value in file: ',name_in(:istrln(name_in))
      call camxerr()
c
 7003 continue
      write(iout,'(//,A)') 'ERROR in NCF_RDBCRT:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)',ERR=9999) 'Cannot find global variable for ',
     &                                       'type of file: CAMx_NAME'
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
c-----------------------------------------------------------------------
c    Format statements:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
c
      return
      end
