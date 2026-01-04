c**** NCF_RDTCDDM.F
c
      subroutine ncf_rdtcddm(toptim,topdate)
      use filunit
      use grid
      use chmstry
      use camxfld
      use camxcom
      use tracer
      implicit none
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routine fills one hour of top concentrations for the DDM
c   species. It then places these concentrations in the  appropriate 
c   place in the gridded array used for tracer concentrations.  This
c   version is for the NetCDF files.
c
c      Copyright 1996 - 2022
c     Ramboll
c
c      Argument description:
c       Outputs:
c       Inputs:
c        toptim            model simulation time (HHMM)
c        topdate           model simulation date (YYJJJ)
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     11/28/19   --gwilson--       Created
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
c    Arguement declarations:
c-----------------------------------------------------------------------
c
      integer topdate
      real    toptim
c
c-----------------------------------------------------------------------
c    External functions:
c-----------------------------------------------------------------------
c
      integer ncf_chkfile
      integer ncf_get_tstep
      integer istrln
c
c-----------------------------------------------------------------------
c    Local variables:
c----------------------------------------------------------------------
c
      character*200 action
      character*10  tcname, this_var, cspec, name_in
      integer       data_start(3), data_count(3), this_varid, i, j
      integer       nedge, ioff, ierr, this_tstep, ispc, n3d, iptr
      integer       itcdate_time_tflag, itcdate_time_etflag, ip, ic
      integer       ig, jerr, iddm, itmp
      logical       lexist, luse
c
      real, allocatable, dimension(:,:) :: tpconc
c
c-----------------------------------------------------------------------
c    Data statements:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c  --- set the number of BC edges --
c
      if( lbndry ) then
        nedge = 5
        ioff = IDXBTP
      else
        nedge = 1
        ioff = 1
      endif
c
c  --- open the TC ---
c
      tcname = 'AIRQUALITY'
      inquire(file=tcfil,exist=lexist)
      if( .NOT. lexist ) goto 7000
      action = 'Opening DDM Top Conocentrations file.'
      ierr = ncf_chkfile(iortc,tcfil,action,tcname)
      if( ierr .NE. ISUCES ) then
         jerr = nf_open(tcfil,NF_NOWRITE,iortc)
         if( jerr .NE. NF_NOERR ) goto 7001
      else
         return
      endif
c
c --- get the type of file to make sure it is a IC file ---
c
      this_var = 'CAMx_NAME'
      tcname = 'TOPCONC   '
      ierr = nf_get_att_text(iortc, NF_GLOBAL, this_var, name_in)
      if( ierr .NE. NF_NOERR ) then
         this_var = 'NAME'
         ierr = nf_get_att_text(iortc, NF_GLOBAL, this_var, name_in)
         if( ierr .NE. NF_NOERR ) goto 7003
       endif
      if( name_in(:7) .NE. tcname(:7) ) goto 7002
c
c --- call routine to make sure grid defintion is consistent ---
c
      call ncf_chk_griddef(iortc,action,1,.TRUE.,.TRUE.,.FALSE.,.FALSE.,itmp)
c
c --- call routine to make sure file spans the episode ---
c
      call ncf_chk_tstep(iortc,action,topdate,toptim,topdate,toptim,.FALSE.)
c
c  --- call routine to zero out the array ---
c
      allocate( tpconc(ncol(1),nrow(1)) )
c
c  --- set the parameters for how to read the this data ---
c
      data_start(1) = 1
      data_count(1) = ncol(1)
      data_start(2) = 1
      data_count(2) = nrow(1)
c
c  --- get the timestep for the first hour ---
c
      this_tstep = ncf_get_tstep(iortc,action,topdate,toptim,
     &          itcdate_time_tflag,itcdate_time_etflag,.FALSE.,.TRUE.)
      data_start(3) = this_tstep
      data_count(3) = 1
c
c  --- read the concentrations for each species ---
c
      do ispc=1,nspec
         this_var = spname(ispc)
         cspec = spname(ispc)
         ierr = nf_inq_varid(ioric,this_var,this_varid)
         if( ierr .NE. NF_NOERR ) cycle
         ierr = nf_get_vara_real(iortc,this_varid,data_start,data_count,tpconc)
         if( ierr .NE. NF_NOERR) then
           write(iout,'(//,a)') 'ERROR in NCF_RDTCDDM:'
           write(iout,*)'Cannot read top concenntration data for species: ',
     &                                   spname(ispc)(:istrln(spname(ispc)))
           write(iout,'(A,I5)') 'NetCDF error code: ',ierr
           call camxerr()
         endif
c
c  --- find this species in the modeled species list ---
c
         do iddm=1,nbcddm
           luse = .FALSE.
           if( bcddmsp(iddm) .EQ. cspec ) luse = .TRUE. 
           if( bcddmsp(iddm) .EQ. NAMALL  ) luse = .TRUE. 
           if( bcddmsp(iddm) .EQ. NAMVOC 
     &                               .AND. lvocsp(ispc) ) luse = .TRUE. 
           if( bcddmsp(iddm) .EQ. NAMNOX 
     &                               .AND. lnoxsp(ispc) ) luse = .TRUE. 
           if( bcddmsp(iddm) .EQ. NAMHRV 
     &                               .AND. lhrvoc(ispc) ) luse = .TRUE.
c
c  --- if this DDM species matches this modeled species, load it ---
c
           if( luse ) then
              iptr = iptddm(ispc) + (iddm-1)*nedge + nicddm + ioff - 1
              do j=1,nrow(1)
              do i=1,ncol(1)
                 n3d = i + ncol(1)*(j - 1) + ncol(1)*nrow(1)*(iptr - 1)
c
c  --- load the top concentrations ---
c
                 if( tpconc(i,j) .GT. bdnl(ispc) ) then
                     pttop(n3d) = tpconc(i,j)
                 else
                     pttop(n3d) = 0.0
                 endif
              enddo
              enddo
c
c  --- interpolate top con sens to all nested grids
c
              if( ngrid .GT. 1 ) then
                 do ip = 1,ngrid
                 do ic = 1,nchdrn(ip)
                    ig = idchdrn(ic,ip)
                    call interp2d(ncol(ip),nrow(ip),1,
     &                   i1(ig),j1(ig),nmesh(ig),ncol(ig),nrow(ig),
     &                   pttop( ipsa2d(ip) + ncol(ip)*nrow(ip)*(iptr - 1) ),
     &                   pttop( ipsa2d(ig) + ncol(ig)*nrow(ig)*(iptr - 1) ))
                 enddo
                 enddo
              endif
           endif
c
c  --- check the next BC DDM species ---
c
         enddo
c
c  --- read next species worth of data ---
c
      enddo
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue 
      write(iout,'(//,A)') 'ERROR in NCF_RDTCDDM:' 
      write(iout,*) 'ERROR:  TC file for DDM does not exist: ',TRIM(tcfil)
      call camxerr()
c
 7001 continue
      write(iout,'(//,A)') 'ERROR in NCF_RDTCDDM:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Could not open file: ',
     &                                   tcfil(:istrln(tcfil))
      write(iout,'(3A)') 'If this is a NCF file, check that its',
     &                   ' format is consistent with the NCF library',
     &                   ' used to build CAMx.'
      write(iout,'(2A)') 'The CAMx makefile supports netCDF3-Classic,',
     &                   ' and netCDF4 compressed or uncompressed.'
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7002 continue
      write(iout,'(//,A)') 'ERROR in NCF_RDTCDDM:' 
      write(iout,*) 'ERROR:  TC file for DDM is not labelled ',TRIM(tcname)
      call camxerr()
c
 7003 continue
      write(iout,'(//,A)') 'ERROR in NCF_RDTCDDM:'
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
