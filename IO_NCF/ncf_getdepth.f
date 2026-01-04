C**** NCF_GETDEPTH
c
      subroutine ncf_getdepth(nvar,ncols,nrows,nlays,ibgdhp,idtnow,btimhp,
     &                    timnow,iunit,height,depth)
      use filunit
      implicit none
c
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c     GETDEPTH reads the 3D met file until to the current hour,
c     calculates layer depth, and passes back the data for OSAT calculations
c
c
c      Copyright 1996 - 2022
c     Ramboll
c
c      Argument description:
c       Inputs:
c        nvar   I number of variables on 3D met file
c        ncols  I number of columns
c        nrows  I number of rows
c        nlays  I number of layers
c        idtnow I current date (YYJJJ)
c        timnow R current time (HH.HH)
c        iunit  I height/pressure file unit number
c             
c     Output arguments: 
c        ibgdhp I date from height/pressure file (YYJJJ)
c        btimhp R time from height/pressure file (HH.HH) 
c        height R array of gridded layer heights (m)
c        depth  R array of gridded layer depths (m)
c             
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
      include "camx.prm"
      include "netcdf.inc"
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
c
      integer nvar
      integer ncols
      integer nrows
      integer nlays
      integer ibgdhp
      integer idtnow
      real    btimhp
      real    timnow
      integer iunit
      real    height(ncols,nrows,nlays)
      real    depth(ncols,nrows,nlays)
c
c-----------------------------------------------------------------------
c    External functions:
c-----------------------------------------------------------------------
c
      integer istrln
      integer ncf_get_tstep
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      character*200 action
      character*20  this_var
      integer       this_varid, ierr, data_start(4), data_count(4)
      integer       this_date, this_time_tflag, this_time_etflag
      integer       this_tstep, i, j, k
      real          this_time
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
      if( .NOT. is_netcdf_i3dmet(1) ) goto 9999
      action = 'Getting layer heights from 3DMet file for grid 1: '
c
c   --- get the timestep that encompasses this date/time ---
c
      this_date = idtnow
      this_time = timnow
      if( this_time .GE. 24. ) then
        this_time = ANINT(this_time - 24.)
        this_date = this_date + 1
        if( MOD(this_date,1000) .GT. 365 ) then
           if( MOD(INT(this_date/1000),4) .EQ. 0 ) then
              if( MOD(this_date,1000) .EQ. 367 )
     &           this_date = (INT(this_date/1000)+1)*1000 + 1
           else
              this_date = (INT(this_date/1000)+1)*1000 + 1
           endif
        endif
      endif
      this_tstep = ncf_get_tstep(i3dmet(1),action,this_date,this_time*100.,
     &                 this_time_tflag,this_time_etflag,.FALSE.,.FALSE.)
c
c   --- set the indexes for what to read ---
c
      data_start(1) = 1
      data_count(1) = ncols
      data_start(2) = 1
      data_count(2) = nrows
      data_start(3) = 1
      data_count(3) = nlays
      data_start(4) = this_tstep
      data_count(4) = 1
c
c ---- layer heights ---
c
      this_var = "z"
      ierr = nf_inq_varid(i3dmet(1), this_var, this_varid)
      if( ierr .NE. NF_NOERR) goto 7000
      ierr = nf_get_vara_real(i3dmet(1),this_varid,data_start,
     &                                          data_count,height)
      if( ierr .NE. NF_NOERR) goto 7001
c
c  --- Calculate the depths from the heights ---
c
      do k = 1,nlays
        do j = 1,nrows
          do i = 1,ncols
            depth(i,j,k) = height(i,j,k)
            if( k .GT. 1 ) depth(i,j,k) = height(i,j,k) - height(i,j,k-1)
           enddo
         enddo
      enddo
c
c  --- successful completion ---
c
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,a)') 'ERROR in NCF_GETFEPTH:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot find variable id for: ',
     &                                      this_var(:istrln(this_var))
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7001 continue
      write(iout,'(//,a)') 'ERROR in NCF_GETFEPTH:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot read data for variable: ',
     &                                      this_var(:istrln(this_var))
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 9999 continue
      return
      end
