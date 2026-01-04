C**** NCF_METINIT_3D
c
      subroutine ncf_metinit_3d(igrd,ncol,nrow,nlay,height,press,depth,
     &                                          windu,windv,tempk,water)
      use filunit
      implicit none
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routines initalizes the met fields by reading the NetCDF
c   met files and loading the data into the global arrays.
c 
c      Copyright 1996 - 2022
c     Ramboll
c           
c     Input arguments:
c        igrd                grid index
c        ncol                number of columns
c        nrow                number of rows
c        nlay                number of layers
c
c     Output arguments:
c        height              layer interface height field (m)
c        press               layer midpoint pressure field (mb)
c        depth               layer depth (m)
c        windu               u-component wind field (m/s)
c        windv               u-component wind field (m/s)
c        tempk               temperature field (K)
c        water               water vapor field (ppm)
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     02/20/17   --gwilson--    Original development
c     07/23/19   Staggered wind flag is now grid-specific
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'camx.inc'
      include 'flags.inc'
      include 'netcdf.inc'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      integer igrd
      integer ncol
      integer nrow
      integer nlay
      real    height(ncol,nrow,nlay)
      real    press(ncol,nrow,nlay)
      real    depth(ncol,nrow,nlay)
      real    windu(ncol,nrow,nlay)
      real    windv(ncol,nrow,nlay)
      real    tempk(ncol,nrow,nlay)
      real    water(ncol,nrow,nlay)
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
      integer       this_tstep, this_time_tflag, this_time_etflag
      integer       istag, i, j, k
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
      if( .NOT. is_netcdf_i3dmet(igrd) ) goto 9999
c
      write(action,'(A,I3)') 'Reading 3D met file for grid ',igrd
c
c   --- get the timestep that encompasses this date/time ---
c
      this_tstep = ncf_get_tstep(i3dmet(igrd),action,begdate,begtim,
     &                 this_time_tflag,this_time_etflag,.FALSE.,.FALSE.)
c
c   --- set the indexes for what to read ---
c
      data_start(1) = 1
      data_count(1) = ncol
      data_start(2) = 1
      data_count(2) = nrow
      data_start(3) = 1
      data_count(3) = nlay
      data_start(4) = this_tstep
      data_count(4) = 1
c
c ---- layer heights ---
c
      this_var = "z"
      ierr = nf_inq_varid(i3dmet(igrd), this_var, this_varid)
      if( ierr .NE. NF_NOERR) goto 7000
      ierr = nf_get_vara_real(i3dmet(igrd),this_varid,data_start,
     &                                          data_count,height)
      if( ierr .NE. NF_NOERR) goto 7001
c
c ---- pressure ---
c
      this_var = "pressure"
      ierr = nf_inq_varid(i3dmet(igrd), this_var, this_varid)
      if( ierr .NE. NF_NOERR) goto 7000
      ierr = nf_get_vara_real(i3dmet(igrd),this_varid,data_start,
     &                                          data_count,press)
      if( ierr .NE. NF_NOERR) goto 7001
c
c ---- temperature ---
c
      this_var = "temperature"
      ierr = nf_inq_varid(i3dmet(igrd), this_var, this_varid)
      if( ierr .NE. NF_NOERR) goto 7000
      ierr = nf_get_vara_real(i3dmet(igrd),this_varid,data_start,
     &                                          data_count,tempk)
      if( ierr .NE. NF_NOERR) goto 7001
c
c ---- humidity ---
c
      this_var = "humidity"
      ierr = nf_inq_varid(i3dmet(igrd), this_var, this_varid)
      if( ierr .NE. NF_NOERR) goto 7000
      ierr = nf_get_vara_real(i3dmet(igrd),this_varid,data_start,
     &                                          data_count,water)
      if( ierr .NE. NF_NOERR) goto 7001
c
c ---- U-component of wind ---
c
      this_var = "uwind"
      ierr = nf_inq_varid(i3dmet(igrd), this_var, this_varid)
      if( ierr .NE. NF_NOERR) goto 7000
      ierr = nf_get_vara_real(i3dmet(igrd),this_varid,data_start,
     &                                          data_count,windu)
      if( ierr .NE. NF_NOERR) goto 7001
c
c ---- V-component of wind ---
c
      this_var = "vwind"
      ierr = nf_inq_varid(i3dmet(igrd), this_var, this_varid)
      if( ierr .NE. NF_NOERR) goto 7000
      ierr = nf_get_vara_real(i3dmet(igrd),this_varid,data_start,
     &                                          data_count,windv)
      if( ierr .NE. NF_NOERR) goto 7001
c
c-----Calculate layer depth
c
      do k = 1,nlay
        do j = 1,nrow
          do i = 1,ncol
            if (k.eq.1) then
              depth(i,j,k) = height(i,j,k)
            else
              depth(i,j,k) = height(i,j,k) - height(i,j,k-1)
            endif
            if(depth(i,j,k) .LE. 0.0 ) then
              write(iout,'(//,a)') 'ERROR in NCF_METINIT_3D:'
              write(iout,'(a,3i3)')
     &                  'Negative depth in metinit, i,j,k = ',i,j,k
              call camxerr()
            endif
          enddo
        enddo
      enddo
c
c-----Get staggered wind flag and Convert windu and windv 
c     if they are not staggered ---
c
      this_var = 'ISTAG'
      ierr = nf_get_att_int(i3dmet(igrd), NF_GLOBAL, 'ISTAG', istag)
      if( ierr .NE. NF_NOERR ) goto 7002
      lstagw(igrd) = .TRUE.
      if( istag .EQ. 0 ) lstagw(igrd) = .FALSE.
      if( lstagw(igrd) .NEQV. lstagw(1) ) then
        write(iout,'(//,a)') 'WARNING in RDMETHDR:'
        write(iout,'(A)') action(:istrln(action))
        write(iout,*)'Mismatch in staggered wind flag with master grid'
        write(iout,*) '         Master grid stagger is: ',lstagw(1)
        write(iout,*) '         Nested grid stagger is: ',lstagw(igrd)
        write(iout,*)
      endif
      if (.not.lstagw(igrd)) call cvtwind(ncol,nrow,nlay,windu,windv)
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
      write(iout,'(//,a)') 'ERROR in NCF_METINIT_3D:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot find variable id for: ',
     &                                      this_var(:istrln(this_var))
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7001 continue
      write(iout,'(//,a)') 'ERROR in NCF_METINIT_3D:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot read data for variable: ',
     &                                      this_var(:istrln(this_var))
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7002 continue
      write(iout,'(//,a)') 'ERROR in NCF_METINIT_3D:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot find global attribute: ',
     &                                      this_var(:istrln(this_var))
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 9999 continue
      return
      end
