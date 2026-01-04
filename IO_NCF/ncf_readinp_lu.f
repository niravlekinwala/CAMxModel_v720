C**** NCF_READINP_LU
c
      subroutine ncf_readinp_lu(igrd,ncol,nrow,num_luse,fsurf,topo,
     &                                           lai,lrdlai,icdocn)
      use filunit
      implicit none
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routines reads the next hour of data in the landuse file and
c   loads the data into global arrays
c 
c      Copyright 1996 - 2022
c     Ramboll
c           
c     Input arguments:
c        igrd                grid index
c        ncol                number of columns
c        nrow                number of rows
c        num_luse            number of landuse categeries
c        fsurf               fractional landuse field
c        topo                topographic elevation (m MSL, optional)
c        lai                 lead area index field (optional)
c        lrdlai              flag indicating if LAI was read
c        icdocn              ocean cell flag
c
c     Output arguments:
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
      integer num_luse
      real    fsurf(ncol,nrow,num_luse)
      real    topo(ncol,nrow)
      real    lai(ncol,nrow)
      logical lrdlai
      integer icdocn(ncol,nrow)
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
      character*10  namsrf(NLUZ03)
      integer       lumap(NLUZ03), mapvar(NLUZ03), i, j, l, ilu, ifresh
      integer       this_varid, ierr, data_start(4), data_count(4)
      real          areatot
c
      real, allocatable, dimension(:,:) :: array2d
c
c-----------------------------------------------------------------------
c    Data statements:
c-----------------------------------------------------------------------
c
      data namsrf   /'water     ',
     &               'ice       ',
     &               'lake      ',
     &               'eneedl    ',
     &               'ebroad    ',
     &               'dneedl    ',
     &               'dbroad    ',
     &               'tbroad    ',
     &               'ddecid    ',
     &               'eshrub    ',
     &               'dshrub    ',
     &               'tshrub    ',
     &               'sgrass    ',
     &               'lgrass    ',
     &               'crops     ',
     &               'rice      ',
     &               'sugar     ',
     &               'maize     ',
     &               'cotton    ',
     &               'icrops    ',
     &               'urban     ',
     &               'tundra    ',
     &               'swamp     ',
     &               'desert    ',
     &               'mwood     ',
     &               'tforest   '/
      data lumap / 7, 8, 7, 5, 5, 4, 4, 5, 4, 3,
     &             3, 3, 3,10, 2, 2, 2, 2, 2, 2,
     &             1,11, 9, 8, 6, 6/
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c  ---- initialize the data ----
c
      lrdlai = .FALSE.
      topo = 0.
      fsurf = 0.
      lai = 0.
      icdocn = 0
c
c   --- set the indexes for what to read ---
c
      data_start(1) = 1
      data_count(1) = ncol
      data_start(2) = 1
      data_count(2) = nrow
      data_start(3) = 1
      data_count(3) = 1
      data_start(4) = 1
      data_count(4) = 1
c
c  ---- allocate temperary array ---
c
      allocate( array2d(ncol,nrow) )
c
      if( .NOT. is_netcdf_isurf(igrd) ) goto 9999
      write(action,'(A,I3)') 'Reading landuse file for grid ',igrd
c
c  --- loop over all landuse species and check if it is in file ---
c
      do ilu=1,NLUZ03
        this_var = namsrf(ilu)
        ierr = nf_inq_varid(isurf(igrd), this_var, this_varid)
        if( ierr .NE. NF_NOERR) goto 7000
c
c   --- read data for this category ---
c
        ierr = nf_get_vara_real(isurf(igrd),this_varid,data_start,
     &                                            data_count,array2d)
        if( ierr .NE. NF_NOERR) goto 7001
c
c  ---- set ocean flag based on water ---
c
        if( namsrf(ilu) .EQ. 'water     ') then
            ifresh = 0
            do j=1,nrow
              do i=1,ncol
                if( array2d(i,j) .GT. 0. ) icdocn(i,j) = 1
              enddo
            enddo
        endif
c
c  ---- check for fresh water ---
c
        if( namsrf(ilu) .EQ. 'lake      ') then
            ifresh = 0
            do j=1,nrow
              do i=1,ncol
                if( array2d(i,j) .GT. 0. ) ifresh = 1
              enddo
            enddo
        endif
c
c   --- read data and load into global array ---
c
        if( idrydep .EQ. 2) then
            fsurf(:,:,ilu) = array2d(:,:)
        else
            fsurf(:,:,lumap(ilu)) = fsurf(:,:,lumap(ilu)) + array2d(:,:)
        endif
      enddo
c
c  --- make sure data is value for inline emissions ---
c
      if( ifresh .EQ. 0 .AND. lixemis ) then
        write(iout,'(//,A)') 'WARNING in NCF_READINP_LU:'
        write(iout,'(A)') 'No fresh water coverage found in LU file'
        write(iout,'(A,i3,a,2i4)') 'Grid = ', igrd
        write(iout,'(2A)') 'If you are running with in-line halogen',
     &                     ' emissions, halogens will be emitted'
        write(iout,'(A)') 'for all water bodies.'
        write(iout,'(2A)') 'Use the WATERMASK program to convert your LU',
     &                   ' file to differentiate between'
        write(iout,'(A,//)') 'ocean and fresh water bodies.'
        write(*,'(//,A)') 'WARNING in NCF_READINP_LU:'
        write(*,'(A)') 'No fresh water coverage found in LU file'
        write(*,'(A,i3,a,2i4)') 'Grid = ', igrd
        write(*,'(2A)') 'If you are running with in-line halogen',
     &                     ' emissions, halogens will be emitted'
        write(*,'(A)') 'for all water bodies.'
        write(*,'(2A)') 'Use the WATERMASK program to convert your LU',
     &                   ' file to differentiate between'
        write(*,'(A,//)') 'ocean and fresh water bodies.'
      endif
c
c  --- Leaf area index ---
c
      this_var = 'lai'
      ierr = nf_inq_varid(isurf(igrd), this_var, this_varid)
      if( ierr .EQ. NF_NOERR) then
        lrdlai = .TRUE.
        ierr = nf_get_vara_int(isurf(igrd),this_varid,data_start,
     &                                            data_count,lai)
        if( ierr .NE. NF_NOERR) goto 7001
      endif
c
c  --- topographic elevation ---
c
      this_var = 'topo'
      ierr = nf_inq_varid(isurf(igrd), this_var, this_varid)
      if( ierr .EQ. NF_NOERR) then
        write(iout,'(a40,15x,a,i3)') 'Read TOPO field',' grid',igrd
        ierr = nf_get_vara_real(isurf(igrd),this_varid,data_start,
     &                                            data_count,topo)
        if( ierr .NE. NF_NOERR) goto 7001
      endif
c
c --- Check that the LU data are reasonable and adjust out
c     minor inconsistencies ---
c
      do j = 1,nrow
        do i = 1,ncol
           areatot = 0.0
           do l = 1,num_luse
             areatot = areatot + fsurf(i,j,l)
           enddo
           if( areatot .LT. 0.95 .OR. areatot .GT. 1.05) goto 7002
           do l = 1,num_luse
             fsurf(i,j,l) = fsurf(i,j,l)/areatot
           enddo
         enddo
      enddo
c
c  ---- deallocate temperary array ---
c
      deallocate( array2d )
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
      write(iout,'(//,a)') 'ERROR in NCF_READINP_LU:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot find variable id for: ',
     &                                      this_var(:istrln(this_var))
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7001 continue
      write(iout,'(//,a)') 'ERROR in NCF_READINP_LU:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Cannot read data for variable: ',
     &                                      this_var(:istrln(this_var))
      write(iout,'(A,I5)') 'NetCDF error code: ',ierr
      call camxerr()
c
 7002 write(iout,'(//,A)') 'ERROR in NCF_READINP_LU::'
      write(iout,'(/,2A)') 'Sum of landuse fractions differs from 1.0 ',
     &                     'by more than 5%'
      write(iout,'(A,i3,a,2i4)')
     &       'Grid = ', igrd, '  Cell(i,j) = ', i, j
      write(iout,'(A)') 'Table of input landuse data follows:'
      write(iout,'(/,A)') ' Class    Fraction'
      write(iout,'(i6,F10.3)') (l,fsurf(i,j,l),l=1,num_luse)
      write(iout,'(A6,F10.3)') 'Total', areatot
      write(iout,'(/,A)') 'Check your input landuse data file'
      call camxerr()
c
 9999 continue
      return
      end
