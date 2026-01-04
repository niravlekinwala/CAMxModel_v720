      Module filunit
      include 'filunit.inc'
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c        This allocates the dynamic memory arrays in the FILUNIT.COM
c        include file.
c-----------------------------------------------------------------------
c
c    Argument descriptions:
c     Input:  
c     Output:  
c
c      Copyright 1996 - 2022
c     Ramboll
c
c-----------------------------------------------------------------------
c    LOG:
c        01/04/11  Revised for new met input format
c        04/30/13  Added surface model
c        09/02/14  Added subgrid convective model
c-----------------------------------------------------------------------
c
      Contains
c
c-----------------------------------------------------------------------
c   BEGIN SUBROUTINE ALLOC_FILUNIT
c-----------------------------------------------------------------------
c
         subroutine alloc_filunit(numgrds,numemissfiles,numpointfiles)
c
c-----------------------------------------------------------------------
c    Modules used:
c-----------------------------------------------------------------------
c
         use camx_includes
         implicit none
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c    Argument descriptions:
c     Input:  
c        numgrds        I  number of grids
c        numemissfiles  I  number of emissions files
c        numpointfiles  I  number of point source files
c     Output:  
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
         integer :: numgrds
         integer :: numemissfiles(*)
         integer :: numpointfiles
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c         
         integer :: max_emissfiles
         integer :: max_ptsfiles
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c         
         allocate( iavg  (numgrds) )
         allocate( idep  (numgrds) )
         allocate( ismin (numgrds) )
         allocate( ismout(numgrds) )
         max_emissfiles = maxval( numemissfiles(1:numgrds) )
         max_emissfiles = MAX( max_emissfiles,1 )
         max_ptsfiles = MAX( numpointfiles,1 )
         allocate( iarem (numgrds,max_emissfiles) )
         allocate( is_netcdf_iarem (numgrds,max_emissfiles) )
         allocate( buffer_offset_iarem (numgrds,max_emissfiles) )
         if( max_ptsfiles .GT. 0 ) 
     &                  allocate( iptem (max_ptsfiles) )
         if( max_ptsfiles .GT. 0 ) 
     &                  allocate( is_netcdf_iptem (max_ptsfiles) )
         allocate( isurf (numgrds) )
         allocate( is_netcdf_isurf (numgrds) )
         allocate( i3dmet(numgrds) )
         allocate( is_netcdf_i3dmet(numgrds) )
         allocate( i2dmet(numgrds) )
         allocate( is_netcdf_i2dmet(numgrds) )
         allocate( ikv   (numgrds) )
         allocate( is_netcdf_ikv   (numgrds) )
         allocate( icld  (numgrds) )
         allocate( lcig  (numgrds) )
c
         allocate( n3dmet(numgrds) )
         allocate( n2dmet(numgrds) )
         allocate( nkvmet(numgrds) )
         allocate( ncldmet(numgrds) )
         allocate( nsrfvar(numgrds) )
c
         return
         end subroutine
c
c-----------------------------------------------------------------------
c   END SUBROUTINE ALLOC_FILUNIT
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c   BEGIN SUBROUTINE ALLOC_FILUNIT_SAMPLE
c-----------------------------------------------------------------------
c
         subroutine alloc_filunit_sample(numsamples)
c
c-----------------------------------------------------------------------
c    Modules used:
c-----------------------------------------------------------------------
c
         use camx_includes
         implicit none
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c    Argument descriptions:
c     Input:  
c        numsamples   I  number of sampling grids
c     Output:  
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
         integer :: numsamples
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c         
         allocate( isample(numsamples) )
c
         return
         end subroutine
c-----------------------------------------------------------------------
c   END SUBROUTINE ALLOC_FILUNIT_SAMPLE
c-----------------------------------------------------------------------
c
      end Module
