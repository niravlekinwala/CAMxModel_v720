      subroutine sa_icprep(begtim,begdate,endtim,enddate)
      use tracer
      use filunit
      use grid
      implicit none
c 
c----CAMx v7.20 220430
c 
c     SA_ICPREP reads the initial conditions file for Source Appointionment
c     and sets up the data structures to map the species into the 
c     tracer arrays.
c                           
c      Copyright 1996 - 2022
c     Ramboll
c           
c     Modifications: 
c     Input arguments: 
c        begtim              model start time (HHMM)
c        begdate             model start date (YYJJJ)
c        endtim              model end time (HHMM)
c        enddate             model end date (YYJJJ)
c             
c     Output arguments: 
c        none
c             
c     Routines Called: 
c        none
c             
c     Called by: 
c        INITSA 
c
      include 'camx.prm'
      include 'flags.inc'
c
c     Argument delcarations:
c
      integer begdate
      real    begtim
      integer enddate
      real    endtim
c
c     Local variables 
c
      character*10 icname, infil, icspc
      character*4  ifile(10), note(60)
      character*3  first_class
      integer      nseg, idat1, idat2, n, nx, ny, nz, izone
      integer      l
      real         tim1, tim2, orgx, orgy, utmx, utmy, dx, dy
c
      character*4 icspec(10,MXTRSP)
c
      data icname /'AIRQUALITY'/
c
c-----Entry point
c
c-----Read 1st IC header record and check inputs
c
      rewind(ioric)
      read(ioric) ifile,note,nseg,num_ioric,idat1,tim1,idat2,tim2
c
c----check for array overflow ---
c
      if( num_ioric .GT. MXTRSP ) then
         write(iout,'(//,A)') 'ERROR in SA_ICPREP:'
         write(iout,*) 'A parameter in the camx.prm is not ',
     &                                        'sufficiently large.'
         write(iout,*) 'Please change the value for parameter: MXTRSP'
         write(iout,*) 'It should be set to a value of at least: ',num_ioric
         call flush(iout)
         call camxerr()
      endif
      allocate( spc_ioric(num_ioric) )
c
      if(INT(tim2) .EQ. 24) then
        idat2 = idat2 + 1
        tim2 = 0.
        if( MOD(idat2,1000) .GT. 365 ) then
           if( MOD(INT(idat2/1000),4) .EQ. 0 ) then
              if( MOD(idat2,1000) .EQ. 367 )
     &                     idat2 = (INT(idat2/1000)+1)*1000 + 1
           else
              idat2 = (INT(idat2/1000)+1)*1000 + 1
           endif
        endif
      endif
c
      write(infil,'(10a1)') (ifile(n),n=1,10)
      if( infil .NE. icname ) then
        write(iout,'(//,a)') 'ERROR in SA_ICPREP:'
        write(iout,*)'IC input file is not labelled AIRQUALIY'
        call camxerr()
      endif
c
      tim1 = 100.*tim1
      tim2 = 100.*tim2
      if( idat1 .GT. begdate ) then
        write(iout,'(//,a)') 'ERROR in SA_ICPREP:'
        write(iout,*)'Source Apportionment IC start date > simulation start date'
        write(iout,*)'  IC file: ',idat1
        write(iout,*)'Sim start: ',begdate
        call camxerr()
      elseif( idat1 .EQ. begdate .AND. tim1 .gt. begtim) then
        write(iout,'(//,a)') 'ERROR in SA_ICPREP:'
        write(iout,*)'Source Apportionment IC start time > simulation start time'
        write(iout,*)'  IC file: ',tim1
        write(iout,*)'Sim start: ',begtim
        call camxerr()
      endif
c
c-----Read 2nd IC header record and check inputs
c
      read(ioric) orgx,orgy,izone,utmx,utmy,dx,dy,nx,ny,nz
      if(.NOT. llatlon ) then
        dx = dx/1000.
        dy = dy/1000.
      endif
      if( ABS(dx-delx) .GT. 0.001 .OR. ABS(dy-dely) .GT. 0.001 ) then
        write(iout,'(//,a)') 'WARNING in SA_ICPREP:'
        write(iout,*)'Source Apportionment IC cell size not equal to model cell size'
        write(iout,'(a,2f10.4)')'  IC file: ',dx,dy
        write(iout,'(a,2f10.4)')'    model: ',delx,dely
        write(iout,*)
      elseif( nx .NE. ncol(1) .OR. ny .NE. nrow(1) .OR. nz .NE. nlay(1) ) then
        write(iout,'(//,a)') 'ERROR in SA_ICPREP:'
        write(iout,*)'Source Apportionment IC grid size not equal to model grid size'
        write(iout,*)'IC file: ',nx,ny,nz
        write(iout,*)'  model: ',ncol(1),nrow(1),nlay(1)
        write(iout,*)
        call camxerr()
      endif 
c
      ncls_ioric = 0
      read(ioric)
      read(ioric) ((icspec(n,l),n=1,10),l=1,num_ioric)
      write(first_class,'(3A1)') (icspec(n,1),n=1,3)
      do l=1,num_ioric
        write(spc_ioric(l),'(10A1)') (icspec(n,l),n=1,10)
        if( spc_ioric(l)(1:3) .EQ. first_class ) ncls_ioric = ncls_ioric + 1
      enddo
c
      return
      end
