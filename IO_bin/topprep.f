      subroutine topprep(begtim,begdate,endtim,enddate)
      use filunit
      use grid
      use chmstry
      use bndary
c 
c----CAMx v7.20 220430
c 
c     TOPPREP reads the top concentrations file header and sets
c     variables that map the TC species list to the internal CAMx
c     species list
c                           
c      Copyright 1996 - 2022
c     Ramboll
c           
c     Modifications: 
c        None
c 
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
c        STARTUP 
c
      include 'camx.prm'
      include 'flags.inc'
c
      integer begdate
      integer enddate
c
      character*10 tcfil, infil, tcspc
      character*4  ifile(10), note(60)
c
      character*4 tcspec(10,MXSPEC)
      integer     indx(MXCELLS,MXCELLS)
c
      data tcfil /'AIRQUALITY'/
c
c-----Entry point
c
c-----Read 1st TC header record and check inputs
c
      rewind(itc)
      read(itc) ifile,note,nseg,ntcspc,idat1,tim1,idat2,tim2
c
c----check for array overflow ---
c
      if( ntcspc .GT. MXSPEC ) then
         write(iout,'(//,A)') 'ERROR in TOPPREP:'
         write(iout,*) 'A parameter in the camx.prm is not ',
     &                                        'sufficiently large.'
         write(iout,*) 'Please change the value for parameter: MXSPEC'
         write(iout,*) 'It should be set to a value of at least: ',ntcspc
         call flush(iout)
         call camxerr()
      endif
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
      write(infil,'(10a1)') (ifile(n),n=1,10)
      if (infil.ne.tcfil) then
        write(iout,'(//,a)') 'ERROR in TOPPREP 1:'
        write(iout,*)'TC input file is not labelled AIRQUALITY'
        call camxerr()
      endif
      tim1 = 100.*tim1
      tim2 = 100.*tim2
      if (idat1.gt.begdate) then
        write(iout,'(//,a)') 'ERROR in TOPPREP 2:'
        write(iout,*)'TC start date > simulation start date'
        write(iout,*)'  TC file: ',idat1
        write(iout,*)'Sim start: ',begdate
        call camxerr()
      elseif (idat1.eq.begdate .and. tim1.gt.begtim) then
        write(iout,'(//,a)') 'ERROR in TOPPREP 3:'
        write(iout,*)'TC start time > simulation start time'
        write(iout,*)'  TC file: ',tim1
        write(iout,*)'Sim start: ',begtim
        call camxerr()
      elseif (idat2.lt.enddate) then
        write(iout,'(//,a)') 'ERROR in TOPPREP 4:'
        write(iout,*)'TC end date < simulation end date'
        write(iout,*)'TC file: ',idat2
        write(iout,*)'Sim end: ',enddate
        call camxerr()
      elseif (idat2.eq.enddate .and. tim2.lt.endtim) then
        write(iout,'(//,a)') 'ERROR in TOPPREP 5:'
        write(iout,*)'TC end time < simulation end time'
        write(iout,*)'TC file: ',tim2
        write(iout,*)'Sim end: ',endtim
        call camxerr()
      endif
c
c-----Read 2nd TC header record and check inputs
c
      read(itc) orgx,orgy,izone,utmx,utmy,dx,dy,nx,ny,nz
      if (.NOT.llatlon) then
        dx = dx/1000.
        dy = dy/1000.
      endif
      if (abs(dx-delx).gt.0.001 .or. abs(dy-dely).gt.0.001) then
        write(iout,'(//,a)') 'WARNING in TOPPREP:'
        write(iout,*)'TC cell size not equal to model cell size'
        write(iout,'(a,2f10.4)')'  TC file: ',dx,dy
        write(iout,'(a,2f10.4)')'    model: ',delx,dely
        write(iout,*)
      elseif (nx.ne.ncol(1) .or. ny.ne.nrow(1)) then
        write(iout,'(//,a)') 'ERROR in TOPPREP 6:'
        write(iout,*)'TC grid size not equal to model grid size'
        write(iout,*)'TC file: ',nx,ny
        write(iout,*)'  model: ',ncol(1),nrow(1)
        write(iout,*)
        call camxerr()
      endif 
c
c-----Read 3rd & 4th TC header 
c
      read(itc) 
      read(itc) ((tcspec(n,l),n=1,10),l=1,ntcspc)
c
c-----Map TC species to model species
c
      do 20 l = 1,nspec
        ltcmap(l) = 0
        do 15 ltc = 1,ntcspc
          write(tcspc,'(10a1)') (tcspec(n,ltc),n=1,10)
          if (tcspc.eq.'HNO2      ') tcspc = 'HONO      '
          if (tcspc.eq.'HCHO      ' .and. kHCHO.eq.nspec+1)
     &                                        tcspc = 'FORM      '
          if (tcspc.eq.spname(l)) then
            ltcmap(l) = ltc
            write(idiag,'(2(a,i5,2x,a))')'Topcon species   ',ltc,tcspc,
     &                   ' mapped to model species ',l,spname(l)
            goto 20
          endif
 15     continue
        write(idiag,*)'Did not find species: ',spname(l),' on TC file'
 20   continue
c
      return
      end
