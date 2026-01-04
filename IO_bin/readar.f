      subroutine readar(igrd,idxfile,ncol,nrow,nlay_ems,iounit,iout_unit,
     &                                     aremis,nspcar,nspcem,nspcs)
      use chmstry
      use camxcom
      use filunit
c
c----CAMx v7.20 220430
c
c     READAR reads the time-variant records of the area source file for
c     the given grid and cycles through to current time/date to load 
c     area emission rates
c
c      Copyright 1996 - 2022
c     Ramboll
c          
c     Modifications:
c        10/20/00  Added check for negative emission rates
c        05/01/03  Time span of emissions must now match emiss update interval
c        12/07/14  Revised for VBS emissions
c        01/08/16  Updated for Revised VBS emissions
c
c     Input arguments:
c        igrd                grid index
c        idxfile             emissions file number
c        ncol                number of columns
c        nrow                number of rows
c        iounit              file unit for area emissions file
c        iout_unit           file unit number for output file
c        nspcar              number of input gridded emission species
c        nspcem              number of total emissons species
c        nspcs               number of model species
c
c     Output arguments:
c        aremis              area emissions rate (mole/s)
c
c     Routines Called:
c        none
c
c     Called by:
c        CAMx
c
      include 'camx.prm'
      include 'flags.inc'
      include 'vbs.inc'
c
      character*4 arspec(10)
      character*10 arspcname
      real aremis(ncol,nrow,nlay_ems,nspcem)
      real poa_gv_em(ncol,nrow),poa_dv_em(ncol,nrow),
     &     poa_mc_em(ncol,nrow),poa_op_em(ncol,nrow),
     &     poa_bb_em(ncol,nrow)
      integer buffer_offset
      logical lcheck
c
      real argrid(MXCELLS,MXCELLS)
c
c-----Entry point
c
      buffer_offset = buffer_offset_iarem(igrd,idxfile)
c
      if (iounit.eq.0) goto 999
      kount = 1
 100  read(iounit,end=900) idat1,tim1,idat2,tim2 
      if (NINT(1000*tim2) .EQ. 0) then
        tim2 = 24.
        idat2 = idat1
      endif
      ichktm1 = NINT( 1000*(tim1) )
      if( le1day ) then
         ichktm2 = NINT( 1000*(tim2) )
      else
         ichktm2 = NINT( 1000*(tim2)+24000*(idat2-idat1) )
      endif
      if( ichktm2 .EQ. 0 ) ichktm2 = 24000
      ichkems = NINT( 1000*(dtems/60.) )
      if( (ichktm2 - ichktm1) .NE. ichkems ) then
          write(iout_unit,'(//,a)')'ERROR in READAR:'
          write(iout_unit,*) 'Time interval in surface emissions file does'
          write(iout_unit,*)  ' not match emissions update time interval.'
          write(iout_unit,*) '   Beginning Date/Time (Hour): ',idat1,tim1
          write(iout_unit,*) '   Ending Date/Time    (Hour): ',idat2,tim2
          write(iout_unit,*) '   Emiss Input interval(Hour): ',dtems/60.
          call camxerr()
      endif
      tim1 = 100.*aint(tim1) + 60.*amod(tim1,1.)
      tim2 = 100.*aint(tim2) + 60.*amod(tim2,1.)
      if (lvbs .and. LVBSPREPROC) then ! initialize temporary VBS POA emiss array
        poa_gv_em = 0.0
        poa_dv_em = 0.0
        poa_mc_em = 0.0
        poa_op_em = 0.0
        poa_bb_em = 0.0
      endif
      do 50 ll = 1,nspcar
        read(iounit) idum,(arspec(i),i=1,10), ((argrid(i,j),
     &                    i=1+buffer_offset,ncol-buffer_offset),
     &                        j=1+buffer_offset,nrow-buffer_offset)
        write(arspcname,'(10A1)') arspec
C
c-----Check times only if LE1DAY = T, otherwise check both time and date
c
         lcheck = .FALSE.
         if (le1day) then
           if (abs(tim1-time).lt.0.01 .and. tim2.gt.time) lcheck = .TRUE.
         else
           if ((idat1.lt.date .or.
     &         (idat1.eq.date .and. abs(tim1-time).lt.0.01)) .and.
     &         (idat2.gt.date .or.
     &         (idat2.eq.date .and. tim2.gt.time)))
     &       lcheck = .TRUE.
         endif
         if( .NOT. lcheck ) goto 50
c
c-----Put the data into the global array ----
c
        if( lvbs .AND. LVBSPREPROC ) then
           if( arspcname .EQ. 'POA_OP    ' ) then 
              do i=1,ncol
                poa_op_em(i,1) = argrid(i,2)
                poa_op_em(i,nrow) = argrid(i,nrow-1)
              enddo
              do j=1,nrow
                poa_op_em(1,j) = argrid(2,j)
                poa_op_em(ncol,j) = argrid(ncol-1,j)
              enddo
              do j=1+buffer_offset,nrow-buffer_offset
                 do i=1+buffer_offset,ncol-buffer_offset
                   poa_op_em(i,j) = poa_op_em(i,j) + argrid(i,j)/(60.*dtems)
                 enddo
              enddo
           endif
           if( arspcname .EQ. 'POA_GV    ' ) then
              do i=1,ncol
                poa_gv_em(i,1) = argrid(i,2)
                poa_gv_em(i,nrow) = argrid(i,nrow-1)
              enddo
              do j=1,nrow
                poa_gv_em(1,j) = argrid(2,j)
                poa_gv_em(ncol,j) = argrid(ncol-1,j)
              enddo
              do j=1+buffer_offset,nrow-buffer_offset
                 do i=1+buffer_offset,ncol-buffer_offset
                   poa_gv_em(i,j) = poa_gv_em(i,j) + argrid(i,j)/(60.*dtems)
                 enddo
              enddo
           endif
           if( arspcname .EQ. 'POA_DV    ' ) then
              do i=1,ncol
                poa_dv_em(i,1) = argrid(i,2)
                poa_dv_em(i,nrow) = argrid(i,nrow-1)
              enddo
              do j=1,nrow
                poa_dv_em(1,j) = argrid(2,j)
                poa_dv_em(ncol,j) = argrid(ncol-1,j)
              enddo
              do j=1+buffer_offset,nrow-buffer_offset
                 do i=1+buffer_offset,ncol-buffer_offset
                   poa_dv_em(i,j) = poa_dv_em(i,j) + argrid(i,j)/(60.*dtems)
                 enddo
              enddo
           endif
           if( arspcname .EQ. 'POA_MC    ' ) then
              do i=1,ncol
                poa_mc_em(i,1) = argrid(i,2)
                poa_mc_em(i,nrow) = argrid(i,nrow-1)
              enddo
              do j=1,nrow
                poa_mc_em(1,j) = argrid(2,j)
                poa_mc_em(ncol,j) = argrid(ncol-1,j)
              enddo
              do j=1+buffer_offset,nrow-buffer_offset
                 do i=1+buffer_offset,ncol-buffer_offset
                   poa_mc_em(i,j) = poa_mc_em(i,j) + argrid(i,j)/(60.*dtems)
                 enddo
              enddo
           endif
           if( arspcname .EQ. 'POA_BB    ' ) then
              do i=1,ncol
                poa_bb_em(i,1) = argrid(i,2)
                poa_bb_em(i,nrow) = argrid(i,nrow-1)
              enddo
              do j=1,nrow
                poa_bb_em(1,j) = argrid(2,j)
                poa_bb_em(ncol,j) = argrid(ncol-1,j)
              enddo
              do j=1+buffer_offset,nrow-buffer_offset
                 do i=1+buffer_offset,ncol-buffer_offset
                   poa_bb_em(i,j) = poa_bb_em(i,j) + argrid(i,j)/(60.*dtems)
                 enddo
              enddo
           endif
        endif
c
        do ispc=1,nemspc
           if( arspcname .EQ. emspcname(ispc) ) then   
             do i=1,ncol
               aremis(i,1,1,ispc) = argrid(i,2)
               aremis(i,nrow,1,ispc) = argrid(i,nrow-1)
             enddo
             do j=1,nrow
               aremis(1,j,1,ispc) = argrid(2,j)
               aremis(ncol,j,1,ispc) = argrid(ncol-1,j)
             enddo
             do j=1+buffer_offset,nrow-buffer_offset
                 do i=1+buffer_offset,ncol-buffer_offset
                   aremis(i,j,1,ispc) = aremis(i,j,1,ispc) + 
     &                                            argrid(i,j)/(60.*dtems) 
                enddo
             enddo
           endif
        enddo
c
 50   continue
c
c  --- calculate linear combination for VBS species and load it ---
c
      if(lvbs .AND. LVBSPREPROC ) then
         do ibin= 0,NVOLBIN
           do ispc=1,nemspc
              if( emspcname(ispc) .EQ. spname(kpap_c(ibin)) ) then
                 lemmap(kpap_c(ibin)) = ispc
                 aremis(:,:,1,ispc) = aremis(:,:,1,ispc) +
     &                           poa_op_em(:,:) * poa_op_ef(ibin)
     &                          + poa_gv_em(:,:) * poa_gv_ef(ibin)
     &                          + poa_dv_em(:,:) * poa_dv_ef(ibin)
              endif
           enddo
         enddo
         do ibin= 0,NVOLBIN
           do ispc=1,nemspc
              if( emspcname(ispc) .EQ. spname(kpcp_c(ibin)) ) then
                     lemmap(kpcp_c(ibin)) = ispc
                     aremis(:,:,1,ispc) = aremis(:,:,1,ispc) +
     &                           poa_mc_em(:,:) * poa_mc_ef(ibin)
              endif
           enddo
         enddo
         do ibin= 0,NVOLBIN
           do ispc=1,nemspc
              if( emspcname(ispc) .EQ. spname(kpfp_c(ibin)) ) then
                lemmap(kpfp_c(ibin)) = ispc
                     aremis(:,:,1,ispc) = aremis(:,:,1,ispc) +
     &                           poa_bb_em(:,:) * poa_bb_ef(ibin)
              endif
           enddo
         enddo
      endif
      write(iout_unit,'(a40,2(f7.0,i8.5),a,i3)')
     &      'Read area source file at ',tim1,idat1,tim2,idat2,
     &      ' grid',igrd
       call flush(iout_unit)
c
c-----Check times only if LE1DAY = T, otherwise check both time and date
c
      if (le1day) then
        if (abs(tim1-time).lt.0.01 .and. tim2.gt.time) goto 200
        if (tim1-time.ge.0.01) goto 900
      else
        if ((idat1.lt.date .or.
     &      (idat1.eq.date .and. abs(tim1-time).lt.0.01)) .and.
     &      (idat2.gt.date .or.
     &      (idat2.eq.date .and. tim2.gt.time)))
     &    goto 200
      endif
      goto 100
c 
c-----Convert emission rates from moles/(dtems-hours) to moles/s for gas
c     or g/(dtems-hours) to g/s for aero species
c 
 200  continue
      do 10 l = 1,nemspc
        do j = 1,nrow
          do i = 1,ncol
            do k=1,nlay_ems
              if (aremis(i,j,k,l).lt.0.) then
                write(iout_unit,'(//,a)') 'ERROR in READAR:'
                write(iout_unit,'(a,i3)') 'Negative emissions for grid:',igrd
                write(iout_unit,'(a,4i3)') '(i,j,k,l): ',i,j,k,l
                call camxerr()
              endif
            enddo
          enddo 
        enddo
 10   continue 
      goto 999
c
c-----End of file reached; if 1-day emissions requested, rewind and read 
c     through header once more.  Otherwise, report error and exit
c
 900  continue
      if (le1day) then
        if (kount.ge.2) then
          write(iout_unit,'(//,a)') 'ERROR in READAR:'
          write(iout_unit,*)'Cannot match model time with area source time'
          call camxerr()
        endif 
        rewind(iounit)
        read(iounit) idum 
        read(iounit) dum  
        read(iounit) idum  
        read(iounit) idum 
        kount = kount + 1
        goto 100
      else
        write(iout_unit,'(//,a)') 'ERROR in READAR:'
        write(iout_unit,*)'End of area source file reached'
        call camxerr()
      endif
c
 999  continue
c
      return
      end
