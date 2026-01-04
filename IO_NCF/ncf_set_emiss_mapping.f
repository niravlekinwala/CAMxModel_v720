c**** NCF_SET_EMISS_MAPPING
c
      subroutine ncf_set_emiss_mapping(iounit,action)
      use filunit
      use chmstry
      implicit none
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routine checks which species are in the file and adds them to
c   the master emissions species list. It also sets the array for
c   the index of species into this file.
c
c      Copyright 1996 - 2022
c     Ramboll
c      Argument description:
c       Inputs:
c         iounit        I NCF file ID
c         action        C string that describes file being read
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
      include 'vbs.inc'
      include 'netcdf.inc'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      integer       iounit
      character*(*) action
c
c-----------------------------------------------------------------------
c    External functions:
c-----------------------------------------------------------------------
c
      integer istrln
c
c-----------------------------------------------------------------------
c    Local parameters:
c-----------------------------------------------------------------------
c
      integer NUMVBS
      parameter( NUMVBS = 9 )
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      character*20  this_var(NUMVBS)
      integer ispc, ierr, this_varid, i
      logical lfound
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
      this_var(1) = 'IVOG'
      this_var(2) = 'IVOD'
      this_var(3) = 'IVOA'
      this_var(4) = 'IVOB'
      this_var(5) = 'POA_OP'
      this_var(6) = 'POA_GV'
      this_var(7) = 'POA_DV'
      this_var(8) = 'POA_MC'
      this_var(9) = 'POA_BB'
c
c   --- loop ober species in list ----
c
      do ispc=1,nspec
c
c   --- seek this species ---
c
         ierr = nf_inq_varid(iounit,spname(ispc),this_varid)
         if( ierr .EQ. NF_NOERR  ) then
c
c   --- now check if this species is already in master list ---
c
           lfound = .FALSE.
           do i=1,nemspc
             if(spname(ispc) .EQ. emspcname(i) ) then
                lfound = .TRUE.
                lemmap(ispc) = i
             endif
           enddo
           if( .NOT. lfound ) then
             nemspc = nemspc + 1
             emspcname(nemspc) = spname(ispc)
             lemmap(ispc) = nemspc
           endif
         endif
c
c  --- next species ---
c
      enddo
c
c  --- do VBS as a separate case ---
c
      if( lvbs .AND. LVBSPREPROC ) then
         do ispc=1,NUMVBS
            ierr = nf_inq_varid(iounit,this_var(ispc),this_varid)
            if( ierr .NE. NF_NOERR  ) cycle
            lfound = .FALSE.
            do i=1,nemspc
              if(this_var(ispc) .EQ. emspcname(i) ) then
                 lfound = .TRUE.
                 lemmap(ispc) = i
              endif
            enddo
            if( .NOT. lfound ) then
              nemspc = nemspc + 1
              emspcname(nemspc) = this_var(ispc)
              lemmap(ispc) = nemspc
            endif
         enddo
      endif
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
 
