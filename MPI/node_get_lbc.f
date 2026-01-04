      subroutine node_get_lbc(iproc_id,nspc,ngr)                          
c
c-----------------------------------------------------------------------
c    Modules used:
c-----------------------------------------------------------------------
c
      use node_mod
      use grid
      use camxfld
      use filunit
c
      implicit none
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c    Argument descriptions:
c     Input:  
c     Output:  
c
c    Called by:
c       EMISTRNS
c    Subroutines called:
c       PAR_WAIT
c       PAR_ASSOC_BUFF
c       PAR_GET_INT
c       PAR_GET_FLOAT
c       EX_LBC4_BUFF
c
c      Copyright 1996 - 2022
c     Ramboll
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
      integer :: iproc_id
      integer :: nspc
      integer :: ngr
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      integer :: itype
      integer :: nm
      integer :: ibytes
      integer :: msgid
      integer :: ihostnum
      integer :: ii1
      integer :: ii2
      integer :: jj1
      integer :: jj2
      integer :: node_src
      integer :: mtc
      integer :: mtp
      integer (kind=8) :: mtp8
      integer :: nptsxy
      character*20 tmp_string
      integer (kind=8) :: maxbytes
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
      maxbytes = DBLE(2**15)*DBLE(2**16) - 1
      itype=1
c
c  --- make sure sends are all finished and de-allocated ---
c
      do nm=1,nmachs
         if (ipaths(1,itype,ngr,nm) .ne. 0) then
            call par_wait(isend_req(nm),ibytes,msgid,ihostnum)
         endif
      enddo
c
c  --- wait on receives ---
c
      do nm=1,nmachs
         if (iget_paths(itype,ngr,nm) .ne. 0) then
            call par_wait(irecv_req(nm),ibytes,msgid,ihostnum)
         endif
      enddo
c
c  --- unpack it into appropriate space. ---
c
      do nm=1,nmachs
         if (iget_paths(itype,ngr,nm) .ne. 0) then
            call par_assoc_buff(node_buffs(nm)%lbc_recv_buff(1),
     &                          node_buffs(nm)%nrecv            ) 
            call par_get_int(ii1,1)
            call par_get_int(ii2,1)
            call par_get_int(jj1,1)
            call par_get_int(jj2,1)
c
            mtp8=DBLE(4)*DBLE(mzp)*
     &           (DBLE(ii2)-DBLE(ii1)+1)*(DBLE(jj2)-DBLE(jj1)+1)*DBLE(nspc)
            if( mtp8 .GT. maxbytes ) then
              write(*,'(//,a)')'ERROR in NODE_GET_LBC:'
              write(*,'(a)')'Message Passing for tracer speices.'
              write(*,'(a)')'Number of bytes exceeds 2,147,483,648 (2^31).'
              write(*,'(a)')'Reduce your Probing Tool application.'
              write(iout,'(/a)')'ERROR in NODE_GET_LBC:'
              write(iout,'(a)')'Message Passing for tracer speices.'
              write(iout,'(a)')'Number of bytes exceeds 2,147,483,648 (2^31).'
              write(iout,'(a)')'Reduce your Probing Tool application.'
              call camxerr()
            endif
c
            call par_get_int(node_src,1)
            nptsxy=(ii2-ii1+1)*(jj2-jj1+1)
            mtp=mzp * nptsxy * nspc
            call par_get_float(scr1(1),mtp)
            call ex_lbc4_buff(mxp,myp,mzp,nspc,conc(iptr4d(ngr)), 
     &                        scr1(1),ii1-i0,ii2-i0,jj1-j0,jj2-j0,mtc)
         endif
      enddo
c
      return
      end
c
c-----------------------------------------------------------------------
c   END subroutine node_get_lbc
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c   BEGIN subroutine ex_lbc4_buff
c-----------------------------------------------------------------------
c
c    Called by:
c       NODE_GET_LBC
c       NODE_GET_LBC_PT
c
      subroutine ex_lbc4_buff(n1,n2,n3,n4,a,b,il,ir,jb,jt,ind)
c
      implicit none
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      integer :: n1
      integer :: n2
      integer :: n3
      integer :: n4
c
      real    :: a(n1,n2,n3,n4)
      real    :: b(*)
c
      integer :: il
      integer :: ir
      integer :: jb
      integer :: jt
      integer :: ind
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      integer :: i
      integer :: j
      integer :: k
      integer :: m
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
      ind=0
c
      do m=1,n4
         do k=1,n3
            do j=jb,jt
               do i=il,ir
                  ind=ind+1
                  a(i,j,k,m)=b(ind)
               enddo
            enddo
         enddo
      enddo
c
      return
      end
c
c-----------------------------------------------------------------------
c   END subroutine ex_lbc4_buff
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c   BEGIN subroutine ex_lbc4_buff_pt
c-----------------------------------------------------------------------
c
      subroutine ex_lbc4_buff_pt(n1,n2,n3,n4,a,b,il,ir,jb,jt,ind)
c
      implicit none
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      integer :: n1
      integer :: n2
      integer :: n3
      integer :: n4
c
      real    :: a(n1,n2,n3,n4)
      real    :: b(*)
c
      integer :: il
      integer :: ir
      integer :: jb
      integer :: jt
      integer :: ind
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      integer :: i
      integer :: j
      integer :: k
      integer :: m
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
      ind=0
c
      do m=1,n4
         do k=1,n3
            do j=jb,jt
               do i=il,ir
                  ind=ind+1
                  a(i,j,k,m)=b(ind)
               enddo
            enddo
         enddo
      enddo
c
      return
      end
c
c-----------------------------------------------------------------------
c   END subroutine ex_lbc4_buff_pt
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c   BEGIN subroutine node_get_lbc_pt
c-----------------------------------------------------------------------
c
      subroutine node_get_lbc_pt(iproc_id,nspc,ngr)                      
c
c-----------------------------------------------------------------------
c    Modules used:
c-----------------------------------------------------------------------
c
      use node_mod
      use grid
      use camxfld
      use tracer
      use filunit
c
      implicit none
c
c----CAMx v7.20 220430
c
c-----------------------------------------------------------------------
c    Description:
c        probing tool
c-----------------------------------------------------------------------
c
c    Argument descriptions:
c     Input:  
c     Output:  
c
c    Called by:
c       EMISTRNS
c    Subroutines called:
c       PAR_WAIT
c       PAR_ASSOC_BUFF
c       PAR_GET_INT
c       PAR_GET_FLOAT
c       EX_LBC4_BUFF
c
c      Copyright 1996 - 2022
c     Ramboll
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      integer :: iproc_id
      integer :: nspc
      integer :: ngr
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      integer :: itype
      integer :: nm
      integer :: ibytes
      integer :: msgid
      integer :: ihostnum
      integer :: ii1
      integer :: ii2
      integer :: jj1
      integer :: jj2
      integer :: node_src
      integer :: mtc
      integer :: mtp
      integer (kind=8) :: mtp8
      integer :: nptsxy
      integer (kind=8) :: maxbytes
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c  --- make sure sends are all finished and de-allocated ---
c
      itype = 1
      do nm=1,nmachs
         if (ipaths(1,itype,ngr,nm) .ne. 0) then
            call par_wait(isend_req_pt(nm),ibytes,msgid,ihostnum)
         endif
      enddo
c
c  --- wait on receives ---
c
      do nm=1,nmachs
         if (iget_paths(itype,ngr,nm) .ne. 0) then
            call par_wait(irecv_req_pt(nm),ibytes,msgid,ihostnum)
         endif
      enddo
c
c  --- unpack it into appropriate space. ---
c
      do nm=1,nmachs
         if (iget_paths(itype,ngr,nm) .ne. 0) then
            call par_assoc_buff(node_buffs(nm)%lbc_pt_recv_buff(1),  !cbwpt
     &                          node_buffs(nm)%nrecv)
            call par_get_int(ii1,1)
            call par_get_int(ii2,1)
            call par_get_int(jj1,1)
            call par_get_int(jj2,1)
            maxbytes = DBLE(2**15)*DBLE(2**16) - 1
            mtp8=DBLE(4)*DBLE(mzp)*
     &           (DBLE(ii2)-DBLE(ii1)+1)*(DBLE(jj2)-DBLE(jj1)+1)*DBLE(nspc)
            if( mtp8 .GT. maxbytes ) then
              write(*,'(//,a)')'ERROR in NODE_GET_LBC_PT:'
              write(*,'(a)')'Message Passing for tracer speices.'
              write(*,'(a)')'Number of bytes exceeds 2,147,483,648 (2^31).'
              write(*,'(a)')'Reduce your Probing Tool application.'
              write(iout,'(//,a)')'ERROR in NODE_GET_LBC:'
              write(iout,'(a)')'Message Passing for tracer speices.'
              write(iout,'(a)')'Number of bytes exceeds 2,147,483,648 (2^31).'
              write(iout,'(a)')'Reduce your Probing Tool application.'
              call camxerr()
            endif
            call par_get_int(node_src,1)
            nptsxy=(ii2-ii1+1)*(jj2-jj1+1)
            mtp=mzp * nptsxy * nspc
            call par_get_float(scr1_pt(1),mtp)
            call ex_lbc4_buff(mxp,myp,mzp,nspc,ptconc(ipsa3d(ngr)), 
     &                     scr1_pt(1),ii1-i0,ii2-i0,jj1-j0,jj2-j0,mtc)
         endif
      enddo
c
      return
      end
