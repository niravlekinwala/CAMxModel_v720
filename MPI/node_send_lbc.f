      subroutine node_send_lbc(iproc_id,nspc,ngr)                                   
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
c       PAR_GET_NOBLOCK
c       PAR_INIT_PUT
c       PAR_PUT_INT
c       MK_LBC4_BUFF
c       PAR_PUT_FLOAT
c       PAR_SEND_NOBLOCK
c
c      Copyright 1996 - 2022
c     Ramboll
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c    Include files:
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
      integer :: ii1
      integer :: ii2
      integer :: jj1
      integer :: jj2
      integer :: mtp
      integer (kind=8) :: mtp8
      integer (kind=8) :: maxbytes
      integer :: ierr
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
      maxbytes = DBLE(2**15)*DBLE(2**16) - 1
c
      itype=1
c  --- before send anything, post the receives. ---
c
      do nm=1,nmachs
         if (iget_paths(itype,ngr,nm) .ne. 0) then
            mtp8 = DBLE(node_buffs(nm)%nrecv)*DBLE(4)
            if( mtp8 .GT. maxbytes ) then
               write(*,'(//,a)')'ERROR in NODE_SEND_LBC:'
               write(*,'(a)')'Message Passing for tracer speices.'
               write(*,'(a)')'Number of bytes exceeds 2,147,483,648 (2^31).'
               write(*,'(2a)')'Reduce your Probing Tool application or',
     &                        ' use a different number of cores for MPI.'
               write(iout,'(//,a)')'ERROR in NODE_SEND_LBC:'
               write(iout,'(a)')'Message Passing for tracer speices.'
               write(iout,'(a)')'Number of bytes exceeds 2,147,483,648 (2^31).'
               write(iout,'(2a)')'Reduce your Probing Tool application or',
     &                        ' use a different number of cores for MPI.'
               call camxerr()
            endif
            call par_get_noblock(node_buffs(nm)%lbc_recv_buff(1),
     &                           node_buffs(nm)%nrecv,20000+ngr,machs(nm),
     &                           irecv_req(nm)                            )
         endif
      enddo
c     
c  --- sending the stuff ---
c
      do nm=1,nmachs
         isend_req(nm) = nm
         if (ipaths(1,itype,ngr,nm) .ne. 0) then
            ii1=ipaths(1,itype,ngr,nm)
            ii2=ipaths(2,itype,ngr,nm)
            jj1=ipaths(3,itype,ngr,nm)
            jj2=ipaths(4,itype,ngr,nm)
            mtp8 = DBLE(4)* DBLE(mzp)*DBLE(nspc)*DBLE(ii2-ii1+1)*DBLE(jj2-jj1+1)
            if( mtp8 .GT. maxbytes ) then
               write(*,'(//,a)')'ERROR in NODE_SEND_LBC:'
               write(*,'(a)')'Message Passing for tracer speices.'
               write(*,'(a)')'Number of bytes exceeds 2,147,483,648 (2^31).'
               write(*,'(2a)')'Reduce your Probing Tool application or',
     &                        ' use a different number of cores for MPI.'
               write(iout,'(//,a)')'ERROR in NODE_SEND_LBC:'
               write(iout,'(a)')'Message Passing for tracer speices.'
               write(iout,'(a)')'Number of bytes exceeds 2,147,483,648 (2^31).'
               write(iout,'(2a)')'Reduce your Probing Tool application or',
     &                        ' use a different number of cores for MPI.'
               call camxerr()
            endif
            call par_init_put(node_buffs(nm)%lbc_send_buff(1),
     &                        node_buffs(nm)%nsend            )
            call par_put_int(ii1,1)
            call par_put_int(ii2,1)
            call par_put_int(jj1,1)
            call par_put_int(jj2,1)
            call par_put_int(mynum,1)
            call mk_lbc4_buff(mxp,myp,mzp,nspc,conc(iptr4d(ngr)),
     &                        scr1(1),ii1-i0,ii2-i0,jj1-j0,jj2-j0,mtp)
            call par_put_float(scr1(1),mtp)  
            call par_send_noblock(ipaths(5,itype,ngr,nm),
     &                            20000+ngr,isend_req(nm))
         endif
      enddo
c
      return
      end
c
c-----------------------------------------------------------------------
c    END subroutine node_send_lbc:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c    BEGIN subroutine mk_lbc4_buff:
c-----------------------------------------------------------------------
c
c    Called by:
c       NODE_SEND_LBC
c       NODE_SEND_LBC_PT
c
      subroutine mk_lbc4_buff(n1,n2,n3,n4,a,b,il,ir,jb,jt,ind)
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
                  b(ind)=a(i,j,k,m)
               enddo
            enddo
         enddo
      enddo
c
      return
      end
c
c-----------------------------------------------------------------------
c    END subroutine mk_lbc4_buff:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c    BEGIN subroutine mk_lbc4_buff_pt:
c-----------------------------------------------------------------------
c
      subroutine mk_lbc4_buff_pt(n1,n2,n3,n4,a,b,il,ir,jb,jt,ind)
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
                  b(ind)=a(i,j,k,m)
               enddo
            enddo
         enddo
      enddo
c
      return
      end
c
c-----------------------------------------------------------------------
c    END subroutine mk_lbc4_buff_pt:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c    BEGIN subroutine node_send_lbc_pt:
c-----------------------------------------------------------------------
c
      subroutine node_send_lbc_pt(iproc_id,nspc,ngr)    
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
c-----------------------------------------------------------------------
c
c    Argument descriptions:
c     Input:  
c     Output:  
c
c    Called by:
c       EMISTRNS
c    Subroutines called:
c       PAR_GET_NOBLOCK
c       PAR_INIT_PUT
c       PAR_PUT_INT
c       MK_LBC4_BUFF
c       PAR_PUT_FLOAT
c       PAR_SEND_NOBLOCK
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
      integer :: ii1
      integer :: ii2
      integer :: jj1
      integer :: jj2
      integer :: mtp
      integer (kind=8) :: mtp8
      integer (kind=8) :: maxbytes
      integer :: ierr
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
      maxbytes = DBLE(2**15)*DBLE(2**16) - 1
      itype=1
c      
c  --- before send anything, post the receives. ---
c
      do nm=1,nmachs
         if (iget_paths(itype,ngr,nm) .ne. 0) then
            mtp8 = DBLE(node_buffs(nm)%nrecv)*DBLE(4)
            if( mtp8 .GT. maxbytes ) then
               write(*,'(//,a)')'ERROR in NODE_SEND_LBC_PT:'
               write(*,'(a)')'Message Passing for tracer speices.'
               write(*,'(a)')'Number of bytes exceeds 2,147,483,648 (2^31).'
               write(*,'(2a)')'Reduce your Probing Tool application or',
     &                        ' use a different number of cores for MPI.'
               write(iout,'(//,a)')'ERROR in NODE_SEND_LBC_PT:'
               write(iout,'(a)')'Message Passing for tracer speices.'
               write(iout,'(a)')'Number of bytes exceeds 2,147,483,648 (2^31).'
               write(iout,'(2a)')'Reduce your Probing Tool application or',
     &                        ' use a different number of cores for MPI.'
               call camxerr()
            endif
            call par_get_noblock(node_buffs(nm)%lbc_pt_recv_buff(1),  !cbwpt
     &                           node_buffs(nm)%nrecv,20000+ngr+pt_identifier,
     &                           machs(nm),irecv_req_pt(nm)                   )
         endif
      enddo
c     
c  --- sending the stuff ---
c
      do nm=1,nmachs
         isend_req_pt(nm) = 900+nm
         if (ipaths(1,itype,ngr,nm) .ne. 0) then
            ii1=ipaths(1,itype,ngr,nm)
            ii2=ipaths(2,itype,ngr,nm)
            jj1=ipaths(3,itype,ngr,nm)
            jj2=ipaths(4,itype,ngr,nm)
            call par_init_put(node_buffs(nm)%lbc_pt_send_buff(1), !cbwpt
     &                        node_buffs(nm)%nsend               )
            call par_put_int(ii1,  1)
            call par_put_int(ii2,  1)
            call par_put_int(jj1,  1)
            call par_put_int(jj2,  1)
            call par_put_int(mynum,1)
            mtp8 = DBLE(4)* DBLE(mzp)*DBLE(nspc)*DBLE(ii2-ii1+1)*DBLE(jj2-jj1+1)
            if( mtp8 .GT. maxbytes ) then
               write(*,'(//,a)')'ERROR in NODE_SEND_LBC_PT:'
               write(*,'(a)')'Message Passing for tracer speices.'
               write(*,'(a)')'Number of bytes exceeds 2,147,483,648 (2^31).'
               write(*,'(2a)')'Reduce your Probing Tool application or',
     &                        ' use a different number of cores for MPI.'
               write(iout,'(//,a)')'ERROR in NODE_SEND_LBC_PT:'
               write(iout,'(a)')'Message Passing for tracer speices.'
               write(iout,'(a)')'Number of bytes exceeds 2,147,483,648 (2^31).'
               write(iout,'(2a)')'Reduce your Probing Tool application or',
     &                        ' use a different number of cores for MPI.'
               call camxerr()
            endif
            call mk_lbc4_buff_pt(mxp,myp,mzp,nspc,ptconc(ipsa3d(ngr)),
     &                           scr1_pt(1),ii1-i0,ii2-i0,jj1-j0,jj2-j0,mtp)
            call par_put_float(scr1_pt(1),mtp) 
            call par_send_noblock(ipaths(5,itype,ngr,nm),
     &                            20000+ngr+pt_identifier,isend_req_pt(nm))
         endif
      enddo
c
      return
      end
