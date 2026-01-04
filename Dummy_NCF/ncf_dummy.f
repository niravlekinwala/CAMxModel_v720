      function nf_create( path, cmode, ncid )
      integer nf_create
c
      character*(*) path
      integer       cmode
      integer       ncdid

      write(*,'(//,a)') 'ERROR in NF_CREATE:'
      write(*,'(2A)') 'You have requested NetCDF output, but the ',
     &               'model has not been'
      write(*,'(A)') 'built with NetCDF support.'
      write(*,'(A,//)') 'Please recompile the model with the NCF option.'
      call camxerr()
c
      nf_create = 0
      end
c
      function nf_inq_varid( ncid, name, varid )
      integer nf_inq_varid
c
      integer      ncid
      character(*) name
      integer      varid
c
      write(*,'(//,a)') 'ERROR in NF_INQ_VARID:'
      write(*,'(2A)') 'You have requested NetCDF output, but the ',
     &               'model has not been'
      write(*,'(A)') 'built with NetCDF support.'
      write(*,'(A,//)') 'Please recompile the model with the NCF option.'
      call camxerr()
c
      nf_inq_varid = 0
      end
c
      function nf_put_var_int( ncid, varid, ivals )
      integer nf_put_var_int
c
      integer ncid
      integer varid
      integer ivals(*)
c
      write(*,'(//,a)') 'ERROR in NF_INQ_VAR_INT:'
      write(*,'(2A)') 'You have requested NetCDF output, but the ',
     &               'model has not been'
      write(*,'(A)') 'built with NetCDF support.'
      write(*,'(A,//)') 'Please recompile the model with the NCF option.'
      call camxerr()
c
      nf_put_var_int = 0
      end

      function nf_put_vara_int( ncid, varid, start, count, ivals )
      integer nf_put_vara_int
c
      integer ncid
      integer varid
      integer start(*)
      integer count(*)
      integer ivals(*)
c
      write(*,'(//,a)') 'ERROR in NF_INQ_VARA_INT:'
      write(*,'(2A)') 'You have requested NetCDF output, but the ',
     &               'model has not been'
      write(*,'(A)') 'built with NetCDF support.'
      write(*,'(A,//)') 'Please recompile the model with the NCF option.'
      call camxerr()
c
      nf_put_vara_int = 0
      end

      function nf_put_vara_double( ncid, varid, start, count, dvals )
      integer nf_put_vara_double
c
      integer ncid
      integer varid
      integer start(*)
      integer count(*)
      real*8  dvals(*)
c
      write(*,'(//,a)') 'ERROR in NF_PUT_VARA_DOUBLE:'
      write(*,'(2A)') 'You have requested NetCDF output, but the ',
     &               'model has not been'
      write(*,'(A)') 'built with NetCDF support.'
      write(*,'(A,//)') 'Please recompile the model with the NCF option.'
      call camxerr()
c
       nf_put_vara_double = 0
      end

      function nf_put_vara_real( ncid, varid, start, count, rvals )
      integer nf_put_vara_real
c
      integer ncid
      integer varid
      integer start(*)
      integer count(*)
      real    rvals(*)
c
      write(*,'(//,a)') 'ERROR in NF_PUT_VARA_REAL:'
      write(*,'(2A)') 'You have requested NetCDF output, but the ',
     &               'model has not been'
      write(*,'(A)') 'built with NetCDF support.'
      write(*,'(A,//)') 'Please recompile the model with the NCF option.'
      call camxerr()
c
       nf_put_vara_real = 0
      end

      function nf_def_var( ncid, name, datatype, ndims, dimids, varid )
      integer      nf_def_var
c
      integer      ncid
      character(*) name
      integer      datatype
      integer      ndims
      integer      dimids(*)
      integer      varid
c
      write(*,'(//,a)') 'ERROR in NF_DEF_VAR:'
      write(*,'(2A)') 'You have requested NetCDF output, but the ',
     &               'model has not been'
      write(*,'(A)') 'built with NetCDF support.'
      write(*,'(A,//)') 'Please recompile the model with the NCF option.'
      call camxerr()
c
      nf_def_var = 0
      end
c
      function nf_put_var_real( ncid, varid, rvals )
      integer nf_put_var_real
c
      integer ncid
      integer varid
      real*4  rvals(*)
c
      write(*,'(//,a)') 'ERROR in NF_PUT_VAR_REAL:'
      write(*,'(2A)') 'You have requested NetCDF output, but the ',
     &               'model has not been'
      write(*,'(A)') 'built with NetCDF support.'
      write(*,'(A,//)') 'Please recompile the model with the NCF option.'
      call camxerr()
c
      nf_put_var_real = 0
      end
c
      function nf_put_var_double( ncid, varid, dvals )
      integer nf_put_var_double
c
      integer ncid
      integer varid
      real*8  dvals(*)
c
      write(*,'(//,a)') 'ERROR in NF_PUT_VAR_DOUBLE:'
      write(*,'(2A)') 'You have requested NetCDF output, but the ',
     &               'model has not been'
      write(*,'(A)') 'built with NetCDF support.'
      write(*,'(A,//)') 'Please recompile the model with the NCF option.'
      call camxerr()
c
      nf_put_var_double = 0
      end
c
      function nf_put_att_text( ncid, varid, name, len, text )
      integer nf_put_att_text
c
      integer       ncid
      integer       varid
      character*(*) name
      integer       len
      character*(*) text
c
      write(*,'(//,a)') 'ERROR in NF_PUT_ATT_TEXT:'
      write(*,'(2A)') 'You have requested NetCDF output, but the ',
     &               'model has not been'
      write(*,'(A)') 'built with NetCDF support.'
      write(*,'(A,//)') 'Please recompile the model with the NCF option.'
      call camxerr()
c
      nf_put_att_text = 0
      end
c
      function nf_close( ncid )
      integer nf_close
c
      integer ncid
c
      write(*,'(//,a)') 'ERROR in NF_CLOSE:'
      write(*,'(2A)') 'You have requested NetCDF output, but the ',
     &               'model has not been'
      write(*,'(A)') 'built with NetCDF support.'
      write(*,'(A,//)') 'Please recompile the model with the NCF option.'
      call exit(1)
c
      nf_close = 0
      end
c
      function nf_set_default_format( format, old_format )
      integer nf_set_default_format
c
      integer format
      integer old_format
c
      write(*,'(//,a)') 'ERROR in NF_SET_DEFAULT_FORMAT:'
      write(*,'(2A)') 'You have requested NetCDF output, but the ',
     &               'model has not been'
      write(*,'(A)') 'built with NetCDF support.'
      write(*,'(A,//)') 'Please recompile the model with the NCF option.'
      call camxerr()
c
      nf_set_default_format = 0
      end
c
      function nf_enddef( ncid )
      integer nf_enddef
c
      integer ncid
c
      write(*,'(//,a)') 'ERROR in NF_ENDDEF:'
      write(*,'(2A)') 'You have requested NetCDF output, but the ',
     &               'model has not been'
      write(*,'(A)') 'built with NetCDF support.'
      write(*,'(A,//)') 'Please recompile the model with the NCF option.'
      call camxerr()
c
      nf_enddef = 0
      end
c
      function nf_put_att_real( ncid, varid, name, xtype, len, rvals )
      integer nf_put_att_real
c
      integer       ncid
      integer       varid
      character*(*) name
      integer       xtype
      integer       len
      real          rvals(*)
c
      write(*,'(//,a)') 'ERROR in NF_PUT_ATT_REAL:'
      write(*,'(2A)') 'You have requested NetCDF output, but the ',
     &               'model has not been'
      write(*,'(A)') 'built with NetCDF support.'
      write(*,'(A,//)') 'Please recompile the model with the NCF option.'
      call camxerr()
c
      nf_put_att_real = 0
      end
c
      function nf_put_att_double(ncid,varid,name,xtype,len,dvals)
      integer nf_put_att_double
      integer      ncid
      integer      varid
      character(*) name
      integer      xtype
      integer      len
      real*8       dvals(*)
c
      write(*,'(//,a)') 'ERROR in NF_PUT_ATT_DOUBLE:'
      write(*,'(2A)') 'You have requested NetCDF output, but the ',
     &               'model has not been'
      write(*,'(A)') 'built with NetCDF support.'
      write(*,'(A,//)') 'Please recompile the model with the NCF option.'
      call camxerr()
c
      nf_put_att_double = 0
      end
c
      function nf_put_att_int( ncid, varid, name, xtype, len, ivals )
      integer nf_put_att_int
c
      integer       ncid
      integer       varid
      character*(*) name
      integer       xtype
      integer       len
      integer       ivals(*)
c
      write(*,'(//,a)') 'ERROR in NF_PUT_ATT_INT:'
      write(*,'(2A)') 'You have requested NetCDF output, but the ',
     &               'model has not been'
      write(*,'(A)') 'built with NetCDF support.'
      write(*,'(A,//)') 'Please recompile the model with the NCF option.'
      call camxerr()
c
      nf_put_att_int = 0
      end
c
      function nf_def_dim( ncid, name, len, dimid )
      integer nf_def_dim
c
      integer      ncid
      integer      len
      integer      dimid
c
      write(*,'(//,a)') 'ERROR in NF_DEF_DIM:'
      write(*,'(2A)') 'You have requested NetCDF output, but the ',
     &               'model has not been'
      write(*,'(A)') 'built with NetCDF support.'
      write(*,'(A,//)') 'Please recompile the model with the NCF option.'
      call camxerr()
c
      nf_def_dim = 0
      end
c
      function nf_def_var_chunking(ncid, varid, storage, chunksizes)
c
      integer ncid
      integer varid
      integer storage
      integer chunksizes
c
      write(*,'(//,a)') 'ERROR in NF_DEF_VAR_CHUNKING:'
      write(*,'(2A)') 'You have requested NetCDF output, but the ',
     &               'model has not been'
      write(*,'(A)') 'built with NetCDF support.'
      write(*,'(A,//)') 'Please recompile the model with the NCF option.'
      call camxerr()
c
      nf_def_var_chunking = 0
      end
c
      function nf_def_var_deflate(ncid, varid, shuffle, deflate, deflate_level)
c
      integer ncid
      integer varid
      integer shuffle
      integer deflatE
      integer deflate_level
c
      write(*,'(//,a)') 'ERROR in NF_DEF_VAR_DEFLATE:'
      write(*,'(2A)') 'You have requested NetCDF output, but the ',
     &               'model has not been'
      write(*,'(A)') 'built with NetCDF support.'
      write(*,'(A,//)') 'Please recompile the model with the NCF option.'
      call camxerr()
c
      nf_def_var_deflate = 0
      end
c
      function nf_set_chunk_cache(size, nelems, preemption);
c
      integer size
      integer nelems
      integer preemption
c
      write(*,'(//,a)') 'ERROR in NF_SET_CHUNK_CACHE:'
      write(*,'(2A)') 'You have requested NetCDF output, but the ',
     &               'model has not been'
      write(*,'(A)') 'built with NetCDF support.'
      write(*,'(A,//)') 'Please recompile the model with the NCF option.'
      call camxerr()
c
      nf_set_chunk_cache = 0
      end
c
      function nf_open(path,mode,ncid)
      integer nf_open
      character*(*) path
      integer       mode
      integer       ncid
c
      write(*,'(//,a)') 'ERROR in NF_OPEN:'
      write(*,'(2A)') 'You have requested NetCDF output, but the ',
     &               'model has not been'
      write(*,'(A)') 'built with NetCDF support.'
      write(*,'(A,//)') 'Please recompile the model with the NCF option.'
      call camxerr()
c
      nf_open = 0
      end

      function nf_get_att_int(ncid,varid,name,ivals)
      integer nf_get_att_int
      integer       ncid
      integer       varid
      character*(*) name
      integer       ivals(1)
c
      write(*,'(//,a)') 'ERROR in NF_GET_ATT_INT:'
      write(*,'(2A)') 'You have requested NetCDF output, but the ',
     &               'model has not been'
      write(*,'(A)') 'built with NetCDF support.'
      write(*,'(A,//)') 'Please recompile the model with the NCF option.'
      call camxerr()
c
      nf_get_att_int = 0
      end

      function nf_get_att_double(ncid,varid,name,dvals)
      integer nf_get_att_double
      integer       ncid
      integer       varid
      character*(*) name
      real*8        dvals(1)
c
      write(*,'(//,a)') 'ERROR in NF_GET_ATT_DOUBLE:'
      write(*,'(2A)') 'You have requested NetCDF output, but the ',
     &               'model has not been'
      write(*,'(A)') 'built with NetCDF support.'
      write(*,'(A,//)') 'Please recompile the model with the NCF option.'
      call camxerr()
c
      nf_get_att_double = 0
      end

      function nf_get_vara_int(ncid,varid,start,count,ivals)
      integer nf_get_vara_int
      integer ncid
      integer varid
      integer start(1)
      integer count(1)
      integer ivals(1) 
      write(*,'(//,a)') 'ERROR in NF_GET_VARA_INT:'
      write(*,'(2A)') 'You have requested NetCDF output, but the ',
     &               'model has not been'
      write(*,'(A)') 'built with NetCDF support.'
      write(*,'(A,//)') 'Please recompile the model with the NCF option.'
      call camxerr()
c
      nf_get_vara_int = 0
      end
c
      function nf_get_att_text(ncid,varid,name,text)
      integer nf_get_att_text
      integer       ncid
      integer       varid
      character*(*) name
      character*(*) text
c
      write(*,'(//,a)') 'ERROR in NF_GET_ATT_TEXT:'
      write(*,'(2A)') 'You have requested NetCDF output, but the ',
     &               'model has not been'
      write(*,'(A)') 'built with NetCDF support.'
      write(*,'(A,//)') 'Please recompile the model with the NCF option.'
      call camxerr()
c
      nf_get_att_text = 0
      end
c
      function nf_get_vara_real(ncid,varid,start,count,rvals)
      integer nf_get_vara_real
c
      integer ncid
      integer varid
      integer start
      integer count
      real    rvals
c
      write(*,'(//,a)') 'ERROR in NF_GET_VARA_REAL:'
      write(*,'(2A)') 'You have requested NetCDF output, but the ',
     &               'model has not been'
      write(*,'(A)') 'built with NetCDF support.'
      write(*,'(A,//)') 'Please recompile the model with the NCF option.'
      call camxerr()
c
      nf_get_vara_real = 0
      end

      function nf_inq_dimid(ncid,name,dimid)
      integer nf_inq_dimid
c
      integer       ncid
      character*(*) name
      integer       dimid
c
      write(*,'(//,a)') 'ERROR in NF_INQ_DIMID:'
      write(*,'(2A)') 'You have requested NetCDF output, but the ',
     &               'model has not been'
      write(*,'(A)') 'built with NetCDF support.'
      write(*,'(A,//)') 'Please recompile the model with the NCF option.'
      call camxerr()
c
      nf_inq_dimid = 0
      end

      function nf_get_var_real(ncid,varid,ivals)
      integer nf_get_var_real
c
      integer       ncid
      integer       varid
      real          rvals(1)
c
      write(*,'(//,a)') 'ERROR in NF_GET_VAR_REAL:'
      write(*,'(2A)') 'You have requested NetCDF output, but the ',
     &               'model has not been'
      write(*,'(A)') 'built with NetCDF support.'
      write(*,'(A,//)') 'Please recompile the model with the NCF option.'
      call camxerr()
c
      nf_get_var_real = 0
      end

      function nf_get_var_int(ncid,varid,ivals)
      integer nf_get_var_int
c
      integer       ncid
      integer       varid
      integer       ivals(1)
c
      write(*,'(//,a)') 'ERROR in NF_GET_VAR_INT:'
      write(*,'(2A)') 'You have requested NetCDF output, but the ',
     &               'model has not been'
      write(*,'(A)') 'built with NetCDF support.'
      write(*,'(A,//)') 'Please recompile the model with the NCF option.'
      call camxerr()
c
      nf_get_var_int = 0
      end

      function nf_inq_dimlen(ncid,dimid,len)
      integer nf_inq_dimlen
c
      integer       ncid
      integer       dimid
      integer       len
c
      write(*,'(//,a)') 'ERROR in NF_INQ_DIMLEN:'
      write(*,'(2A)') 'You have requested NetCDF output, but the ',
     &               'model has not been'
      write(*,'(A)') 'built with NetCDF support.'
      write(*,'(A,//)') 'Please recompile the model with the NCF option.'
      call camxerr()
c
      nf_inq_dimlen = 0
      end
