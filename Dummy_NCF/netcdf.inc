      integer NF_NOWRITE
      parameter(NF_NOWRITE = 0)
      integer NF_NOERR
      parameter(NF_NOERR = 0)
      integer NF_FLOAT
      parameter(NF_FLOAT = 5)
      integer NF_REAL
      parameter(NF_REAL = 5)
      integer NF_DOUBLE
      parameter(NF_DOUBLE = 6)
      integer NF_INT
      parameter(NF_INT = 4)
      integer NF_GLOBAL
      parameter (NF_GLOBAL = 0)
      integer NF_CHUNKED
      parameter (NF_CHUNKED = 0)
 
      integer nf_get_vara_int
      external nf_get_vara_int
 
      integer nf_get_att_text
      external nf_get_att_text
 
      integer nf_get_att_int
      external nf_get_att_int
 
      integer nf_get_att_double
      external nf_get_att_double
 
      integer nf_inq_varid
      external nf_inq_varid
 
      integer nf_put_var_real
      external nf_put_var_real
 
      integer nf_put_var_double
      external nf_put_var_double
 
      integer nf_put_var_int
      external nf_put_var_int

      integer nf_put_vara_int
      external nf_put_vara_int
 
      integer nf_put_vara_real
      external nf_put_vara_real
 
      integer nf_get_vara_real
      external nf_get_vara_real
 
      integer nf_put_vara_double
      external nf_put_vara_double
 
      integer nf_def_var
      external nf_def_var
 
      integer nf_put_att_text
      external  nf_put_att_text
 
      integer nf_close
      external nf_close
 
      integer nf_set_default_format
      external nf_set_default_format
 
      integer nf_enddef
      external nf_enddef
 
      integer nf_put_att_real
      external nf_put_att_real
 
      integer nf_put_att_double
      external nf_put_att_double
 
      integer nf_put_att_int
      external nf_put_att_int
 
      integer nf_def_dim
      external nf_def_dim
 
      integer nf_def_var_chunking
      external nf_def_var_chunking
 
      integer nf_def_var_deflate
      external nf_def_var_deflate
 
      integer nf_set_chunk_cache
      external nf_set_chunk_cache
 
      integer nf_open
      external nf_open

      integer nf_inq_dimid
      external nf_inq_dimid

      integer nf_inq_dimlen
      external nf_inq_dimlen

      integer nf_get_var_real
      external nf_get_var_real

      integer nf_get_var_int
      external nf_get_var_int
