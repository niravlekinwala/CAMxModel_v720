subroutine ncf_set_compress_flag()
!
!----CAMx v7.20 220430
!
!     This routine just sets the compression flag for NetCDF based on the
!     the compiler directive. It needs to be a F90 routine so the preprocssor
!     directives can be used.
!
!      Copyright 1996 - 2022
!     Ramboll
!
!     Modifications:
!        none
!
!     Input arguments:
!        none
!
!     Output arguments:
!        none
!
!     Routines called:
!        none
!
!
!     Called by:
!        READNML
!
use camx_includes
implicit none

#ifndef CHUNK
ncf_compress = .FALSE.
#endif

return
end

