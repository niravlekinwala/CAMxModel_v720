#include "util.h"

extern int par_retrieve_tag(MPI_Request msgtag);
extern void zero_tag_array();

void par_wait_(int *fmsgtag,int *ibytes,int *msgtype,int *ihostnum)
{
   MPI_Request msgtag;
   MPI_Status status;
   int ierr=0;

   if (flag_msgtag == 0) zero_tag_array();
   
   msgtag=par_retrieve_tag(*fmsgtag);
   
   ierr=MPI_Wait(&msgtag,&status);
  
   if(ierr < 0) 
      printf("Error in par_wait\n");
    
   MPI_Get_count(&status,MPI_PACKED,ibytes);
   *msgtype=status.MPI_TAG;
   *ihostnum=status.MPI_SOURCE;
   
}
