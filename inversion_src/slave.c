
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <gc.h>
#include "../common_src/prototypes.h"

#define QUIT 0

static double *recv_buffer;

/******************************************************************
INPUTS:  (IN)  int my_rank  (this slave node's id)
         (IN)  FILE *log_file  (this node's log_file handle)
RETURN:  none
 *****************************************************************/
void slave(int my_rank, FILE *log_file) {

  int ret;
  MPI_Status status;
  
  fprintf(log_file, "Slave[%d] here, creating buffer size = %d\n",my_rank, NUM_OF_PARAMS);
  
  recv_buffer = (double *)GC_MALLOC((size_t)NUM_OF_PARAMS * sizeof(double));
  if (recv_buffer == NULL) {
    fprintf(log_file, "No room for receive buffer. Exiting/n");
      return;
  }
  for (;;) {
    ret = MPI_Recv( (void *)recv_buffer, NUM_OF_PARAMS, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status );
    if ( recv_buffer[0] == QUIT ) { 
      fprintf(log_file, "\t[%d]received QUIT . . .\n[%d] ",my_rank, ret);
      break;
    }
    ret = minimizing_func(recv_buffer);
  }
  /* free(recv_buffer); */
  /*fprintf(log_file, "Slave exiting.\n");*/
}
