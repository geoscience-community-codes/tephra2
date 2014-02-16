
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <time.h>
#include <mpi.h>
#include "../common_src/prototypes.h"

/*********************************************************************
INPUTS:  none
RETURN:  double (the best CHI value)
 ********************************************************************/
double master(FILE* log_file) {

  int vert, param;
  
  /* A table containing values for parameters:
   * (surf_to_top, south_edge, north_edge, east_edge, west_edge, etc.) 
   * for each vertex of the simplex. 
   * There must be one more vertex than the number of parameters. */
  double optimal_param[NUM_OF_VERTICES][NUM_OF_PARAMS];
  
  /* An array of values returned by the minimizing function, 
     for each vertex of the simplex. */
  double minimizing_func_value[NUM_OF_VERTICES];

  int num_evals; /* the number of function evaluations taken */

  /* the set of parameters we are trying to optimize */
  double param_val[NUM_OF_PARAMS]; 

#ifdef DEBUG
fprintf(log_file, "ENTER[master]\n");
#endif
 
    /* initial parameter guesses : optimal_parameter[vertex][parameter]*/
 
    init_optimal_params(optimal_param); 
    
#ifdef DEBUG2
    for (vert = 0; vert < ( NUM_OF_PARAMS + 1); vert++)
      for ( param=0; param < NUM_OF_PARAMS; param++) 
	fprintf(log_file, "\t[%d][%d]: %f ", vert, param, optimal_param[vert][param]);
    fprintf(log_file, "\n");
#endif
 /* fprintf(log_file, "%d Params initialized ....\n", NUM_OF_PARAMS);   */
    /* the number of vertices equals one more than the number of parameters being optimized */
    for (vert=0; vert < NUM_OF_VERTICES; vert++) { 
      
      for (param=0; param < NUM_OF_PARAMS; param++) 
							param_val[param] = optimal_param[vert][param];
      
      minimizing_func_value[vert] = minimizing_func( param_val );
      
#ifdef DEBUG
      fprintf(log_file, "[%d]chi=%f\n", vert, minimizing_func_value[vert]);
#endif
    }

    /* the dimension of the simplex equals the number of parameters being optimized */
    /* fprintf(stderr, "CHI TOLERANCE = %e\n", (double)TOLERANCE); */
    fprintf(log_file, "Master optimizing params .....\n");
    optimize_params(optimal_param, 
		    minimizing_func_value,  
		    TOLERANCE, 
		    minimizing_func, 
		    &num_evals);
    /*
    for ( vert=0; vert < NUM_OF_VERTICES; vert++ ) {
      
      fprintf(stderr,"[%d]chi=%f\n", vert, minimizing_func_value[vert]);
    }
    rmse = sqrt(minimizing_func_value[0]);
    fprintf(stderr, "MEAN Squared ERROR = %.2f\n", rmse); 
    for ( param=0; param < NUM_OF_PARAMS; param++) 
	  fprintf(stderr, "\tPrism[%d]: %f\n", param, optimal_param[i][param]); */
    
    printout_points();

    for (param=0; param < NUM_OF_PARAMS; param++)
      param_val[param] =  optimal_param[0][param];

    assign_new_params( param_val );
    printout_model(param_val);
#ifdef DEBUG
    fprintf(log_file, "EXIT[master]\n");
#endif
    return  minimizing_func_value[0];
}










