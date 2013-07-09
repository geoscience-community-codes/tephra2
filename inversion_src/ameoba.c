
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <gc.h>
#include "../common_src/prototypes.h"

#define TINY 1.0e-10

#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}
/* static double **MODEL_GRID; */
/************************************************************************************
 * INPUTS:
 * double op[][NUM_OF_PARAMS]  :  (in/out) a 2-D array of optimal parameters
 * double mfv[]    :  (in) an array of minimizing function return values
 * double psum[]   :  (in/out) an array of values to be optimized
 * double (*funk)  :  (in) pointer to the minimizing function
 * int worst       :  (in) vertex with the highest value
 * double extrapolation_factor  :  (in)

 * RETURN:  double :  try - an optimized value from the minimizing function 
 ***************************************************************************************/ 

double evaluate(double op[][NUM_OF_PARAMS], double mfv[], double psum[], 
		double (*funk)(double []), 
		int worst, double extrapolation_factor) {
  int param, back=2;
  double fac1,fac2,try,*ptry;
  int tephra_param; 
  
  /* if (DEBUG == 2) fprintf(stderr, "ENTER[evaluate]: worst=%d\n", worst); */
  ptry = (double *)GC_MALLOC((size_t)NUM_OF_PARAMS * sizeof(double));

  fac1 = (1.0 - extrapolation_factor) / NUM_OF_PARAMS;
  fac2 = fac1 - extrapolation_factor;
  
  for (param = 0; param < NUM_OF_PARAMS; param++){
    ptry[param] = psum[param] * fac1 - op[worst][param] * fac2;
    tephra_param = param;
    if (tephra_param >= LAST_PARAM) {
      tephra_param = LAST_PARAM - back;
      back = (back == 2) ? 1: 2;
    }
    test_bounds(tephra_param, &ptry[param], ptry[0]);
  }
  /* smooth_model(&ptry[SURF_TO_BOT]); */
 
  /* Evaluate the function at the trial vertex. */
  try = (*funk)(ptry);

  /* If <try> value is better than the worse, move the worse vertex. */
  if (try < mfv[worst]) {
    /* if (DEBUG) fprintf(stderr, "TRY=%.0f  ",try); */
    mfv[worst] = try; 
    for (param = 0; param < NUM_OF_PARAMS; param++) {
      psum[param] += ptry[param] - op[worst][param];
      op[worst][param] = ptry[param];
    }
  }
  /* if (DEBUG) fprintf(stderr, "EXIT[evaluate]: try=%f mfv[worst]=%f \n",try, mfv[worst]); */
  /* free(ptry); */
  return try;
}


/***********************************************************************************************
 * INPUTS:
 * double op[][NUM_OF_PARAMS]  :  a 2-D array,
 *                                each row represents a vertice of the simplex,
 *                                each column represents a parameter to be optimized (0..n-1)
 *                                (i.e. the coordinates of a vertex in n-dimensional space)
 *
 * double mfv[]  :  array of minimizing function returns for each vertex
 * double tol    :  tolerance
 * double (*funk):  minimizing_function

 * RETURN: none
 ************************************************************************************************/ 
void optimize_params(double op[][NUM_OF_PARAMS], double mfv[], double tol,
		     double (*funk)(double []), int *num_evals) {
  int param, vert;
  int worst; /* vertex with the highest value */
  int better; /* vertex with the next-highest value */
  int best; /* vertex with the lowest value */
  double rtol, sum, swap, save, try, *psum, fit;

  
  fprintf(stderr, "ENTER[optimize_params] ...\n");

  psum=(double *)GC_MALLOC((size_t)NUM_OF_PARAMS*sizeof(double));
  if (psum == NULL) {
    fprintf(stderr, "\t[optimize_params]Cannot malloc memory for psum:[%s]\n",
	    strerror(errno));
    return;
  }
  *num_evals = 0;

  /*
    MODEL_GRID = (double **)malloc(ROWS * sizeof(double));
    if (MODEL_GRID == NULL) {
    fprintf(stderr, "Cannot malloc memory for MODEL_GRID rows:[%s]\n",
    strerror(errno));
    return;
    } 
    else {
    for (i=0; i < ROWS; i++) {
    MODEL_GRID[i] = (double *)malloc(COLS * sizeof(double));
    if (MODEL_GRID[i] == NULL) {
    fprintf(stderr, "Cannot malloc memory for grid row %d:[%s]\n",
    i, strerror(errno));
    return;
    }
    }
    } 
  */

  /* GET PSUM (i.e. sum up each column of parameter values) */
  for (param = 0; param < NUM_OF_PARAMS; param++) {
    for (sum = 0.0, vert = 0; vert < NUM_OF_VERTICES; vert++) 
      sum += op[vert][param];
    psum[param] = sum;

  } 

  for (;;) {
   
    best = 0;
    worst = (mfv[0] > mfv[1]) ? (better = 1,0) : (better = 0,1);
    for (vert=0; vert < NUM_OF_VERTICES; vert++) {
      if (mfv[vert] <= mfv[best]) best = vert;
      if (mfv[vert] > mfv[worst]) {
								better = worst;
								worst = vert;
      } else if (mfv[vert] > mfv[better] && vert != worst) better = vert;
    }
    
    rtol = 2.0 * fabs(mfv[worst] - mfv[best]) / (fabs(mfv[worst]) + fabs(mfv[best]) + TINY);
    if (rtol < tol) {
      SWAP(mfv[0], mfv[best])
	for (param = 0; param < NUM_OF_PARAMS; param++) 
	  SWAP(op[0][param], op[best][param]) 
	    break;
    }
    
    if (*num_evals >= NMAX) {

      fprintf(stderr, "\t[optimize_params]NMAX[%d] exceeded\n",NMAX);
      SWAP(mfv[0], mfv[best])
	for (param = 0; param < NUM_OF_PARAMS; param++) 
	  SWAP(op[0][param], op[best][param]) 
	    break;
    }
    
    *num_evals += 2;
 
    /* Begin a new iteration 
       First extrapolate by a factor of -1. 
    */

    try = evaluate(op, mfv, psum, funk, worst, -1.0);
    
    /* If <try> gives a result better than the best,
       then try an extrapolation by a factor of 2.
    */
    if (try <= mfv[best]) {
      fprintf(stderr, "%d[%d]FIT=%.2f  ", *num_evals, best, try);
      try = evaluate(op, mfv, psum, funk, worst, 2.0);
    
    /* If <try> is worse than the 'better' ,
       look for an intermediate 'better' . */
    } else if (try >= mfv[better]) {
      save = mfv[worst];
      try = evaluate(op, mfv, psum, funk, worst, 0.5);
      
      /* If <try> is still worse than the worst,
	 contract around the best vertex. */
      if (try >= save) {
	
	for (vert = 0; vert < NUM_OF_VERTICES; vert++) {
	  if (vert != best) {
	    for (param = 0; param < NUM_OF_PARAMS; param++)
	      op[vert][param] = psum[param] = 0.5 *(op[vert][param] + op[best][param]);
	    mfv[vert] = (*funk)(psum);
	  }
	  else
	    fprintf(stderr, "%d-CTV[%d] ", *num_evals, vert);
	}
	*num_evals += NUM_OF_PARAMS;
	
	/* GET PSUM (i.e. sum up each column of parameters) */
	for (param = 0; param < NUM_OF_PARAMS; param++) {
	  for (sum = 0.0, vert = 0; vert < NUM_OF_VERTICES; vert++) sum += op[vert][param];
	  psum[param] = sum;
	}
      }
      
    } else --(*num_evals);
    /* if (DEBUG) fprintf(stderr, "\t[optimize_params]NUM_EVAL=%d VERT=%d CHI=%f\n", *num_evals, best, try); */
    /*
      if (!(*num_evals % 1000)) {
      fprintf(stderr, "Printing out data ");
      printout_model();
      printout_points();
      }
    */
  }
  /*free(psum);*/
  fit=mfv[best];
  fprintf(stderr,"\nEXIT[optimize_params]: \nNUM_EVAL=%d FIT=%.2f\n", *num_evals, fit);
}

