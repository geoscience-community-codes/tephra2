#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include "../common_src/prototypes.h"
/***************************************************************
FUNCTION: variance
DESCRIPTION: This function calculates the variance in the
observed mass at each tephra location and sets the variance
value to be used in calculating chi.
INPUTS: none
RETURN : none
**************************************************************/
double set_variance(FILE *log_file, int num_pts, POINT *p_all) {

  double sum = 0.0;
  int j;
  /* because this varialbe is static, the declaration value is set only once */
  static double variance=0.0;
  double average=0.0;
  double ep=0.0;

#ifdef DEBUG
  fprintf(log_file, "    ENTER[variance] num_pts=%d \n", num_pts);
#endif

if (variance) return variance;

  for (j = 0; j < num_pts; j++)
    sum += (p_all+j)->observed_mass;
  average = sum/(double) num_pts;


  for (j = 0; j < num_pts; j++) {
    ep = (p_all+j)->observed_mass - average;
    sum += ep*ep;

  }
  variance = sum/((double)num_pts-1.0);
  fprintf(stderr, "SUM=%f VARIANCE = %f\n", sum, variance);

#ifdef DEBUG
  fprintf(log_file,
	  "    EXIT[variance]\tAVG=%f\tVAR=%f\tSTD DEV=%f\n",
	  average, variance, sqrt(variance));
#endif
return variance;
}

/********************************************************************************
FUNCTION: chi_squared
DESCRIPTION: The chi-squared test is used as the goodness-of-fit test.
INPUTS: none
OUTPUT:  double, the chi-squared value,
         based on the calulated mass and the observed mass, at each location.
********************************************************************************/
double chi_squared(FILE *log_file, int num_pts, POINT *p_all) {

  int i;
  double chi=0.0, error;

#ifdef DEBUG
 fprintf(log_file,"   ENTER[chi_squared] ...\n");
#endif
 
  for (i=0; i < num_pts; i++) {
#ifdef DEBUG
    fprintf (stderr, "[%d]%f %f/ ", i, (p_all+i)->calculated_mass, (p_all+i)->observed_mass);
#endif

   /* error = log( (p_all+i)->calculated_mass ) - log( (p_all+i)->observed_mass);  */
    error = (p_all+i)->calculated_mass - (p_all+i)->observed_mass; 
    chi += (error*error)/(p_all+i)->observed_mass;

  }

#ifdef DEBUG
  fprintf(stderr, "\n");
#endif

#ifdef DEBUG
  fprintf(log_file,"   EXIT[chi-squared] [ret=%f]\n\n", chi);
#endif

  return chi;
}

/********************************************************************************
FUNCTION: rmse
DESCRIPTION: The root mean squared error test is used as the goodness-of-fit test.
INPUTS: none
OUTPUT:  double, the goodness-of-fit value

         Based on the difference between the calulated mass 
         and the observed mass, at each location.
********************************************************************************/
double rmse(FILE *log_file, int num_pts, POINT *p_all) {

  int i;
  double rmse=0.0, error;

#ifdef DEBUG
  fprintf(log_file,"   ENTER[rmse] ...\n");
#endif

  for (i=0; i < num_pts; i++) {
#ifdef DEBUG
    fprintf (stderr, "[%d]%f %f/ ", i, (p_all+i)->calculated_mass, (p_all+i)->observed_mass);
#endif

    error = (p_all+i)->calculated_mass - (p_all+i)->observed_mass; 
    rmse += (error*error);
  }
  rmse /= (num_pts);
  rmse = sqrt(rmse);

#ifdef DEBUG
  fprintf(stderr, "\n");
#endif

#ifdef DEBUG
  fprintf(log_file,"   EXIT[rmse] [ret=%f]\n\n", rmse);
#endif

  return rmse;
}

/********************************************************************************
FUNCTION: log_test
DESCRIPTION: The Tokyo log test is used as the goodness-of-fit test.
INPUTS: none
OUTPUT:  double, the fit value,
         based on the calulated mass and the observed mass, at each location.
         sum[0..i] ([log_base10(calculated_mass[i]/observed_mass[i])]^2)
********************************************************************************/
double log_test(FILE *log_file, int num_pts, POINT *p_all) {

	int i;
  double fit = 0.0, ratio, log_ratio;

#ifdef DEBUG
  fprintf(log_file,"   ENTER[log_test] ...\n");
#endif

  for (i=0; i < num_pts; i++) {
  
#ifdef DEBUG
    fprintf (stderr, "[%d]%f %f/ ", i, (p_all+i)->calculated_mass, (p_all+i)->observed_mass);
#endif 
		ratio = (p_all+i)->calculated_mass / (p_all+i)->observed_mass;
		if (ratio <=0) log_ratio = 0;
		else log_ratio = log10(ratio);
		fit += (log_ratio * log_ratio);
	}
	
#ifdef DEBUG
  fprintf(stderr, "\n");
#endif

#ifdef DEBUG
  fprintf(log_file,"   EXIT[log_test] [ret=%f]\n\n", fit);
#endif

	return fit;
}
