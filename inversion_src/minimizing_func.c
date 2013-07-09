#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <gc.h>
#include "../common_src/prototypes.h"

#define README "parameters.README"

int COL_STEPS;
int PART_STEPS;
int FIT_TEST;
static unsigned int SEED;
/*
double WIND_INTERVAL;
int WIND_LEVELS;
int WIND_DAYS = 1;
int WIND_COLUMNS = 3;
double LITHIC_DENSITY = 2000.0;
double PUMICE_DENSITY = 400.0;
int PLUME_MODEL = 2;

static double variance = 0.0;
*/

/*define the following data structures for the code
  full descriptions are found in common_structures.h */
static ERUPTION e;
static WIND **W;
/*static STATS stats; */

static int num_pts = 0; /*total number of points used in the analysis */
static POINT *p_all;

/* local node varialbles */
static POINT *pt;
static int local_n; /* my number of points */
static int procs;
static int my_rank;
static int my_count; /* byte count of node's POINT array */
static int *displ, *recv_ct; /* pointers to arrays of integers based on total number of nodes */

static FILE *log_file;

/****************************************************************
FUNCTION: test_bounds
DESCRIPTION: This function bounds a value if it is
greater than or less than its preset boundaries.
INPUTS: (IN) int param  (the index of the parameter being tested)
        (IN/OUT) double *try  (the actual value being tested)
        (IN) double bound  (a variable boundary value)
OUTPUT: none
*****************************************************************/
void test_bounds(int param, double *try, double bound) {

#ifdef DEBUG
  if (param == MAX_COL_HT) fprintf(log_file, "\n");
  fprintf(log_file, "Param=%d ", param);
#endif

  if (*try < LO_PARAM(param)) {
    *try = LO_PARAM(param);

  } else if (*try > HI_PARAM(param)) {
    *try = HI_PARAM(param);

  }

}

/****************************************************************
FUNCTION: init_globals
DESCRIPTION: This function read a configuration file
and sets some global variables based on the MODEL global parameter.
INPUTS:  (IN) char *  (complete path to the configuation file)
OUTPUTS: int 1=error, 0=no error
****************************************************************/
int init_globals(char *config_file) {

  FILE *in_config;
  char buf[1][30], **ptr1;
  char line[MAX_LINE];
  char space[4] = "\n\t ";
  char *token;
  int i, WIND_DAYS=1, WIND_LEVELS;

  /* Find out how many processes are being used */
  MPI_Comm_size(MPI_COMM_WORLD, &procs);

  /* Get my process rank, an integer */
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  in_config = fopen(config_file, "r");
  if (in_config == NULL) {
    fprintf(stderr,
	    "[%d-of-%d]\tCannot open configuration file=[%s]:[%s]. Exiting.\n",
	    my_rank, procs, config_file, strerror(errno));
    return 1;
  }

  ptr1 = (char **)&buf[0];
  while (fgets(line, MAX_LINE, in_config) != NULL) {
    /* fprintf(stderr, "%s", line); */
    if (line[0] == '#' || line[0] == '\n' || line[0] == ' ') continue;

    token = strtok_r(line, space, ptr1);
    if (!strncmp(token, "TOLERANCE", strlen("TOLERANCE"))) {
      token = strtok_r(NULL, space, ptr1);
      TOLERANCE = strtod(token, NULL);
      fprintf(log_file, "TOLERANCE = %f\n", TOLERANCE);
      if (!TOLERANCE) return 1;
    }
    else if (!strncmp(token, "SEED", strlen("SEED"))) {
      token = strtok_r(NULL, space, ptr1);
      SEED = (unsigned int)atoi(token);
      fprintf(log_file, "SEED = %u\n", SEED);
    }
    else if (!strncmp(token, "COLUMN_INTEGRATION_STEPS", strlen("COLUMN_INTEGRATION_STEPS"))) {
      token = strtok_r(NULL, space, ptr1);
      COL_STEPS = (int)atoi(token);
      fprintf(log_file, "COLUMN_INTEGRATION_STEPS = %d\n", COL_STEPS);
    }
    else if (!strncmp(token, "PARTICLE_INTEGRATION_STEPS", strlen("PARTICLE_INTEGRATION_STEPS"))) {
      token = strtok_r(NULL, space, ptr1);
      PART_STEPS = (int)atoi(token);
      fprintf(log_file, "PARTICLE_INTEGRATION_STEPS = %d\n", PART_STEPS);
    }
    else if (!strncmp(token, "PLUME_MODEL", strlen("PLUME_MODEL"))) {
      token = strtok_r(NULL,space,ptr1);
      e.plume_model = (int)atoi(token);

      if (e.plume_model == 0) {
        pdf = plume_pdf0;
        fprintf(log_file, "PLUME_MODEL=[%d]%s\n", e.plume_model, "Uniform Distribution with threshold");
      }
      else if (e.plume_model == 1) {
        pdf = plume_pdf1;
        fprintf(log_file, "PLUME_MODEL=[%d]%s\n",e.plume_model, "log-normal Distribution using beta");
      }
      else if (e.plume_model == 2) {
        pdf = plume_pdf2;
        fprintf(log_file, "PLUME_MODEL=[%d]%s\n", e.plume_model, "Beta distribution with alpha and beta parameters");
      }
    }
    else if (!strncmp(token, "FIT_TEST", strlen("FIT_TEST"))) {
      token = strtok_r(NULL, space,ptr1);
      FIT_TEST = (int)atoi(token);
      if (FIT_TEST == CHI2) {
        fit = chi_squared;
        fprintf(log_file, "Goodness-of-fit test=[%d]%s\n", FIT_TEST, "chi-squared test");
      }
      else if (FIT_TEST == RMSE) {
        fit = rmse;
        fprintf(log_file, "Goodness-of-fit test=[%d]%s\n", FIT_TEST, "root-mean-squared-error test");
      }
    }
    else if (!strncmp(token, "DIFFUSION_COEFFICIENT", strlen("DIFFUSION_COEFFICIENT"))) {
      token = strtok_r(NULL,space,ptr1);
       _LO[DIFF_COEF] = strtod(token, NULL);
       if (_LO[DIFF_COEF] < 1.0) _LO[DIFF_COEF] = 1.0;
      token = strtok_r(NULL,space,ptr1);
      _HI[DIFF_COEF] = strtod(token, NULL);
      if (_HI[DIFF_COEF] < 1.0) _HI[DIFF_COEF] = 1.0;
      e.diffusion_coefficient = _LO[DIFF_COEF];
      fprintf(log_file, "DIFFUSION_COEFFICIENT=%.0f to %.0f\n", _LO[DIFF_COEF], _HI[DIFF_COEF]);
    }
    else if (!strncmp(token, "EDDY_CONSTANT", strlen("EDDY_CONSTANT"))) {
      token = strtok_r(NULL,space,ptr1);
      _LO[EDDY_COEF] = strtod(token, NULL);
      token = strtok_r(NULL,space,ptr1);
      _HI[EDDY_COEF] = strtod(token, NULL);
      e.eddy_constant = _HI[EDDY_COEF];
      fprintf(log_file, "EDDY_CONSTANT=%.2f to %.2f\n", _LO[EDDY_COEF], _HI[EDDY_COEF]);
    }
    else if (!strncmp(token, "FALL_TIME_THRESHOLD", strlen("FALL_TIME_THRESHOLD"))) {
      token = strtok_r(NULL,space,ptr1);
      _LO[FALLTIME_THRESH] = strtod(token, NULL);
      token = strtok_r(NULL,space,ptr1);
      _HI[FALLTIME_THRESH] = strtod(token, NULL);
      e.fall_time_threshold = _LO[FALLTIME_THRESH];
      fprintf(log_file, "FALL_TIME_THRESHOLD=%.0f to %.0f\n", _LO[FALLTIME_THRESH], _HI[FALLTIME_THRESH]);
    }
    else if (!strncmp(token, "LITHIC_DENSITY", strlen("LITHIC_DENSITY"))) {
      token = strtok_r(NULL, space, ptr1);
      e.lithic_density = strtod(token, NULL);
      fprintf(log_file, "LITHIC_DENSITY = %.0f\n", e.lithic_density);
    }
    else if (!strncmp(token, "PUMICE_DENSITY", strlen("PUMICE_DENSITY"))) {
      token = strtok_r(NULL, space, ptr1);
      e.pumice_density = strtod(token, NULL);
      fprintf(log_file, "PUMICE_DENSITY = %.0f\n", e.pumice_density);
    }
    else if (!strncmp(token, "VENT_NORTHING", strlen("VENT_NORTHING"))) {
      token = strtok_r(NULL, space, ptr1);
      e.vent_northing = strtod(token, NULL);
      fprintf(log_file, "Volcano Location[northing]: = %.0f\n", e.vent_northing);
    }
    else if (!strncmp(token, "VENT_EASTING", strlen("VENT_EASTING"))) {
      token = strtok_r(NULL, space, ptr1);
      e.vent_easting = strtod(token, NULL);
      fprintf(log_file, "Volcano Location[easting]: = %.0f\n", e.vent_easting);
    }
    else if (!strncmp(token, "PHI_RANGE", strlen("PHI_RANGE"))) {
      token = strtok_r(NULL,space,ptr1);
      e.min_phi = strtod(token, NULL);
      token = strtok_r(NULL,space,ptr1);
      e.max_phi = strtod(token, NULL);
      fprintf(log_file, "MIN_PHI = %g\nMAX_PHI = %g \n", e.min_phi, e.max_phi);
    }
    else if (!strncmp(token, "VENT_ELEVATION", strlen("VENT_ELEVATION"))) {
      token = strtok_r(NULL,space,ptr1);
      e.vent_elevation = strtod(token, NULL);
      fprintf(log_file, "VENT_ELEVATION=%.0f\n", e.vent_elevation);
    }
    else if (!strncmp(token, "MAX_PLUME_ELEVATION", strlen("MAX_PLUME_ELEVATION"))) {
      token = strtok_r(NULL,space,ptr1);
       _LO[MAX_COL_HT] = strtod(token, NULL);
      token = strtok_r(NULL,space,ptr1);
      _HI[MAX_COL_HT] = strtod(token, NULL);
      e.max_plume_elevation = _HI[MAX_COL_HT];
      fprintf(log_file, "MAX_PLUME_ELEVATION=%.0f to %.0f\n", _LO[MAX_COL_HT], _HI[MAX_COL_HT]);
    }
    else if (!strncmp(token, "PLUME_RATIO", strlen("PLUME_RATIO"))) {
      token = strtok_r(NULL,space,ptr1);
       _LO[COL_RATIO] = strtod(token, NULL);
      token = strtok_r(NULL,space,ptr1);
      _HI[COL_RATIO] = strtod(token, NULL);
      e.plume_ratio = _LO[COL_RATIO];
      fprintf(log_file, "PLUME_RATIO=%.2f to %.2f\n", _LO[COL_RATIO], _HI[COL_RATIO]);
    }
    else if (!strncmp(token, "ALPHA", strlen("ALPHA"))) {
       token = strtok_r(NULL,space,ptr1);
      _LO[ALPHAP] = strtod(token, NULL);
      token = strtok_r(NULL,space,ptr1);
      _HI[ALPHAP] = strtod(token, NULL);
      e.alpha = 1.0;
      fprintf(log_file, "ALPHA=%g to %g\n", _LO[ALPHAP], _HI[ALPHAP]);
    }
    else if (!strncmp(token, "BETA", strlen("BETA"))) {
       token = strtok_r(NULL,space,ptr1);
      _LO[BETAP] = strtod(token, NULL);
      token = strtok_r(NULL,space,ptr1);
      _HI[BETAP] = strtod(token, NULL);
      e.beta = 1.0;
      fprintf(log_file, "BETA=%g to %g\n", _LO[BETAP], _HI[BETAP]);
    }
    else if (!strncmp(token, "TOTAL_MASS", strlen("TOTAL_MASS"))) {
       token = strtok_r(NULL,space,ptr1);
      _LO[TOTAL_MASS] = strtod(token, NULL);
      token = strtok_r(NULL,space,ptr1);
      _HI[TOTAL_MASS] = strtod(token, NULL);
      e.total_ash_mass = 0.0;
      fprintf(log_file, "TOTAL_MASS=%g to %g\n", _LO[TOTAL_MASS], _HI[TOTAL_MASS]);
    }
    else if (!strncmp(token, "MEDIAN_PHI", strlen("MEDIAN_PHI"))) {
      token = strtok_r(NULL,space, ptr1);
      _LO[MED_PHI] = strtod(token, NULL);
      token = strtok_r(NULL,space,ptr1);
      _HI[MED_PHI] = strtod(token, NULL);
      e.mean_phi = 0.0;
      fprintf(log_file, "MEDIAN_PHI=%g to %g\n", _LO[MED_PHI], _HI[MED_PHI]);
    }
    else if (!strncmp(token, "SIGMA_PHI", strlen("SIGMA_PHI"))) {
      token = strtok_r(NULL,space,ptr1);
      _LO[STD_DEV_PHI] = strtod(token, NULL);
      /* if (_LO[STD_DEV_PHI] < 1.0) _LO[STD_DEV_PHI] = 1.0; */
      token = strtok_r(NULL,space,ptr1);
      _HI[STD_DEV_PHI] = strtod(token, NULL);
      /* if (_HI[STD_DEV_PHI] < 1.0) _HI[STD_DEV_PHI] = 1.0; */
      e.sigma_phi = 1.0;
      fprintf(log_file, "SIGMA_PHI=%g to %g\n", _LO[STD_DEV_PHI], _HI[STD_DEV_PHI]);
    }
    else if (!strncmp(token, "WIND_SPEED", strlen("WIND_SPEED"))) {
      token = strtok_r(NULL,space,ptr1);
      _LO[WIND_SPEED] = strtod(token, NULL);
      token = strtok_r(NULL,space,ptr1);
      _HI[WIND_SPEED] = strtod(token, NULL);
      fprintf(log_file, "WIND_SPEED=%.0f to %.0f\n", _LO[WIND_SPEED], _HI[WIND_SPEED]);
    }
    else if (!strncmp(token, "WIND_DIRECTION", strlen("WIND_DIRECTION"))) {
      token = strtok_r(NULL,space,ptr1);
      _LO[WIND_DIRECTION] = strtod(token, NULL);
      _LO[WIND_DIRECTION] *= DEG2RAD;
      token = strtok_r(NULL,space,ptr1);
      _HI[WIND_DIRECTION] = strtod(token, NULL);
      _HI[WIND_DIRECTION] *= DEG2RAD;
      fprintf(log_file, "WIND_DIRECTION=%.0f to %.0f\n", _LO[WIND_DIRECTION], _HI[WIND_DIRECTION]);
    }
    else continue;
  }
	fflush(log_file);
	WIND_DAYS = 1;
	WIND_LEVELS = COL_STEPS + 1; 
  /* WIND_LEVELS = (_HI[MAX_COL_HT]/1000) + 1; */
  if (!my_rank) fprintf(stderr, "WIND_LEVELS=%d\n", WIND_LEVELS);
  
  /*WIND_INTERVAL = (_HI[MAX_COL_HT] - e.vent_elevation)/(double)COL_STEPS;*/
 /* WIND_INTERVAL = (_HI[MAX_COL_HT] - e.vent_elevation) /1000; */


  set_global_values(log_file);

  W = (WIND**)GC_MALLOC((size_t)WIND_DAYS * sizeof(WIND));
  if (W == NULL) {
    fprintf(stderr, "[%d-of-%d]\tCannot malloc memory for columns of wind data:[%s]\n",
	    my_rank, procs, strerror(errno));
    return -1;
  }
  else {
    for (i=0; i < WIND_DAYS; i++) {
      W[i] = (WIND *)GC_MALLOC((size_t)WIND_LEVELS * sizeof(WIND));
      if (W[i] == NULL) {
				fprintf(stderr,
				"[%d-of-%d]\tCannot malloc memory for singel row of wind data %d:[%s]\n",
				my_rank, procs, i, strerror(errno));
				return -1;
      }
    }
  }
  /* Set the Wind heights as they remain constant throughtout the program.
   * The heights just refer to the integration steps, plus one more level,
   * The level between the vent height and the ground.
   *
    for (i=0; i < WIND_DAYS; i++) {
        W[i][0].wind_height =  e.vent_elevation;
	    if (!my_rank) fprintf (log_file, "WIND HEIGHTS: %.1f ", W[i][0].wind_height);
        for (j=1; j < WIND_LEVELS; j++) {
            W[i][j].wind_height = W[i][j-1].wind_height + WIND_INTERVAL;
		    if (!my_rank) fprintf (log_file, "%.1f ", W[i][j].wind_height);
	    }
    }
  if (!my_rank) fprintf (stderr, "\n");

  for (i=0; i < WIND_DAYS; i++)
    W[i][0].wind_height =  e.vent_elevation;
    for (j=1; j < WIND_LEVELS; j++)
      W[i][j].wind_height = (double)(j * WIND_INTERVAL);
*/


  NUM_OF_PARAMS = (int) (2 * WIND_DAYS * WIND_LEVELS);
  NUM_OF_PARAMS += (LAST_PARAM - 2);
  fprintf(log_file, "NUM_OF_PARAMS=%d\n", NUM_OF_PARAMS);
  NUM_OF_VERTICES = NUM_OF_PARAMS + 1;

  (void) fclose(in_config);
  return 0;
}



/*****************************************************************
FUNCTION:  get_points
DESCRIPTION:  This function reads easting northing and grain size
distribution  into a POINTS array.

   easting northing mass/kg^2 median sigma)

The total number of points read are divided up between beowulf
nodes so that each node can calulate the grainsize distribution
for its portion of the points read.

INPUTS: (IN) FILE *in  (file handle from which to read)
OUTPUTS: int -1=error, 0=no error
****************************************************************/
int get_points(FILE *in) {

  char line[MAX_LINE]; /*maximum line read */
  int i, ret;

  int extra_pts = 0; /* remaining points to calculate if total does not divide evenly amount nodes */
  int my_start; /* starting line in points file (local) */
  int pts_read = 0; /* number of points read so far (local) */

#ifdef DEBUG
  fprintf(log_file, "ENTER[get_points]\n");
#endif

  while (fgets(line, MAX_LINE, in) != NULL)  {
    if (line[0] == '#' || line[0] == '\n') continue;
    num_pts++;
  }
  rewind(in);

  extra_pts = num_pts % procs;
  if (extra_pts) {
  	local_n = (int)(num_pts/procs);
  	my_start = my_rank*local_n;
	if (my_rank == procs -1) local_n += extra_pts;
  }
  else {
  local_n = (int)(num_pts/procs);
  my_start = my_rank*local_n;
  }
  /*fprintf(stderr, "[%d]my local_n=%d\n", my_rank, local_n); */
  if (!my_rank) {
    displ = (int*)GC_MALLOC((size_t)procs*sizeof(int));
    if (displ == NULL) {
      fprintf(stderr, "[%d-of-%d]\tCannot malloc memory for displ array:[%s]\n",
	      my_rank, procs, strerror(errno));
      return -1;
    }
    recv_ct = (int*)GC_MALLOC((size_t)procs*sizeof(int));
    if (recv_ct == NULL) {
      fprintf(stderr, "[%d-of-%d]\tCannot malloc memory for recv_ct array:[%s]\n",
	      my_rank, procs, strerror(errno));
      return -1;
    }

    displ[0] = 0;
    for (i=0; i<procs; i++){
      recv_ct[i] = local_n*sizeof(POINT);
      if (extra_pts)
      	if (i == procs - 1) recv_ct[i] = (extra_pts+local_n)*sizeof(POINT);

      if (i>0) displ[i] = displ[i-1] + recv_ct[i-1];

    }

    p_all = (POINT*)GC_MALLOC((size_t)num_pts*sizeof(POINT));
    if (p_all == NULL) {
      fprintf(stderr, "[%d-of-%d]\tCannot malloc memory for p_all:[%s]\n",
	      my_rank, procs, strerror(errno));
      return -1;
    }
  }

  my_count = local_n*sizeof(POINT);
  /*fprintf(stderr, "[%d]my count=%d\n", my_rank, my_count); */
  pt = (POINT*)GC_MALLOC((size_t)(local_n)* sizeof(POINT));
  if (pt == NULL) {
    fprintf(stderr, "[%d-of-%d]\tCannot malloc memory for eruptions:[%s]\n",
            my_rank, procs, strerror(errno));
    return -1;
  }

  /* for ( i=0; i < num_pts; i++) { */
  i=0;
  while (i < num_pts) {
    fgets(line, MAX_LINE, in);
    if (line[0] == '#' || line[0] == '\n') continue;
    else {
      while (ret = sscanf(line,
			  "%lf %lf %lf %lf",
			  &(pt+pts_read)->easting,
			  &(pt+pts_read)->northing,
			  &(pt+pts_read)->elevation,
			  &(pt+pts_read)->observed_mass),
	     ret != 4) {

	if (ret == EOF && errno == EINTR) continue;
	fprintf(stderr, "[%d-of-%d]\t[line=%d,ret=%d] Did not read in location data:[%s]\n",
		my_rank, procs, i+1,ret, strerror(errno));
	return -1;
      }
      fprintf(log_file, "%lf %lf %lf %lf\n",
                        (pt+pts_read)->easting,
                        (pt+pts_read)->northing,
                        (pt+pts_read)->elevation,
                        (pt+pts_read)->observed_mass);
      if(i>=my_start){
	pts_read++;
	if (pts_read == local_n) break;
      }
    }
    i++;
  }


  fprintf(log_file,"EXIT[get_points]:[%d-of-%d]Read %d points.\n",
	  my_rank, procs, pts_read);
  fflush(log_file);
  return 0;
}




/***************************************************************
FUNCTION: variance
DESCRIPTION: This function calculates the variance in the
observed mass at each tephra location and sets the variance
value to be used in calculating chi.
INPUTS: none
RETURN : none
**************************************************************/
double set_variance(void) {

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
double chi_squared(void) {

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
    chi += (error*error)/(p_all+i)->calculated_mass;

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
double rmse(void) {

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

/*****************************************************************
FUNCTION: minimizing_func
DESCRIPTION: this is where the nodes assign new parameter values
to the prisms, calculate new magnetic values and compare the
calculated value with the observed value using the chi-squared test.
INPUTS: (IN)  double param[]  (an array of new prism parameters)
RETURN:  double, the result of the chi-squared test
*****************************************************************/
double minimizing_func(double param[]) {

  int i, ret;
  STATS stats;
  double chi_calc;

#ifdef DEBUG
  fprintf(log_file, "  ENTER[minimizing_func]node=%d\n", my_rank);
#endif

  if ( !my_rank ) {
    /* Send the updated parameters to the slave nodes */
    for (i = 0; i < procs; i++) {
      ret = MPI_Send((void *)param, NUM_OF_PARAMS, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
#ifdef DEBUG
      fprintf(log_file, "  \tParameters sent to node %d, MPIret=%d\n", i, ret);
#endif
    }
  }

  /* Every node assisgn the new parameters to their copy of the array of PRISM's */
  assign_new_params( param );
  
#ifdef DEBUG
  fprintf(log_file, "[%d] Calculating new value.\n", my_rank);
#endif

  /* Every node can now calculate a tephra blanket for their subset of POINTs */
  set_eruption_values(&e, W[0]);
  for (i=0; i < local_n; i++)
    tephra_calc(&e, pt+i, W[0], &stats);
#ifdef DEBUG
  fprintf(log_file, "Calculated mass...\n");
  for (i = 0;  i < local_n;  i++) {
    fprintf(log_file, "[%d] [%.0f][%.0f]%g ", i, (pt+i)->easting, (pt+i)->northing, (pt+i)->calculated_mass) ;
  fprintf(log_file, "\n"); 
  }
  fprintf(log_file, "[%d]Finished.\n", my_rank);
  fflush(log_file);
    
#endif

/* fprintf(log_file,"\nMin Particle Fall Time = %gs\nMax Particle Fall Time = %gs\n",stats.min_falltime, stats.max_falltime); */

  /* Gather all of the calculated values from each node into a single POINT array.
     This MPI function orders the bytes of data from each node in rank order (0, 1, etc).

     fprintf(stderr, "[NODE_%d] my_count=%d, my_data=%d bytes p_all=%d bytes\n",
     my_rank, my_count, (int)sizeof recv_ct, (int)sizeof p_all);
  */
  if (ret = MPI_Gatherv(pt, /* address of each proc's BYTE array of data */
			my_count, /* number of BYTES of data from each node */
			MPI_BYTE, /* data_type of each unit of sent of data */
			p_all, /* address of global point array */
			recv_ct, /* address of array of received data */
			displ, /* array of received displacements */
			MPI_BYTE, /* received datatype */
			0, /* root proc */
			MPI_COMM_WORLD), !ret)	{
    if ( !my_rank ) {
      /* Only the master node calculates a new chi_value */
      chi_calc = (*fit)();
      /* fprintf(stderr,"CHI=%f\n", chi_calc); */
    } else
      chi_calc = 0.0;
  }
  else {
    fprintf(stderr, "ERROR: ret=%d\n", ret);
    return 0.0;
  }
#ifdef DEBUG
  fprintf(log_file, "  EXIT[minimizing_func]\t[%d-of-%d] ret=%f\n\n",
	  my_rank, procs, chi_calc);
#endif
  return chi_calc;
}

/******************************************************************
FUNCTION:  assign_new_params
The function assigns updated parameter values to the anomaly being modeled.
This happens before each calulation of the magnetic value.
INPUTS: (IN)  double param[]  (an array of new prism parameters to be tested)
RETURN:  none
*******************************************************************/
void assign_new_params( double param[]) {

  int parm, wind_level;
	double wind_interval;
	
#ifdef DEBUG
  fprintf(log_file, "ENTER[assign_new_params] num_param=%d\n",NUM_OF_PARAMS);
#endif

    e.max_plume_elevation = param[MAX_COL_HT];
    e.alpha = param[ALPHAP];
    e.beta = param[BETAP];
    e.plume_ratio  = param[COL_RATIO];
    e.diffusion_coefficient = param[DIFF_COEF];
    e.total_ash_mass = param[TOTAL_MASS];
    e.mean_phi = param[MED_PHI];
    e.sigma_phi = param[STD_DEV_PHI];
    e.fall_time_threshold = param[FALLTIME_THRESH];
    e.eddy_constant = param[EDDY_COEF];



    #ifdef DEBUG
fprintf(log_file, "%.0f %g %g %.2f %.0f %g %g %g %.0f %.2f\n",
     e.max_plume_elevation,
    e.alpha,
    e.beta,
    e.plume_ratio,
    e.diffusion_coefficient,
    e.total_ash_mass,
    e.median_phi,
    e.sigma_phi,
    e.fall_time_threshold,
    e.eddy_constant);
#endif

wind_interval = (e.max_plume_elevation - e.vent_elevation)/(double)COL_STEPS;

    for (wind_level=0, parm=WIND_SPEED; parm < NUM_OF_PARAMS; parm +=2, wind_level++) {
        W[0][wind_level].wind_height = e.vent_elevation +((double)wind_level * wind_interval); 
        W[0][wind_level].windspeed = param[parm];
        W[0][wind_level].wind_direction = param[parm+1];
        
#ifdef DEBUG
        if (!my_rank)   {
            fprintf(log_file, "[\n%.0f %.0f %.0f]",
            W[0][wind_level].wind_height, 
            W[0][wind_level].windspeed,
            W[0][wind_level].wind_direction/DEG2RAD);
        }
#endif

    }
#ifdef DEBUG
fprintf(log_file, "\nEXIT[assign_new_params].\n");
#endif
}

void set_LOG(FILE *log  ) {
  log_file = log;
}

void close_logfile(void) {
  fclose(log_file);
}

/****************************************************************************
FUNCTION: init_optimal_params
This function initially sets values for all sets of parameters. The values
are randomly chosen between the minimum and maximum values specified for that
parameter.
INPUTS: (IN/OUT) double op[][NUM_OF_PARAMS]  the 2-D array of parameter sets
RETURN:  none
******************************************************************************/
void init_optimal_params(double op[][NUM_OF_PARAMS]) { /* init_optimal_params */

  int vert=0, parm=0;

  fprintf(stderr, "ENTER[init_optimal_params]: NUM_OF_PARAMS=%d \n",
	  NUM_OF_PARAMS);
  srand(SEED);

    for (vert=0; vert < NUM_OF_VERTICES; vert++) { /* for loop */

        /* Maintaining the order of these parameters is critical.
        The first parameter is the column_height parameter of an eruption.
        For each set of possible parameters randomly select an initial column_height value for
        the eruption. This random value should fall within the LO_PARAM - HI_PARAM range (meters).
        */
        op[vert][MAX_COL_HT] =
            (double)LO_PARAM(MAX_COL_HT) +
            ((double)(HI_PARAM(MAX_COL_HT) - LO_PARAM(MAX_COL_HT)) * (double)rand()/(RAND_MAX+1.0));
#ifdef DEBUG2
fprintf(stderr, "  param[%d][%d]=%f ", vert, MAX_COL_HT, op[vert][MAX_COL_HT]);
#endif

        /* The next parameter is the alpha parameter of the beta distribution.
        For each set of possible parameters randomly select an initial alpha value for
        the eruption. This random value should fall within the LO_PARAM - HI_PARAM range (greater than zero).
        */
        op[vert][ALPHAP] =
            (double)LO_PARAM(ALPHAP) +
            ((double)(HI_PARAM(ALPHAP) - LO_PARAM(ALPHAP)) * (double)rand()/(RAND_MAX+1.0));
#ifdef DEBUG2
fprintf(stderr, "  param[%d][%d]=%f ", vert, ALPHAP, op[vert][ALPHAP]);
#endif

        /* The next parameter is the beta parameter of the beta distribution.
        For each set of possible parameters randomly select an initial beta value for
        the eruption. This random value should fall within the LO_PARAM - HI_PARAM range (greater than zero).
        */
        op[vert][BETAP] =
            (double)LO_PARAM(BETAP) +
            ((double)(HI_PARAM(BETAP) - LO_PARAM(BETAP)) * (double)rand()/(RAND_MAX+1.0));
#ifdef DEBUG2
fprintf(stderr, "  param[%d][%d]=%f ", vert, BETAP, op[vert][BETAP]);
#endif

      /* The next parameter is the plume_ratio.
	 For each set of eruption parameters randomly select an initial beta value
	 within the LO_PARAM - HI_PARAM range (no units).
      */
      op[vert][COL_RATIO] =
	(double)LO_PARAM(COL_RATIO) +
	((double)(HI_PARAM(COL_RATIO)-LO_PARAM(COL_RATIO)) * (double)rand()/(RAND_MAX+1.0));
#ifdef DEBUG2
      fprintf(stderr, "  param[%d][%d]=%f ",
			      vert, COL_RATIO, op[vert][COL_RATIO]);
#endif

      /* The next parameter is the diffusion coeffient.
	 For each set of eruption parameters ramdomly select a diffusion coeffient (m^2/sec).
      */
    do {op[vert][DIFF_COEF] =
        (double)LO_PARAM(DIFF_COEF) +
        ((double)(HI_PARAM(DIFF_COEF)-LO_PARAM(DIFF_COEF)) * (double)rand()/(RAND_MAX+1.0));
    } while (  !op[vert][DIFF_COEF]);
#ifdef DEBUG2
      fprintf(stderr, "  param[%d][%d]=%f ", vert, DIFF_COEF, op[vert][DIFF_COEF]);
#endif

      /* The next parameter is the total_mass_erupted.
	 For each set of eruption parameters ramdomly select a total_mass_ejected (kg).
      */
      op[vert][TOTAL_MASS] =
	(double)LO_PARAM(TOTAL_MASS) +
	((double)(HI_PARAM(TOTAL_MASS)-LO_PARAM(TOTAL_MASS)) * (double)rand()/(RAND_MAX+1.0));
#ifdef DEBUG2
      fprintf(stderr, "  param[%d][%d]=%f ",
			      vert, TOTAL_MASS, op[vert][TOTAL_MASS]);
#endif

      /* The next parameter is the median_particle_size in phi units.
	 For each set of eruption parameters ramdomly select a median_phi value).
      */
      op[vert][MED_PHI] =
	(double)LO_PARAM(MED_PHI) +
	((double)(HI_PARAM(MED_PHI)-LO_PARAM(MED_PHI)) * (double)rand()/(RAND_MAX+1.0));
#ifdef DEBUG2
      fprintf(stderr, "  param[%d][%d]=%f ",
	      vert, MED_PHI, op[vert][MED_PHI]);
#endif

      /* The next parameter is the particle size std. deviation (sigma) in phi units.
	 For each set of eruption parameters ramdomly select a sigma_phi value).
      
    do { op[vert][STD_DEV_PHI] =
            (double)LO_PARAM(STD_DEV_PHI) +
            ((double)(HI_PARAM(STD_DEV_PHI)-LO_PARAM(STD_DEV_PHI)) * (double)rand()/(RAND_MAX+1.0)); 
        } while (op[vert][STD_DEV_PHI] < 1.0); */
     op[vert][STD_DEV_PHI] =
            (double)LO_PARAM(STD_DEV_PHI) +
            ((double)(HI_PARAM(STD_DEV_PHI)-LO_PARAM(STD_DEV_PHI)) * (double)rand()/(RAND_MAX+1.0));    
#ifdef DEBUG2
      fprintf(stderr, "  param[%d][%d]=%f ",
	      vert,STD_DEV_PHI , op[vert][STD_DEV_PHI]);
#endif

      /* The next parameter is the fall_time_threshold.
	 For each set of eruption parameters ramdomly select a fall_time_threshold (sec)).
      */
      op[vert][FALLTIME_THRESH] =
	(double)LO_PARAM(FALLTIME_THRESH) +
	((double)(HI_PARAM(FALLTIME_THRESH)-LO_PARAM(FALLTIME_THRESH)) * (double)rand()/(RAND_MAX+1.0));
#ifdef DEBUG2
      fprintf(stderr, "  param[%d][%d]=%f ",
	      vert, FALLTIME_THRESH, op[vert][FALLTIME_THRESH]);
#endif

       /* The next parameter is the eddy constant.
	 For each set of eruption parameters ramdomly select an eddy_constant (sec)).
      */
      op[vert][EDDY_COEF] =
	(double)LO_PARAM(EDDY_COEF) +
	((double)(HI_PARAM(EDDY_COEF)-LO_PARAM(EDDY_COEF)) * (double)rand()/(RAND_MAX+1.0));
#ifdef DEBUG2
      fprintf(stderr, "  param[%d][%d]=%f ",
	      vert, EDDY_COEF, op[vert][EDDY_COEF]);
#endif

      /* The remaining parameters are the wind speed and wind direction for each wind level.
	 Each eruption set will use these wind parameters.
	 initally gets a random value within the LO_PARAM - HI_PARAM
	 range.  This value must not be greater than the value selected for the
	 surface-to-bottom value for the current parameter set.
      */
      for (parm = WIND_SPEED; parm < NUM_OF_PARAMS; parm += 2) { /* for loop */
				 /*rand()/(RAND_MAX+1.0);*/
        op[vert][parm] =LO_PARAM(WIND_SPEED) + ((HI_PARAM(WIND_SPEED) - LO_PARAM(WIND_SPEED)) * rand()/(RAND_MAX+1.0));

/* op[vert][parm] = LO_PARAM(WIND_DIRECTION); */
/*	op[vert][parm] =
	  LO_PARAM(WIND_SPEED) + (HI_PARAM(WIND_SPEED) - LO_PARAM(WIND_SPEED)) * rand()/(RAND_MAX+1.0));
*/
				op[vert][parm+1] =
  			LO_PARAM(WIND_DIRECTION) + ((HI_PARAM(WIND_DIRECTION) - LO_PARAM(WIND_DIRECTION)) * rand()/(RAND_MAX+1.0));

#ifdef DEBUG2
	fprintf(log_file, "  param[%d][%d]=%f param[%d][%d]=%f",
		vert, parm, op[vert][parm], vert, parm+1, op[vert][parm+1]);
#endif
      } /* END for loop */
#ifdef DEBUG2
    fprintf(log_file, "\n");
#endif
/*fprintf(log_file, "\n"); */
    } /* END for loop */


  /*
    for ( param=0; param < NUM_OF_PARAMS; param++)
    printf("%f ", optimal_param[vert][param]);
    }
  */
  fprintf(stderr, "\nEXIT[init_optimal_params].\n");
}

/*************************************************************************
FUNCTION:  printout_points
DESCRIPTION:  This function prints out the northing and easting
coordinates along with the stored calculated magnetic value to the
file "points.out".
INPUTS:  none
OUTPUTS:  none
 ************************************************************************/
void printout_points(void) {

  int i, /* phi_bins, */ bin;
  double val, bin_width;
  FILE *out_pt;
  FILE *out;

  out_pt = fopen(POINTS_OUT, "w");
  if (out_pt == NULL) {
    fprintf(stderr, "Cannot open POINTS file=[%s]:[%s]. Printing to STDOUT.\n",
	    POINTS_OUT, strerror(errno));
    out = stdout;
  } else
    out = out_pt;

  /* phi_bins = (int)(e.max_phi - e.min_phi);
  if (phi_bins < 1) phi_bins = 1; */
  bin_width = (e.max_phi - e.min_phi)/PART_STEPS;
  fprintf(out, "#EAST NORTH ELEV Kg/m^2 ");
  for (bin = 0; bin < PART_STEPS; bin++)
    fprintf(out, "[%g->%g) ",
	   e.min_phi + bin_width*bin, e.min_phi + bin_width*(bin + 1));
  fprintf(out, "(percent)\n");
  
  fprintf(stderr, "PART_STEPS=%d bin_width=%g\n", PART_STEPS, bin_width);

  
  for (i=0; i < num_pts; i++) {
    fprintf(out, "%.0f %.0f %.0f %g ",
	    (p_all+i)->easting,
	    (p_all+i)->northing,
	    (p_all+i)->elevation,
	    (p_all+i)->calculated_mass);
	    /* Print out percent of total*/
    for (bin=0; bin < PART_STEPS; bin++){
	    val = ((p_all+i)->calculated_phi[bin] / (p_all+i)->calculated_mass)*100.0;
      fprintf(out, "%g ", val);
    }
    fprintf(out, "\n");
  }
  if (out == out_pt) fclose(out);
}

/*************************************************************************
FUNCTION:   printout_model
DESCRIPTION:  This function prints out the prism locations and each prism's
set of parameters to the file "prisms.out".
INPUTS:  none
OUTPUTS:  none
 ************************************************************************/
void printout_model(void) {

  int i, j,  WIND_DAYS = 1;
  double wind_levels;
  FILE *out1, *out2;

  out1 = fopen(GRAIN_SZ_MODEL, "w");
  if (out1 == NULL) {
    fprintf(stderr,
	    "Cannot output model to file:[%s]. Printing to STDOUT.\n", strerror(errno));
    out1 = stdout;
  }
  fprintf(out1, 
"Volcano Location:\n\t\
Northing: %.0f\n\t\
Easting: %.0f\n\n\
Total Grainsize Distribution:\n\t\
Max Particle Size %g\n\t\
Min Particle Size: %g\n\t\
Median Diameter: %g\n\t\
Standard Deviation: %g\n\t\
Beta Distribution Param-Alpha: %g\n\t\
Beta Distribution Param-Beta: %g\n\n\
Diffusion Coefficient: %.0f\n\
Eddy Constant: %.2f\n\
Total Mass Erupted: %g\n\
Eruption Plume Height: %.0f\n\
Vent Height: %.0f\n\
Fall Time Threshold: %.0f\n",
        e.vent_northing,
        e.vent_easting,
        e.min_phi,
        e.max_phi,
        e.mean_phi,
        e.sigma_phi,
        /* e.plume_ratio, */
        e.alpha,
        e.beta,
        e.diffusion_coefficient,
       e.eddy_constant,
        e.total_ash_mass,
        e.max_plume_elevation,
        e.vent_elevation,
        e.fall_time_threshold);
    (void) fclose(out1);

 
  wind_levels = COL_STEPS + 1;
  out2 = fopen(WIND_MODEL, "w");
  if (out2 == NULL) {
    fprintf(stderr,
	    "Cannot output model to file:[%s]. Printing to STDOUT.\n", strerror(errno));
    out2 = stdout;
  }
  fprintf(out2, "#HEIGHT\tSPEED\tDIRECTION\n");
  for (i=0; i < WIND_DAYS; i++)
    for (j=0; j < wind_levels; j++) {
            fprintf(out2, "%.1f %.1f %.1f\n",
	      W[i][j].wind_height,
	      W[i][j].windspeed,
	      W[i][j].wind_direction/DEG2RAD);
	      
	  }
  (void) fclose(out2);

}

/*************************************************************************
FUNCTION:   printout_parameters
DESCRIPTION:  This function prints out a time/date stamp and the
parameters used and modified during the modeling process.
INPUTS:  double chi  The best chi value.
OUTPUTS:  none
 ************************************************************************/
void printout_parameters(double chi) {

  FILE *out;
  time_t mytime;
  double rmse;

  out = fopen(README, "w");
  if (out == NULL) {
    fprintf(stderr,
	    "Cannot output parameters to file:[%s]. Printing to STDOUT.\n",
	    strerror(errno));
    out = stdout;
  }
  mytime = time(&mytime);
  rmse = sqrt(chi);
  fprintf(out, "%s\n\
FIT = %.2f\n\
Modeled Values:\n\tMax Column Height: %.0f (m)\n\t\
Alpha Param: %g\n\t\
Beta Param: %g\n\t\
Diffusion Coefficient: %.0f (m^2/s)\n\t\
Fall Time Threshold: %.0f (s)\n\t\
Eddy Constant: %.2f (m^2/s)\n\n\
Total Mass Ejected: %g (kg)\n\n\
Grainsize Distribution\n\tMax Particle Size: %g (phi)\n\t\
Min Particle Size: %g (phi)\n\t\
Median Size: %g (phi)\n\t\
Std. Dev.: %g (phi)\n\n\
Number of Parameters = %d\n\
Number of Vertices = %d\n\
Tolerance = %g\n\n\
Volcano Location:\n\t%.0f (Easting)\n\t\
%.0f (Northing)\n\n\
Parameter Ranges:\n\tMaximum Column Height: %.0f  %.0f (m)\n\t\
Alpha Param: %g  %g\n\t\
Beta Param: %g  %g\n\t\
Diffusion Coefficient: %.0f  %.0f (m^2/s)\n\t\
Eddy Constant: %.2f %.2f (m^2/s)\n\t\
Total Mass Ejected: %g  %g (kg)\n\t\
Median Grain Size: %g  %g (phi)\n\t\
Std. Dev. in particle diameter: %g  %g (phi)\n\t\
Fall Time Threshold: %.0f  %.0f (s)\n\t\
Wind Speed: %.1f  %.1f (m/s)\n\t\
Wind Direction: %.0f  %.0f (+/- degrees from N)\n",
	  asctime(localtime(&mytime)),
	  rmse,
    e.max_plume_elevation,
	 /* e.plume_ratio, */
	  e.alpha,
	  e.beta,
	  e.diffusion_coefficient,
	  e.fall_time_threshold,
     e.eddy_constant, 
	  e.total_ash_mass,
	  e.min_phi,
	  e.max_phi,
	  e.mean_phi,
	  e.sigma_phi,
	  NUM_OF_PARAMS,
	  NUM_OF_VERTICES,
	  TOLERANCE,
	  e.vent_easting, e.vent_northing,
	  LO_PARAM(MAX_COL_HT), HI_PARAM(MAX_COL_HT),
	  LO_PARAM(ALPHAP), HI_PARAM(ALPHAP),
	  LO_PARAM(BETAP), HI_PARAM(BETAP),
	  /*LO_PARAM(COL_RATIO), HI_PARAM(COL_RATIO), */
	  LO_PARAM(DIFF_COEF), HI_PARAM(DIFF_COEF),
	  LO_PARAM(EDDY_COEF), HI_PARAM(EDDY_COEF),
	  LO_PARAM(TOTAL_MASS), HI_PARAM(TOTAL_MASS),
	  LO_PARAM(MED_PHI), HI_PARAM(MED_PHI),
	  LO_PARAM(STD_DEV_PHI), HI_PARAM(STD_DEV_PHI),
	  LO_PARAM(FALLTIME_THRESH), HI_PARAM(FALLTIME_THRESH),
	  LO_PARAM(WIND_SPEED), HI_PARAM(WIND_SPEED),
	  LO_PARAM(WIND_DIRECTION)/DEG2RAD, HI_PARAM(WIND_DIRECTION)/DEG2RAD);

  if (out != stdout) fclose(out);

}

void print_for_stats(double chi) {

  FILE *out;
  double rmse;

  out = fopen("plume.dat", "w");
  if (out == NULL) {
    fprintf(stderr,
	    "Cannot output parameters to file:[%s]. Printing to STDOUT.\n",
	    strerror(errno));
    out = stdout;
  }
  rmse = sqrt(chi);
  /* ME | Max Column Height | Total Mass Ejected | Median Size | Std. Dev. in Distribution | Column Ratio */
  /*fprintf(out, "%.2f %.0f %g %.2f %.2f %.2f\n",
	  rmse, e.max_plume_elevation, e.total_ash_mass, e.median_phi, e.sigma_phi, e.plume_ratio);*/
	  fprintf (out, "%.0f %.0f  %d %g %g", e.max_plume_elevation, e.vent_elevation, (int) COL_STEPS, e.alpha, e.beta);

  if (out != stdout) fclose(out);

}
