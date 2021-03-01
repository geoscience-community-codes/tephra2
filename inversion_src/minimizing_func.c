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
int WIND_LEVELS;
int INV_WIND_LEVELS;
static unsigned int SEED;

/*define the following data structures for the code
  full descriptions are found in common_structures.h */
static ERUPTION e;
static WIND *W;
static WIND *FW;

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
//double (*fit)(FILE *, int, POINT *);

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

  if (!FIXED_WIND) { 
  	if (*try < LO_PARAM(param)) *try = LO_PARAM(param);
  	else if (*try > HI_PARAM(param)) *try = HI_PARAM(param);
  }
  else {
    if (param < WIND_SPEED) {
    	if (*try < LO_PARAM(param)) *try = LO_PARAM(param);
  		else if (*try > HI_PARAM(param)) *try = HI_PARAM(param);
  	}
  }

 /*if (*try < LO_PARAM(param)) *try = LO_PARAM(param);
 else if (*try > HI_PARAM(param)) *try = HI_PARAM(param); 
*/
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
  int WIND_LEVELS;

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
       e.pdf = &plume_pdf0;
        fprintf(log_file, "PLUME_MODEL=[%d]%s\n", e.plume_model, "Uniform Distribution with threshold");
      }
      else if (e.plume_model == 1) {
        e.pdf = &plume_pdf1;
        fprintf(log_file, "PLUME_MODEL=[%d]%s\n",e.plume_model, "log-normal Distribution using beta");
      }
      else if (e.plume_model == 2) {
       e.pdf = &plume_pdf2;
        fprintf(log_file, "PLUME_MODEL=[%d]%s\n", e.plume_model, "Beta distribution with alpha and beta parameters");
      }
    }
    else if (!strncmp(token, "FIT_TEST", strlen("FIT_TEST"))) {
      token = strtok_r(NULL, space,ptr1);
      FIT_TEST = (int)atoi(token);
      if (FIT_TEST == CHI2) {
       e.fit = &chi_squared;
        fprintf(log_file, "Goodness-of-fit test=[%d]%s\n", FIT_TEST, "chi-squared test");
      }
      else if (FIT_TEST == RMSE) {
        e.fit = &rmse;
        fprintf(log_file, "Goodness-of-fit test=[%d]%s\n", FIT_TEST, "root-mean-squared-error test");
      }
      else if (FIT_TEST == LOG_TEST) {
      	e.fit = &log_test;
      	fprintf(log_file, "Goodness-of-fit test=[%d]%s\n", FIT_TEST, "Tokyo log test");
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
      token = strtok_r(NULL,space,ptr1);
      _HI[WIND_DIRECTION] = strtod(token, NULL);
      _HI[WIND_DIRECTION] *= DEG2RAD;
      _LO[WIND_DIRECTION] *= DEG2RAD;
      fprintf(log_file, "WIND_DIRECTION=%.0f to %.0f\n", _LO[WIND_DIRECTION], _HI[WIND_DIRECTION]);
    }
    else if (!strncmp(token, "WIND_LEVELS", strlen("WIND_LEVELS"))) {
      token = strtok_r(NULL, space, ptr1);
      INV_WIND_LEVELS = (int)atoi(token);
      fprintf(log_file, "INV_WIND_LEVELS = %d\n", INV_WIND_LEVELS);
      NUM_OF_PARAMS = (int) (2 * INV_WIND_LEVELS); 
      NUM_OF_PARAMS += (LAST_PARAM - 2);
      fprintf(log_file, "NUM_OF_PARAMS=%d\n", NUM_OF_PARAMS);
      NUM_OF_VERTICES = NUM_OF_PARAMS + 1;
    }
    else if (!strncmp(token, "FIXED_WIND", strlen("FIXED_WIND"))) {
      token = strtok_r(NULL, space, ptr1);
      FIXED_WIND = (int)atoi(token);
    }
    else continue;
  }
	fflush(log_file);
	WIND_LEVELS = COL_STEPS + 1; 
	fprintf(log_file, "FIXED_WIND = %s\n", (FIXED_WIND) ? "TRUE" : "FALSE");
  	if (!my_rank) fprintf(stderr, "WIND_LEVELS=%d\n", WIND_LEVELS);

  set_global_values(log_file);
  
  W = (WIND*)GC_MALLOC((size_t)WIND_LEVELS*sizeof(WIND));
  if (W == NULL) {
      fprintf(stderr, "[%d-of-%d]\tCannot malloc memory for W (wind levels):[%s]\n",
	      my_rank, procs, strerror(errno));
      return ERROR;
  }
  (void) fclose(in_config);
  return FALSE;
}



/*****************************************************************
FUNCTION:  get_points
DESCRIPTION:  This function opens the grid file and reads
easting northing elevation mass into a POINTS array.

   easting northing mass/kg^2)

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
      return ERROR;
    }
    recv_ct = (int*)GC_MALLOC((size_t)procs*sizeof(int));
    if (recv_ct == NULL) {
      fprintf(stderr, "[%d-of-%d]\tCannot malloc memory for recv_ct array:[%s]\n",
	      my_rank, procs, strerror(errno));
      return ERROR;
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
      return ERROR;
    }
  }

  my_count = local_n*sizeof(POINT);
  /*fprintf(stderr, "[%d]my count=%d\n", my_rank, my_count); */
  pt = (POINT*)GC_MALLOC((size_t)(local_n)* sizeof(POINT));
  if (pt == NULL) {
    fprintf(stderr, "[%d-of-%d]\tCannot malloc memory for eruptions:[%s]\n",
            my_rank, procs, strerror(errno));
    return ERROR;
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
	return ERROR;
      }
      if (!(pt+pts_read)->observed_mass) (pt+pts_read)->observed_mass = .000000001;
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
  return FALSE;
}

/**************************************************************
FUNCTION:  get_wind
DESCRIPTION:  This function reads wind data into the
FW WIND array. Each node stores all of the wind data. It is only called
if FIXED_WIND = 1 (TRUE)

INPUTS: (IN) FILE *in  (file handle from which to read)
OUTPUTS: int -1=error, 0=no error
***************************************************************/
int get_wind(FILE *in) {

  int i=0, levels=0, ret;
  char line[MAX_LINE];
  
  fprintf(log_file,"ENTER[get_wind].\n");
  
  /* Count the number of wind levels */
  while (fgets(line, MAX_LINE, in) != NULL)  {
    if (line[0] == '#' || line[0] == '\n') continue;
    levels++;
  }
  rewind(in);
           
  FW = (WIND*)GC_MALLOC((size_t) levels*sizeof(WIND));
  if (FW == NULL) {
    fprintf(stderr, "Cannot malloc memory for wind columns:[%s]\n", strerror(errno));
    return ERROR;
  } 
  
  while (NULL != fgets(line, MAX_LINE, in)) { 
  	if (line[0] == '#') continue;
	else {
		while (ret = sscanf(line, "%lf %lf %lf", 
	            &FW[i].wind_height,
	            &FW[i].windspeed,
	            &FW[i].wind_direction), ret != 3) { 
	     	if (ret == EOF && errno == EINTR) continue;
	        fprintf(stderr, 
	                "[line=%d,ret=%d] Did not read in 3 parameters:[%s]\n", 
	                i+1,ret, strerror(errno));
            return ERROR;
         } /* END while */
         
	     fprintf(log_file, "%lf %lf %lf\n",
                   FW[i].wind_height,
                   FW[i].windspeed,
                   FW[i].wind_direction);
          FW[i].wind_direction *= DEG2RAD; /* change to radians */
          i++;
    }
  } 
	 /* Since the wind levels are the 'last param' the total
	     number of params is determined when the number
	     of wind levels is determined. When FIXED_WIND is
	     TRUE, then the number of parameters and vertices is set here.
	 */
	INV_WIND_LEVELS = levels; 
	 fprintf(log_file, "INV_WIND_LEVELS=%d\n", INV_WIND_LEVELS);   
	/* NUM_OF_PARAMS = (int) (2 * levels); */
   	NUM_OF_PARAMS = (LAST_PARAM - 2);
   	fprintf(log_file, "NUM_OF_PARAMS=%d\n", NUM_OF_PARAMS);
   	NUM_OF_VERTICES = NUM_OF_PARAMS + 1;      
   
  	fprintf(log_file, "\tRead %d wind levels.\n", levels);
  	fprintf(log_file, "EXIT[get_wind].\n");	  	  
  	return FALSE;
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

  /* Every node assisgn the new eruption parameters to their copy of the model*/
  assign_new_params( param );
  
#ifdef DEBUG
  fprintf(log_file, "[%d] Calculating new value.\n", my_rank);
#endif

  /* Every node can now calculate a tephra blanket for their subset of POINTs */
    set_eruption_values(&e, W);
  for (i=0; i < local_n; i++)
    tephra_calc(&e, pt+i, W, &stats);
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
      //chi_calc = (*fit)(log_file, num_pts, p_all);
      chi_calc = e.fit(log_file, num_pts, p_all);
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

  int parm, wind_level=0, j;
	double wind_interval,wind_inversion_interval;
	double height, direction, speed, level, dir0, ht0, sp0;
	
#ifdef DEBUG
  	fprintf(log_file, "ENTER[assign_new_params] num_param=%d\n",NUM_OF_PARAMS);
#endif

    e.max_plume_elevation = param[MAX_COL_HT];
    e.alpha = param[ALPHAP];
    e.beta = param[BETAP];
    e.diffusion_coefficient = param[DIFF_COEF];
    e.total_ash_mass = param[TOTAL_MASS];
    e.mean_phi = param[MED_PHI];
    e.sigma_phi = param[STD_DEV_PHI];
    e.fall_time_threshold = param[FALLTIME_THRESH];
    e.eddy_constant = param[EDDY_COEF];

/*	fprintf(log_file, "%.0f %g %g %.0f %g %g %g %.0f %.2f\n",
    e.max_plume_elevation,
    e.alpha,
    e.beta,
    e.diffusion_coefficient,
    e.total_ash_mass,
    e.mean_phi,
    e.sigma_phi,
    e.fall_time_threshold,
    e.eddy_constant);
*/
	/*wind_interval is the column_step interval (mapped to)*/
	wind_interval = (e.max_plume_elevation - e.vent_elevation) / (double)(COL_STEPS);
	/* wind_inversion_interval is the inversion wind interval (mapped from)*/
	wind_inversion_interval = (e.max_plume_elevation - e.vent_elevation)/(double)INV_WIND_LEVELS;

	/* start at the vent */ 
	level = e.vent_elevation; /* level of column wind levels */ 
		
	for (j=0; j <= COL_STEPS; j++) { 
	  ht0 = 0.0;
		dir0 = 0.0;
		sp0 = 0.0; 
		height = e.vent_elevation; /* level of assigner */ 
		W[j].wind_height = 0; /* no wind level is assigned yet */
		
		if (FIXED_WIND) {
			for (wind_level=0; wind_level < INV_WIND_LEVELS; wind_level++) {
				height = FW[wind_level].wind_height;
				speed = FW[wind_level].windspeed;
				direction = FW[wind_level].wind_direction;
				/* This is the case where the first height is equal to
	         * or greater than the level that we are assigning.
	         */
        if (height >= level) {
        	if(height == level) {
        		W[j].wind_direction = direction;
	          W[j].windspeed = speed;
          } 
          else { /* interpolate */
	          W[j].wind_direction = 
	          	((direction - dir0) * (level - ht0) / (height - ht0)) + dir0;
            W[j].windspeed = 
	          	((speed - sp0) * (level - ht0) / (height - ht0)) + sp0;
	        }
	        W[j].wind_height = level;
	        break; /* ready to rescan the file for a match for the next level */
	      }
	      else { /* height is less than the level being assigned. */
	      	/* Maintain the real wind values for possible interpolation 
	        * at the next level.
	        */
	      	ht0 = height;
	        dir0 = direction;
	        sp0 = speed;
	      }
      } /* END for (wind_level */
	  } /* END if FIXED_WIND */
	  
   	else { /* We are using inversion to determine the wind */
   		for (wind_level=0, parm=WIND_SPEED;  parm < NUM_OF_PARAMS; parm +=2, wind_level++) { 
   				/* height += wind_inversion_interval; */
   				height = e.vent_elevation +((double)wind_level * wind_inversion_interval);
   				speed = param[parm];
        	direction = param[parm+1];
        	    	
        /* This is the case where the first height is equal to
	       * or greater than the level that we are assigning.
	       */
        if (height >= level) {
	      	if(height == level) {
	        	W[j].wind_direction = direction;
	          W[j].windspeed = speed;
          } 
          else { /* interpolate */
	          W[j].wind_direction = 
	          	((direction - dir0) * (level - ht0) / (height - ht0)) + dir0;
            W[j].windspeed = 
	            ((speed - sp0) * (level - ht0) / (height - ht0)) + sp0;
	        }
	        W[j].wind_height = level;
	        break; /* ready to rescan the file for a match for the next level */
	      }
	      else { /* height is less than the level being assigned. */
	        /* Maintain the real wind values for possible interpolation 
	        * at the next level.
	        */
	        ht0 = height;
	        dir0 = direction;
	        sp0 = speed;
	      }
      } /* END for (wind_level */
    } /* END else */
    
     /* If we finish scanning the inversion winds and all heights are below the level we are
      * currently assigning, then just use the direction and speed
	    * at the upper-most height.
	    */
	  if (!W[j].wind_height) {
	   	W[j].wind_height = level;
	    W[j].windspeed = sp0;
	    W[j].wind_direction = dir0;
	  }
/*	  fprintf(log_file, "[%d] %.0f %.2f %.2f\n", j, W[j].wind_height, W[j].windspeed, W[j].wind_direction);
*/
	  level += wind_interval;   
		/* Go to the next column height */  	
	} /* END for j=0 */
}
/******************************/
void set_LOG(FILE *log  ) {
  log_file = log;
}
/*****************************/
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
	 I  FIXED_WIND = 0, then initally these values get a random value within the LO_PARAM - HI_PARAM
	 range.  If FIXED_WIND = 1 then these values are set from the FW array (fixed wind array);
      */
      if (!FIXED_WIND) {
      	for (parm = WIND_SPEED; parm < NUM_OF_PARAMS; parm += 2) { /* for loop */
        	op[vert][parm] =
        	LO_PARAM(WIND_SPEED) + ((HI_PARAM(WIND_SPEED) - LO_PARAM(WIND_SPEED)) * rand()/(RAND_MAX+1.0));
			op[vert][parm+1] =
  			LO_PARAM(WIND_DIRECTION) + ((HI_PARAM(WIND_DIRECTION) - LO_PARAM(WIND_DIRECTION)) * rand()/(RAND_MAX+1.0));
      	} /* END for loop */
      } /* END if !FIXED_WIND */
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
INPUTS:  The array of inversion parameters, param[]
OUTPUTS:  none
 ************************************************************************/
void printout_model(double param[]) {

	int j = 0, parm;
  double height, speed, direction, wind_inversion_interval;
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
        e.alpha,
        e.beta,
        e.diffusion_coefficient,
       e.eddy_constant,
        e.total_ash_mass,
        e.max_plume_elevation,
        e.vent_elevation,
        e.fall_time_threshold);
    (void) fclose(out1);

  out2 = fopen(WIND_MODEL, "w");
  if (out2 == NULL) {
    fprintf(stderr,
	    "Cannot output model to file:[%s]. Printing to STDOUT.\n", strerror(errno));
    out2 = stdout;
  }
  
  fprintf(out2, "#HEIGHT\tSPEED\tDIRECTION\n");
  if (FIXED_WIND)  {
	for (j=0; j < INV_WIND_LEVELS; j++) { 
		height = FW[j].wind_height;
		direction = FW[j].wind_direction/DEG2RAD;
		speed = FW[j].windspeed;
		fprintf(out2, "%.0f %.2f %.2f\n", height, speed, direction);
		fprintf(stderr, "[%d] %.0f %.2f %.2f\n", j, height, speed, direction);
	}
 }
 else {
 	height = e.vent_elevation;
 	wind_inversion_interval = (e.max_plume_elevation - e.vent_elevation)/(double)INV_WIND_LEVELS;
 	for (parm = WIND_SPEED; parm < NUM_OF_PARAMS; parm += 2) { /* for loop */
   		speed = param[parm];
        direction = param[parm+1]/DEG2RAD;
        fprintf(out2, "%.0f %.2f %.2f\n", height, speed, direction);
		fprintf(stderr, "[%d] %.0f %.2f %.2f\n", j, height, speed, direction);
		height += wind_inversion_interval;
 	}
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
void printout_parameters(double fit) {

  FILE *out;
  time_t mytime;

  out = fopen(README, "w");
  if (out == NULL) {
    fprintf(stderr,
	    "Cannot output parameters to file:[%s]. Printing to STDOUT.\n",
	    strerror(errno));
    out = stdout;
  }
  mytime = time(&mytime);
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
	  fit,
    e.max_plume_elevation,
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
  /*double rmse; */

  out = fopen("plume.dat", "w");
  if (out == NULL) {
    fprintf(stderr,
	    "Cannot output parameters to file:[%s]. Printing to STDOUT.\n",
	    strerror(errno));
    out = stdout;
  }
  /*rmse = sqrt(chi);*/
	  fprintf (out, "%.0f %.0f  %d %g %g", e.max_plume_elevation, e.vent_elevation, (int) COL_STEPS, e.alpha, e.beta);

  if (out != stdout) fclose(out);

}
