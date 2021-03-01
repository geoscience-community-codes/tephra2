/*
tephra2

Copyright (C) 2003  C. Bonadonna, C.B. Connor, L.J. Connor, T. Hincks
By: C. Bonadonna, C.B. Connor, L.J. Connor, T. Hincks

This file, new_tephra.c, is part of tephra2.

tephra2 is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.
tephra2 is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with tephra2.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <gc.h>
#include "../common_src/prototypes.h"

/* The following Global Variables are assigned some default values 

double DIFFUSION_COEFFICIENT = 200.0;
double FALL_TIME_THRESHOLD = 180.0;
double EDDY_CONST = .04;
double LITHIC_DENSITY = 2350.0;
double PUMICE_DENSITY = 1000.0;
int PLUME_MODEL = 2;
double PLUME_RATIO = 0.1;
double ALPHA = 1.0;
double BETA = 1.0;
double WIND_INTERVAL;
int WIND_DAYS = 1;
int WIND_COLUMNS = 3;
double PLUME_HEIGHT = 5.0;
double ERUPTION_MASS = 1e10;
double MAX_GRAINSIZE = -7.0;
double MIN_GRAINSIZE = 7.0;
double MEDIAN_GRAINSIZE = 1.5;
double STD_GRAINSIZE = 2.0;
double VENT_EASTING = 0.0;
double VENT_NORTHING = 0.0;
double VENT_ELEVATION = 0.0;
*/

int COL_STEPS;
int PART_STEPS;

/*define the following data structures for the code
  full descriptions are found in common_structures.h */
static ERUPTION e;
static WIND **W;
static POINT *pt;

/*define the following variables passed among functions in this file */
static int num_pts = 0; /*total number of points used in the analysis */
/* static int num_eruptions = 1; total number of eruptions used in the analysis */
/* static int num_wind_data = 0; number of wind data points in the analysis */

/*static int local_n; */

/* FILE *in_eruptions; */
FILE *in_points;
FILE *in_wind;
FILE *log_file;

int cmp(double *x, double *y) {
  if (*x < *y) return -1;
  if (*x > *y) return 1;
  return 0;
}

void exit_now(int e) {
  (void) fclose(in_points);
  (void) fclose(in_wind);

#ifdef DEBUG
  (void) fclose(log_file);
#endif
  exit(e);
}

int main(int argc, char *argv[]) { /* MAIN CODE STARTS HERE */
  
  int i, /* j, */ bin;
  double val, bin_width;
  STATS stats;
  char log_name[25];


  /* Check for correct number of comand line arguments */
  if (argc != 4) {
      fprintf(stderr, 
	      "Missing comand line arguments,\nUSAGE: <program name> <config file> <points file> <wind file>\n\n");
    exit(1);
  }
  

  /* Each node opens a file for logging */
  sprintf(log_name, "%s", LOG_FILE);
  fprintf(stderr, "%s\n", log_name);
  log_file  = fopen(log_name, "w+");
  if (log_file == NULL) {
    fprintf(stderr, "Cannot open LOG file=[%s]:[%s]. Exiting.\n", 
	    log_name, strerror(errno));
    exit(1);
  }

  
  /* Initialize the global variables (see top of file) with inputs from the configuration file. */
  if ( init_globals(argv[1]) ) {
    exit(1);
  }
  
#ifdef _PRINT
  fflush(log_file); 
#endif

  /*make sure the points file exists*/
  in_points= fopen(argv[2], "r");
  if (in_points == NULL) {
    fprintf(stderr, "Cannot open points  file=[%s]:[%s]. Exiting.\n", 
	    argv[2], strerror(errno));
    exit_now(1);
  }
  
  /* Input the data points from a file using the
     get_points function. 
  */
  
  if (get_points(in_points) ) {
    exit_now(1);
  }

#ifdef _PRINT
  fflush(log_file); 
#endif
  
  /*make sure the wind file exists*/
  in_wind= fopen(argv[3], "r");
  if (in_wind == NULL) {
    fprintf(stderr, "Cannot open wind file=[%s]:[%s]. Exiting.\n", 
	    argv[3], strerror(errno));
    exit_now(1);
  }
  
  if (get_wind(in_wind) ) {
    exit_now(1);
  }
  
#ifdef _PRINT
  fflush(log_file); 
#endif
  
  set_global_values(log_file);
  
#ifdef _PRINT
  fflush(log_file); 
#endif
 

  set_eruption_values(&e, W[0]);
  for (i = 0; i < num_pts; i++)   /* For each location */
     tephra_calc(&e, pt+i, W[0], &stats);
 #ifdef _PRINT  
  fprintf(log_file, "Calculated mass...\n"); 
  for (i = 0;  i < num_pts;  i++) {
    fprintf(log_file, "[%.0f][%.0f] %g ", (pt+i)->easting, (pt+i)->northing, (pt+i)->calculated_mass) ;
    fprintf(log_file, "\n"); 
}
  fprintf(log_file, "Finished.\n");
  fflush(log_file);
    
  #endif
    
  fprintf(stderr,"\nMin Particle Fall Time = %gs\nMax Particle Fall Time = %gs\n",stats.min_falltime, stats.max_falltime);
	
  /*phi_bins = (int)((erupt+j)->max_phi - (erupt+j)->min_phi); 
  if (phi_bins < 1) phi_bins = 1; */
  bin_width = (e.max_phi - e.min_phi)/PART_STEPS;
  printf("#EAST NORTH ELEV Kg/m^2 "); 
  for (bin = 0; bin < PART_STEPS; bin++) 
    printf("[%g->%g) ", 
	e.min_phi + bin_width*bin, 
	e.min_phi + bin_width*(bin+1));
  
  printf("(percent)\n");
	
  fprintf(stderr, "PART_STEPS=%d bin_width=%g\n", PART_STEPS, bin_width);
  for (i = 0; i < num_pts; i++) {
    printf("%.0f %.0f %.0f %g ", 
	(pt+i)->easting, 
	(pt+i)->northing,
	(pt+i)->elevation, 
	(pt+i)->calculated_mass);
	/* Print out percent of total */
    for (bin = 0; bin < PART_STEPS; bin++) {
      val = ((pt+i)->calculated_phi[bin]/(pt+i)->calculated_mass) * 100.0;
      /*fprintf(stderr, "bin = %g mass = %g ", (pt+i)->phi[bin], (pt+i)->mass ); */
      printf("%g ", val);
    }
    printf("\n");
  }
  print_for_stats(0); 
  fprintf(log_file, "Finished.\n");
  exit_now(0);
  return 1;
}

void print_for_stats(double val) {

  FILE *out;

  out = fopen("plume2.dat", "w");
  if (out == NULL) {
    fprintf(stderr,
	    "Cannot output parameters to file:[%s]. Printing to STDOUT.\n",
	    strerror(errno));
    out = stdout;
  }
  fprintf (out, "%.0f %.0f  %d %g %g", e.max_plume_elevation, e.vent_elevation, (int) COL_STEPS, e.alpha, e.beta);
  if (out != stdout) fclose(out);
}


/**************************************************************
FUNCTION:  get_points
DESCRIPTION:  This function reads eruption data into the
ERUPTION array. Each node stores 
all of the eruption  parameters which are varied and then 
used in calculating the mass loading value at each point. 
INPUTS: (IN) FILE *in  (file handle from which to read)
OUTPUTS: int -1=error, 0=no error
***************************************************************/
int get_points(FILE *in) {
  
  char line[MAX_LINE];
  int i, j, ret, /* my_start=0,*/ pts_read=0;
  
#ifdef _PRINT
  fprintf(log_file,"ENTER[get_points]\n");
#endif
  while (fgets(line, MAX_LINE, in) != NULL)  {
    if (line[0] == '#' || line[0] == '\n') continue;
    num_pts++;
  }
  rewind(in);
  

 /* local_n = num_pts;*/
  
#ifdef _PRINT  
  fprintf(log_file, "Total locations: %d.\n", num_pts);
#endif    
  
  pt = (POINT *)GC_MALLOC((size_t)num_pts * sizeof(POINT));
  if (pt == NULL) {
    fprintf(stderr, "Cannot malloc memory for my points:[%s]\n", strerror(errno));
    return -1;
  } 
  
#ifdef _PRINT  
  fprintf(log_file, "\tReading in %d locations, starting at line %d.\n",
	  local_n, my_start);
#endif    
   
  i = 0;
  while (i < num_pts) {
    fgets(line, MAX_LINE, in);
    if (line[0] == '#' || line[0] == '\n') continue;
    else {
      while (ret = sscanf(line,
			  "%lf %lf %lf",
			  &(pt+pts_read)->easting,
			  &(pt+pts_read)->northing,
			  &(pt+pts_read)->elevation),
	     ret != 3) { 
	
        if (ret == EOF && errno == EINTR) continue;
        fprintf(stderr, "[line=%d,ret=%d] Did not read in 3 coordinates:[%s]\n", i+1,ret, strerror(errno));
        return -1;
      }
      /* Initialize some values in the point structure */
      /*(pt+pts_read)->cum_mass = 0.0; */
	    
      (pt+pts_read)->calculated_mass = 0.0;
      
      /*(pt+pts_read)->calculated_phi = (double *)GC_MALLOC((size_t)PART_STEPS * sizeof(double));
      if ((pt+pts_read)->calculated_phi == NULL) {
        fprintf(stderr, "Cannot malloc memory for phi_bins:[%s]\n", strerror(errno));
        return -1;
      } */
	    
      for (j = 0; j < PART_STEPS; j++)
        (pt+pts_read)->calculated_phi[j] = 0.0;
	      
      pts_read++;
      /* if (i >= my_start) {
	      pts_read++;
	      if (pts_read == local_n) break;
      } */
    }
    i++;
  }
  fprintf(log_file,"EXIT[get_points]:Read %d points.\n",
	  pts_read);
  
#ifdef _PRINT  
  fprintf(log_file, "EXIT[get_points].\n");
#endif		  
  return 0;
}

/**************************************************************
FUNCTION:  get_wind
DESCRIPTION:  This function reads wind data into the
WIND array. Each node stores all of the wind data. 

INPUTS: (IN) FILE *in  (file handle from which to read)
OUTPUTS: int -1=error, 0=no error
***************************************************************/
int get_wind(FILE *in) {
  
  int i=0, j=0, ret, WIND_DAYS = 1;
  char line[MAX_LINE];
  double wind_height, wind_dir, windspeed, dir0, ht0, sp0;
  double level, WIND_INTERVAL;
  
  fprintf(log_file,"ENTER[get_wind].\n");
  WIND_INTERVAL = (e.max_plume_elevation - e.vent_elevation)/(double)COL_STEPS;
  
  W = (WIND**)GC_MALLOC((size_t)WIND_DAYS * sizeof(WIND *));
  if (W == NULL) {
    fprintf(stderr, "Cannot malloc memory for wind columns:[%s]\n", strerror(errno));
    return -1;
  } else {
    for (i=0; i < WIND_DAYS; i++) {
      W[i] = (WIND *)GC_MALLOC((size_t)(COL_STEPS+1) * sizeof(WIND));
      if (W[i] == NULL) {
	fprintf(stderr, "Cannot malloc memory for wind rows %d:[%s]\n", i, strerror(errno));
	return -1;
      }
    }
  }
  /* Assume one wind day */
  i = 0;
  /* start at the vent */ 
  level = e.vent_elevation;
    
  /* Do for each column step */
  /* j = 0 is for the interval between the vent and the ground.
   * Here we set the wind speed and direction to be at the level of the vent;
   * The values used in the calculations change for each location and are 
   * set in the tephra_calc routine when the point elevation is known. 
   * The last interval ends at the top of the column. 
   */
   W[i][0].wind_height = e.vent_elevation;
    /* fprintf(log_file, "[%d]%f %f %f\n",  0, W[i][j].wind_height, W[i][j].windspeed, W[i][j].wind_direction); */
	
   for (j=0; j <= COL_STEPS; j++) { 
     ht0 = 0.0;
     dir0 = 0.0;
     sp0 = 0.0;
    
     /* Find wind elevation just greater than current level
           Start scanning the wind file for the best match.
           Each new level starts scanning the file from the beginning.
      */
      while (NULL != fgets(line, MAX_LINE, in)) {
        if (line[0] == '#') continue;
        else {
          while (ret = sscanf(line, "%lf %lf %lf", 
	        &wind_height,
	        &windspeed,
	        &wind_dir), ret != 3) { 
            if (ret == EOF && errno == EINTR) continue;
            fprintf(stderr, 
	           "[line=%d,ret=%d] Did not read in 3 parameters:[%s]\n", 
	            i+1,ret, strerror(errno));
            return -1;
       }
     }
     /* This is the case where we find the first height that is equal to
      * or greater that the level that we are assigning.
      */
      if (wind_height >= level) {
        if (wind_height == level) {
          W[i][j].wind_direction = wind_dir;
          W[i][j].windspeed = windspeed;
        } else { /* interpolate */
          W[i][j].wind_direction = 
	         ((wind_dir - dir0) * (level - ht0) / (wind_height - ht0)) + dir0;
          W[i][j].windspeed = 
	         ((windspeed - sp0) * (level - ht0) / (wind_height - ht0)) + sp0;
        }
        W[i][j].wind_height = level;
        break; /* ready to rescan the file for a match for the next level */
      }
      /* This is the case where the scanned height is less than the level
       * we are assigning.
       */
      else {
	/* Maintain the scanned values for possible interpolation 
	 * at the next level.
	 */
        ht0 = wind_height;
        dir0 = wind_dir;
        sp0 = windspeed;
      }  
    }
    /* If we finish scanning the file and all heights are below the level we are
    * currently assigning, then just use the direction and speed
    * at the upper-most height.
    */
    if (!W[i][j].wind_height) {
      W[i][j].wind_height = level;
      W[i][j].windspeed = sp0;
      W[i][j].wind_direction = dir0;
    }
    /* Go to the next column height */
    fprintf(log_file, "[%d] %f %f %f\n", 
	   j, W[i][j].wind_height, W[i][j].windspeed, W[i][j].wind_direction);
	   
    W[i][j].wind_direction *= DEG2RAD; /* change to radians */
    rewind(in); 
    level += WIND_INTERVAL; 
  } 

  fprintf(log_file, "\tRead %d wind days with %d wind levels per day.\n", i+1, j);
  fprintf(log_file, "EXIT[get_wind].\n"); 	  
  return 0;
}

/**************************************************************
FUNCTION:  get_config_data
DESCRIPTION:  This function reads the configuration file,
and sets some global variables.

INPUTS: (IN) FILE *in  (file handle from which to read)
OUTPUTS: int -1=error, 0=no error
***************************************************************/

int init_globals(char *config_file) {

  FILE *in_config;
  char buf[1][30], **ptr1;
  char line[MAX_LINE];
  char space[4] = "\n\t ";
  char *token;
 
#ifdef _PRINT   
  fprintf(log_file, "ENTER[init_globals].\n");
#endif
 
  in_config = fopen(config_file, "r");
  if (in_config == NULL) {
    fprintf(stderr, 
	    "Cannot open configuration file=[%s]:[%s]. Exiting.\n", config_file, strerror(errno)); 
    return 1;
  }
  
  ptr1 = (char **)&buf[0];
  while (fgets(line, MAX_LINE, in_config) != NULL) {
    /* fprintf(stderr, "%s\n", line); */
    if (line[0] == '#' || line[0] == '\n' || line[0] == ' ')  continue;
    token = strtok_r(line, space, ptr1);
    if (!strncmp(token, "DIFFUSION_COEFFICIENT", strlen("DIFFUSION_COEFFICIENT"))) {
      token = strtok_r(NULL,space,ptr1);
      e.diffusion_coefficient = strtod(token, NULL);
      /* DIFFUSION_COEFFICIENT can never be 0 as it is used in divisions */
      if (e.diffusion_coefficient < 1.0) e.diffusion_coefficient = 1.0;
      fprintf(stderr, "DIFFUSION_COEFFICIENT=%.0f\n", e.diffusion_coefficient);
    } 
    else if (!strncmp(token, "EDDY_CONST", strlen("EDDY_CONST"))) {
      token = strtok_r(NULL,space,ptr1);
      e.eddy_constant = strtod(token, NULL);
      fprintf(stderr, "EDDY_CONST=%.2f\n", e.eddy_constant);
    }
    else if (!strncmp(token, "FALL_TIME_THRESHOLD", strlen("FALL_TIME_THRESHOLD"))) {
      token = strtok_r(NULL,space,ptr1);
      e.fall_time_threshold = strtod(token, NULL);
      fprintf(stderr, "FALL_TIME_THRESHOLD=%.0f\n", e.fall_time_threshold);
    }
    else if (!strncmp(token, "LITHIC_DENSITY", strlen("LITHIC_DENSITY"))) {
      token = strtok_r(NULL,space,ptr1);
      e.lithic_density = strtod(token, NULL);
      fprintf(stderr, "LITHIC_DENSITY=%.0f\n", e.lithic_density);
    }
    else if (!strncmp(token, "PUMICE_DENSITY", strlen("PUMICE_DENSITY"))) {
      token = strtok_r(NULL,space,ptr1);
      e.pumice_density = strtod(token, NULL);
      fprintf(stderr, "PUMICE_DENSITY=%.0f\n", e.pumice_density);
    }
    else if (!strncmp(token, "PART_STEPS", strlen("PART_STEPS"))) {
      token = strtok_r(NULL, space, ptr1);
      PART_STEPS = (int)atoi(token);
      fprintf(stderr, "PART_STEPS = %d\n", PART_STEPS);
    }
    else if (!strncmp(token, "COL_STEPS", strlen("COL_STEPS"))) {
      token = strtok_r(NULL, space, ptr1);
      COL_STEPS = (int)atoi(token);
      fprintf(stderr, "COL_STEPS = %d\n", COL_STEPS);
    }
    else if (!strncmp(token, "PLUME_MODEL", strlen("PLUME_MODEL"))) {
      token = strtok_r(NULL,space,ptr1);
      e.plume_model = (int)atoi(token);
      
      if (e.plume_model == 0) {
	      e.pdf = &plume_pdf0;
	      fprintf(stderr, "PLUME_MODEL=[%d]%s\n", e.plume_model, "Uniform Distribution with threshold");
      }
      else if (e.plume_model == 1) {
	      e.pdf = &plume_pdf1;
	      fprintf(stderr, "PLUME_MODEL=[%d]%s\n", e.plume_model, "log-normal Distribution using beta");
      }
      else if (e.plume_model == 2) {
	     e.pdf = &plume_pdf2;
	      fprintf(log_file, "PLUME_MODEL=[%d]%s\n", e.plume_model, "Beta Distribution with ALPHA and BETA parameteres");
      }
    }
    else if (!strncmp(token, "PLUME_RATIO", strlen("PLUME_RATIO"))) {
      token = strtok_r(NULL, space, ptr1);
      e.plume_ratio = strtod(token, NULL);
      if (!e.plume_model) fprintf(stderr, "PLUME_RATIO = %.2f\n", e.plume_ratio);
    }
    else if (!strncmp(token, "ALPHA", strlen("ALPHA"))) {
      token = strtok_r(NULL, space, ptr1);
      e.alpha = strtod(token, NULL);
      if (e.plume_model==2) fprintf(stderr, "ALPHA = %g\n", e.alpha);
    }
    else if (!strncmp(token, "BETA", strlen("BETA"))) {
      token = strtok_r(NULL, space, ptr1);
      e.beta = strtod(token, NULL);
      if (e.plume_model==2) fprintf(stderr, "BETA = %g\n", e.beta);
    }
    else if (!strncmp(token, "PLUME_HEIGHT", strlen("PLUME_HEIGHT"))) {
      token = strtok_r(NULL, space, ptr1);
      e.max_plume_elevation = strtod(token, NULL);
      fprintf(stderr, "PLUME_HEIGHT = %.0f\n", e.max_plume_elevation);
    }
    else if (!strncmp(token, "ERUPTION_MASS", strlen("ERUPTION_MASS"))) {
      token = strtok_r(NULL, space, ptr1);
      e.total_ash_mass = strtod(token, NULL);
      fprintf(stderr, "ERUPTION_MASS = %g\n", e.total_ash_mass);
    }
    else if (!strncmp(token, "MAX_GRAINSIZE", strlen("MAX_GRAINSIZE"))) {
      token = strtok_r(NULL, space, ptr1);
      e.min_phi = strtod(token, NULL);
      fprintf(stderr, "MAX_GRAINSIZE = %g (phi)\n", e.min_phi);
    }
    else if (!strncmp(token, "MIN_GRAINSIZE", strlen("MIN_GRAINSIZE"))) {
      token = strtok_r(NULL, space, ptr1);
      e.max_phi = strtod(token, NULL);
      fprintf(stderr, "MIN_GRAINSIZE = %g (phi)\n", e.max_phi);
    }
    else if (!strncmp(token, "MEDIAN_GRAINSIZE", strlen("MEDIAN_GRAINSIZE"))) {
      token = strtok_r(NULL, space, ptr1);
      e.mean_phi = strtod(token, NULL);
      fprintf(stderr, "MEDIAN_GRAINSIZE = %g (phi)\n", e.mean_phi);
    }
    else if (!strncmp(token, "STD_GRAINSIZE", strlen("STD_GRAINSIZE"))) {
      token = strtok_r(NULL, space, ptr1);
      e.sigma_phi = strtod(token, NULL);
      /* if (STD_GRAINSIZE < 1.0) STD_GRAINSIZE = 1.0; */
      fprintf(stderr, "STD_GRAINSIZE = %g (phi)\n", e.sigma_phi);
    }
    else if (!strncmp(token, "VENT_EASTING", strlen("VENT_EASTING"))) {
      token = strtok_r(NULL, space, ptr1);
      e.vent_easting = strtod(token, NULL);
      fprintf(stderr, "VENT_EASTING = %.0f\n", e.vent_easting);
    }
    else if (!strncmp(token, "VENT_NORTHING", strlen("VENT_NORTHING"))) {
      token = strtok_r(NULL, space, ptr1);
      e.vent_northing = strtod(token, NULL);
      fprintf(stderr, "VENT_NORTHING = %.0f\n", e.vent_northing);
    }
    else if (!strncmp(token, "VENT_ELEVATION", strlen("VENT_ELEVATION"))) {
      token = strtok_r(NULL, space, ptr1);
      e.vent_elevation = strtod(token, NULL);
      fprintf(stderr, "VENT_ELEVATION = %.0f\n", e.vent_elevation);
    }
    else continue;
  }
  (void) fclose(in_config);

#ifdef _PRINT
  fprintf(log_file, "EXIT[init_globals].\n");
#endif

  return 0;
}
