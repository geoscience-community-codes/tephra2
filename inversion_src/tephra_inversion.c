#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <errno.h>
#include <time.h>
#include <mpi.h>
#include "../common_src/prototypes.h"



/* The following are Global Variables defined in parameters.h*/
int FIXED_WIND;
int NUM_OF_PARAMS;
int NUM_OF_VERTICES;
double TOLERANCE;

double _LO[LAST_PARAM] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,0.1};
double _HI[LAST_PARAM] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,1.0};

/* Input file handles */								
FILE *in_pts;  
FILE *log_file;
FILE *wind_file;
/******************************************************************************************** 
Program: tephra_inversion.c
Authors: Laura Connor, Chuck Connor, Costanza Bonadonna
Date: March 2004
Language: Ansi C

Purpose: This program will perform an inversion on collected tephra samples to discover the
         eruption dynamics of the eruption producing the sampled tephra blanket. 

*********************************************************************************************/

/* 
	 Subroutine to shutdown program and close open files
*/

void exit_now(int level, int e) {
	
	int i;
	for (i = 0; i < level; i++) {
		if (i==0) (void) fclose(log_file);
  		if (i==1)(void) fclose(in_pts);
  		if (i==2) (void) fclose(wind_file);
  }
  MPI_Finalize();
  exit(e);
}

/*
	Main program
*/

int main(int argc, char *argv[]) {

  char log_name[25];
  double quit = 0.0;
  int i;
  int my_rank; /* process rank of each node (local) */
  int procs; /* number of nodes used for processing */
  double chi;

  /* Start up MPI */
  MPI_Init(&argc, &argv);

   /* Get my process rank, an integer */
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	/* Check for correct number of comand line arguments */
  if (argc < 3) {
    if (!my_rank)
      fprintf(stderr, 
	      "Missing comand line arguments,\nUSAGE: <program name> <config file> <points file> <wind file (optional)>\n\n");
    MPI_Finalize();
    exit(1);
  }

  /* Find out how many processes are being used */
  MPI_Comm_size(MPI_COMM_WORLD, &procs);
  
	/* Each node opens a file for logging */
  sprintf(log_name, "%s%d", LOG_FILE, my_rank);
  fprintf(stderr, "[%d]logging to %s....\n",my_rank, log_name);

  log_file  = fopen(log_name, "w+");
  if (log_file == NULL) {
    fprintf(stderr, "Cannot open LOG file=[%s]:[%s]. Exiting.\n", 
	    log_name, strerror(errno));
    exit_now(1,1);
  }
  
  set_LOG(log_file);
  
#ifdef DEBUG
  fprintf(stderr, "Opened log file.\n");
#endif

	/* Initialize the global variables with inputs from the configuration file. */
	if ( init_globals(argv[1])) {
		fflush(stderr); 
		exit_now(1,1);
	}

#ifdef DEBUG
  fprintf(stderr, "Initialized global variables.\n");
#endif

  /* Only the Master node (node_0) will print out some values to the console */
  if (!my_rank) 
    fprintf(stderr, 
"Using fixed wind? %s\n\
NUM_OF_PARAMS = %d\n\
NUM_OF_VERTICES = %d\n\
TOLERANCE = %g\n\
Wind Speed: min=%.0f max=%.0f\n\
Wind Direction: %.0f to %.0f\n\
Plume Height: min=%.0f max=%.0f\n\
Alpha Param: min=%g max=%g\n\
Beta Param: min=%g max=%g\n\
Fall Time Threshhold: min=%.0f max=%.0f\n\
Diffusion Coefficient: low=%.0f high=%.0f\n\
Eddy Constant: low=%.2f high=%.2f\n\
Total Mass Ejected: min=%g max=%g\n\
Grain Size:\n\
\tMedian PHI: min=%g max=%g\n\
\tStd. Deviation PHI: min=%g max=%g\n",
      (FIXED_WIND) ? "TRUE" : "FALSE",
	    NUM_OF_PARAMS, 
	    NUM_OF_VERTICES, 
	    TOLERANCE,
	    LO_PARAM(WIND_SPEED), HI_PARAM(WIND_SPEED),
	    LO_PARAM(WIND_DIRECTION)/DEG2RAD, HI_PARAM(WIND_DIRECTION)/DEG2RAD,
	    LO_PARAM(MAX_COL_HT), HI_PARAM(MAX_COL_HT),
	    LO_PARAM(ALPHAP), HI_PARAM(ALPHAP),
	    LO_PARAM(BETAP), HI_PARAM(BETAP),
	    LO_PARAM(FALLTIME_THRESH), HI_PARAM(FALLTIME_THRESH),
	    LO_PARAM(DIFF_COEF), HI_PARAM(DIFF_COEF),
        LO_PARAM(EDDY_COEF), HI_PARAM(EDDY_COEF),
	    LO_PARAM(TOTAL_MASS), HI_PARAM(TOTAL_MASS),
	    LO_PARAM(MED_PHI), HI_PARAM(MED_PHI),
	    LO_PARAM(STD_DEV_PHI), HI_PARAM(STD_DEV_PHI));
 
  /* Read in the data values and locations */ 
  in_pts = fopen(argv[2], "r");
  if (in_pts == NULL) {
    fprintf(stderr, "Cannot open POINTS file=[%s]:[%s]. Exiting.\n", 
	    argv[2], strerror(errno));
    exit_now(1,1);
  }
  if (get_points(in_pts) ) 
    exit_now(2,1);
  
  if (FIXED_WIND) {
    if (argc == 4) {
  	  /*make sure the wind file exists*/
  	  wind_file = fopen(argv[3], "r");
  	  if (wind_file == NULL) {
    	  fprintf(stderr, "Cannot open wind file=[%s]:[%s]. Exiting.\n", 
	      argv[3], strerror(errno));
	  	  exit_now(2,1);  
  	  }
    }
    if (get_wind(wind_file) ) {
      exit_now(3,1);
    }
  }
  
  fprintf(stderr, "[%d] Finished with input - running the optimization\n", my_rank);

  if ( my_rank ) { /* slave */
    /* the slave nodes now wait to be called
     * upon by the master to calculate their portion of the magnetic values
     */
     fprintf(log_file, "[%d ] Slave ready .....", my_rank);
    slave(my_rank, log_file); 
  } else { /* master */
    fprintf(log_file, "Master ready .....\n");
    chi = master(log_file);
    
    /* Send all slaves a quitin' time signal (i.e. a single zero value */
    for ( i = 1; i < procs; i++ )
      MPI_Send((void *)&quit, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
    
    /* The Master node prints out a README file listing some input parameters and changed values */
    printout_parameters(chi);
    print_for_stats(chi);
    
  } /* end master code */
  
  /* _free();*/ 
  
  exit_now(1,0);
  return(0);
}
