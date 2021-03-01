#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <gc.h>
#include "prototypes.h"


static double PLUME_THRESHOLD;
static double AIR_VISCOSITY_x_18;
static double LITHIC_DENSITY_minus_PUMICE_DENSITY;
static double PUMICE_DIAMETER_THRESHOLD_minus_LITHIC_DIAMETER_THRESHOLD;
static double ONE_THIRD;
static double AIR_VISCOSITY_x_225;
static double GRAV_SQRD_x_4;
static double SQRT_TWO_PI; /* add new line, 2-22-2011 */
static double BETA_x_SQRT_TWO_PI;
static double TWO_BETA_SQRD;
static double PDF_GRAINSIZE_DEMON1;
static double TWO_x_PART_SIGMA_SIZE;
static TABLE **T;
static FILE *log_file;

/*
tephra2
Copyright (C) 2003  C. Bonadonna, C.B. Connor, L.J. Connor, T. Hincks
By: C. Bonadonna, C.B. Connor, L.J. Connor, T. Hincks

This file, tephra_calc.c, is part of tephra2.

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
/* ---------------------------------------------------------------------
 * FUNCTION:tephra_calc.c 
 *
 * Purpose: This function calculates and returns the expected accumulation
 * of volcanic ash (kg/m2) at a specific geogrphic location (x,y)
 * due to an eruption with specific input parameters. 
 * These points may be random or on  a UTM grid (m)
 * 
 * This implementation accounts for variation in wind velocity with height.
 * The model is discretized w.r.t. height and particle size. 

 * This function is called for each point (x,y,) If more than one eruption is
 * involved, for example in a probabilistic analysis, the function is called for
 * each set of eruption parameters. 
 
 * INPUTS:
 * ERUPTION *erupt: pointer to array of eruption parameters
 * POINT *pt: pointer to an array of location specific parameters, 
 * WIND *day: pointer to a day of wind data :
              height asl in m; 
              wind speed in ms-1
              wind direction in degrees N
 *	     
 * OUTPUTs:
 *   the value of the mass accumulated at the input location (northing, easting) in kg/m2 
 *
 *   a distribution of particle sizes
 *   the exact number of binss (i.e. sizes) and phi size used per bin is an integer and 
 *   is determined by (erupt->max_phi - erupt->min_phi)
 *   each bin accumulates phi sizes up to its integer size
 *   ex. bin[0] holds grainsizes [min_phi to min_phi+1)  
 ***************************************************************************/

void tephra_calc(ERUPTION *erupt, POINT *pt, WIND *day, STATS *stats) { /* tephra_calc starts ... */
 
   /**********************************************************************************
   * WIND structure:
   * day->windspeed: wind speed in m/s
   * day->wind_direction: wind direction in +/- degrees from north
   * day->wind_height: meters above sea level
   *
   * See common_structures_strat.h for structure definitions.
   **********************************************************************************/
  
  
  int i, j=0; /*, bin = -1; */
  double new_xspace, new_yspace, cos_wind = 0.0, sin_wind = 0.0, windspeed = 0.0;
  double sigma, demon2=0, demon3=0, ash_fall = 0.0, layer, fall_time_adj = 0.0, total_fall_time=0.0;
  double average_windspeed_x, average_windspeed_y, average_wind_direction, average_windspeed =0.0;
  double x_adj = 0.0, y_adj = 0.0;
  static double min=10e6, max=0.0;
  
#ifdef _PRINT
  fprintf(log_file, "IN tephra_calc ...");
#endif

  /* Initialize mass to zero */
  pt->calculated_mass = 0.0;
	
  /* Transform the volcano location coordinate to 0,0 */
  new_xspace = pt->northing - erupt->vent_northing;
  new_yspace = pt->easting - erupt->vent_easting;
/* fprintf(log_file, "point[x][y]=%.0f, %.0f\n",new_xspace, new_yspace); */

  /* do the double integration over grainsize and column height */
  
  /* Interpolate to fine the wind speed and direction below the height of the vent.
     Find one average speed and direction between vent and grid elevation point. 
   * The first values in the wind array give the wind speed and direction at the
   * vent height.
   */
  layer = erupt->vent_elevation - pt->elevation;
  windspeed = (day[0].windspeed * pt->elevation) / erupt->vent_elevation;
  cos_wind = cos(day[0].wind_direction) * windspeed;
	sin_wind = sin(day[0].wind_direction) * windspeed;

  for (i = 0; i < PART_STEPS; i++) { /* PART_STEPS_LOOP */
    fall_time_adj = 0.0;
    /* Accumulate the particle sizes into bins of whole numbered phi sizes 
    if (!(i % 10)) {
      bin++;
	Initiialize new phi accumulator to zero */
	pt->calculated_phi[i] = 0.0;
#ifdef _PRINT
    fprintf(log_file, "PART_STEP=%d phi[%d] = %g\n", i, i, pt->calculated_phi[i]);
#endif
   /* } */
  
    /* Adjust the total fall time of each particle size (i) 
       by the time it takes to descend from vent height to the grid cell (pt) 
     */           
     if (layer > 0) {
     	 fall_time_adj = 
     	 part_fall_time(erupt->vent_elevation, layer, T[i][0].ashdiam, T[i][0].part_density);
     
    }

     for (j = 0; j < COL_STEPS; j++) { /* COL_STEPS_LOOP */
     
    	total_fall_time = T[i][j].total_fall_time + fall_time_adj;
	   
	   /* Calc. x-component and y-component windspeed adjustments  
    	 * for each particle size  falling from each level.
    	 */
     /* removed 2 lines, 2-22-2011
   		x_adj = cos_wind * fall_time_adj * windspeed;		
	  	y_adj = sin_wind * fall_time_adj * windspeed;
	  	*/
	  	/* change 2 lines, 2-22-2011 */
	  	x_adj = cos_wind * fall_time_adj;
	  	y_adj = sin_wind * fall_time_adj;
	  	
			/* Now add the adjustments to the already summed
	     * wind components and
    	 * find the average windspeed in the x and y directions 
    	 * over the total fall time for each particle size at each level.
			*/
			average_windspeed_x = 
			(T[i][j].wind_sum_x + x_adj)/total_fall_time;
			
			average_windspeed_y = 
			(T[i][j].wind_sum_y + y_adj)/total_fall_time;
   	    	
    	/* If zero, make windspeed a very small value (cannot divide by zero in next step) */
      if (!average_windspeed_x) average_windspeed_x = .001;
      if (!average_windspeed_y) average_windspeed_y = .001;
      
      /* Find the average wind direction (direction of the velocity vector) */
    	if (average_windspeed_x < 0) {
    	  average_wind_direction =
    	  atan(average_windspeed_y/average_windspeed_x ) + pi;
    	}
    	else 
    	  average_wind_direction = atan(average_windspeed_y/average_windspeed_x);
    	
    	/* Find the average wind speed ( magnitude of the velocity vector) */
      average_windspeed =
      sqrt(average_windspeed_x * average_windspeed_x + 
      average_windspeed_y * average_windspeed_y);

    	if (total_fall_time > max) max = total_fall_time;
    	if (total_fall_time < min) min = total_fall_time;
    	
    	/* calculate the value of sigma (dispersion) based on total_fall_time  
    	    to acct for the change in the shape of the column with ht - increasing radius
          ht_above_vent = T[i][j].particle_ht - erupt->vent_elevation; 
       */
       
      /* falltime for fine particles */
      if (total_fall_time >= erupt->fall_time_threshold) {      
         sigma = 
         erupt->eddy_constant * pow((total_fall_time + T[i][j].plume_diffusion_fine_particle), 2.5);
      } 
      
      /* falltime for coarse particles */  
      else {
         sigma =
         4.0 * erupt->diffusion_coefficient * (total_fall_time + T[i][j].plume_diffusion_coarse_particle);
      }     
    
      demon2 =  pi * sigma;

      /* Modify fall time by the variation of wind velocity with height */
      demon3 = strat_average( average_wind_direction, 
                       average_windspeed,             
			                 new_xspace, new_yspace, 
			                 total_fall_time,
			                 sigma); 
			         
			 ash_fall = (T[i][j].demon1 / demon2) * demon3;
			 pt->calculated_mass += ash_fall;
			 pt->calculated_phi[i] += ash_fall;
		}  /* COL_STEPS_LOOP */  
   #ifdef _PRINT   
fprintf(log_file, "bin[%g] mass[%g] part[%g]", pt->calculated_phi[i], pt->calculated_mass, ash_fall); 
			fprintf(log_file, "\n"); 
			 

  #endif
  } /* PART_STEPS_LOOP */

/*fprintf(log_file, "COL_STEPS=%d PART_STEP=%d phi[%d] = %g mass=%g ash_fall=%g\n", j, i, bin, pt->calculated_phi[bin], pt->calculated_mass, ash_fall); */
#ifdef _PRINT
  fprintf(log_file, "OUT\n");
#endif

  stats->min_falltime = min;
  stats->max_falltime = max;
}


/* ----------------- New Function Starts Here -------------------- */
/* function phi2m converts the ash diameter from 
   units of phi to m
*/

 double phi2m(double xx) {
   double cms;
   cms = 0.001 * pow(2, -xx);
   return cms;
 }
/* ----------------- New Function Starts Here -------------------- */
/* function particle_density calculates varying particle density
 * based on their grain size diamete using a linear correlation
 * between pumice_threshold (PHI) and lithic_threshold (PHI)
*/

 double particle_density (double phi_slice, double lithic_density, double pumice_density) {

  double mean_density = 0;

  if (phi_slice >= LITHIC_DIAMETER_THRESHOLD) mean_density = lithic_density;
  else if (phi_slice <= PUMICE_DIAMETER_THRESHOLD) mean_density = pumice_density;
  else if (phi_slice < LITHIC_DIAMETER_THRESHOLD && phi_slice > PUMICE_DIAMETER_THRESHOLD) 
    mean_density =     
      lithic_density - LITHIC_DENSITY_minus_PUMICE_DENSITY * 
      (phi_slice - LITHIC_DIAMETER_THRESHOLD) / PUMICE_DIAMETER_THRESHOLD_minus_LITHIC_DIAMETER_THRESHOLD;

   return mean_density;
 }


/* ----------------- New Function Starts Here -------------------- */
/* function part_diff_time determines the particle diffusion time const
   in the atmosphere, where 
   col_ht = height of the particle (m) in the column w.r.t. vent. 
   THis is added to the particle fall time to account for the width of the column,
   which changes as a function of height

   c = eddy diffusivity in the atmosphere
   given in si units, e.g., 0.04 m^2 s^-1
   
   returns the particle diffusion time
*/


/* ----------------- New Function Starts Here -------------------- */
/* function part_fall_time determines the time of particle fall within each falling step
   falling steps are here:
   set = particle rising steps = ht_step_width
   
   returns the particle fall time within each falling step

This function follows the approach outlined in Bonadonna et al. (1998) 
Briefly, particle fall time is calculated based on terminal velocities in
layers that are 1000 m thick. The terminal velocity is a function of the
particle Reynolds number, which varies with grainsize, air properties.

The thickness of the first layer (closest to the ground) is equal to the vent 
height. The vent_height is in meters above sea level. The area the ash falls on
is considered to be at sea level. This leads to some assumptions (!) near the
volcano...
*/


double part_fall_time(double particle_ht, double layer, double ashdiam, double part_density) {
  
	double rho, hz, temp0, temp1;
   	double vtl, vti, vtt;
   	double reynolds_number;
   	double particle_term_vel;
   	double particle_fall_time;
 
  	particle_fall_time = 0.0;
  	hz = particle_ht;  /* height of the particle above sea level */
    
  	/*rho is the density of air (kg/m^3) at the elevation of the current particle*/
  	temp0 = -hz / 8200.0;
  	rho = AIR_DENSITY * exp(temp0);
  
	/*
   	(friction due to the air) :
    	vtl is terminal velocity (m/s) in laminar regime RE<6 
    	vti is terminal velocity (m/s) in intermediate regime 6<RE<500
    	vtt is terminal velocity (m/s) in turbulent regime RE>500
  	*/
  	vtl = part_density * GRAVITY * ashdiam * ashdiam / AIR_VISCOSITY_x_18; /* 18.0 * AIR_VISCOSITY */
  
  	/*
    	vti = ashdiam * 
    	pow(((4.0*GRAVITY*GRAVITY*erupt->part_mean_density *erupt->part_mean_density )/		(225.0*AIR_VISCOSITY*rho)),(1.0/3.0));
    	vtt=sqrt(3.1*erupt->part_mean_density *GRAVITY*ashdiam/rho);
  	*/
  
  	/*
    	RE is calculated using vtl (RE is Reynolds Number)
  	*/
  	reynolds_number = ashdiam * rho * vtl / AIR_VISCOSITY;
  	particle_term_vel = vtl;
  	temp0 = ashdiam * rho;


  	/*
    	c...if laminar RE>6 (intermediate regime), RE is calculated again considering vti
  	*/
  
	if (reynolds_number >= 6.0) {

    		/*4.0 * GRAVITY * GRAVITY * part_density * part_density / AIR_VISCOSITY * 225.0 * rho */
    		/*temp1 = GRAV_SQRD_x_4 * part_density * part_density / AIR_VISCOSITY_x_225 * rho; */
    		/*Added parentheses around deniminator - 04/23/2012 */
    		temp1 = GRAV_SQRD_x_4 * part_density * part_density / (AIR_VISCOSITY_x_225 * rho); 
    		vti = ashdiam * pow(temp1, ONE_THIRD); /* ONE_THIRD = 1.0/3.0 */    
    		reynolds_number = temp0 * vti / AIR_VISCOSITY;
    		particle_term_vel = vti;
    		/*
    		c...if intermediate RE>500 (turbulent regime), RE is calculated again considering vtt 
  		*/ 
  		if (reynolds_number >= 500.0) {
    			vtt = sqrt( 3.1 * part_density * GRAVITY * ashdiam / rho);
    			reynolds_number =  temp0 * vtt / AIR_VISCOSITY; 
    			particle_term_vel = vtt;
  		}  
  	}
/* Calculate the time it takes this particle to fall through this distance=layer */
  particle_fall_time = layer / particle_term_vel;
  
/* particle fall time is in sec   */
  
  //printf("i= %d, layer = %f, hz = %f, particle_term_vel = %f, diam=%f, reynolds = %f\n", i,layer_thickness, hz, particle_term_vel, a//shdiam, reynolds_number);
  
  return particle_fall_time;
}



/* ----------------- New Function Starts Here -------------------- */
/* this function calculates the expected fraction of particles
   in a given grainsize class (part_size_slice) assuming a normal 
   distribution in phi units about the mean, dmean,
   with standard deviation sigma.
*/

double pdf_grainsize(double part_mean_size, double part_size_slice, double part_step_width) {
  
  double func_rho, temp;
  double demon3, demon2;
  
  /* PDF_GRAINSIZE_DEMON1 = 1.0 / 2.506628 * erupt->part_sigma_size; */
  demon3   = part_size_slice - part_mean_size;
  temp = -demon3 * demon3 / TWO_x_PART_SIGMA_SIZE; 
  /* 2.0 * erupt->part_sigma_size * erupt->part_sigma_size */
  demon2   = exp(temp);
  func_rho = PDF_GRAINSIZE_DEMON1 * demon2 * part_step_width; 
/*  if (func_rho <= 0.0 || isnan(func_rho)) 
    fprintf(stderr,
    "[method pdf_grainsize] error: func_rho=%f]: demon2=%f step_width=%f\n", func_rho, demon2, part_step_width); */
  return func_rho;
} 
/* ----------------- New Function Starts Here -------------------- */
/* Function strat_average accounts for the variation in wind velocity 
   with height by using the average velocity value
   
		exp[ -5{ (x'-ut)^2 + y'^2} / {8*pi*C(t+td)} ]
		
   over the path of the particle as it falls from 
   its column release height to the ground.
    
   u = wind velocity (m/s), varies with height
   t = particle fall time 
   td = particle diffusion time
   
   The Suzuki equation has been formulated s.t. the wind is in the x direction
   We therefore need to transform the coordinates (xspace and yspace) to xprime
   and yprime, with xprime increasing in the downwind direction:
   x' = x cos a + y sin a
   y' = y cos a - x sin a

   
*/
double strat_average(
    double average_wind_direction, 
    double average_windspeed,
    double xspace, double yspace,
    double total_fall_time,
    double sigma) {

    double temp0, temp1, xprime, yprime, demon1, demon3;

    temp0 = cos(average_wind_direction);
    temp1 = sin(average_wind_direction);
    
    xprime = xspace * temp0 + yspace * temp1;
    yprime = yspace * temp0 - xspace * temp1;
    
    temp0 = xprime - average_windspeed * total_fall_time;
    demon1 = temp0 * temp0 + yprime * yprime;
    demon3 = exp(-demon1/sigma); /* where sigma is calculated for the total fall time */
    return demon3;
}

/***************************************************************
   inputs:
   x: height of a particle within the plume, relative to vent height
   slice: integration step (index)
   ht_section_width: the width of an integration step
   none: not used
   
   output: the probability that a given grainsize will be released from a given height
********************************************************************************/
double plume_pdf0(double x_norm, double step, double none0, double none1, double sum_prob) {

    double probability = 0.0;
    static double prob = 0;
    int num_slices_left= 0;
    int i;
    double x;
    /*fprintf(stderr, "x_norm=%g step=%g sum_prob=%g\n", x_norm, step, sum_prob);*/
    x = x_norm;
    if (!sum_prob) {
        for (i=0; i < COL_STEPS; i++) {
            x += step;
            if (x >= PLUME_THRESHOLD) { /* height at which tephra begins to be released */
                if (!num_slices_left) 
                    num_slices_left = COL_STEPS - i;
                prob = 1.0 / (double)num_slices_left;
                /* fprintf(stderr, "slices left = %d\n ", num_slices_left); */
                probability += prob;
            }
        }
    }
    else {
        if (x >= PLUME_THRESHOLD) /* height at which tephra begins to be released */
            /* Just calculate the probability for one column step */
            probability = prob;
    }
    return probability; 
}
/****************************************
    double probability;
    static int num_slices_left = 0;
    static double plume_slice = 0.0;
    
    probability = 0.0;
    if (x > PLUME_THRESHOLD) { height at which tephra begins to be released 
        if (!num_slices_left) {
            num_slices_left = COL_STEPS - (INT)slice;
            plume_slice = 1.0 / (double)num_slices_left;
            fprintf(stderr, "slices left = %d\n ", num_slices_left); 
        }
        probability = plume_slice;
    }
    return probability; 
}
************************************************/

/*************************************************************** 
   inputs:
   x: height of a particle within the plume, relative to vent height
   slice: integration step
   beta: the column beta parameter
   none: not used
   
   output: the probability that a given grainsize will be released from a given height
********************************************************************************/
double plume_pdf1(double x, double slice, double plume, double total, double none) {

  double col_ht, probability, beta_limit;
  double temp1, temp0, demon1, demon2;

  beta_limit = plume;
  col_ht = beta_limit - beta_limit * x / total;
  if (col_ht <= 0.0) col_ht = 1e-9; 
  temp1 = log(col_ht);
  temp1 *= temp1;
  temp0 = -temp1/TWO_BETA_SQRD; /* 2.0 * beta * beta */
  demon1 = exp(temp0);
  demon2 = col_ht * BETA_x_SQRT_TWO_PI; /* beta * sqrt(2.0 * PI) */
  

  probability = demon1 / demon2;
/*  if (isnan(probability)) {
    fprintf(stderr, "ht=%g  demon1=%g demon2=%g temp0=%g\n", 
	    col_ht, demon1, demon2, temp0); 
    exit(-1);
  }*/
  
/*  if (probability < 0.0) 
    fprintf(stderr, "col_ht=%f demon1=%f demon2=%f prob=%f\n", col_ht, demon1, demon2, probability); 
*/
  return probability;
}

/******************************************************************
This plume model uses the beta distribution, parameterized by 2 shape
parameters, alpha and beta. These parameters alter the shape of the 
distribution curve. The probability density function uses values between 
0 and 1, so the actual value will be normalized to within the range [0...1].
The alpha and and beta parameters are positive values greater than zero.
When alpha and beta equal 1, we have the standard normal distribution.
That is, for each slice of the column, a given grainsize has equal probability
of being released.

Inputs:
   x: height of a particle within the plume, relative to vent height
   (normalized to a value between 0 and 1)
   slice: integration step (normalized to a value btween 0 and 1)
   alpha: shape parameter (positive value)
   betat: shape parameter (positive value)
   
   output: the probability that a given grainsize will be released from a given height
********************************************************************************/
double plume_pdf2(double x_norm, double step, double alpha, double beta, double sum_prob) {

    double probability = 0.0;
    double prob;
    double x, x1, a1, b1, x2;
    int i=0;
    
    x = x_norm;
    a1 = alpha - 1.0;
    b1 = beta -1.0;
    
    if (!sum_prob) {
        for (i=0; i < COL_STEPS; i++) {
            /* step is the small slice of the column as a fraction of the whol */
            x += step;
            if (x <= (double)0) {
                x1 = x + 0.001;
                x2 = 1.0 - x1;                 
                prob = pow(x1, a1) * pow(x2, b1);
            }
            else if (x >= (double)1) {
               x1 =  x - 0.001;
               x2 = 1.0 - x1;
                prob = pow(x1, a1) * pow(x2, b1);
            }
            else {
                x1 = 1.0 - x;
                prob = pow(x, a1) * pow(x1, b1);
             }
            /* fprintf(log_file, "[sum_prob=0][%d] x=%g step=%g prob=%g\n", i, x, step, prob );  */
             probability += prob;
            
        }
    }
    /* Just calculate the probability for one column step */
    else {
        x1 = 1.0 - x;
        probability = pow(x, a1 ) * pow(x1, b1);
        /* fprintf(log_file, "[%d] x=%g step=%g prob=%g\n", i, x, step, probability ); */
    }
        
    if (isnan(probability)) {
	fprintf(stderr,"ERROR: x_norm=%g step=%g sum_prob=%g prob=%g\n", x_norm, step, sum_prob, probability);
	probability = 0;
}    
    return probability;
}

/******************************************
 * Set the values for globally used parameters
 ******************************************/
void set_global_values(FILE *log) {

#ifdef _PRINT  
fprintf(log_file, "IN set_global_values ...");
#endif
    /* Set values for global static variables */
    log_file = log;
  
#ifdef _PRINT  
fprintf(log_file, "OUT");
#endif
}

void set_eruption_values(ERUPTION *erupt, WIND *wind) { /* set_eruption_values */

/************************************************************
   * The following parameters are the properties of a eruption
   * each eruption must have all of these parameters defined:
   *
   * erupt->total_ash_mass is the total amount of ash erupted by
   * the volcano over the course of the entire eruption or calculation period
   * erupt->max_part_size is the maximum particle diameter considered
   * in the calculation. This is input in phi units (so it will likely be
   * a negative number like -5 and appear to be less than min_part_size)
   * erupt->min_part_size is the minimum particle diameter condsidered in the
   * calculation. This input is in phi units.
   *
   * Note: erupt->max/min_part_size are used to set the limits of integration
   * on the calculation. Particles outside this range are not considered at all.
   *
   * erupt->part_mean_size is the mean particle diameter erupted in phi units
   * erupt->part_sigma_size is the standard deviation in particle diameter in phi units
   * erupt-> vent_height is the elevation of the vent m.a.s.l. in meters
   * erupt->max_column_height is the eruption column height m.a.s.l. 
   * (not used) erupt->column_beta is the shape factor governing 
   * the particle size distribution in the eruption column. 
   * A large vlaue of beta (1) places most of the particles
   * high in the eruption column, a low value of beta (0.01) spreads the particle density
   * lower in the column.
   *********************************************************************************/

    int i, j;
    double y, x;
    double total_P_col=0.0, total_P_part=0.0, cum_prob_part=0.0, cum_prob_col=0.0, total_P=0.0;
    double cum_fall_time = 0.0, wind_x, wind_y, ht_above_vent, temp;
    double prob =0.0, col_prob, part_prob;

    double ht_section_width;
    double part_section_width;
    double ht_step_width;
    double part_step_width;
    double x_norm, step_norm;

    double pmin=10e6, pmax=0.0;
		
#ifdef _PRINT
fprintf(log_file, "IN set_eruption_values ... ");
#endif

   /* PART_STEPS = (erupt->max_phi - erupt->min_phi) * 10; 
    
#ifdef _PRINT
printf(log_file, "PART_STEPS=%d\n", PART_STEPS);
#endif
  */
    /* PLUME_THRESHOLD = erupt->plume_ratio * (erupt->max_plume_elevation - erupt->vent_elevation); replaced 2-22-2011 */
    PLUME_THRESHOLD = erupt->vent_elevation + (erupt->plume_ratio * (erupt->max_plume_elevation - erupt->vent_elevation)); /*new line 2-22-2011 */
    SQRT_TWO_PI = sqrt(2.0 * pi); /* new line,2-22-2011 */
    BETA_x_SQRT_TWO_PI = erupt->column_beta * SQRT_TWO_PI;
    TWO_BETA_SQRD = 2.0 * erupt->column_beta * erupt->column_beta;
    /* PDF_GRAINSIZE_DEMON1 = 1.0 / 2.506628 * erupt->sigma_phi; line removed, 2-22-2011 */
    PDF_GRAINSIZE_DEMON1 = 1.0 / (2.506628 * erupt->sigma_phi); /* line added 2-22-2011 */
    TWO_x_PART_SIGMA_SIZE = 2.0 * erupt->sigma_phi * erupt->sigma_phi;
    erupt->eddy_constant  = erupt->eddy_constant * 8.0 / 5.0;
    AIR_VISCOSITY_x_18 = 18.0 * AIR_VISCOSITY;
    LITHIC_DENSITY_minus_PUMICE_DENSITY = erupt->lithic_density - erupt->pumice_density;
    PUMICE_DIAMETER_THRESHOLD_minus_LITHIC_DIAMETER_THRESHOLD = PUMICE_DIAMETER_THRESHOLD - LITHIC_DIAMETER_THRESHOLD;
    ONE_THIRD = 1.0 / 3.0;
    AIR_VISCOSITY_x_225 = AIR_VISCOSITY * 225.0;
    GRAV_SQRD_x_4 = 4.0 * GRAVITY * GRAVITY;
     T = NULL;
  
    /*define the limits of integration */ 
    part_section_width = erupt->max_phi - erupt->min_phi;
    part_step_width = part_section_width / (double)PART_STEPS;
    ht_section_width = erupt->max_plume_elevation - erupt->vent_elevation; 
    ht_step_width = ht_section_width / (double)COL_STEPS; 
    PLUME_THRESHOLD = erupt->plume_ratio;
    step_norm = ht_step_width/ht_section_width;
    x_norm = 0.0;
/******************* steps for nomalization of probabilities **********************************/
/*    x = erupt->vent_elevation;
    for (i=0; i < COL_STEPS; i++) {
        x += ht_step_width;
        prob = (*pdf)(x, (double)i, ht_section_width, erupt->max_plume_elevation);
        cum_prob_col += prob;
    }
    */
    cum_prob_col = erupt->pdf (x_norm, step_norm, erupt->alpha, erupt->beta, total_P_col);
    total_P_col = cum_prob_col;

    y = (erupt)->min_phi;
    for (i=0; i < PART_STEPS; i++) {
        prob = pdf_grainsize(erupt->mean_phi, y, part_step_width);
        cum_prob_part += prob;
        y += part_step_width;
    }
    total_P_part = cum_prob_part;

    /* Normalization constant */
    total_P = (total_P_col * total_P_part);
    /*fprintf(log_file, "Total_PC=%g  Total_P=%g  Total_C=%g\n",total_P, total_P_col, total_P_part);*/
/************************End of normalization steps  ******************************************/

    /* Dynamically allocated table for storing integration data.
     Used in the double integration steps below for each point considered.
    */
    if (T == NULL) {
        T = (TABLE **)GC_MALLOC((size_t)PART_STEPS * sizeof(TABLE *));
        if (T == NULL) {
            fprintf(log_file, 
            "Cannot malloc memory for Integration Table:[%s]\n", strerror(errno));
            exit(1);
        }
        for (i=0; i<PART_STEPS; i++) {
            T[i] = (TABLE *)GC_MALLOC((size_t)COL_STEPS * sizeof(TABLE));
            if (T[i] == NULL) {
                fprintf(log_file, 
                "Cannot malloc memory for Integration Table[%d]:[%s]\n", i, strerror(errno));
                exit(1);
            }
        }
    } 
    else {
        T = (TABLE **)GC_REALLOC(T, (size_t)PART_STEPS * sizeof(TABLE *));
        if (T == NULL) {
            fprintf(log_file, 
            "Cannot malloc memory for Integration Table:[%s]\n", strerror(errno));
            exit(1);
        }
        for (i=0; i<PART_STEPS; i++) {
            T[i] = (TABLE *)GC_REALLOC(T[i], (size_t)COL_STEPS * sizeof(TABLE));
            if (T[i] == NULL) {
                fprintf(log_file, 
                "Cannot malloc memory for Integration Table[%d]:[%s]\n", i, strerror(errno));
                exit(1);
            }
        }
    }

    /* Start with the maximum particle size */
    y = (erupt)->min_phi;
    for (i = 0; i < PART_STEPS; i++) { /* PART_STEPS_LOOP */
        /* y += part_step_width; changed to end of loop 2-22-2011 */
        T[i][0].part_density  =  particle_density(y, erupt->lithic_density, erupt->pumice_density);    
        T[i][0].ashdiam = phi2m(y);
        /* the expected fraction of particles of this size based on given mean and std deviation */
        part_prob = pdf_grainsize(erupt->mean_phi, y, part_step_width);
        cum_fall_time = 0.0;
        wind_x = 0.0;
        wind_y = 0.0;
      
        /* Start at the height of the vent */
        x = erupt->vent_elevation;
        x_norm = 0.0;
        for (j = 0; j < COL_STEPS; j++) { /* COL_STEPS_LOOP */
            /* define the small slice dz */
            x += ht_step_width;
            x_norm += step_norm;
            
            /* Calculate the time it takes a particle to fall from its release point
                in the column to the next column release point.
            */
            T[i][j].fall_time = 
            part_fall_time(x, ht_step_width, T[i][0].ashdiam, T[i][0].part_density); 
            
            /* Particle diffusion time (seconds) */
            ht_above_vent = x - erupt->vent_elevation;	     
            temp = 0.2 * ht_above_vent * ht_above_vent; 
            T[i][j].plume_diffusion_fine_particle = pow(temp, 0.4); /* 0.4 = 2.0/5.0 */
            
            T[i][j].plume_diffusion_coarse_particle =
            0.0032 *  (ht_above_vent *  ht_above_vent) / erupt->diffusion_coefficient; 
            /* Sum the windspeed and wind_direction for a each particle size
            * falling from each level. 
            * In the wind array, the first wind level
            * gives windspeed and direction at the vent height. 
            * Start with the next wind level, 
            * so that we are using the windspeed and direction 
            * starting from one step above the vent. 
            */
            /*fprintf(log_file, "[%d][%d] ht=%.2f wspeed=%.2f wdir=%.2f\n",i,j,x, wind[j+1].windspeed, wind[j+1].wind_direction);*/
            wind_x += 
                T[i][j].fall_time * wind[j+1].windspeed * cos(wind[j+1].wind_direction);
                
            wind_y += 
                T[i][j].fall_time * wind[j+1].windspeed * sin(wind[j+1].wind_direction);
                
            T[i][j].wind_sum_x = wind_x;
            T[i][j].wind_sum_y = wind_y;
            
            /* Accumulate the time it takes each particle size to descend
            * from its release point in the column  down
            * to its final resting place.This part of the code just 
            * calculates the fall_time from the release point in the eruption
            * column to the height of the vent.
            * The time it takes a particle to fall from the vent height 
            * to a grid cell will be calculated later. 
            */
            cum_fall_time += T[i][j].fall_time;
            T[i][j].total_fall_time = cum_fall_time;
            /*fprintf(log_file, "[%d][%d]Total fall time=%g\n",i,j,T[i][j].total_fall_time);*/
            if (T[i][j].total_fall_time > pmax) pmax = T[i][j].total_fall_time;
            if (T[i][j].total_fall_time < pmin) pmin = T[i][j].total_fall_time;
            if (x_norm <= (double)0)
                col_prob =
                   // (*pdf)((x_norm+0.001), step_norm, erupt->alpha, erupt->beta, total_P_col);
                  erupt-> pdf((x_norm+0.001), step_norm, erupt->alpha, erupt->beta, total_P_col);
            else if (x_norm >= (double)1)
                col_prob =
                   // (*pdf)((x_norm-0.001), step_norm, erupt->alpha, erupt->beta, total_P_col);
                   erupt->pdf((x_norm-0.001), step_norm, erupt->alpha, erupt->beta, total_P_col);
            else
                col_prob =
                  //  (*pdf)(x_norm, step_norm, erupt->alpha, erupt->beta, total_P_col);
                  erupt->pdf(x_norm, step_norm, erupt->alpha, erupt->beta, total_P_col);
                    
 /*fprintf(log_file, "[%d][%d] %.2f %.2f %g %g %g\n", i,j, x, x_norm, col_prob, part_prob, total_P); */
 
            /* Normalization is now done here */
              T[i][j].demon1 = (erupt->total_ash_mass * col_prob  * part_prob)/total_P;
              
/*fprintf(log_file, "[%d][%d]mass=%g, col_prob=%g, part_prob=%g, total_P=%g\n", i,j, erupt->total_ash_mass, col_prob, part_prob, total_P); */
              
            /*if (T[i][j].demon1 < 1 &&  T[i][j].demon1 >= 0)
                fprintf(log_file, "[%d][%d]mass=%g, col_prob=%g, part_prob=%g, total_P=%g\n", i,j, erupt->total_ash_mass, col_prob, part_prob, total_P);*/
                
           /* if (isnan(T[i][j].demon1)) {
                fprintf(stderr, "[%d][%d]total ash mass=%g, col_prob=%g part_prob=%g total_P=%g x_norm=%f\n", i, j, erupt->total_ash_mass, col_prob, part_prob, total_P, x_norm);
                exit(0);
            }*/
            
            T[i][j].particle_ht = x;
/*            
            fprintf(log_file,
	      "[%d][%d]%g %g %g %g %g %g %g %g %g %g\n",
	      i,j,
	      T[i][j].particle_ht, 
	      T[i][j].ashdiam, 
	      T[i][j].part_density, 
	      T[i][j].fall_time,
	      T[i][j].plume_diffusion_fine_particle,
	      T[i][j].plume_diffusion_coarse_particle,
	      T[i][j].total_fall_time,
	      T[i][j].wind_sum_x,
	      T[i][j].wind_sum_y,
	      T[i][j].demon1);
*/
        } /* END COL_STEPS_LOOP */ 
        
        y += part_step_width;   /* moved from beg of loop 2-22-2011 */
        
    } /* END PART_STEPS_LOOP */

	/*fprintf(log_file, "MIN particle fall time = %.0f\n", pmin);
	fprintf(log_file, "MAX particle fall time = %.0f\n", pmax);
*/
}
