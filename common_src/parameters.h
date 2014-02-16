/* parameters.h */
/* Define the physical constants used in the new_tephra.c code */

#define LOG_FILE "node_"
#define ERROR -1
/* #define DEBUG 1 */
/* #define _PRINT 1 */
#define R 0.61803399
#define C (1.0-R)
#define NMAX 100000
#define pi 3.141592654
#define GRAVITY 9.81
#define LITHIC_DIAMETER_THRESHOLD 7.0
#define PUMICE_DIAMETER_THRESHOLD -1.0

/* Value to multiply by degrees to get radians 
	 #define DEG2RAD (double)(M_PI/180.0)*/
#define DEG2RAD   .01745329251994329576

/* air density at sea level in kg/m3 */
#define AIR_DENSITY 1.293
   
/* dynamic viscosity of air */
#define AIR_VISCOSITY 0.000018325

/* The maximum line length */
#define MAX_LINE 200

#define POINTS_OUT "tephra.out"
#define GRAIN_SZ_MODEL "model.out"
#define WIND_MODEL "wind_levels.out"
#define LO_PARAM(p) (double)_LO[(p)]
#define HI_PARAM(p) (double)_HI[(p)]

extern int NUM_OF_PARAMS;
extern int NUM_OF_VERTICES;
extern double _LO[];
extern double _HI[];

/* This variable determines whether to perform the inversion
    for wind speed and direction or to use a supplied wind
    file */
extern int FIXED_WIND;

/* These are now defined in config file. 
extern int PLUME_MODEL; */
extern int PART_STEPS;
extern int COL_STEPS;
extern double TOLERANCE;

/*density model for the pyroclasts 
extern double LITHIC_DENSITY;
extern double PUMICE_DENSITY;
*/

/* wind data read in at 0.5km intervals - this can be changed to suit available data
 * make sure input file intervals agree 
extern double WIND_INTERVAL;
extern int WIND_LEVELS;
extern int WIND_DAYS;
*/
/* mixed diffusion model*/
/* eddy diff for small particles in m2/s (400 cm2/s) */
/* extern double EDDY_CONST; */

/* diffusion coeff for large particles (m2/s) */
/* extern double DIFFUSION_COEFFICIENT; */ 

/* threshold for change in diffusion (seconds fall time) */
/* extern double FALL_TIME_THRESHOLD; */

/* extern double PLUME_RATIO; Hb/Ht[area_of_release] of the laterally spreading cloud */
/* the parameters being modeled 

   1) column height [1000 - 30000 meters]
   2) plume ratio (area_of_release) [Hb/Ht of the laterally spreading cloud]
   3) diffusion coef. [10 to 50,000 m^2/sec]
   4) total mass ejected [10^6 to 10^12 kg]
   5)  median particle size [-3 to 5 phi]
   6) std. deviation in particle size (sigma) [1 to 5 phi]
   7) fall time threshhold [18 sec to 18,000 sec]
   8) eddy constant 
   9 to ..) wind speed [0 to 50 m/sec]
   .. to last param) wind direction [0 to 360 degrees]
*/
enum{FALSE, TRUE};
enum{DISTR_1, DISTR_2};
enum {CHI2, RMSE, LOG_TEST};
enum {MAX_COL_HT, ALPHAP, BETAP, DIFF_COEF, TOTAL_MASS, MED_PHI, STD_DEV_PHI, FALLTIME_THRESH, EDDY_COEF, WIND_SPEED, WIND_DIRECTION, LAST_PARAM}; 



