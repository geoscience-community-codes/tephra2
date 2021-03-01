/* common_structures.h */

typedef struct {
  /* geographic location of a point where ash mass in estimated */
  double easting; /*utm coordinate in meters */
  double northing; /*utm coordinate in meters */
  double elevation; /*elevation at the point in masl */
  double observed_mass;  /*mass of material accumulated at the point (gm) */
  float confidence; /* value [0 - 1] (i.e. [least - greatest]) indicating confidence of measurement */
  double calculated_mass; /*calculated mass of material at the point (gm) */
  double observed_phi[100]; /* pointer to the grainsize distribution, 100 possible bins, from collected samples */
  double calculated_phi[100]; /*pointer to a grainsize distribution of 100 possible bins or grain sizes */
  /* double std_dev_phi; the graphic standard deviation in grainsize diameter */
  /* double skewness_phi; the graphic skewness in grainsize diameter */
  /* double median_phi; the mean diameter in the gainsize distribution */
} POINT;

/* These values are calculated ahead of time and used in the double integration [PART_STEPS][COL_STEPS]*/
typedef struct {

  double particle_ht;
  double ashdiam;
  double part_density;
  double fall_time;
  double plume_diffusion_fine_particle;
  double plume_diffusion_coarse_particle;
  double total_fall_time;
  double wind_sum_x;
  double wind_sum_y;
  double demon1;
} TABLE;

typedef struct {
  /* eruption parameters */
  int plume_model;
  double vent_northing; /*volcano location in UTM north (meters) */
  double vent_easting; /*volcano location in UTM east (meters) */
  double total_ash_mass; /* is the total amount of ash erupted (kg) */
  double min_phi;  /*the maximum particle diameter considered (phi)*/
  double max_phi; /* is the minimum particle diameter considered (phi) */
   
  /*Note: erupt->max/min_part_size are used to set the limits of integration
   * on the calculation. Particles outside this range are not considered at all */

  /*Note: phi uits are such that max will appear to be less than min, this is
    accounted for in the conversion to cm, which is internal */

  double mean_phi;  /*the mean particle diameter erupted in phi units */
  double sigma_phi; /*standard deviation in particle diameter in phi units */
  double vent_elevation;  /* elevation of the vent (meters above sea level)*/
  double max_plume_elevation;  /*eruption column height (meters above sea level)*/
 
  double plume_ratio; /* parameter governing the particle size distribution in column */
  double pumice_density;
  double lithic_density;
  /* Note: A large value  of beta (1) places most of the particles
   * high in the eruption column, a low value of beta (0.01) spreads the particle density
   * lower in the column. Particle release models based on "corner" models etc strongly
   * suggest that a larger value for beta should be used. */
  double column_beta;
    /* These are two parameters used for calculating the PDF for the beta distribution 
    * Their can be any value greater than zero. When alpha and beta both equal 1.0, the
    * resulting distribution is the uniform normal distribution.
    */
    double alpha;
    double beta;
  double diffusion_coefficient; /* diffusion coeff for large particles (m2/s)  */
  double fall_time_threshold; /* threshold for change in diffusion (seconds fall time) */
  double eddy_constant;

  /* double duration; => duration of the sustained phase of eruption; this is calculated from the empirical power-law eq. 
                      total_erupted_mass = density of fallout deposit * ( plume_height/1670)^4 * duration */ 
  double (*pdf)(double, double, double, double, double);
  double (*fit)(FILE *, int, POINT *);
} ERUPTION;

typedef struct {
    double wind_height; /* height a.s.l. in km */
    double windspeed; 	/* the average windspeed in m/s */
    double wind_direction;  	/* average wind direction in +/- degrees from north */
} WIND;

typedef double (*PFI)(double, double, double);

typedef struct {
  double min_falltime;
  double max_falltime;
} STATS;
