/* prototypes.h */

#include "parameters.h"
#include "common_structures.h"

// extern  double (*fit)(FILE *, int, POINT *);
// extern double (*pdf)(double, double, double, double, double);
void optimize_params(double [][NUM_OF_PARAMS], double [], double, double (*funk)(double []), int *);
double minimizing_func(double []);
void test_bounds(int, double *, double);
void init_optimal_params( double [][NUM_OF_PARAMS]);
void assign_new_params( double []);
int init_globals(char *);
void set_global_values(FILE *);
/*void set_eruption_values(ERUPTION *erupt);*/
void set_eruption_values(ERUPTION *, WIND *);
int get_points(FILE *);
void printout_model(double []);
void printout_points(void);
void printout_parameters(double);
void print_for_stats(double);
double set_variance(FILE *, int, POINT *);
double chi_squared(FILE *, int, POINT *);
double rmse(FILE *, int, POINT *);
double log_test(FILE *, int, POINT *);
// double (*fit)(FILE *, int, POINT *);

double phi2m( double);
double particle_density (double, double, double);
double upward_part_vel( ERUPTION *, double);
double pdf_ash_diff (ERUPTION *, double, double, double);
double part_diff_time( double);
double part_fall_time(double, double, double, double);
double part_term_vel(double, double);
//double (*pdf)(double, double, double, double, double);
double plume_pdf0(double, double, double, double, double);
double plume_pdf1(double, double, double, double, double);
double plume_pdf2(double, double, double, double, double);
double pdf_grainsize(double, double, double);
double ash_accum(ERUPTION *, WIND *, double, double, double, double, double, double, double);
/*double strat_average(WIND *meteor, double col_ht, double xspace, double yspace, double t1, double t2, double sigma);

void  tephra_calc(ERUPTION *erupt, POINT *pt, WIND *level);*/
double strat_average( double, double, double, double, double,double);
void tephra_calc(ERUPTION *, POINT *, WIND *, STATS *);
double part_fall_time_vent2ground (double, double, double, double);

void slave(int my_rank, FILE *log_file);
double master(FILE *log_file);
void set_LOG(FILE *log_file);
void close_logfile(void);
void _free(void);

int get_eruptions(void);
int get_wind(FILE *);
