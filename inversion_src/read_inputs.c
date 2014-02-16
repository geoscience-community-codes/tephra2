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
  double wind_height, wind_dir, windspeed;

	fprintf(log_file,"ENTER[get_wind].\n");
	
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
  i=0;
  /* start at the vent */ 
    level = e.vent_elevation;
  
  fprintf(log_file, "\tRead %d wind levels.\n", i+1, j);
  fprintf(log_file, "EXIT[get_wind].\n");	  	  
  return 0;
}
