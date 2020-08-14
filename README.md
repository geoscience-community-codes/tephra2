

### OVERVIEW
tephra2 (VERSION 2.0) is a tephra dispersion simulation tool, written in C, used to estimate the mass of tephra that will accumulate at a site or over a region, given explosive eruption conditions. The tephra2 model simplifies eruption conditions. The code calculates an expected tephra accumulation for a given total eruption mass, eruption column height, and particle size distribution. Wind velocity varies as a function of height but remains constant within a given stratum. Wind velocities are sampled from the NOAA REANALYSIS database at the location of the volcano. Tephra2 outputs the mass per unit area (kg per meter squared) and the weight percent of individual grain sizes at specified locations. The output is presented as a text table of data that can be used to create an isomass map, to visualize expected accumulation over a region.

Two executables can be compiled, tephra2_2020 (the forward model) and tephra2-inversion_2020 (the inversion model). The inversion model requires MPI (Message Passing Interface) libraries and should be run on a compute cluster with multiple compute nodes. 

### COMPILING
tephra2 is written in C, a compiled language, and must be compiled before it can be executed. 

### INSTALL THESE DEPENDENCIES BEFORE COMPILING
- gcc
- openmpi
- gc, gc-devel, libatomic_ops(opensuse) or libgc-dev, libgc1c2(ubuntu)

#### WIND DATA
[How to download NOAA REANALYSIS data](plotting_scripts/readme.wind)

### USAGE
Both the forward model and the inversion model require a configuration file. The configuration file is arranged as a list of keywords followed by their corresponding value(s). The program user can change any of the values; the keywords must not be changed.

To run the forward model, at the command line, type:
```
tephra2_2020 tephra2.conf grid_file wind_file > tephra2.out
```
where,
- **tephra2_2020** is the name of the executable
- **tephra2.conf** is the name of the file of configuration parameters
- **grid_file** is is a text file of 3 columns separated by spaces following the format: 
    ```
    Easting(m)  Northing(m)  Grid-Elevation(m)
    ```  
- **wind file** is a 3-column text file of wind data following the format:
   ```
   Height(masl)  Wind-Speed(m/s)  Wind-Direction(wind vector azimuth in degrees)
   ```
- **tephra2.out** is the output file name where the tephra accumulation values will be written, following the format:
    ```
    Easting(m)  Northing(m)  Grid-Elevation(m)  Mass(kg/m^2) [weight percent of modeled phi fractions]
    ```


### ADDITIONAL DEPENDENCIES


#### <I>perl</I> DEPENDENCIES


#### <i>gmt</i> PLOTTING DEPENDENCIES











=========================================================================
###### versions --
01-27-2018
Formally versioning this version of tephra2 as Version 2.0, to distinguish it from other, earlier versions.

03-12-2014
Initial upload of the tephra2 probability scripts (perl scripts) which use the tephra2 model of ash dispersion to create survivor curves for specific locations around a volcano.

07-09-2013
Initial upload to GitHub. tephra2 is a tephra dispersion model based on an analytical (closed-form) solution of the advection-diffusion equation. This initial version uses a plume model based on the beta function and incorporates a 2-D wind model that varies in speed and direction with changes in elevation. 
