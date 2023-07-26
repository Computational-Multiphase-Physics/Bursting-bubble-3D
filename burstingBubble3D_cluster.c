/*
First version: Anil 
Second (Current) version: Vatsal (contact@vatsalsanjay.com)

Tags:
Date: 26th July 2023, 9:00h AMS time.
*/

#include "grid/octree.h"
#include "navier-stokes/centered.h"
#define FILTERED
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "distance.h"
#include "adapt_wavelet_limited_v2.h"


#define fErr (1e-3) // error tolerance in VOF
#define KErr (1e-6) // error tolerance in curvature 
#define VelErr (1e-3) // error tolerances in velocity

#define tsnap (0.01) // timestep used to save simulation snapshot

char nameOut[80], dumpFile[80];

int MAXlevel;
double Ldomain, tmax, ShiftDist;

uf.n[bottom] = 0.0;
uf.n[left] = 0.0;

u.t[bottom] = neumann(0.0);
u.t[left] = neumann(0.0);

u.r[bottom] = neumann(0.0);
u.r[left] = neumann(0.0);

int main(int argc, char *argv[])
{
    MAXlevel = atoi(argv[1]);
    Ldomain = atof(argv[2]);
    ShiftDist = atof(argv[3]);

    tmax = 1.0;
    
    f.sigma = 1.0; // this is the liquid-gas surface tension coefficient

    rho1 = 1.; // this is the density of the liquid
    rho2 = 0.001; // this is the gas-liquid density ratio
    mu1 = 0.01; // this is the liquid Ohnesorge number
    mu2 = 1e-4; // this is the gas Ohnesorge number

    size(Ldomain);
    origin(0.0, 0.0, -(Ldomain-3.0));
    init_grid(256);

    char comm[80];
    sprintf (comm, "mkdir -p intermediate");
    system(comm);
    // Name of the restart file. See writingFiles event.
    sprintf (dumpFile, "RunningRestartFile");

    run();
}

// For region based refinement in adapt_wavelet_limited
int refRegion(double x, double y, double z){
  return (sq(x)+sq(y)<sq(0.1)? MAXlevel+1: sq(x)+sq(y)<sq(1.1)? MAXlevel: sq(x)+sq(y)<sq(2.5)? MAXlevel-1: MAXlevel-2);
}

event init(t = 0){
  if(!restore (file = dumpFile)){
    char filename[60];
    sprintf(filename,"dump-%4.3f", ShiftDist);
    FILE * fp = fopen (filename, "r");
    if (fp == NULL){
      fprintf(ferr, "There is no file named %s\n", filename);
      return 1;
    }
    restore (file = filename);
  }
}

// adaptive meshing
scalar KAPPA[];
event adapt(i++)
{
  curvature(f, KAPPA);
  adapt_wavelet_limited ((scalar *){f, u.x, u.y, KAPPA},
     (double[]){fErr, VelErr, VelErr, KErr},
     refRegion);
}

// for dump output, it will save files in the intermediate folder with name snapshot-<timestamp>
event writingFiles (t = 0; t += tsnap; t <= tmax) {
  dump (file = dumpFile);
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump(file=nameOut);
}

// end of simulation
event end (t = end) {
  fprintf(ferr, "Done: Level %d\n", MAXlevel);
}

// writing log files every timestep
event logfile2(i++)
{
    double ke = 0.; // good metric to track if things blow up
    foreach (reduction(+:ke)){
        ke += (2*pi*y)*(0.5*(f[])*(sq(u.x[]) + sq(u.y[])))*sq(Delta);
    }
    if (pid() == 0)
    {
        FILE *ftr = (i == 0) ? fopen("track.out", "w") : fopen("track.out", "a");
        fprintf(ftr, "\tt: %-20g dt: %-20g ke: %-20g\n", t, dt, ke);
        fprintf(ferr, "\tt: %-20g dt: %-20g ke: %-20g \n", t, dt, ke);
        fclose(ftr);
    }

    if (ke > 1e3 || ke < 1e-6){
        if (i > 1e2){
            fprintf(ferr, "KE is too high or too low. Exiting...\n");
            return 1;
        }
    }
}

