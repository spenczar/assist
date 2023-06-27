/**
 * A simple benchmark to measure the speed of the integration (no interpolation is done).
 */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include "rebound.h"
#include "assist.h"

double runtime_holman(int const num_bodies, int const num_timesteps, double const timestep){
    struct assist_ephem* ephem = assist_ephem_create(
	    "../../data/linux_p1550p2650.440",
            "../../data/sb441-n16.bsp");
    if (ephem == NULL){
        fprintf(stderr,"Error initializing assist_ephem.\n");
        exit(1);
    }
    struct reb_simulation* r = reb_create_simulation();
    struct assist_extras* ax = assist_attach(r, ephem);

    r->t = 8416.5;

    // Asteroid Holman with slightly different initial conditions
    for (int i=0; i<num_bodies; i++){
        reb_add_fmt(r, "x y z vx vy vz",
                -2.724183384883979E+00+i/1e-10, -3.523994546329214E-02, 9.036596202793466E-02, 
                -1.374545432301129E-04, -1.027075301472321E-02, -4.195690627695180E-03); 
    }
   
    struct timeval time_beginning;
    struct timeval time_end;

    
    gettimeofday(&time_beginning,NULL);
    for (int i=0; i<num_timesteps; i++){
      enum REB_STATUS status = reb_integrate(r, r->t + timestep);
      assert(status == REB_EXIT_SUCCESS);
    }
    gettimeofday(&time_end,NULL);
    
    assist_free(ax);
    assist_ephem_free(ephem);
    reb_free_simulation(r);

    return (time_end.tv_sec-time_beginning.tv_sec)*1e6+(time_end.tv_usec-time_beginning.tv_usec);
}

void remove_outliers(double* const values, const int n, double* const result, int *const result_n) {
  // Remove values which are at least 5 sigma away from the mean
  double mean = 0;
  double M2 = 0;
  for (int i=0; i<n; i++){
    double r = values[i];
    double delta = r - mean;
    mean += delta/(i+1);
    M2 += delta*(r-mean);
  }
  double std = sqrt(M2 / (n -1));

  int j = 0;
  for (int i=0; i<n; i++){
    if (fabs(values[i]-mean) < 5*std){
      result[j] = values[i];
      j++;
    }
  }
  *result_n = j;
}

void mean_runtime(int const num_bodies, int const num_timesteps, double const timestep, int const N_samples, double* mean_ptr, double* std_ptr){
    double* runtimes = malloc(sizeof(double)*N_samples);
    for (int i=0; i<N_samples; i++){
      runtimes[i] = runtime_holman(num_bodies, num_timesteps, timestep);
    }

    double *runtimes_trimmed = malloc(sizeof(double)*(N_samples));
    int N_samples_trimmed = N_samples;
    remove_outliers(runtimes, N_samples, runtimes_trimmed, &N_samples_trimmed);
    free(runtimes);

    // Calculate mean and std
    double mean = 0;
    double M2 = 0;
    for (int i=0; i<N_samples-2; i++){
        double r = runtimes_trimmed[i];
        double delta = r - mean;
        mean += delta/(i+1);
        M2 += delta*(r-mean);
    }
    double std = sqrt(M2 / (N_samples -1));
    free(runtimes_trimmed);
    *mean_ptr = mean;
    *std_ptr = std;
}

int main(int argc, char* argv[]) {
  double mean, std;

  printf("one-particle");
  mean_runtime(1, 1, 1, 100, &mean, &std);
  printf("\t%.0f\t%.0f\n", mean, std);

  printf("many-particles");
  mean_runtime(1000, 1, 1, 100, &mean, &std);  
  printf("\t%.0f\t%.0f\n", mean, std);

  printf("many-small-timesteps");
  mean_runtime(1, 8640, 1.0/8640.0, 10, &mean, &std);
  printf("\t%.0f\t%.0f\n", mean, std);

  printf("few-big-timesteps");
  mean_runtime(1, 1000, 1.0, 100, &mean, &std);
  printf("\t%.0f\t%.0f\n", mean, std);  
}
