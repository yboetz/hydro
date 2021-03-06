/*
  A simple 2D hydro code
  (C) Romain Teyssier : CEA/IRFU           -- original F90 code
  (C) Pierre-Francois Lavallee : IDRIS      -- original F90 code
  (C) Guillaume Colin de Verdiere : CEA/DAM -- for the C version
*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <omp.h>
#include <unistd.h>

#include "parametres.h"
#include "hydro_funcs.h"
#include "vtkfile.h"
#include "compute_deltat.h"
#include "hydro_godunov.h"
#include "mpi_sync.h"
#include "utils.h"

hydroparam_t H;
hydrovar_t Hv;                  // nvar
hydrovarwork_t Hvw;             // nvar
hydrowork_t Hw;
unsigned long flops = 0;


int
main(int argc, char **argv)
{
  // Initialize MPI
  MPI_Init(NULL,NULL);

  // Get number of processes
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // Get rank of processes
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);


  //int nb_th=1;
  double dt = 0;
  long nvtk = 0;
  char outnum[80];
  long time_output = 0;
  
  // double output_time = 0.0;
  double next_output_time = 0;
  double start_time = 0, end_time = 0;
  double start_iter = 0, end_iter = 0;
  double start_dt = 0, end_dt = 0;
  double start_sync = 0, end_sync = 0;
  double start_god = 0, end_god = 0;
  double elaps_dt = 0, elaps_sync = 0, elaps_god = 0, elaps_iter = 0;
  double elaps = 0;

  start_time = cclock();
  process_args(argc, argv, &H);
  hydro_init(&H, &Hv);
  PRINTUOLD(H, &Hv);
  
  if(rank == 0)
  {
    #if defined(_OPENMP)
      // Get number of threads per process
      int nthreads = omp_get_max_threads();
      printf("Hydro starts - MPI-OpenMP version.\n");
      printf("Working on %d nodes, with %d threads each.\n\n",world_size,nthreads);
    #else
      printf("Hydro starts - MPI version.\n");
      printf("Working on %d nodes.\n\n",world_size);
    #endif
  }
  
  // vtkfile(nvtk, H, &Hv, 0b0001);
  if (H.dtoutput > 0)
    {	
      // outputs are in physical time not in time steps
      time_output = 1;
      next_output_time = next_output_time + H.dtoutput;
    }
  //usleep(50000); // This time has to be added to start_time
  MPI_Barrier(MPI_COMM_WORLD);
  

  while ((H.t < H.tend) && (H.nstep < H.nstepmax)) 
    {	
      start_iter = cclock();
      outnum[0] = 0;
      flops = 0;

      start_sync = cclock();
      mpi_sync(H, &Hv); // Sync boundaries of processes
      end_sync = cclock();

      start_dt = cclock();
      if ((H.nstep % 2) == 0) 
      {
        compute_deltat(&dt, H, &Hv); // Computed in parallel
        if (H.nstep == 0) dt = dt / 2.0;
        // Get global minima of all dt and broadcast it
        MPI_Allreduce(&dt, &dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      }
      end_dt = cclock();

      start_god = cclock();
// Start omp parallel region
#pragma omp parallel private(Hw, Hvw) // Each thread needs his own workspace
{
      if ((H.nstep % 2) == 0) 
      {
        hydro_godunov(1, dt, H, &Hv, &Hw, &Hvw);
        hydro_godunov(2, dt, H, &Hv, &Hw, &Hvw); 
      } 
      else 
      {
        hydro_godunov(2, dt, H, &Hv, &Hw, &Hvw);
        hydro_godunov(1, dt, H, &Hv, &Hw, &Hvw); 
      }
} // End omp parallel
      end_god = cclock();

      end_iter = cclock();
      H.nstep++;
      H.t += dt;

      elaps_dt += (double)(end_dt - start_dt);
      elaps_god += (double)(end_god - start_god);
      elaps_sync += (double)(end_sync - start_sync);
      MPI_Allreduce(&elaps_dt, &elaps_dt, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      MPI_Allreduce(&elaps_god, &elaps_god, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      MPI_Allreduce(&elaps_sync, &elaps_sync, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

      if(flops > 0)
      {
        double iter_time = (double) (end_iter - start_iter);
        // Get global max of iter_time and broadcast it
        MPI_Allreduce(&iter_time, &iter_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        if(iter_time > 1.e-9)
        {
          double mflops = (double) flops / (double) 1.e+6 / iter_time;
          sprintf(outnum, "%s {%.3f Mflops} (%.3fs)", outnum, mflops, iter_time);
        }
      } 
      else
      {
        double iter_time = (double) (end_iter - start_iter);
        // Get global max of iter_time and broadcast it
        MPI_Allreduce(&iter_time, &iter_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        sprintf(outnum, "%s (%.3fs)", outnum, iter_time);
      }
      if(time_output == 0)
      {
        if ((H.nstep % H.noutput) == 0)
        {
          vtkfile(++nvtk, H, &Hv, 0b0001);
          sprintf(outnum, "%s [%04ld]", outnum, nvtk);
        }
      }
      else
      {
        if(H.t >= next_output_time)
        {
          vtkfile(++nvtk, H, &Hv, 0b0001);
          next_output_time = next_output_time + H.dtoutput;
          sprintf(outnum, "%s [%04ld]", outnum, nvtk);
        }
      }
      if(rank == 0) fprintf(stdout, "--> step=%-4ld %12.5e, %10.5e %s\n", H.nstep, H.t, dt, outnum);
    } // End while loop


  hydro_finish(H, &Hv);
  end_time = cclock();
  elaps_iter = elaps_dt + elaps_god + elaps_sync;
  elaps = (double) (end_time - start_time);
  timeToString(outnum, elaps); 
  
  if(rank == 0)
  {
    fprintf(stdout, "dt = %.3lfs (%.4lf%%); mpi = %.3lfs (%.4lf%%); god = %.3lfs (%.4lf%%); total = %.3lfs.\n",
            elaps_dt, elaps_dt / elaps_iter*100, elaps_sync, elaps_sync / elaps_iter*100, elaps_god, elaps_god / elaps_iter*100, elaps_iter);
    fprintf(stdout, "Hydro ends in %ss (%.3lf).\n", outnum, elaps);
  }

  MPI_Finalize();

return 0;
}
    
