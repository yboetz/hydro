/*
  A simple 2D hydro code
  (C) Romain Teyssier : CEA/IRFU           -- original F90 code
  (C) Pierre-Francois Lavallee : IDRIS      -- original F90 code
  (C) Guillaume Colin de Verdiere : CEA/DAM -- for the C version
*/
#include <stdio.h>
#include <time.h>
#include <mpi.h>
#include <unistd.h>

#include "parametres.h"
#include "hydro_funcs.h"
#include "vtkfile.h"
#include "compute_deltat.h"
#include "hydro_godunov.h"
#include "utils.h"

hydroparam_t H;
hydrovar_t Hv;                  // nvar
hydrovarwork_t Hvw;             // nvar
hydrowork_t Hw;
unsigned long flops = 0;

int
main(int argc, char **argv)
{

  MPI_Init(NULL,NULL);

  // Get number of processes
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // Get rank of processes
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Get name of the processor
  char hostname[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(hostname, &name_len);

  // Print off a hello world message
  printf("Hello world from host %s, rank %d"
          " out of %d processors\n",
          hostname, rank, world_size);
  

  int value;
  int tag = 100;
  MPI_Status status;

  /*
  // Simple communication
  if(rank == 0)
  {
    value = 1000;
    MPI_Send(&value, 1, MPI_INT, 3, tag, MPI_COMM_WORLD);
  }
  else if(rank == 3)
  {
    MPI_Recv(&value, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
    printf("Rank %d received value %d from rank %d\n", 3,value,0);
  }
  */

  
  // Send value around processes
  int n = 0;
  value = 0;
  while(n < 100)
  {
    if(rank == (n+1) % 4)
    {
      MPI_Recv(&value, 1, MPI_INT, n % 4, tag, MPI_COMM_WORLD, &status);
      printf("Rank %d got value %d from rank %d.\n", rank, value, n % 4);
      value += 1;
      usleep(100000);
    }
    else if(rank == n % 4)
    {
      MPI_Send(&value, 1, MPI_INT, (n+1)%4, tag, MPI_COMM_WORLD);
      printf("Rank %d sent value %d to rank %d.\n", rank, value, (n+1)%4);
    }
    n++;
  }
  
  if(rank == 5) // Only run hydro code on one process for now
  {

  int nb_th=1;
  double dt = 0;
  long nvtk = 0;
  char outnum[80];
  long time_output = 0;
  
  // double output_time = 0.0;
  double next_output_time = 0;
  double start_time = 0, end_time = 0;
  double start_iter = 0, end_iter = 0;
  double elaps = 0;

  start_time = cclock();
  process_args(argc, argv, &H);
  hydro_init(&H, &Hv);
  PRINTUOLD(H, &Hv);
  
  printf("Hydro starts - mpi version \n");

  // vtkfile(nvtk, H, &Hv);
  if (H.dtoutput > 0) 
    {	
      // outputs are in physical time not in time steps
      time_output = 1;
      next_output_time = next_output_time + H.dtoutput;
    }

  while ((H.t < H.tend) && (H.nstep < H.nstepmax)) 
    {	
      start_iter = cclock();
      outnum[0] = 0;
      flops = 0;
      if ((H.nstep % 2) == 0) 
	{
	  compute_deltat(&dt, H, &Hw, &Hv, &Hvw);
	  if (H.nstep == 0) {
	    dt = dt / 2.0;
	  }
	}
      
      if ((H.nstep % 2) == 0) {
	hydro_godunov(1, dt, H, &Hv, &Hw, &Hvw);
	hydro_godunov(2, dt, H, &Hv, &Hw, &Hvw); 
      } else {
	hydro_godunov(2, dt, H, &Hv, &Hw, &Hvw);
	hydro_godunov(1, dt, H, &Hv, &Hw, &Hvw); 
      }
      
      end_iter = cclock();
      H.nstep++;
      H.t += dt;
      
      if (flops > 0) {
	double iter_time = (double) (end_iter - start_iter);
	if (iter_time > 1.e-9) {
	  double mflops = (double) flops / (double) 1.e+6 / iter_time;
	  sprintf(outnum, "%s {%.3f Mflops} (%.3fs)", outnum, mflops, iter_time);
	}
      } else {
	double iter_time = (double) (end_iter - start_iter);
	sprintf(outnum, "%s (%.3fs)", outnum, iter_time);
      }
      if (time_output == 0) {
	if ((H.nstep % H.noutput) == 0) {
	  vtkfile(++nvtk, H, &Hv);
	  sprintf(outnum, "%s [%04ld]", outnum, nvtk);
	}
      } else {
	if (H.t >= next_output_time) {
	  vtkfile(++nvtk, H, &Hv);
	  next_output_time = next_output_time + H.dtoutput;
	  sprintf(outnum, "%s [%04ld]", outnum, nvtk);
	}
      }
	fprintf(stdout, "--> step=%-4ld %12.5e, %10.5e %s\n", H.nstep, H.t, dt, outnum);
    }   // end while loop
  hydro_finish(H, &Hv);
  end_time = cclock();
  elaps = (double) (end_time - start_time);
  timeToString(outnum, elaps); 
  
  fprintf(stdout, "Hydro ends in %ss (%.3lf).\n", outnum, elaps);

  } // End of hydro code

  MPI_Finalize();

return 0;
}
    
