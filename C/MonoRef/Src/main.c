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
  // Initialize MPI
  MPI_Init(NULL,NULL);

  // Get number of processes
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // Get rank of processes
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  /*
  // Print off a hello world message
  printf("Hello world from rank %d"
          " out of %d processors\n",
          rank, world_size);
  
  MPI_Barrier(MPI_COMM_WORLD);
  

  int value;
  int tag = 100;
  MPI_Status status;

  // MPI test
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
  
  // Send value around processes
  int n = 0;
  value = 0;
  while(n < 8)
  {
    if(rank == (n+1) % world_size)
    {
      MPI_Recv(&value, 1, MPI_INT, n % world_size, tag, MPI_COMM_WORLD, &status);
      printf("Rank %d got value %d from rank %d.\n\n", rank, value, n % world_size);
      value += 1;
      usleep(50000);
    }
    else if(rank == n % world_size)
    {
      printf("Rank %d sent value %d to rank %d.\n", rank, value, (n+1) % world_size);
      MPI_Send(&value, 1, MPI_INT, (n+1) % world_size, tag, MPI_COMM_WORLD);
    }
    n++;
  }
  */

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
  
  if(rank == 0) printf("Hydro starts - MPI version \n");
  
  // vtkfile(nvtk, H, &Hv);
  if (H.dtoutput > 0) 
    {	
      // outputs are in physical time not in time steps
      time_output = 1;
      next_output_time = next_output_time + H.dtoutput;
    }

  MPI_Barrier(MPI_COMM_WORLD);
  
  while ((H.t < H.tend) && (H.nstep < H.nstepmax)) 
    {	
      start_iter = cclock();
      outnum[0] = 0;
      flops = 0;

      // Start exchange of boundary conditions
      int dest;       // Destination rank for boundary exchange
      int tag = 111;  // Random tag

      int buffer_size = 2 * H.ny * H.nvar;
      double* buffer_a = (double*)malloc(sizeof(double) * buffer_size); // Buffer for boundary data on left
      double* buffer_b = (double*)malloc(sizeof(double) * buffer_size); // and right site

      /* Initialize buffers with correct values. For even ranks 
      buffer_a is left boundary and buffer_b right, for odd ranks
      vice versa */
      for(int nv = 0; nv < H.nvar; nv++)
      {
        #define IHv(i, j, v) ((i) + (H.nxt * (H.nyt * (v)+ (j))))
        
        for(int j = H.jmin + ExtraLayer; j < H.jmax - ExtraLayer; j++)
          {
            // Left boundary
            for(int i = H.imin + ExtraLayer; i < H.imin + 2*ExtraLayer; i++)
            {
              #define IHbuff(i, j, v) ((i-H.imin-ExtraLayer) + (2 * (H.ny * (v)+ (j-ExtraLayer))))
              if(rank % 2 == 0) 
                buffer_a[IHbuff(i,j,nv)] = Hv.uold[IHv(i,j,nv)];
              else
                buffer_b[IHbuff(i,j,nv)] = Hv.uold[IHv(i,j,nv)];
              #undef IHbuff
            }
            // Right boundary
            for(int i = H.imax - 2*ExtraLayer; i < H.imax - ExtraLayer; i++)
            {
              #define IHbuff(i, j, v) ((i-H.imax+2*ExtraLayer) + (2 * (H.ny * (v)+ (j-ExtraLayer))))
              if(rank % 2 == 1) 
                buffer_a[IHbuff(i,j,nv)] = Hv.uold[IHv(i,j,nv)];
              else
                buffer_b[IHbuff(i,j,nv)] = Hv.uold[IHv(i,j,nv)];
              #undef IHbuff
            }
          } // End for j
          #undef IHv
      } // End for nv


      // Send buffer_b between 0 & 1, 2 & 3, ...
      if(rank < world_size - world_size % 2)
      {
        if(rank % 2 == 0) dest = rank + 1;
        if(rank % 2 == 1) dest = rank - 1;
        
        MPI_Sendrecv_replace(buffer_b, buffer_size, MPI_DOUBLE, dest, tag, dest, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //printf("rank %d got buffer_b of size %d from rank %d\n",rank,buffer_size,dest);
      }
      MPI_Barrier(MPI_COMM_WORLD);

      // Send buffer_a between 1 & 2, 3 & 4, ... 
      if(0 < rank && rank < world_size - 1 + world_size % 2)
      {
        if(rank % 2 == 0) dest = rank - 1;
        if(rank % 2 == 1) dest = rank + 1;

        MPI_Sendrecv_replace(buffer_a, buffer_size, MPI_DOUBLE, dest, tag, dest, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //printf("rank %d got buffer_a of size %d from rank %d\n",rank,buffer_size,dest);
      }
      MPI_Barrier(MPI_COMM_WORLD);


      /* Write received data back to grid.*/
      for(int nv = 0; nv < H.nvar; nv++)
      {
        #define IHv(i, j, v) ((i) + (H.nxt * (H.nyt * (v)+ (j))))

        for(int j = H.jmin + ExtraLayer; j < H.jmax - ExtraLayer; j++)
          {
            if(rank > 0) // Don't write left boundary on rank 0
            {
              // Left boundary
              for(int i = H.imin; i < H.imin + ExtraLayer; i++)
              {
                #define IHbuff(i, j, v) ((i-H.imin) + (2 * (H.ny * (v)+ (j-ExtraLayer))))
                if(rank % 2 == 0) 
                  Hv.uold[IHv(i,j,nv)] = buffer_a[IHbuff(i,j,nv)];
                else
                  Hv.uold[IHv(i,j,nv)] = buffer_b[IHbuff(i,j,nv)];
                #undef IHbuff
              }
            }

            if(rank < world_size-1) // Don't write right boundary on rank world_size-1
            {
              // Right boundary
              for(int i = H.imax - ExtraLayer; i < H.imax; i++)
              {
                #define IHbuff(i, j, v) ((i-H.imax+ExtraLayer) + (2 * (H.ny * (v)+ (j-ExtraLayer))))
                if(rank % 2 == 1) 
                  Hv.uold[IHv(i,j,nv)] = buffer_a[IHbuff(i,j,nv)];
                else
                  Hv.uold[IHv(i,j,nv)] = buffer_b[IHbuff(i,j,nv)];
                #undef IHbuff
              }
            }

          } // End for j

          #undef IHv
      } // End for nv
      

      free(buffer_a); // Free up workspace
      free(buffer_b);
      MPI_Barrier(MPI_COMM_WORLD);
      //break; // Remove at some point


      if ((H.nstep % 2) == 0) 
      {
        compute_deltat(&dt, H, &Hw, &Hv, &Hvw);
        if (H.nstep == 0) dt = dt / 2.0;
        // Get global minima of all dt and broadcast it
        MPI_Allreduce(&dt, &dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      }
      
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
      
      end_iter = cclock();
      H.nstep++;
      H.t += dt;

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
          vtkfile(++nvtk, H, &Hv);
          sprintf(outnum, "%s [%04ld]", outnum, nvtk);
        }
      }
      else
      {
        if(H.t >= next_output_time)
        {
          vtkfile(++nvtk, H, &Hv);
          next_output_time = next_output_time + H.dtoutput;
          sprintf(outnum, "%s [%04ld]", outnum, nvtk);
        }
      }
      if(rank == 0) fprintf(stdout, "--> step=%-4ld %12.5e, %10.5e %s\n", H.nstep, H.t, dt, outnum);

    } // End while loop

  hydro_finish(H, &Hv);
  end_time = cclock();
  elaps = (double) (end_time - start_time);
  timeToString(outnum, elaps); 
  
  if(rank == 0) fprintf(stdout, "Hydro ends in %ss (%.3lf).\n", outnum, elaps);

  MPI_Finalize();

return 0;
}
    
