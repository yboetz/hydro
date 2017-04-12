/* Function to sync boundaries accross processes */

#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <mpi.h>

#include "parametres.h"
#include "mpi_sync.h"
void
mpi_sync(const hydroparam_t H, hydrovar_t * Hv)
{
    int j, nv;

    // Get number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get rank of processes
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    
    // Start exchange of boundary conditions
    int dest, tag; // Destination rank for boundary exchange

    int buffer_size = 2 * H.ny * H.nvar;
    double* buffer_a = (double*)malloc(sizeof(double) * buffer_size); // Buffer for boundary data
    double* buffer_b = (double*)malloc(sizeof(double) * buffer_size);

    /* Initialize buffers with correct values. For even ranks 
    buffer_a is left boundary and buffer_b right, for odd ranks
    vice versa */
    for(nv = 0; nv < H.nvar; nv++)
    {
        #define IHv(i, j, v) ((i) + (H.nxt * (H.nyt * (v)+ (j))))
        #define IHbuff(i, j, v) ((i) + (2 * (H.ny * (v)+ (j-ExtraLayer))))
        for(j = H.jmin + ExtraLayer; j < H.jmax - ExtraLayer; j++)
            {
                if(rank % 2 == 0)
                {   // Left boundary
                    buffer_a[IHbuff(0,j,nv)] = Hv->uold[IHv(2,j,nv)];
                    buffer_a[IHbuff(1,j,nv)] = Hv->uold[IHv(3,j,nv)];
                    // Right boundary
                    buffer_b[IHbuff(0,j,nv)] = Hv->uold[IHv(H.imax-4,j,nv)];
                    buffer_b[IHbuff(1,j,nv)] = Hv->uold[IHv(H.imax-3,j,nv)];
                }
                else
                {
                    // Left boundary
                    buffer_b[IHbuff(0,j,nv)] = Hv->uold[IHv(2,j,nv)];
                    buffer_b[IHbuff(1,j,nv)] = Hv->uold[IHv(3,j,nv)];
                    // Right boundary
                    buffer_a[IHbuff(0,j,nv)] = Hv->uold[IHv(H.imax-4,j,nv)];
                    buffer_a[IHbuff(1,j,nv)] = Hv->uold[IHv(H.imax-3,j,nv)];
                }
            } // End for j
        #undef IHbuff
        #undef IHv
    } // End for nv


    // Send buffer_b between 0 & 1, 2 & 3, ...
    if(rank < world_size - world_size % 2)
    {
        tag = rank - rank % 2;
        if(rank % 2 == 0)   dest = rank + 1;
        else                dest = rank - 1;
        MPI_Sendrecv_replace(buffer_b, buffer_size, MPI_DOUBLE, dest, tag, dest, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    //MPI_Barrier(MPI_COMM_WORLD);

    // Send buffer_a between 1 & 2, 3 & 4, ... 
    if(0 < rank && rank < world_size - 1 + world_size % 2)
    {
        tag = rank - (rank+1) % 2;
        if(rank % 2 == 0)   dest = rank - 1;
        else                dest = rank + 1;
        MPI_Sendrecv_replace(buffer_a, buffer_size, MPI_DOUBLE, dest, tag, dest, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    //MPI_Barrier(MPI_COMM_WORLD);


    /* Write received data back to grid.*/
    for(nv = 0; nv < H.nvar; nv++)
    {
        #define IHv(i, j, v) ((i) + (H.nxt * (H.nyt * (v)+ (j))))
        #define IHbuff(i, j, v) ((i) + (2 * (H.ny * (v)+ (j-ExtraLayer))))
        for(j = H.jmin + ExtraLayer; j < H.jmax - ExtraLayer; j++)
            {
                if(rank % 2 == 0)
                    {
                    // Left boundary
                    Hv->uold[IHv(0,j,nv)] = buffer_a[IHbuff(0,j,nv)];
                    Hv->uold[IHv(1,j,nv)] = buffer_a[IHbuff(1,j,nv)];
                    // Right boundary
                    Hv->uold[IHv(H.imax-2,j,nv)] = buffer_b[IHbuff(0,j,nv)];
                    Hv->uold[IHv(H.imax-1,j,nv)] = buffer_b[IHbuff(1,j,nv)];
                    }
                else
                    {
                    // Left boundary
                    Hv->uold[IHv(0,j,nv)] = buffer_b[IHbuff(0,j,nv)];
                    Hv->uold[IHv(1,j,nv)] = buffer_b[IHbuff(1,j,nv)];
                    // Right boundary
                    Hv->uold[IHv(H.imax-2,j,nv)] = buffer_a[IHbuff(0,j,nv)];
                    Hv->uold[IHv(H.imax-1,j,nv)] = buffer_a[IHbuff(1,j,nv)];
                    }
            } // End for j
        #undef IHbuff
        #undef IHv
    } // End for nv


    free(buffer_a); // Free up workspace
    free(buffer_b);
    MPI_Barrier(MPI_COMM_WORLD);
}