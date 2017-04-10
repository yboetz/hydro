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
    int i, j, nv;

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
    
    for(j = H.jmin + ExtraLayer; j < H.jmax - ExtraLayer; j++)
        {
            // Left boundary
            for(i = H.imin + ExtraLayer; i < H.imin + 2*ExtraLayer; i++)
            {
                #define IHbuff(i, j, v) ((i-H.imin-ExtraLayer) + (2 * (H.ny * (v)+ (j-ExtraLayer))))
                if(rank % 2 == 0) 
                    buffer_a[IHbuff(i,j,nv)] = Hv->uold[IHv(i,j,nv)];
                else
                    buffer_b[IHbuff(i,j,nv)] = Hv->uold[IHv(i,j,nv)];
                #undef IHbuff
            }
            // Right boundary
            for(i = H.imax - 2*ExtraLayer; i < H.imax - ExtraLayer; i++)
            {
                #define IHbuff(i, j, v) ((i-H.imax+2*ExtraLayer) + (2 * (H.ny * (v)+ (j-ExtraLayer))))
                if(rank % 2 == 1) 
                    buffer_a[IHbuff(i,j,nv)] = Hv->uold[IHv(i,j,nv)];
                else
                    buffer_b[IHbuff(i,j,nv)] = Hv->uold[IHv(i,j,nv)];
                #undef IHbuff
            }
        } // End for j
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

        for(j = H.jmin + ExtraLayer; j < H.jmax - ExtraLayer; j++)
            {
                if(rank > 0) // Don't write left boundary on rank 0
                {
                    // Left boundary
                    for(i = H.imin; i < H.imin + ExtraLayer; i++)
                    {
                        #define IHbuff(i, j, v) ((i-H.imin) + (2 * (H.ny * (v)+ (j-ExtraLayer))))
                        if(rank % 2 == 0) 
                            Hv->uold[IHv(i,j,nv)] = buffer_a[IHbuff(i,j,nv)];
                        else
                            Hv->uold[IHv(i,j,nv)] = buffer_b[IHbuff(i,j,nv)];
                        #undef IHbuff
                    }
                }

                if(rank < world_size-1) // Don't write right boundary on rank world_size-1
                {
                    // Right boundary
                    for(i = H.imax - ExtraLayer; i < H.imax; i++)
                    {
                        #define IHbuff(i, j, v) ((i-H.imax+ExtraLayer) + (2 * (H.ny * (v)+ (j-ExtraLayer))))
                        if(rank % 2 == 1) 
                            Hv->uold[IHv(i,j,nv)] = buffer_a[IHbuff(i,j,nv)];
                        else
                            Hv->uold[IHv(i,j,nv)] = buffer_b[IHbuff(i,j,nv)];
                        #undef IHbuff
                    }
                }

            } // End for j
        #undef IHv
    } // End for nv


    free(buffer_a); // Free up workspace
    free(buffer_b);
    MPI_Barrier(MPI_COMM_WORLD);
}