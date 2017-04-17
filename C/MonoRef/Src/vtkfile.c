/*
  A simple 2D hydro code
  (C) Romain Teyssier : CEA/IRFU           -- original F90 code
  (C) Pierre-Francois Lavallee : IDRIS      -- original F90 code
  (C) Guillaume Colin de Verdiere : CEA/DAM -- for the C version
*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <strings.h>
#include <mpi.h>

#include "parametres.h"
#include "utils.h"
#include "vtkfile.h"
void
vtkfile(long step, const hydroparam_t H, hydrovar_t * Hv, int mask)
{
    // Get number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    // Get rank of processes
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Get number of output variables
    int ifID = (mask / 1) % 2;
    int ifIU = (mask / 2) % 2;
    int ifIV = (mask / 4) % 2;
    int ifIP = (mask / 8) % 2;
    int nvArr[4] = {ifID, ifIU, ifIV, ifIP};
    //int num = ifID + ifIU + ifIV + ifIP;
    
    if(rank == 0)
    {
        char name[160];
        FILE *fic;
        int i, j, k, nv;

        WHERE("vtkfile");
        sprintf(name, "outputvtk_%05ld.vts", step);
        fic = fopen(name, "w");
        if (fic == NULL) {
            fprintf(stderr, "Ouverture du fichier %s impossible\n", name);
            exit(1);
        }
        fprintf(fic, "<?xml version=\"1.0\"?>\n");
        fprintf(fic, "<VTKFile type=\"StructuredGrid\">\n");
        fprintf(fic, "<StructuredGrid WholeExtent=\" %ld %ld %ld %ld %ld %ld\">\n", (long)0,
                H.nxg, (long)0, H.nyg, (long)0, (long)0);
        fprintf(fic, "<Piece Extent=\" %ld %ld %ld %ld %ld %ld\">\n", (long)0, H.nxg, (long)0, H.nyg, (long)0, (long)0);
        fprintf(fic, "<Points>\n");
        fprintf(fic,
                "<DataArray type=\"Float32\" format=\"ascii\" NumberOfComponents=\"3\">\n");
        for (j = 0; j < H.nyg + 1; j++) {
            for (i = 0; i < H.nxg + 1; i++) {
                fprintf(fic, "%.2f %.2f %.2f\n", i * H.dx, j * H.dx, 0.0);
            }
        }
        fprintf(fic, "</DataArray>\n");
        fprintf(fic, "</Points>\n");
        name[0] = 0;
        for (nv = 0; nv < IP; nv++) {
            if(nvArr[nv] == 0)
                continue;
            if (nv == ID)
                sprintf(name, "%s varID", name);
            if (nv == IU)
                sprintf(name, "%s varIU", name);
            if (nv == IV)
                sprintf(name, "%s varIV", name);
            if (nv == IP)
                sprintf(name, "%s varIP", name);
        }


        // Create temporary buffer of size world_size*H->uold of rank 0 (rank 0 is always the largest)
        int usize = (H.nxt)*(H.nyt)*(H.nvar);
        double* buffer = (double*)malloc(sizeof(double) * world_size*usize);
        // Write data from rank 0 to buffer
        for(k = 0; k < usize; k++)
        {
            buffer[k] = Hv->uold[k];
        }

        // Rank 0 collects data from every process and writes it to temporary buffer
        for(k = 1; k < world_size; k++)
        {
            int size, nxl, nyl;
            nxl = (H.nxg / world_size) + ExtraLayerTot;
            nyl = H.nyt;
            if(k < H.nxg % world_size) nxl += 1; // Depending on rank nx is 1 entry longer
            size = nxl * nyl * H.nvar;

            // Get data from process k and write it to temporary buffer at position k*usize           
            MPI_Recv(buffer+k*usize, size, MPI_DOUBLE, k, k, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }


        // declaration of the variable list
        fprintf(fic, "<CellData Scalars=\"%s\">\n", name);
        name[0] = 0;
        for (nv = 0; nv <= IP; nv++) {
            if(nvArr[nv] == 0)
                continue;
            if (nv == ID)
                sprintf(name, "varID");
            if (nv == IU)
                sprintf(name, "varIU");
            if (nv == IV)
                sprintf(name, "varIV");
            if (nv == IP)
                sprintf(name, "varIP");

            //Definition of the cell values
            fprintf(fic, "<DataArray type=\"Float32\" Name=\"%s\" format=\"ascii\">\n", name);

            // The image is the interior of the computed domain
            for(j = H.jmin + ExtraLayer; j < H.jmax - ExtraLayer; j++)
            {
                // Go over all processes and append one after the other
                for(k = 0; k < world_size; k++)
                {
                    int nxl, nyl;
                    nxl = (H.nxg / world_size) + ExtraLayerTot;
                    nyl = H.nyt;
                    if(k < H.nxg % world_size) nxl += 1; // Depending on rank nx is 1 entry longer

                    for(i = 0 + ExtraLayer; i < nxl - ExtraLayer; i++)
                    {
                        #define IHv_loc(i, j, v) ((i) + (nxl * (nyl * (v)+ (j))))
                        fprintf(fic, "%.3e ", *(buffer + k*usize + IHv_loc(i, j, nv)));
                        #undef IHv_loc
                    }
                }
                fprintf(fic, "\n");
            } // End for j

            fprintf(fic, "</DataArray>\n");
        } // End of declaration of the variable list

        // Clean up temporary buffer
        free(buffer);

        fprintf(fic, "</CellData>\n");
        fprintf(fic, "</Piece>\n");
        fprintf(fic, "</StructuredGrid>\n");
        fprintf(fic, "</VTKFile>\n");
        fclose(fic);
    } // End execution of rank 0
    else
    {
        MPI_Send(Hv->uold, (H.nxt)*(H.nyt)*(H.nvar), MPI_DOUBLE, 0, rank, MPI_COMM_WORLD);
    }
    //printf("rank %d before barrier\n",rank);
    MPI_Barrier(MPI_COMM_WORLD);
    //printf("rank %d after barrier\n",rank);
}
