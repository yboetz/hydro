## 2D hydro code ##

Code for the lecture 'High Performance Computing' in spring semester of 2017.

Based on a simple 2D hydro code by

  + (C) Romain Teyssier : CEA/IRFU (original F90 code)
  + (C) Pierre-Francois Lavallee : IDRIS (original F90 code)
  + (C) Guillaume Colin de Verdiere : CEA/DAM (for the C version)

MPI/OpenMP functionallity added by Yannick Boetzel.


### How to run it ###

Got to `Src` folder and compile the code:

    cd C/MonoRef/Src
    make clean
    make {mpi|mpiomp|omp}

Choose either one of the options in `{mpi|mpiomp|omp}` depending if you want to run the MPI, OpenMP
or hybrid MPI/OpenMP version. `mpi` and `mpiomp` can both be used on the `mpiomp` branch, `omp` can
only be used on the `omp` branch.

The makefile will compile everything and create a file `hydro_{mpi|mpiomp|omp}`. You can then run it
by first specifying the number of OpenMP threads

    export OMP_NUM_THREADS = 2

then run it either by

    mpirun -np 4 --map-by core --bind-to core hydro_mpiomp -i input.nml

for the MPI or hybrid MPI/OpenMP version (4 cores, binds threads to cores), or by

  ./hydro_omp -i input.nml

for the pure OpenMP version. The input files are located in the `Input` folder.

You can also use one of the two python scripts `weak.py` and `strong.py` to test weak and strong
scaling of the code. Usage here is

    python weak.py -g 500 500 -s 100 -n 2 4 8 12 16
    python strong.py -g 500 500 -s 100 -n 2 4 8 12 16

This will automatically run `hydro_mpiomp` on a grid of 500x500 for 100 steps using 2, 4, 8, 12
and then 16 cpu cores. If you want to enable OpenMP just run it with the flag `--omp`

    python weak.py -g 500 500 -s 100 -n 2 4 8 12 16 --omp

An overview of all options is given with

    python weak.py --help
