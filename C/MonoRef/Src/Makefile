# Link make.inc with the targeted architecture in Arch directory
# For example : ln -s make_sp6 make.inc
#--------------------------------------
include ../../../Arch/make.inc

SRC=cmpflx.c conservar.c equation_of_state.c hydro_godunov.c main.c parametres.c riemann.c trace.c vtkfile.c compute_deltat.c constoprim.c hydro_funcs.c make_boundary.c qleftright.c slope.c utils.c mpi_sync.c

OBJ = $(patsubst %.c, %.o, ${SRC})

default: mono

mono: CCOMP = $(CC)
mono: CCOMPFLAGS = $(CFLAGS)
mono: LINKFLAGS = $(LDFLAGS)
mono: $(OBJ)
	$(CCOMP) $(LINKFLAGS) $(OBJ) -o hydro_mono -lm

mpi: CCOMP = $(MPICC)
mpi: CCOMPFLAGS = $(CFLAGS) $(CFLAGS_MPI)
mpi: LINKFLAGS = $(LDFLAGS) $(LDFLAGS_MPI)
mpi: $(OBJ)
	$(CCOMP) $(LINKFLAGS) $(OBJ) -o hydro_mpi -lm

omp: CCOMP = $(CC)
omp: CCOMPFLAGS = $(CFLAGS) $(CFLAGS_OMP)
omp: LINKFLAGS = $(LDFLAGS) $(LDFLAGS_OMP)
omp: $(OBJ)
	$(CCOMP) $(LINKFLAGS) $(OBJ) -o hydro_omp -lm

mpiomp: CCOMP = $(MPICC)
mpiomp: CCOMPFLAGS = $(CFLAGS) $(CFLAGS_MPI) $(CFLAGS_OMP)
mpiomp: LINKFLAGS = $(LDFLAGS) $(LDFLAGS_MPI) $(LDFLAGS_OMP)
mpiomp: $(OBJ)
	$(CCOMP) $(LINKFLAGS) $(OBJ) -o hydro_mpiomp -lm


clean   :
	rm -f *.o core* hydro_mono hydro_mpi hydro_omp hydro_mpiomp

cleanall   : clean
	rm -f *~ *.vts *.out

.SUFFIXES:  .o .c

.c.o    :
	${CCOMP} ${CCOMPFLAGS} -c $<

# DO NOT DELETE

cmpflx.o: parametres.h utils.h cmpflx.h
compute_deltat.o: parametres.h compute_deltat.h utils.h equation_of_state.h
conservar.o: parametres.h utils.h conservar.h
constoprim.o: parametres.h constoprim.h utils.h
equation_of_state.o: equation_of_state.h
equation_of_state.o: utils.h parametres.h
hydro_funcs.o: utils.h parametres.h hydro_funcs.h
hydro_godunov.o: parametres.h hydro_godunov.h
hydro_godunov.o: hydro_funcs.h utils.h make_boundary.h cmpflx.h conservar.h
hydro_godunov.o: equation_of_state.h qleftright.h constoprim.h riemann.h
hydro_godunov.o: trace.h slope.h
main.o: parametres.h hydro_funcs.h vtkfile.h
main.o: compute_deltat.h hydro_godunov.h utils.h
make_boundary.o: parametres.h make_boundary.h
make_boundary.o: utils.h
parametres.o: parametres.h
qleftright.o: parametres.h utils.h qleftright.h
riemann.o: parametres.h utils.h riemann.h
slope.o: parametres.h utils.h slope.h
trace.o: parametres.h utils.h trace.h
utils.o: utils.h parametres.h
vtkfile.o: parametres.h utils.h
vtkfile.o: vtkfile.h
cmpflx.o: utils.h parametres.h
constoprim.o: utils.h parametres.h
equation_of_state.o: utils.h parametres.h
hydro_funcs.o: parametres.h
utils.o: parametres.h
