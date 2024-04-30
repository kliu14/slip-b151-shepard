#---------------------------------------------------------------------#
# makefile for superhydrobphoic channel code 
# Jongmin Seo & Sanjeeb T. Bose 
# Modification by Danah Park
# date created : 05/07/2008 
# date modified : 04/23/2017
#---------------------------------------------------------------------#

MKL = /opt/ohpc/pub/compiler/intel-18/compilers_and_libraries_2018.2.199/linux/mkl/lib/intel64
MPI = /opt/ohpc/pub/compiler/intel-18/compilers_and_libraries_2018.2.199/linux/mpi

PFFT = /home/mani/kimliu/libraries/fft_sandia

LAPACK_LIB = -L/usr/lib64 -lblas -llapack

opt = -r8  -c -O3
opt3 = -c -O3

FC=	mpiifort

FFLAGS=  -I$(MPI)/include64 -O3  -r8 -i4 -fpp 
EFLAGS=  -I$(MPI)/includei64 -O3  -r8 -i4

LLOC=   -L$(MKL) -L$(PFFT)
LIBS=    $(LAPACK_LIB) -lz $(PFFT)/libpfft.a 

OBJ= 	global.o grid.o comm_routines.o  numerics.o numerics_tee.o output.o
MOBJ=	main.o 
POBJ=   post_snapshot.o 
PVOBJ=   vorticity_post.o
PPOBJ = poisson_2d_check.o
TOBJ = test_restart_ab.o
BIN=	../bin

#---------------------------------------------------------------------# 
# main compilation rules 
#---------------------------------------------------------------------#

channel : $(OBJ)  main.o 
	$(FC) $(EFLAGS) $(LLOC)  -o $(BIN)/channel -Xlinker --start-group main.o $(OBJ) $(LIBS)  -Xlinker --end-group

%.o : %.f90 
	$(FC) $(opt3) $(FFLAGS)  $<  

clean : 
	rm -f $(OBJ) $(OBJ:.o=.mod) $(MOBJ)

clean_data: 
	rm -vf ../output/*.hdfviz ../output/*.h5mesh 
	rm -vf ../output/collapsed_field*
	rm -vf ../output/mean_*.tec 
	rm -vf ../output/diagnostics.dat 
	rm -vf ../output/clip_stats.dat 
