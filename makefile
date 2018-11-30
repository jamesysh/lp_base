CC     = g++
DEBUG  = -g
OMP    = -fopenmp
#OMP =
LAPACK_DIR=/gpfs/home/shyyuan/local/lapack-3.8.0
MAIN_DIR:=${CURDIR}
BOUNDARY_DIR=$(MAIN_DIR)/boundary/
STATE_DIR=$(MAIN_DIR)/state/
GEOMETRY_DIR=$(MAIN_DIR)/geometry/
INCS   = -I $(BOUNDARY_DIR) -I $(STATE_DIR) -I $(GEOMETRY_DIR) -I $(MAIN_DIR)
LIBS   = -L $(LAPACK_DIR) 
CFLAGS = -Wall -c -std=c++11 $(DEBUG) $(OMP) $(INCS)
LFLAGS = -Wall  $(DEBUG) $(INCS) $(LIBS) $(OMP)
vpath %.h $(GEOMETRY_DIR) $(STATE_DIR) $(BOUNDARY_DIR)
MAIN_OBJS   = eos.o hexagonal_packing.o initializer.o lp_main.o lp_solver.o ls_solver.o\
         neighbour_searcher.o octree.o particle_data.o\
		 particle_viewer.o registrar.o time_controller.o


BOUNDARY_OBJS   = boundary.o boundary_solid_shocktube.o boundary_solid_shocktube3d.o boundary_solid_tpshocktube.o boundary_solid_gresho.o \
	boundary_rayleightaylor.o boundary_rayleightaylor_periodic.o boundary_rayleightaylor3d.o boundary_kelvinhelmholtz.o \
	boundary_dambreak.o boundary_powder_target.o boundary_powder_target_3d.o boundary_nozzle.o boundary_pellet.o 

GEOMETRY_OBJS = geometry.o geometry_1d.o geometry_nozzle.o geometry_random.o\
		 geometry_ballexp.o geometry_collision.o geometry_gresho.o geometry_powder_target.o geometry_powder_target_3d.o\
		 geometry_jet.o geometry_shocktube.o geometry_shocktube3d.o geometry_pellet.o

STATE_OBJS =  state.o state_1d.o state_ballexp.o state_collision.o state_gresho.o state_powder_target.o state_powder_target_3d.o state_nozzle.o\
	     state_jet.o state_shocktube.o state_pellet.o

B_OBJS := $(foreach OBJ,$(BOUNDARY_OBJS),$(addprefix $(BOUNDARY_DIR),$(OBJ)))
S_OBJS := $(foreach OBJ,$(STATE_OBJS),$(addprefix $(STATE_DIR),$(OBJ)))
G_OBJS := $(foreach OBJ,$(GEOMETRY_OBJS),$(addprefix $(GEOMETRY_DIR),$(OBJ)))

OBJS = $(B_OBJS) $(S_OBJS) $(G_OBJS) $(MAIN_OBJS)
all: lp 
 

lp: build  $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o lp -lgomp -llapacke -llapack -lgfortran -lrefblas

build:
	cd $(BOUNDARY_DIR)&&make;
	cd $(STATE_DIR)&&make;
	cd $(GEOMETRY_DIR)&&make; 
	
eos.o: eos.h eos.cpp
	$(CC) $(CFLAGS) eos.cpp


hexagonal_packing.o: hexagonal_packing.h hexagonal_packing.cpp
	$(CC) $(CFLAGS) hexagonal_packing.cpp

initializer.o: initializer.h initializer.cpp geometry.h state.h eos.h hexagonal_packing.h\
				neighbour_searcher.h
	$(CC) $(CFLAGS) initializer.cpp

lp_main.o: lp_main.cpp initializer.h neighbour_searcher.h particle_data.h particle_viewer.h\
           lp_solver.h time_controller.h
	$(CC) $(CFLAGS) lp_main.cpp

lp_solver.o: lp_solver.h lp_solver.cpp neighbour_searcher.h eos.h particle_data.h\
           initializer.h ls_solver.h hexagonal_packing.h
	$(CC) $(CFLAGS) lp_solver.cpp

ls_solver.o: ls_solver.h ls_solver.cpp
	$(CC) $(CFLAGS) ls_solver.cpp

neighbour_searcher.o: neighbour_searcher.h neighbour_searcher.cpp octree.h initializer.h
	$(CC) $(CFLAGS) neighbour_searcher.cpp

octree.o: octree.h octree.cpp
	$(CC) $(CFLAGS) octree.cpp

particle_data.o: particle_data.h particle_data.cpp initializer.h
	$(CC) $(CFLAGS) particle_data.cpp

particle_viewer.o: particle_viewer.h particle_viewer.cpp particle_data.h
	$(CC) $(CFLAGS) particle_viewer.cpp

registrar.o: registrar.h registrar.cpp geometry.h state.h geometry_collision.h state_collision.h
	$(CC) $(CFLAGS) registrar.cpp

time_controller.o: time_controller.h time_controller.cpp lp_solver.h particle_viewer.h initializer.h
	$(CC) $(CFLAGS) time_controller.cpp


clean:
		cd $(BOUNDARY_DIR)&&make clean;
		cd $(STATE_DIR)&&make clean;
		cd $(GEOMETRY_DIR)&&make clean; 
		rm *.o *~ debug log save_init_param lp -r out

tar:
	tar cvzf lp_backup.tar.gz *.h *.cpp makefile* input* README* dox* run* boundary/* state/* geometry/*

