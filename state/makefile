CC     = g++
DEBUG  = -g
OMP    = -fopenmp
#OMP =
LAPACK_DIR=/gpfs/home/shyyuan/local/lapack-3.8.0
BOUNDARY=/gpfs/home/shyyuan/lp_base/boundary
STATE=/gpfs/home/shyyuan/lp_base/state
GEOMETRY=/gpfs/home/shyyuan/lp_base/geometry
MAIN_DIR=/gpfs/home/shyyuan/lp_base
INCS   = -I /gpfs/home/nnaitlho/lapack-3.8.0/LAPACKE/include -I $(BOUNDARY) -I $(STATE) -I $(GEOMETRY)
LIBS   = -L $(LAPACK_DIR) 
CFLAGS = -Wall -c -std=c++11 $(DEBUG) $(OMP)
LFLAGS = -Wall  $(DEBUG) $(INCS) $(LIBS) $(OMP)
SRC_FILES := $(wildcard *.cpp) 
OBJS := $(SRC_FILES:.cpp=.o)
all: $(OBJS)

state.o: state.h state.cpp
	$(CC) $(CFLAGS) state.cpp

state_1d.o: state.h state_1d.h state_1d.cpp
	$(CC) $(CFLAGS) state_1d.cpp

state_ballexp.o: state.h state_ballexp.h state_ballexp.cpp
	$(CC) $(CFLAGS) state_ballexp.cpp

state_collision.o: state.h state_collision.h state_collision.cpp
	$(CC) $(CFLAGS) state_collision.cpp

state_gresho.o: state.h state_gresho.h state_gresho.cpp
	$(CC) $(CFLAGS) state_gresho.cpp

state_jet.o: state.h state_jet.h state_jet.cpp
	$(CC) $(CFLAGS) state_jet.cpp

state_shocktube.o: state.h state_shocktube.h state_shocktube.cpp
	$(CC) $(CFLAGS) state_shocktube.cpp

state_powder_target.o: state.h state_powder_target.h state_powder_target.cpp
	$(CC) $(CFLAGS) state_powder_target.cpp

state_powder_target_3d.o: state.h state_powder_target_3d.h state_powder_target_3d.cpp
	$(CC) $(CFLAGS) state_powder_target_3d.cpp

state_pellet.o: state.h state_pellet.h state_pellet.cpp
	$(CC) $(CFLAGS) state_pellet.cpp

time_controller.o: time_controller.h time_controller.cpp lp_solver.h particle_viewer.h initializer.h
	$(CC) $(CFLAGS) time_controller.cpp

state_nozzle.o: state.h state_nozzle.h state_nozzle.cpp
	$(CC) $(CFLAGS) state_nozzle.cpp


clean:
	rm *.o



