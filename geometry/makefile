CC     = g++
DEBUG  = -g
OMP    = -fopenmp
#OMP =
LAPACK_DIR=/gpfs/home/shyyuan/local/lapack-3.8.0
BOUNDARY=/gpfs/home/shyyuan/lp_base/boundary
STATE=/gpfs/home/shyyuan/lp_base/state
GEOMETRY=/gpfs/home/shyyuan/lp_base/geometry
MAIN_DIR=/gpfs/home/shyyuan/lp_base
INCS   = -I $(BOUNDARY) -I $(STATE) -I $(GEOMETRY)
LIBS   = -L $(LAPACK_DIR) 
CFLAGS = -Wall -c -std=c++11 $(DEBUG) $(OMP)
LFLAGS = -Wall  $(DEBUG) $(INCS) $(LIBS) $(OMP)
SRC_FILES := $(wildcard *.cpp) 
OBJS := $(SRC_FILES:.cpp=.o)
all: $(OBJS)

geometry.o: geometry.h geometry.cpp
	$(CC) $(CFLAGS) geometry.cpp

geometry_1d.o: geometry.h geometry_1d.h geometry_1d.cpp
	$(CC) $(CFLAGS) geometry_1d.cpp

geometry_ballexp.o: geometry.h geometry_ballexp.h geometry_ballexp.cpp
	$(CC) $(CFLAGS) geometry_ballexp.cpp

geometry_collision.o: geometry.h geometry_collision.h geometry_collision.cpp
	$(CC) $(CFLAGS) geometry_collision.cpp

geometry_gresho.o: geometry.h geometry_gresho.h geometry_gresho.cpp
	$(CC) $(CFLAGS) geometry_gresho.cpp

geometry_jet.o: geometry.h geometry_jet.h geometry_jet.cpp
	$(CC) $(CFLAGS) geometry_jet.cpp

geometry_shocktube.o: geometry.h geometry_shocktube.h geometry_shocktube.cpp
	$(CC) $(CFLAGS) geometry_shocktube.cpp

geometry_shocktube3d.o: geometry.h geometry_shocktube3d.h geometry_shocktube3d.cpp
	$(CC) $(CFLAGS) geometry_shocktube3d.cpp

geometry_powder_target.o: geometry.h geometry_powder_target.h geometry_powder_target.cpp
	$(CC) $(CFLAGS) geometry_powder_target.cpp

geometry_powder_target_3d.o: geometry.h geometry_powder_target_3d.h geometry_powder_target_3d.cpp
	$(CC) $(CFLAGS) geometry_powder_target_3d.cpp

geometry_nozzle.o: geometry.h geometry_nozzle.h geometry_nozzle.cpp
	$(CC) $(CFLAGS) geometry_nozzle.cpp

geometry_random.o: geometry.h geometry_random.h geometry_random.cpp
	$(CC) $(CFLAGS) geometry_random.cpp

geometry_pellet.o: geometry.h geometry_pellet.h geometry_pellet.cpp
	${CC} ${CFLAGS} geometry_pellet.cpp



clean:
	rm *.o


