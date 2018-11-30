CC     = g++
DEBUG  = -g
OMP    = -fopenmp
#OMP =
LAPACK_DIR=/gpfs/home/shyyuan/local/lapack-3.8.0
MAIN_DIR:=${CURDIR}/../
BOUNDARY_DIR=$(MAIN_DIR)/boundary/
STATE_DIR=$(MAIN_DIR)/state/
GEOMETRY_DIR=$(MAIN_DIR)/geometry/
INCS   = -I $(BOUNDARY_DIR) -I $(STATE_DIR) -I $(GEOMETRY_DIR) -I $(MAIN_DIR)
LIBS   = -L $(LAPACK_DIR) 
CFLAGS = -Wall -std=c++11 $(DEBUG) $(OMP) $(INCS) -c
LFLAGS = -Wall  $(DEBUG) $(INCS) $(LIBS) $(OMP)
SRC_FILES := $(wildcard *.cpp) 
OBJS := $(SRC_FILES:.cpp=.o)
all: $(OBJS)


boundary.o: boundary.h boundary.cpp
	$(CC) $(CFLAGS) boundary.cpp

boundary_solid_gresho.o: boundary.h boundary_solid_gresho.h boundary_solid_gresho.cpp
	$(CC) $(CFLAGS) boundary_solid_gresho.cpp

boundary_solid_shocktube.o: boundary.h boundary_solid_shocktube.h boundary_solid_shocktube.cpp
	$(CC) $(CFLAGS) boundary_solid_shocktube.cpp

boundary_solid_tpshocktube.o: boundary.h boundary_solid_tpshocktube.h boundary_solid_tpshocktube.cpp
	$(CC) $(CFLAGS) boundary_solid_tpshocktube.cpp

boundary_solid_shocktube3d.o: boundary.h boundary_solid_shocktube3d.h boundary_solid_shocktube3d.cpp
	$(CC) $(CFLAGS) boundary_solid_shocktube3d.cpp

boundary_rayleightaylor.o: boundary.h boundary_rayleightaylor.h boundary_rayleightaylor.cpp
	$(CC) $(CFLAGS) boundary_rayleightaylor.cpp

boundary_rayleightaylor_periodic.o: boundary.h boundary_rayleightaylor_periodic.h boundary_rayleightaylor_periodic.cpp
	$(CC) $(CFLAGS) boundary_rayleightaylor_periodic.cpp

boundary_rayleightaylor3d.o: boundary.h boundary_rayleightaylor3d.h boundary_rayleightaylor3d.cpp
	$(CC) $(CFLAGS) boundary_rayleightaylor3d.cpp

boundary_kelvinhelmholtz.o: boundary.h boundary_kelvinhelmholtz.h boundary_kelvinhelmholtz.cpp
	$(CC) $(CFLAGS) boundary_kelvinhelmholtz.cpp

boundary_dambreak.o: boundary.h boundary_dambreak.h boundary_dambreak.cpp
	$(CC) $(CFLAGS) boundary_dambreak.cpp

boundary_powder_target.o: boundary.h boundary_powder_target.h boundary_powder_target.cpp
	$(CC) $(CFLAGS) boundary_powder_target.cpp

boundary_powder_target_3d.o: boundary.h boundary_powder_target_3d.h boundary_powder_target_3d.cpp
	$(CC) $(CFLAGS) boundary_powder_target_3d.cpp

boundary_nozzle.o: boundary.h boundary_nozzle.h boundary_nozzle.cpp
	$(CC) $(CFLAGS) boundary_nozzle.cpp

boundary_pellet.o: boundary.h boundary_pellet.h boundary_pellet.cpp
	$(CC) $(CFLAGS) boundary_pellet.cpp



clean:
	rm *.o


