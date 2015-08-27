CC=g++
CFLAGS= -I/usr/include/bullet -lBulletDynamics -lBulletCollision -lLinearMath -lCGAL -std=c++11 -pthread -g

all: main.o utilities.o output.o terrain.o fluid.o lp_grid.o sph.o
	$(CC) main.o utilities.o output.o terrain.o fluid.o lp_grid.o sph.o\
	  -o main $(CFLAGS)

main.o: main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)

utilities.o: utilities.cpp utilities.h
	$(CC) -c utilities.cpp -o utilities.o $(CFLAGS)

output.o: output.cpp output.h
	$(CC) -c output.cpp -o output.o $(CFLAGS)

terrain.o: terrain.cpp terrain.h
	$(CC) -c terrain.cpp -o terrain.o $(CFLAGS)

fluid.o: fluid.cpp fluid.h
	$(CC) -c fluid.cpp -o fluid.o $(CFLAGS)

lp_grid.o: lp_grid.cpp lp_grid.h
	$(CC) -c lp_grid.cpp -o lp_grid.o $(CFLAGS)

sph.o: sph.cpp sph.h
	$(CC) -c sph.cpp -o sph.o $(CFLAGS)

clean:
	rm -rf *o main

reset:
	rm frames/frame*
