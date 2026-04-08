main: main.cpp polynom.h polynom_optimized.h lattice.h tests.cpp
	g++ -Ofast -march=native -mtune=native main.cpp tests.cpp -o main