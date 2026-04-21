main: main.cpp polynom.h polynom_optimized.h lattice.h lattice_optimized.h tests.cpp spins.h
	g++ -Ofast -march=native -mtune=native -fopenmp main.cpp tests.cpp -o main