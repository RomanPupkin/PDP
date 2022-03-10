all:
	g++ -O3 -march=native Vector.cpp energy.cpp main.cpp -o main
openmp:
	g++ -fopenmp Vector.cpp energy.cpp main.cpp -o main
numThreads:
	export OMP_NUM_THREADS=1
clean:
	rm main