all: vtSerial vtOmp vtCuda vtCudaOmp

vtSerial: vtSerial.c
	gcc -O3 $< -o $@.out -lm -g utilFunctions.c
	
vtOmp: vtOmp.c
	gcc -O3 $< -o $@.out -fopenmp -lm -g utilFunctions.c
	
vtCuda: vtCuda.cu
	nvcc $< -o $@.out -O3 -g -lm utilFunctions.c
	
vtCudaOmp: vtCudaOmp.cu
	nvcc $< -o $@.out utilFunctions.c -Xcompiler -fopenmp -lgomp -O3 -lm 
	


	

