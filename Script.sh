#!/bin/bash

make

./vtSerial.out 65536 16 0

./vtOmp.out 65536 16 0 100000

./vtCuda.out 65536 16 0 100000 256 2

./vtCudaOmp.out 65536 16 0 10000 500000 256 2


./vtSerial.out 65536 32 0

./vtOmp.out 65536 32 0 100000

./vtCuda.out 65536 32 0 100000 256 2

./vtCudaOmp.out 65536 32 0 10000 500000 256 2


./vtSerial.out 65536 64 0

./vtOmp.out 65536 64 0 10000

./vtCuda.out 65536 64 0 100000 256 2

./vtCudaOmp.out 65536 64 0 10000 500000 256 2



./vtSerial.out 131072 16 0

./vtOmp.out 131072 16 0 100000

./vtCuda.out 131072 16 0 100000 256 2

./vtCudaOmp.out 131072 16 0 10000 500000 256 2



./vtSerial.out 131072 32 0

./vtOmp.out 131072 32 0 100000

./vtCuda.out 131072 32 0 100000 256 2

./vtCudaOmp.out 131072 32 0 10000 1000000 256 2




./vtSerial.out 131072 64 0

./vtOmp.out 131072 64 0 100000

./vtCuda.out 131072 64 0 100000 256 2

./vtCudaOmp.out 131072 64 0 10000 2000000 256 2




./vtSerial.out 262144 16 0

./vtOmp.out 262144 16 0 100000

./vtCuda.out 262144 16 0 100000 256 2

./vtCudaOmp.out 262144 16 0 10000 1000000 256 2



./vtSerial.out 262144 32 0

./vtOmp.out 262144 32 0 100000

./vtCuda.out 262144 32 0 100000 256 2

./vtCudaOmp.out 262144 32 0 10000 1000000 256 2


./vtSerial.out 262144 64 0

./vtOmp.out 262144 64 0 100000

./vtCuda.out 262144 64 0 100000 256 2

./vtCudaOmp.out 262144 64 0 10000 3000000 256 2


