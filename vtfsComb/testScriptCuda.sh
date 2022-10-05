#!/bin/bash

make

./vtSerial.out 65536 16 0

./vtCuda.out 65536 16 0 1000000 256 4

./vtCuda.out 65536 16 0 100000 256 4

./vtCuda.out 65536 16 0 10000 256 4

./vtCuda.out 65536 16 0 1000 256 4

./vtSerial.out 65536 32 0

./vtCuda.out 65536 32 0 1000000 256 2

./vtCuda.out 65536 32 0 100000 256 2

./vtCuda.out 65536 32 0 10000 256 2

./vtCuda.out 65536 32 0 1000 256 2

./vtSerial.out 65536 64 0

./vtCuda.out 65536 64 0 1000000 256 1

./vtCuda.out 65536 64 0 100000 256 1

./vtCuda.out 65536 64 0 10000 256 1

./vtCuda.out 65536 64 0 1000 256 1


./vtSerial.out 131072 16 0


./vtCuda.out 131072 16 0 1000000 256 4

./vtCuda.out 131072 16 0 100000 256 4

./vtCuda.out 131072 16 0 10000 256 4

./vtCuda.out 131072 16 0 1000 256 4


./vtSerial.out 131072 32 0


./vtCuda.out 131072 32 0 1000000 256 2

./vtCuda.out 131072 32 0 100000 256 2

./vtCuda.out 131072 32 0 10000 256 2

./vtCuda.out 131072 32 0 1000 256 2

./vtSerial.out 131072 64 0

./vtCuda.out 131072 64 0 1000000 256 1

./vtCuda.out 131072 64 0 100000 256 1

./vtCuda.out 131072 64 0 10000 256 1

./vtCuda.out 131072 64 0 1000 256 1


./vtSerial.out 262144 16 0

./vtCuda.out 262144 16 0 1000000 256 4

./vtCuda.out 262144 16 0 100000 256 4

./vtCuda.out 262144 16 0 10000 256 4

./vtCuda.out 262144 16 0 1000 256 4


./vtSerial.out 262144 32 0


./vtCuda.out 262144 32 0 1000000 256 2

./vtCuda.out 262144 32 0 100000 256 2

./vtCuda.out 262144 32 0 10000 256 2

./vtCuda.out 262144 32 0 1000 256 2

./vtSerial.out 262144 64 0

./vtCuda.out 262144 64 0 1000000 256 1

./vtCuda.out 262144 64 0 100000 256 1

./vtCuda.out 262144 64 0 10000 256 1

./vtCuda.out 262144 64 0 1000 256 1

