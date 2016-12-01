#! /usr/bin/env make

CUDA_INC_PATH=/usr/local/cuda50/include
CUDA_LIB_PATH=/usr/local/cuda50/lib64

CC=nvcc
CFLAGS= -O3 -I$(CUDA_INC_PATH) --ptxas-options=-v -arch sm_13
LDFLAGS=-L$(CUDA_LIB_PATH) -lcudart -lcublas 

all: saxpy

saxpy:
	$(CC) $(CFLAGS) $(LDFLAGS) tools.cu Timer.cpp saxpy.cu -o saxpy

clean:
	\rm -f *.o
	\rm -f saxpy
	\rm -f *~
	\rm -f *.linkinfo *.log









