#
# TP CUDA
# Maison de la Simulation
# Saclay, December 2011
#

CC=gcc
CXX=g++

DEBUG=no

# PATH to CUDA toolkit
CUDAROOT=/usr/local/cuda-6.5

NVCC=$(CUDAROOT)/bin/nvcc

#
# nvcc flags. 
# --device-emulation to generate device emulation mode.
# --keep             to keep all intermediate files. 
# --keep --clean     to remove all intermediate files.
# -ptxas-options=-v
# see also nvcc man page
#
ifeq ($(DEBUG),yes)
NVCFLAGS=-m64 -G -g -O0 --compiler-options -O0,-Wall,-fPIC,-g 
CFLAGS=-m64 -g3 -O0
else
NVCFLAGS=-m64 -O3 --compiler-options -fno-strict-aliasing 
CFLAGS=-m64 -O3
endif

NVCFLAGS+= -gencode arch=compute_13,code=sm_13 \
	--ptxas-options=-v --use_fast_math \
	-I./cuda_helper

# usefull for profiling
NVCFLAGS+= -lineinfo

CFLAGS+=-Wall -Wextra -fPIC \
        -L $(CUDAROOT)/lib64 -Wl,-rpath=$(CUDAROOT)/lib64

# C++ flags
CXXFLAGS=$(CFLAGS)

EXTRA_LDFLAGS=
# visualization with MGL
# uncomment the 3 following lines to enable MGL (http://mathgl.sourceforge.net/)
# install package libmgl-dev under Ubuntu
CFLAGS   += -DUSE_MGL
CXXFLAGS += -DUSE_MGL
EXTRA_LDFLAGS= -lmgl

#### HDF5 output (uncomment the following 8 lines to have hdf5 output)
#### HDF5 output is usefull when you want to find easily where CPU and GPU output files differ (h5diff -r -d 1e-7 output1.h5 output2.h5)
# H5CC = h5c++
# CFLAGS         += $(shell $(H5CC) -showconfig | grep '\bCPPFLAGS:' | awk -F: '{print $$2}') -DUSE_HDF5
# CXXFLAGS       += $(shell $(H5CC) -showconfig | grep '\bCPPFLAGS:' | awk -F: '{print $$2}') -DUSE_HDF5
# EXTRA_LDFLAGS  += $(shell $(H5CC) -showconfig | grep '\bLDFLAGS:'  | awk -F: '{print $$2}')
# EXTRA_LDFLAGS += -lhdf5_hl -lhdf5
# # on Ubuntu when using package libhdf5-openmpi
# # one needs to add ldflags coming from openmpi
# EXTRA_LDFLAGS += -pthread -L/usr/lib/openmpi/lib -lmpi_cxx -lmpi -lopen-rte -lopen-pal -ldl -Wl,--export-dynamic -lnsl -lutil -lm -ldl

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

%.o: %.cpp
#	@echo "########"
	@echo "OBJECT : $@"
	$(CXX) $(CXXFLAGS) -c -o $@ $<

%.o: %.cu
#	@echo "########"
	@echo "OBJECT : $@"
	$(NVCC) $(NVCFLAGS) -c -o $@ $<

