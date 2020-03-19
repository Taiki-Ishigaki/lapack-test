# Makefile

CC = gcc

CLA_DIR = $(HOME)/src/CLAPACK-3.2.1
CLA_LIB_DIR = $(CLA_DIR)
F2C_LIB_DIR = $(CLA_DIR)/F2CLIBS
CLA_PLAT = _LINUX
BLAS_LIB = $(CLA_LIB_DIR)/blas$(CLA_PLAT).a
LAPACK_LIB = $(CLA_LIB_DIR)/lapack$(CLA_PLAT).a
FLIBS = $(F2C_LIB_DIR)/libf2c.a
CLA_ALL_LIB = \
	$(LAPACK_LIB) \
	$(BLAS_LIB) \
	$(FLIBS) \

CFLAGS = -Wall -O2 

all: eig-clapack.c
	$(CC) $(CFLAGS) -o a.out eig-clapack.c ${CLA_ALL_LIB}

clean:
	rm -f a.out