#load makefile from arpack
#include /home/sgeraedt/sources/arpack++/Makefile.inc

CC=icpc -g
HOME=/home/geraedts


EIGEN_DIR=$(HOME)/sources/eigen-eigen-bdd17ee3b1b3/
MYCDIR=$(HOME)/myClibrary/

ARPACKPP_DIR = $(HOME)/sources/arpack++/include/
MYDIR=$(HOME)/sources/lapack-3.5.0/lapacke/include/
SUPER_LU = $(HOME)/sources/SuperLU/libsuperlu_3.0.a
SUPER_LUH = $(HOME)/sources/SuperLU/SRC/

CFLAGS=-c -O3 -Wall -I$(MYDIR) -I$(ARPACKPP_DIR) -I$(MYCDIR) -I$(EIGEN_DIR) -I$(SUPER_LUH)
#LDFLAGS=-I$(LIBDIR) -I$(MYDIR) 
OBJECTS= utils.o superLU_complex.o superLU_real.o 

all: $(OBJECTS)

test: test.o utils.o
	g++ test.o utils.o -o test
.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
	
clean:
	rm *.o
