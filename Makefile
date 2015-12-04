#load makefile from arpack
#include /home/sgeraedt/sources/arpack++/Makefile.inc

CC=g++
HOME=/home/geraedts
EIGEN_DIR=$(HOME)/sources/eigen-eigen-bdd17ee3b1b3/
MYDIR=$(HOME)/sources/lapack-3.5.0/lapacke/include/
MYCDIR=$(HOME)/myClibrary/
ARPACKPP_DIR = $(HOME)/sources/arpack++/include/
CFLAGS=-c -O3 -Wall -I$(MYDIR) -I$(ARPACKPP_DIR) -I$(MYCDIR) -I$(EIGEN_DIR)
#LDFLAGS=-I$(LIBDIR) -I$(MYDIR) 
OBJECTS= utils.o

all: $(OBJECTS)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
	
clean:
	rm *.o
