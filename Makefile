#load makefile from arpack
#include /home/sgeraedt/sources/arpack++/Makefile.inc

CC=g++
HOME=/home/sgeraedt
MYDIR=$(HOME)/sources/lapack-3.5.0/lapacke/include/
MYCDIR=$(HOME)/myClibrary/
ARPACKPP_DIR = $(HOME)/sources/arpack++/include/
CFLAGS=-c -O3 -Wall -I$(MYDIR) -I$(ARPACKPP_DIR) -I$(MYCDIR)
#LDFLAGS=-I$(LIBDIR) -I$(MYDIR) 
OBJECTS= utils.o

all: $(OBJECTS)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
	
clean:
	rm *.o
