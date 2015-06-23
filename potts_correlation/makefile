#------------------------------------------------------
# Makefile for C/C++ Program
#------------------------------------------------------

TARGET=potts2d

LDFLAGS=

#------------------------------------------------------
# Default Parameters
#------------------------------------------------------

CC=g++
CPPFLAGS=-O3
LDFLAGS=-lfftw3

#------------------------------------------------------
# Compile Option
#------------------------------------------------------

-include makefile.opt

#------------------------------------------------------
# Definition
#------------------------------------------------------

.SUFFIXES: .c .cc .h. .o

#---
# Source Files
#---

SRC=$(shell ls *.cc)
HED=$(shell ls *.h)
OBJ=$(SRC:.cc=.o)

#------------------------------------------------------
# rules
#------------------------------------------------------

all: $(TARGET)
$(TARGET): $(OBJ)
	$(CC) $(CPPFLAGS) $(LDFLAGS) -o $(TARGET) $(OBJ)

.cc.o:
	$(CC) $(CPPFLAGS) -c $< 

.PHONY: clean dep

dep:
	g++ -MM -MG $(SRC) >makefile.depend

makefile.depend: 
	g++ -MM -MG $(SRC) >makefile.depend

clean:
	rm -f $(TARGET) $(OBJ) gmon.*.out gmon.out

tar:
	tar cvzf $(TARGET).tar.gz *.cfg $(SRC) $(HED) makefile

#--------------------------------------------------
-include makefile.depend
