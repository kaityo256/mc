#------------------------------------------------------
# Makefile for C/C++ Program
#------------------------------------------------------

TARGET=a.out
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

.SUFFIXES: .cpp .h. .o

#---
# Source Files
#---

SRC=$(shell ls *.cpp)
HED=$(shell ls *.h)
OBJ=$(SRC:.cpp=.o)

#------------------------------------------------------
# rules
#------------------------------------------------------

all: $(TARGET)
$(TARGET): $(OBJ)
	$(CC) $(CPPFLAGS) $(LDFLAGS) -o $(TARGET) $(OBJ)

.cpp.o:
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
