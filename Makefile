CFITSIO = $(FITSIOROOT)
CPP = g++
CC = gcc
CFLAGS = -g -Wall -I$(CFITSIO) $(shell root-config --cflags) -O3
LIBS = -L$(CFITSIO) -lcfitsio -lm $(shell root-config --glibs) -lTreeViewer
GLIBS =
GLIBS +=
OBJECTS = main.o realdata.o simulation.o aux_functions.o
HEADERS =
DEBUG_FLAG = -g

ALL : main.exe
	@echo File has been successfully compiled $(NEWLINE)

main.exe : $(OBJECTS)
	$(CPP) $(OBJECTS) $(DEBUG_FLAG) -o main.exe $(LIBS) $(GLIBS) $(CFLAGS)

main.o : main.cpp $(HEADERS)
	$(CPP) -c main.cpp $(DEBUG_FLAG) -o main.o $(CFLAGS)

realdata.o : realdata.cpp $(HEADERS)
	$(CPP) -c realdata.cpp $(DEBUG_FLAG) -o realdata.o $(CFLAGS)

simulation.o : simulation.cpp $(HEADERS)
	$(CPP) -c simulation.cpp $(DEBUG_FLAG) -o simulation.o $(CFLAGS)

aux_functions.o : aux_functions.cpp $(HEADERS)
	$(CPP) -c aux_functions.cpp $(DEBUG_FLAG) -o aux_functions.o $(CFLAGS)

clean:
	rm -f *~ *.o *.exe
