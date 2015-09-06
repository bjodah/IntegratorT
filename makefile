CPP = g++

CFLAGS = -g -c

LFLAGS = -lm

DEPENDPATH = -I. 

objects = main.o \
	IntegratorT.o \
	StiffIntegratorT.o \
	NonStiffIntegratorT.o \
	decsol.o

integrate: $(objects)
	$(CPP) -o integrate $(objects) $(LFLAGS)
	
main.o: main.cpp IntegratorT.o
	$(CPP) $(CFLAGS) $(DEPENDPATH) main.cpp

IntegratorT.o: IntegratorT.cpp IntegratorT.h StiffIntegratorT.o \
				NonStiffIntegratorT.o
	$(CPP) $(CFLAGS) $(DEPENDPATH) IntegratorT.cpp
	
StiffIntegratorT.o: StiffIntegratorT.cpp StiffIntegratorT.h decsol.o
	$(CPP) $(CFLAGS) $(DEPENDPATH) StiffIntegratorT.cpp

NonStiffIntegratorT.o: NonStiffIntegratorT.cpp NonStiffIntegratorT.h
	$(CPP) $(CFLAGS) $(DEPENDPATH) NonStiffIntegratorT.cpp

decsol.o: decsol.cpp decsol.h
	$(CPP) $(CFLAGS) decsol.cpp
	
clean:
	rm integrate $(objects)
