#CC=icpc
CC=g++

IDIR=${INCLUDE}
#LDIR=${MKLROOT}/lib/intel64/
#CFLAGS=-O3  -Wall -openmp -lpthread
#CFLAGS=-O3  -Wall -fopenmp -lpthread -std=c++11 -pg
CFLAGS=-O3  -Wall -fopenmp -lpthread -std=c++11

objects = main.cpp
toclean = MoCHA

MoCHA: $(objects) src/headers.h
	$(CC) -I$(IDIR) -L$(LDIR) $(LIBS) $(CFLAGS) $(objects) -o $@

#MoCHA_CL: $(objects) src/headers.h
#	$(CC) -I$(IDIR) -L$(LDIR) $(LIBS) $(CFLAGS) $(objects) -o $@

%.o: %.c
	$(CC) -I$(IDIR) -L$(LDIR) $(LIBS) $(CFLAGS) -c $<

.PHONY: clean
clean:
	rm -rf $(toclean)
