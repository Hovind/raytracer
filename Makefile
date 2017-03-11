CC=mpicc
CXX=mpic++
ifeq ($(DEBUG),yes)
	CXXFLAGS=-std=c++11 -g -fopenmp -pthread -Wall -pedantic
else
	CXXFLAGS=-std=c++11 -O3 -march=native -fopenmp -pthread -Wall -pedantic
endif

CFLAGS=-std=c99 -fopenmp -pthread -Wall -Wextra -pedantic #-Ofast -march=native 
LDFLAGS= -pthread -lm -pg

EXEC= raytracer.exe
SRC= raytracer.cpp
OBJ= $(SRC:.cpp=.cpp.o)

BIN = raytracer
CSRC = raytracer.c
COBJ = $(CSRC:.c=.c.o)

all: $(EXEC)
ifeq ($(DEBUG),yes)
	@echo "Génération en mode debug"
else
	@echo "Génération en mode production"
endif

$(BIN): $(COBJ)
	@$(CC) -o $@ $^ $(LDFLAGS)

$(COBJ): $(CSRC)
	@$(CC) -o $@ -c $< $(CFLAGS)

$(EXEC): $(OBJ)
	@$(CXX) -o $@ $^ $(LDFLAGS)

$(OBJ): $(SRC)
	@$(CXX) -o $@ -c $< $(CXXFLAGS)

.PHONY: clean cleanall

clean:
	@rm -fr *.o *~

cleanall: clean
	@rm -fr $(EXEC)
	
