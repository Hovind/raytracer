CXX=g++
ifeq ($(DEBUG),yes)
	CFLAGS=-std=c++11 -g -fopenmp -pthread -Wall -pedantic
	LDFLAGS= -pthread -lm
else
	CFLAGS=-std=c++11 -O3 -mnative -fopenmp -pthread -Wall -pedantic
	LDFLAGS= -pthread -lm
endif
EXEC= raytracer.exe
SRC= raytracer.cpp
OBJ= $(SRC: .cpp=.o)

all: $(EXEC)
ifeq ($(DEBUG),yes)
	@echo "Génération en mode debug"
else
	@echo "Génération en mode production"
endif

raytracer.exe: $(OBJ)
	@$(CXX) -o $@ $^ $(LDFLAGS)

%.o : %.cpp
	@$(CXX) -o $@ -c $< $(CFLAGS)

.PHONY: clean cleanall

clean:
	@rm -fr *.o *~

cleanall: clean
	@rm -fr $(EXEC)
	
