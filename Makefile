CXX=mpic++
ifeq ($(DEBUG),yes)
	CXXFLAGS=-std=c++11 -g -fopenmp -pthread -Wall -pedantic
	LDFLAGS= -pthread -lm
else
	CXXFLAGS=-std=c++11 -O3 -march=native -fopenmp -pthread -Wall -pedantic
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
	@$(CXX) -o $@ -c $< $(CXXFLAGS)

.PHONY: clean cleanall

clean:
	@rm -fr *.o *~

cleanall: clean
	@rm -fr $(EXEC)
	
