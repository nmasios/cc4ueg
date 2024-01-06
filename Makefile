
CXX=g++
#ARMA_INCDIR   = /home/fs71035/masios5/code/armadillo-9.700.2/include
ARMA_LIBDIR   = /gpfs/opt/sw/spack-0.17.1/opt/spack/linux-almalinux8-zen3/gcc-11.2.0/armadillo-10.5.0-zzssso6lwzgjpsuubriirjj67cf2rin6/lib64
#ARMA_INC_FLAGS = -I${ARMA_INCDIR}
ARMA_LIB_FLAGS = -L${ARMA_LIBDIR} -larmadillo
OPTIMIZATION = -O3 -g

CXXFLAGS := $(OPTIMIZATION) $(ARMA_LIB_FLAGS) -funroll-loops -std=c++17 -frename-registers -march=native -fopenmp
#CXXFLAGS := -O3 -funroll-loops -std=c++17 -frename-registers -larmadillo -march=native -fopenmp

SRC := $(wildcard *.cpp)
OBJS := $(patsubst %.cpp,%.o,$(SRC))

EXE := ueg

.PHONY: all


all: clean $(EXE)


$(EXE): $(OBJS)
		$(CXX) $? $(CXXFLAGS) -o $@



clean:
	$(RM) $(OBJS)
