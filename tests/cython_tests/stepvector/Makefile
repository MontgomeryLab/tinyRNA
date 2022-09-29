CXX        := clang++

PY_VERS    := python3.9
PY_CFLAGS  := $(shell ${PY_VERS}-config --cflags --libs)
PY_LDFLAGS := $(shell ${PY_VERS}-config --ldflags --embed)

INCLUDE    := -I$(realpath ../../../tiny/rna/counter/stepvector/src)
LDFLAGS    := ${PY_LDFLAGS} -L${CONDA_PREFIX}/lib -rpath ${CONDA_PREFIX}/lib
CXXFLAGS   := -std=c++17 -stdlib=libc++ ${INCLUDE} ${PY_CFLAGS}

test_core : test_core.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ -o $@

clean :
	@rm -f test_core