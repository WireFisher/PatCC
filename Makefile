CXX := mpiicpc
MPI_PATH := /opt/intel/impi/3.2.0.011
#PDLN_USE_OPENCV := true
#PDLN_USE_NETCDF := true
#PDLN_TIMING := true
NETCDF_PATH := /opt/netCDF-gcc4.4.7
OPENCV_PATH := /home/yanghy/opt/opencv
CXXFLAGS :=

SRCDIR := src
OBJDIR := obj
TESTDIR := unittest

INC := -isystem $(MPI_PATH)/include64
INC += -isystem dependency/googletest/include
INC += -isystem dependency/googlemock/include
INC += -I $(SRCDIR)

LIB :=
LIB += -Ldependency/googletest
LIB += -Ldependency/googlemock

LIBS := -lgmock -lgtest


SOURCES  := $(wildcard $(SRCDIR)/*.cxx)
SOURCES  := $(filter-out $(SRCDIR)/opencv_utils.cxx, $(SOURCES))
SOURCES  := $(filter-out $(SRCDIR)/netcdf_utils.cxx, $(SOURCES))

test_objs = obj/testmain.o \
			obj/FullProcess.o \
			obj/ProcessingResourceTest.o \
			obj/DelaunayVoronoi2D.o
			#obj/GridDecomposition.o \

COMMON_FLAGS := -Wall -g -fopenmp -pthread

ifeq ($(PDLN_TIMING),true)
	COMMON_FLAGS += -DTIME_PERF
endif

ifeq ($(PDLN_DEBUG),true)
	COMMON_FLAGS += -DDEBUG -g
else
	COMMON_FLAGS += -O3
endif

ifeq ($(PDLN_USE_NETCDF),true)
	COMMON_FLAGS += -DNETCDF
	INC += -isystem $(NETCDF_PATH)/include
	LIB += -L$(NETCDF_PATH)/lib 
	LIBS += -lnetcdf
	SOURCES += $(SRCDIR)/netcdf_utils.cxx
endif

ifeq ($(PDLN_USE_OPENCV),true)
	COMMON_FLAGS += -DOPENCV
	INC += -isystem $(OPENCV_PATH)/include
	LIB += -L$(OPENCV_PATH)/lib64
	LIBS += -lopencv_core -lopencv_imgproc -lopencv_highgui
	SOURCES += $(SRCDIR)/opencv_utils.cxx
endif

INCLUDES := $(wildcard $(SRCDIR)/*.h)
OBJECTS  := $(SOURCES:$(SRCDIR)/%.cxx=$(OBJDIR)/%.o)
OBJECTS  := $(filter-out $(OBJDIR)/main.o, $(OBJECTS))

COMMON_FLAGS += $(INC)
COMMON_FLAGS += $(LIB)
COMMON_FLAGS += $(LIBS)

VPATH = ./ ./unittest


.PHONY : main
main : $(OBJECTS) $(OBJDIR)/main.o
	$(CXX) $(CXXFLAGS) $(COMMON_FLAGS) $(OBJECTS) $(OBJDIR)/main.o -o patcc

$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.cxx
	$(CXX) $(CXXFLAGS) $(COMMON_FLAGS) -c $< -o $@

$(OBJDIR)/main.o: $(SRCDIR)/main.cxx
	$(CXX) $(CXXFLAGS) $(COMMON_FLAGS) -c $< -o $@


.PHONY : test
test : $(test_objs) $(OBJECTS) $(OBJDIR)/testmain.o
	$(CXX) -DUNITTEST $(OBJECTS) $(test_objs) $(COMMON_FLAGS) -o run_all_test

$(test_objs): $(OBJDIR)/%.o : $(TESTDIR)/%.cxx
	$(CXX) $(CXXFLAGS) $(COMMON_FLAGS) -c $< -o $@


.PHONY : all
all : main test


.PHONY : clean
clean :
	-rm run_all_test patcc obj/* 2>/dev/null
