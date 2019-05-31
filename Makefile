############
# Required #
############
CXX := mpicxx
#MPI_PATH := /opt/intel/impi/3.2.0.011
MPI_PATH := /opt/mpich-3.2-gcc4.4.7/
PAT_OPENCV := true
PAT_NETCDF := true
PAT_TIMING := false
PAT_DEBUG := true
PAT_MUTE := true
NETCDF_PATH := /opt/netCDF-gcc4.4.7
#NETCDF_PATH := /opt/netCDF-intel13-without-hdf5
OPENCV_PATH := /home/yanghy/opt/opencv
CXXFLAGS :=

############
# Optional #
############

SRCDIR := src
OBJDIR := obj
TESTDIR := unittest

INC := -isystem $(MPI_PATH)/include64
INC += -isystem dependency/googletest/include
INC += -isystem dependency/googlemock/include
INC += -I./$(SRCDIR)

LIB :=
LIBS:=

TestLib :=
TestLib += -Ldependency/googletest
TestLib += -Ldependency/googlemock
TestLibs := -lgmock -lgtest

SOURCES  := $(wildcard $(SRCDIR)/*.cxx)
SOURCES  := $(filter-out $(SRCDIR)/opencv_utils.cxx, $(SOURCES))
SOURCES  := $(filter-out $(SRCDIR)/netcdf_utils.cxx, $(SOURCES))

test_objs = obj/testmain.o \
			obj/FullProcess.o \
			obj/ProcessingResourceTest.o \
			obj/DelaunayVoronoi2D.o
			#obj/GridDecomposition.o \

ADDED_FLAGS := -O3
COMMON_FLAGS := -Wall -fopenmp -pthread 


ifeq ($(PAT_TIMING),true)
	COMMON_FLAGS += -DTIME_PERF
endif

ifeq ($(PAT_DEBUG),true)
	COMMON_FLAGS += -DDEBUG -g
else
	COMMON_FLAGS += -O3
endif

ifeq ($(PAT_MUTE),true)
	COMMON_FLAGS += -DDEFAULT_LOGLEVEL=LOG_ERROR
endif

ifeq ($(PAT_NETCDF),true)
	COMMON_FLAGS += -DNETCDF
	INC += -I$(NETCDF_PATH)/include
	LIB += -L$(NETCDF_PATH)/lib 
	LIBS += -lnetcdf
	SOURCES += $(SRCDIR)/netcdf_utils.cxx
endif

ifeq ($(PAT_OPENCV),true)
	COMMON_FLAGS += -DOPENCV
	INC += -isystem $(OPENCV_PATH)/include
	LIB += -L$(OPENCV_PATH)/lib64
	LIBS += -lopencv_core -lopencv_imgproc -lopencv_highgui
	SOURCES += $(SRCDIR)/opencv_utils.cxx
endif

COMMON_FLAGS += $(ADDED_FLAGS)
INCLUDES := $(wildcard $(SRCDIR)/*.h)
OBJECTS  := $(SOURCES:$(SRCDIR)/%.cxx=$(OBJDIR)/%.o)
OBJECTS  := $(filter-out $(OBJDIR)/main.o, $(OBJECTS))

COMMON_FLAGS += $(INC)
COMMON_FLAGS += $(LIB)
COMMON_FLAGS += $(LIBS)

VPATH = ./ ./unittest


.PHONY : main
main : $(OBJDIR) $(OBJECTS) $(OBJDIR)/main.o
	$(CXX) $(CXXFLAGS) $(COMMON_FLAGS) $(OBJECTS) $(OBJDIR)/main.o -o patcc

$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.cxx
	$(CXX) $(CXXFLAGS) $(COMMON_FLAGS) -c $< -o $@

$(OBJDIR)/main.o: $(SRCDIR)/main.cxx
	$(CXX) $(CXXFLAGS) $(COMMON_FLAGS) -c $< -o $@

$(OBJDIR):
	@mkdir $(OBJDIR)


.PHONY : test
test : $(OBJDIR) $(test_objs) $(OBJECTS) $(OBJDIR)/testmain.o
	$(CXX) -DUNITTEST $(OBJECTS) $(test_objs) $(COMMON_FLAGS) $(TestLib) $(TestLibs) -o run_all_test

$(test_objs): $(OBJDIR)/%.o : $(TESTDIR)/%.cxx
	$(CXX) $(CXXFLAGS) $(COMMON_FLAGS) $(TestLib) $(TestLibs) -c $< -o $@


.PHONY : all
all : main test


.PHONY : clean
clean :
	-rm run_all_test patcc obj/* 2>/dev/null
