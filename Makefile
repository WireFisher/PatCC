CXX := mpicxx
MPI_PATH := /opt/intel/impi/3.2.0.011
#PDLN_USE_OPENCV := true
#PDLN_USE_NETCDF := true
NETCDF_PATH := /opt/netCDF-gcc4.4.7
OPENCV_PATH := /home/yanghy/opt/opencv
CXXFLAGS :=

INC := -isystem $(MPI_PATH)/include64
INC += -isystem dependency/googletest/include
INC += -isystem dependency/googlemock/include

LIB :=
LIB += -Ldependency/googletest
LIB += -Ldependency/googlemock

LIBS := -lgmock -lgtest

core_objs = ccpl_utils.o \
            component.o \
			processing_unit_mgt.o \
			delaunay_grid_decomposition_mgt.o \
			delaunay_voronoi_2D.o

test_objs = DelaunayVoronoi2D.o \
			FullProcess.o \
			GridDecomposition.o \
			ProcessingResourceTest.o

COMMON_FLAGS := -Wall -g -fopenmp -pthread #-O3 -DTIME_PERF

ifeq ($(PDLN_USE_NETCDF),true)
	COMMON_FLAGS += -DNETCDF
	INC += -isystem $(NETCDF_PATH)/include
	LIB += -L$(NETCDF_PATH)/lib 
	LIBS += -lnetcdf
	core_objs += netcdf_utils.o
endif

ifeq ($(PDLN_USE_OPENCV),true)
	COMMON_FLAGS += -DOPENCV
	INC += -isystem $(OPENCV_PATH)/include
	LIB += -L$(OPENCV_PATH)/lib64
	LIBS += -lopencv_core -lopencv_imgproc -lopencv_highgui
	core_objs += opencv_utils.o
endif

ifeq ($(PDLN_DEBUG),true)
	COMMON_FLAGS += -DDEBUG -g
endif

COMMON_FLAGS += $(INC)
COMMON_FLAGS += $(LIB)
COMMON_FLAGS += $(LIBS)

VPATH = ./ ./unittest

component.o: delaunay_grid_decomposition_mgt.o processing_unit_mgt.o
delaunay_grid_decomposition_mgt.o: processing_unit_mgt.o delaunay_voronoi_2D.o ccpl_utils.o #netcdf_utils.o opencv_utils.o
delaunay_voronoi_2D.o: merge_sort.h #opencv_utils.o

%.o: %.cxx %.h
	$(CXX) -c -o $@ $(CXXFLAGS) $< $(COMMON_FLAGS)
%.o: %.cxx
	$(CXX) -c -o $@ $(CXXFLAGS) $< $(COMMON_FLAGS)

.PHONY : main
main : main.o $(core_objs)
	$(CXX) -o pDelaunay $(CXXFLAGS) main.o $(core_objs) $(COMMON_FLAGS)

.PHONY : test
test : testmain.o $(test_objs) $(core_objs)
	$(CXX) -o run_all_test -DUNITTEST testmain.o $(core_objs) $(test_objs) $(COMMON_FLAGS)

.PHONY : all
all :

.PHONY : clean
clean :
	-rm run_all_test pDelaunay *.o
