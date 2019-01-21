CXX := mpicxx
MPI_PATH := /opt/intel/impi/3.2.0.011
#PDLN_USE_OPENCV := true
#PDLN_USE_NETCDF := true
#PDLN_TIMING := true
NETCDF_PATH := /opt/netCDF-gcc4.4.7
OPENCV_PATH := /opt/opencv
CXXFLAGS :=

INC := -isystem $(MPI_PATH)/include64

LIB :=

LIBS :=

core_objs = ccpl_utils.o \
            component.o \
			processing_unit_mgt.o \
			delaunay_grid_decomposition_mgt.o \
			delaunay_voronoi_2D.o \
			memory_pool.o \
			coordinate_hash.o

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
	core_objs += netcdf_utils.o
endif

ifeq ($(PDLN_USE_OPENCV),true)
	COMMON_FLAGS += -DOPENCV
	INC += -isystem $(OPENCV_PATH)/include
	LIB += -L$(OPENCV_PATH)/lib64
	LIBS += -lopencv_core -lopencv_imgproc -lopencv_highgui
	core_objs += opencv_utils.o
endif

COMMON_FLAGS += $(INC)
COMMON_FLAGS += $(LIB)
COMMON_FLAGS += $(LIBS)

VPATH = ./

.PHONY : main
main : main.o $(core_objs)
	$(CXX) -o patcc $(CXXFLAGS) main.o $(core_objs) $(COMMON_FLAGS)

component.o: delaunay_grid_decomposition_mgt.o processing_unit_mgt.o
delaunay_grid_decomposition_mgt.o: processing_unit_mgt.o delaunay_voronoi_2D.o ccpl_utils.o #netcdf_utils.o opencv_utils.o
delaunay_voronoi_2D.o: merge_sort.h #opencv_utils.o

%.o: %.cxx %.h
	$(CXX) -c -o $@ $(CXXFLAGS) $< $(COMMON_FLAGS)
%.o: %.cxx
	$(CXX) -c -o $@ $(CXXFLAGS) $< $(COMMON_FLAGS)

.PHONY : all
all :

.PHONY : clean
clean :
	-rm patcc *.o 2>/dev/null
