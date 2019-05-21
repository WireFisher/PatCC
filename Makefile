CXX := mpiicpc
MPI_PATH := /opt/intel/impi/3.2.0.011
#PDLN_USE_OPENCV := true
#PDLN_USE_NETCDF := true
#PDLN_TIMING := true
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

core_objs = projection.o \
            patcc.o \
			processing_unit_mgt.o \
			grid_decomposition.o \
			triangulation.o \
			memory_pool.o \
			coordinate_hash.o \
			logger.o

test_objs = \
			FullProcess.o \
			ProcessingResourceTest.o \
			DelaunayVoronoi2D.o
			#GridDecomposition.o \

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

VPATH = ./ ./unittest

.PHONY : main
main : main.o $(core_objs)
	$(CXX) -o patcc $(CXXFLAGS) main.o $(core_objs) $(COMMON_FLAGS)

patcc.o: patcc.cxx grid_decomposition.h patcc.h common_utils.h projection.h timer.h
grid_decomposition.o: grid_decomposition.cxx grid_decomposition.h common_utils.h projection.h netcdf_utils.h opencv_utils.h timer.h
triangulation.o: triangulation.cxx triangulation.h opencv_utils.h common_utils.h merge_sort.h coordinate_hash.h
projection.o: projection.cxx projection.h common_utils.h
processing_unit_mgt.o: processing_unit_mgt.cxx processing_unit_mgt.h common_utils.h timer.h
coordinate_hash.o: coordinate_hash.cxx coordinate_hash.h common_utils.h

%.o: %.cxx %.h
	$(CXX) -c -o $@ $(CXXFLAGS) $< $(COMMON_FLAGS)
%.o: %.cxx
	$(CXX) -c -o $@ $(CXXFLAGS) $< $(COMMON_FLAGS)

.PHONY : test
test : testmain.o $(test_objs) $(core_objs)
	$(CXX) -o run_all_test -DUNITTEST testmain.o $(core_objs) $(test_objs) $(COMMON_FLAGS)

.PHONY : all
all :

.PHONY : clean
clean :
	-rm run_all_test patcc *.o 2>/dev/null
