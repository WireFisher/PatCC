#PDLN_USE_OPENCV := false
#PDLN_USE_NETCDF := false
CXX := mpicxx
MPI_PATH := /opt/intel/impi/3.2.0.011
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
#%.o: %.c
#processing_unit_mgt.cxx delaunay_grid_decomposition_mgt.cxx delaunay_voronoi_2D.cxx ccpl_utils.cxx opencv_utils.cxx netcdf_utils.cxx component.cxx unittest/*.cxx dependency/googlemock/libgmock.a dependency/googletest/libgtest.a -o run_all_test

COMMON_FLAGS := -fopenmp -pthread
COMMON_FLAGS += $(INC)
COMMON_FLAGS += $(LIB)
COMMON_FLAGS += $(LIBS)

core_objs = ccpl_utils.o \
            component.o \
			processing_unit_mgt.o \
			delaunay_grid_decomposition_mgt.o \
			delaunay_voronoi_2D.o

plot_objs = opencv_utils.o \
			netcdf_utils.o

test_objs = DelaunayVoronoi2D.o \
			FullProcess.o \
			GridDecomposition.o \
			GridDecomposition_GlobalManyDifferent.o \
			GridDecomposition_GlobalMiddle.o \
			GridDecomposition_SmallRegion.o \
			ProcessingResourceTest.o

ifeq ($(PDLN_USE_NETCDF),true)
	COMMON_FLAGS += -DOPENCV
	INC += -isystem $(NETCDF_PATH)/include
	LIB += -L$(NETCDF_PATH)/lib 
	LIBS += -lnetcdf
endif

ifeq ($(PDLN_USE_OPENCV),true)
	COMMON_FLAGS += -DOPENCV
	INC += -isystem $(OPENCV_PATH)/include
	LIB += -L$(OPENCV_PATH)/lib64
	LIBS += -lopencv_core -lopencv_imgproc -lopencv_imgcodecs
	core_objs += $(plot_objs)
endif

ifeq ($(PDLN_DEBUG),true)
	COMMON_FLAGS += -DDEBUG -g
endif

VPATH = ./ ./unittest
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
