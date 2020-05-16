CXX=clang++
CXXFLAGS=-std=c++14 -O3

PKG_CXXFLAGS = $(SHLIB_OPENMP_CXXFLAGS) 
PKG_LIBS = $(SHLIB_OPENMP_CXXFLAGS) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)

CPP_SOURCES = $(wildcard mips/*.cpp) $(wildcard branchbound/*.cpp) $(wildcard *.cpp) $(wildcard other/*.cpp)
OBJECTS = $(CPP_SOURCES:.cpp=.o)

R_LIBS = -L/usr/local/lib -F/Library/Frameworks/R.framework/.. -framework R -L/usr/local/lib -lpcre2-8 -llzma -lbz2 -lz -licucore -ldl -lm -liconv
R_INC = -I/Library/Frameworks/R.framework/Resources/include -I"/Library/Frameworks/R.framework/Versions/4.0/Resources/library/RcppArmadillo/include" -I/Library/Frameworks/R.framework/Versions/4.0/Resources/library/Rcpp/include 

all: WHInter.so 
WHInter.so: $(OBJECTS)

%.o: %.cpp $(CPP_SOURCES)
	$(CXX) $(R_INC) $(CXXFLAGS) -c $< -o $@

clean:
	find . -type f -name '*.o' -exec rm {} \;
	#-rm -f *.o
