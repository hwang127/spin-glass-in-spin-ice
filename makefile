CC=g++ -std=c++11 -O3 -fopenmp
CODE_SOURCES = .*cpp
SOURCES= $(CODE_SOURCES)
LAPACK= -L /cm/shared/apps/lapack/open64/64/3.5.0/ -L /usr/lib -L /cm/shared/apps/blas/gcc/current/lib64/ -llapack -lblas -lm -lgomp -fopenmp -ffast-math -lpthread
BOOST= -I /cm/shared/apps/boost/1.55.0/include/boost/ -I /cm/shared/apps/boost/1.55.0/include/

EXECUTABLE=mcfile

CODE_OBJECTS=\
	./main.o	\
	./search_for.o	\

OBJECTS= $(CODE_OBJECTS)

headers1=./*.h 

$(EXECUTABLE): $(OBJECTS) $(headers1)
	$(CC) $(GSL_INC) $(NO_WRITE) $(OBJECTS) $(CFLAGS) $(LAPACK) $(GSL_LIB) -o $@

$(CODE_OBJECTS) : ./%.o : ./%.cpp
	@echo "compiling $<";
	@$(CC) $(GSL_INC) $(BLAS_INC) $(BOOST) -c $< -o 	$@;

clean:
	rm -f $(OBJECTS) $(PROG)
	rm -f $(EXECUTABLE)

depend:
	makedepend $(INCLUDE) $(GSL_INC) $(BOOST) $(BLAS_INC) $(GSL_LIB) $(BLAS_LIB) -- -o $(CFLAGS) -- $(SOURCES) 



# DO NOT DELETE

