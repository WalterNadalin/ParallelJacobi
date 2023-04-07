# Inputs for the executables
dim = 10
itr  = 100
prc = 4

# Compiler
CC = cc
CFLAGS = -O3

# Directories
INCLUDE = -I./include/

# Source files
SRCS = $(wildcard src/*.c)
OBJ = $(SRCS:.c=.o)
MAIN = $(wildcard *.c)
OBJ += $(MAIN:.c=.o)
EXE = $(MAIN:.c=.x)

all: $(EXE)

# Compiling the object files
%.o: %.c 
	$(CC) -c $< -o $@ $(CFLAGS) $(INCLUDE)

# Compiling the executanle
$(EXE): $(OBJ) 
	$(CC) -o $(EXE) $(OBJ) -O3
	rm $(OBJ)

mpi: CC = mpicc
mpi: CFLAGS += -DMPI
mpi: all

run: clean all
	./$(EXE) $(dim) $(itr)
	
mpirun: clean mpi
	mpirun -np $(prc) ./$(EXE) $(dim) $(itr)

clean:
	@rm -f *$(EXE) plot/solution.dat src/*.o *.o

plot:
	@gnuplot -p plot/plot.plt


.PHONY: clean plot all mpi
