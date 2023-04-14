BLUE = \033[0;34m
NC = \033[0m

# Inputs for the executables
dim = 10
iters = 100
prc = 4
delay = 20
frames = 200

# Compiler
CC = cc
CFLAGS = -O3

# Directories
INCLUDE = -I./include/

# Source files
SRCS = $(wildcard src/*.c)
MAIN = $(wildcard *.c)
OBJ = $(SRCS:.c=.o) $(MAIN:.c=.o)
EXE = $(MAIN:.c=.x)

# Libraries
IMPI = -I${SMPI_ROOT}/include
LMPI = -L${SMPI_ROOT}/lib -lmpiprofilesupport -lmpi_ibm

# Make
all: clean $(EXE)

# Compiling the object files
%.o: %.c 
	$(CC) -c $< -o $@ $(CFLAGS) $(INCLUDE)

# Compiling the executanle
$(EXE): $(OBJ) 
	$(CC) -o $(EXE) $^ $(LINK) -O3
	@rm $(OBJ)

mpi: CC = mpicc
mpi: CFLAGS += -DMPI
mpi: clean all

cuda: CC = pgcc
cuda: INCLUDE += $(IMPI)
cuda: CFLAGS += -DMPI -DCUDA
cuda: LINK += $(LMPI)
cuda: clean all

run: all
	./$(EXE) $(dim) $(iters)
	
mpirun: mpi
	mpirun -np $(prc) ./$(EXE) $(dim) $(iters)

cudarun: cuda
	mpirun -np $(prc) ./$(EXE) $(dim) $(iters)

debug: CFLAGS += -DDEBUG
debug: run

mpidebug: CFLAGS += -DDEBUG
mpidebug: mpirun

cudadebug: CFLAGS += -DDEBUG
cudadebug: cudarun

clean:
	@rm -f *$(EXE) src/*.o *.o 
	
flush:
	@rm -f video/*.png plot/*.png video/animation.gif plot/solution.dat 

plot:
	@if [ -e plot/solution.dat ]; then\
		gnuplot -p plot/plot.plt || echo "Please install ${BLUE}gnuplot${NC} to run this command";\
	else\
		echo "Please generate data using the command ${BLUE}make run${NC}";\
	fi

frames: CFLAGS += -DFRAMES=$(frames)
frames: flush run
	
mpiframes: CFLAGS += -DFRAMES=$(frames)
mpiframes: flush mpirun

gif:
	@if [ -e video/00000.png ]; then\
		magick -delay $(delay) video/*.png video/animation.gif || echo "Please install ${BLUE}imagemagick${NC} to run this command";\
	else\
		echo "Please generate frames using the command ${BLUE}make frames${NC}";\
	fi

format: $(SRCS) $(MAIN)
	@clang-format -i $^ -verbose || echo "Please install ${BLUE}clang-format${NC} to run this command"

.PHONY: clean plot all mpi cuda gif frames run mpirun cudarun debug mpidebug cudadebug format
