BLUE = \033[0;34m
YELLOW = \033[1;33m
NC = \033[0m

# Inputs for the executables
dim = 10
iters = 100
prc = 4
delay = 20
frames = 200
mode = all

# Compiler
CC := cc
CFLAGS := -O3 -Wall

# Directories
INCLUDE := -I./include/
IMPI := -I${SMPI_ROOT}/include
LMPI := -L${SMPI_ROOT}/lib -lmpiprofilesupport -lmpi_ibm

# Source files
SRCS := $(wildcard src/*.c)
MAIN := $(wildcard *.c)
OBJ := $(SRCS:.c=.o) $(MAIN:.c=.o)
EXE := $(MAIN:.c=.x)

# Make
all: $(PREFIX)$(EXE)

mpi: CC := mpicc
mpi: CFLAGS += -DMPI
mpi mpirun: EXE := mpi_$(EXE)

openacc: CC := pgcc
openacc: CFLAGS += -DMPI -DOPENACC -Minfo=all -acc -ta=tesla
openacc: INCLUDE += $(IMPI)
openacc: LINK := $(LMPI)
openacc openaccrun: EXE := openacc_$(EXE)

mpi openacc: all

# Compiling the object files
%.o: %.c 
	$(CC) -c $< -o $@ $(CFLAGS) $(INCLUDE)

# Linking the executable
$(EXE): $(OBJ) 
	$(CC) -o $(PREFIX)$(EXE) $^ $(LINK) $(CFLAGS)
	@rm $(OBJ)

run:
	$(SCRIPT) $(COMMAND) ./$(EXE) $(dim) $(iters)

mpirun openaccrun: SCRIPT := mpirun
mpirun openaccrun: COMMAND := -np $(prc)
mpirun openaccrun: run

debug: CFLAGS += -DDEBUG
	
benchmark: CFLAGS += -DBENCHMARK

frames: CFLAGS += -DFRAMES=$(frames)

debug benchmark frames: $(mode)

plot:
	@if [ -e data/solution.dat ]; then\
		gnuplot -p plot/plot.plt || echo "Please install ${YELLOW}gnuplot${NC} to run this command";\
	else\
		echo "Please generate data using the command ${YELLOW}make run${NC}";\
	fi

gif:
	@if [ -e video/00000.png ]; then\
		magick -delay $(delay) video/*.png video/animation.gif || echo "Please install ${YELLOW}imagemagick${NC} to run this command";\
	else\
		echo "Please generate frames using the command ${YELLOW}make frames run${NC}";\
	fi
	
clean:
	@rm -f *$(EXE) src/*.o *.o 
	
flush:
	@rm -f video/*.png plot/*.png video/animation.gif data/*.dat

format: $(SRCS) $(MAIN)
	@clang-format -i $^ -verbose || echo "Please install ${YELLOW}clang-format${NC} to run this command"

.PHONY: debug benchmark clean flush frames plot gif format
