BLUE = \033[0;34m
YELLOW = \033[1;33m
NC = \033[0m

# Inputs for the executables
dim    ?= 10
iters  ?= 100
prc    ?= 4
delay  ?= 20
frames ?= 200
mode   ?= all

# Compiler
CC     := cc
CFLAGS := -O3 -Wall

# Directories
INCLUDE := -I./include/
IMPI    := -I${SMPI_ROOT}/include
LMPI    := -L${SMPI_ROOT}/lib -lmpiprofilesupport -lmpi_ibm

# Files
SRC  := $(wildcard src/*.c)
MAIN := $(wildcard *.c)
OBJ  := $(SRC:.c=.o) $(MAIN:.c=.o)
EXE  := $(MAIN:.c=.x)

# Conditional flag to toggle the debugging, printing result on a file or printing frames
ifeq ($(debug), yes)
        CFLAGS += -DDEBUG
        EXE    := $(EXE:.x=_debug.x)
endif

ifeq ($(plot), yes)
        CFLAGS += -DPLOT
        EXE    := $(EXE:.x=_plot.x)
endif

ifeq ($(video), yes)
        CFLAGS += -DFRAMES=$(frames)
        EXE    := $(EXE:.x=_video.x)
endif

# Make
all: $(EXE)

mpi: CC     := mpicc
mpi: CFLAGS += -DMPI
mpi: mpi$(EXE)

openacc: CC      := pgcc
openacc: CFLAGS  += -DMPI -DOPENACC -Minfo=all -acc -ta=tesla -fast
openacc: INCLUDE += $(IMPI)
openacc: LINK    := $(LMPI)
openacc: openacc$(EXE)

# Compiling the object files
%.o: %.c 
	$(CC) -c $< -o $@ $(CFLAGS) $(INCLUDE)

# Linking the executable
%.x: $(OBJ) 
	$(CC) -o $@ $^ $(LINK) $(CFLAGS)
	@rm $(OBJ)

# Running the executable
%run: %
	mpirun -np $(prc) ./$^$(EXE) $(dim) $(iters)

run: $(EXE)
	./$(EXE) $(dim) $(iters)

# Other commands
analysis:
	@python analysis/bar_plot.py

plot: data/solution.dat
	@gnuplot -p plot/plot.plt || echo "Please install gnuplot";\

data/solution.dat: CFLAGS += -DPLOT
data/solution.dat: $(EXE:.x=_plot.x)
	./$^ $(dim) $(iters)

gif: video/00000.png
	@magick -delay $(delay) video/*.png video/animation.gif;
	
video/00000.png: CFLAGS += -DFRAMES=$(frames)
video/00000.png: $(EXE:.x=_video.x)
	@mkdir -p video
	./$^ $(dim) $(iters)

clean:
	@rm -f *.x *_t1

flush:
	@rm -f video/*.png plot/*.png video/animation.gif data/solution.dat

format: $(SRCS) $(MAIN)
	@clang-format -i $^ -verbose || echo "Please install clang-format"

.PHONY: analysis debug benchmark clean flush frames plot gif format
.INTERMEDIATE: $(OBJ) data/solution.dat
