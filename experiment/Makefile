try=0
choose=0
CC=mpicc
CFLAGS = -O0 -g -Wall -I./ -I/usr/local/include # -I/u/shared/programs/x86_64/openmpi/1.8.3/gnu/4.9.2/torque/include/ 
LIBS = -L/usr/local/lib -lfftw3_mpi -lfftw3 -lm
EXE = main.x
SRC = main.c
OBJ = $(SRC:.c=.o)

prc=4

all: $(EXE)

$(EXE): $(OBJ)
	$(CC) $^ $(LIBS) -o $@

%.o : %.c $(HEAD)
	$(CC) $(CFLAGS) -c $<

run: all
	$(CHECK) mpirun -np $(prc) ./main.x $(choose) $(try)
	
clean: 
	rm -f *.o *.x *~
	
format: $(SRC)
	@clang-format -i $^ -verbose || echo "Please install clang-format to run this command"
