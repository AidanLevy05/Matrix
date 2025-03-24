# Variables
PROG = prog
CC = mpicc 
CFLAGS = -g -Wall -std=c11
OBJS = main.o MatrixC.o

# Default target
all: $(PROG)

# Link objects to create executable
$(PROG): $(OBJS)
	$(CC) $(CFLAGS) -o $(PROG) $(OBJS) -lm

# Compile main.c to main.o
main.o: main.c MatrixC.h
	$(CC) $(CFLAGS) -c main.c

# Compile MatrixC.c to MatrixC.o
MatrixC.o: MatrixC.c MatrixC.h
	$(CC) $(CFLAGS) -c MatrixC.c

# Clean object files and executable
clean:
	rm -f $(OBJS) $(PROG)

# Rebuild the program
rebuild:
	make clean 
	make
