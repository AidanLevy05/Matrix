CC = mpicc
CFLAGS = -g -Wall -std=c11
OBJS = MatrixC.o

all: multiplyTest refTest rrefTest luTest

multiplyTest: multiplyTest.o $(OBJS)
	$(CC) $(CFLAGS) -o multiplyTest multiplyTest.o $(OBJS) -lm

refTest: refTest.o $(OBJS)
	$(CC) $(CFLAGS) -o refTest refTest.o $(OBJS) -lm

rrefTest: rrefTest.o $(OBJS)
	$(CC) $(CFLAGS) -o rrefTest rrefTest.o $(OBJS) -lm

luTest: luTest.o $(OBJS)
	$(CC) $(CFLAGS) -o luTest luTest.o $(OBJS) -lm

multiplyTest.o: multiplyTest.c MatrixC.h
	$(CC) $(CFLAGS) -c multiplyTest.c

refTest.o: refTest.c MatrixC.h
	$(CC) $(CFLAGS) -c refTest.c

rrefTest.o: rrefTest.c MatrixC.h
	$(CC) $(CFLAGS) -c rrefTest.c

luTest.o: luTest.c MatrixC.h
	$(CC) $(CFLAGS) -c luTest.c

MatrixC.o: MatrixC.c MatrixC.h
	$(CC) $(CFLAGS) -c MatrixC.c

clean:
	rm -f *.o multiplyTest refTest rrefTest luTest

rebuild:
	make clean
	make