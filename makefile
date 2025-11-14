CC      = mpicc
CFLAGS  = -Wall -Wextra -O2

SRC_DIR   = src
TEST_DIR  = tests

MATRIX_SRC = $(SRC_DIR)/MatrixC.c
MATRIX_HDR = $(SRC_DIR)/MatrixC.h

TESTS = luTest multiplyTest refTest rrefTest

all: $(TESTS)

luTest: $(MATRIX_SRC) $(MATRIX_HDR) $(TEST_DIR)/luTest.c
	$(CC) $(CFLAGS) -o luTest $(MATRIX_SRC) $(TEST_DIR)/luTest.c

multiplyTest: $(MATRIX_SRC) $(MATRIX_HDR) $(TEST_DIR)/multiplyTest.c
	$(CC) $(CFLAGS) -o multiplyTest $(MATRIX_SRC) $(TEST_DIR)/multiplyTest.c

refTest: $(MATRIX_SRC) $(MATRIX_HDR) $(TEST_DIR)/refTest.c
	$(CC) $(CFLAGS) -o refTest $(MATRIX_SRC) $(TEST_DIR)/refTest.c

rrefTest: $(MATRIX_SRC) $(MATRIX_HDR) $(TEST_DIR)/rrefTest.c
	$(CC) $(CFLAGS) -o rrefTest $(MATRIX_SRC) $(TEST_DIR)/rrefTest.c

.PHONY: clean
clean:
	rm -f $(TESTS)
