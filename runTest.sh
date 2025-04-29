#!/bin/bash

echo "Enter matrix size (N for NxN): "
read N

echo "Run sequential also? (1 = yes, 0 = no): "
read USE_SEQ

echo "Enter number of MPI processes to use: "
read PROCS

echo "Which test to run?"
echo "1. multiplyTest"
echo "2. refTest"
echo "3. rrefTest"
echo "4. luTest"
read CHOICE

# Pick executable
case $CHOICE in
  1) EXEC="./multiplyTest" ;;
  2) EXEC="./refTest" ;;
  3) EXEC="./rrefTest" ;;
  4) EXEC="./luTest" ;;
  *) echo "Invalid choice"; exit 1 ;;
esac

echo
echo "Running $EXEC with:"
echo "- Matrix size: $N"
echo "- Sequential: $USE_SEQ"
echo "- Processes: $PROCS"
echo

# Feed stdin input to the executable using a here string
mpirun --oversubscribe -np $PROCS $EXEC <<< "$N
$USE_SEQ"
