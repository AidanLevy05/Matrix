# Matrix — High-Performance Linear Algebra in C (MPI-Ready)

This project implements core matrix operations in pure C using fixed-size static arrays for deterministic performance, predictable memory access patterns, and easy integration with OpenMPI. The library supports fundamental linear algebra routines including matrix multiplication, REF, RREF, and LU components, and is designed as the computational foundation for distributed parallel algorithms.

The project is structured for clarity and scalability, separating implementation, tests, documentation, timing results, and figures.

---

## Project Structure

```
Matrix/
│
├── src/              # Core library: MatrixC.c, MatrixC.h
├── tests/            # Standalone test drivers (luTest.c, refTest.c, etc.)
├── results/          # Timing tables, speedup data (.ods files)
├── figs/             # Performance plots and figures (.png)
├── docs/             # Project report, slides (.tex, .pdf)
├── makefile          # Build rules for all test programs
└── README.md         # This file
```

---

## Features

### **Matrix Operations**
- Matrix multiplication  
- Row-Echelon Form (REF)  
- Reduced Row-Echelon Form (RREF)  
- Partial LU component routines  
- Deterministic static-array structure (no dynamic allocation)

### **Design Philosophy**
- No C++ or malloc — pure C for MPI portability  
- Predictable memory layout for cache-friendly parallelization  
- Clean separation between library and tests  
- Easy to extend into multi-process MPI execution

---

## Compilation & Usage

### **Build All Tests**
Run:

```bash
make
```

### **Run Specific Test Executables**

```bash
./luTest
./refTest
./rrefTest
./multiplyTest
```

(These executables are generated from the sources in `tests/`.)

### **Clean Build Artifacts**

```bash
make clean
```

---

## Performance Analysis

This project includes real timing and scaling measurements collected across multiple matrix sizes.

### Included Results:
Located in `results/`:

- `multiply.ods` — Multiplication timings & speedups  
- `ref.ods`       — REF benchmarks  
- `rref.ods`      — RREF benchmarks  
- `lu.ods`        — LU decomposition–related timings  

### Included Plots:
Located in `figs/`:

- `matrix.png` — Overall performance summary  
- `ref.png`    — REF timing curve  
- `rref.png`   — RREF timing curve  
- `lu.png`     — LU component timing  

These results were produced on real hardware using optimized GCC builds.

---

## Documentation

All formal documentation is in `docs/`:

- **Matrix_Project_Report.pdf** — Complete writeup  
- **Matrix_Project_Slides.pdf** — Presentation / overview  
- **Matrix_Project_Report.tex** — LaTeX source  

These files describe:
- Mathematical background  
- Algorithm descriptions  
- Experimental design  
- Timing results  
- Speedup/efficiency interpretation  

---

## Source Code Details

The core of the project lives in `src/MatrixC.c` and `src/MatrixC.h`.

### Highlights:

- `Matrix` struct uses static arrays for deterministic size.
- Each operation is isolated and testable.
- Functions are optimized for predictable loops and minimal branching.
- Designed to be extended into distributed operations with MPI:
  - Scatter blocks  
  - Local computations  
  - Gather results  

---

## Test Programs

Each test program in `tests/` focuses on one operation:

- `luTest.c`  
- `refTest.c`  
- `rrefTest.c`  
- `multiplyTest.c`

These programs validate the core matrix operations (LU, REF, RREF, and multiplication).

## How to Run Tests

All testing is handled through the unified script:

```bash
./runTest.sh
```

This script:
- prompts for matrix size
- selects sequential or parallel mode
- asks how many MPI processes to use 
- allows choosing which operation to test 
- launches the selected program with `mpirun`

This is the only recommended way to run and test in this project.

---

## Future Work

- Add OpenMPI parallel versions of all operations  
- Integrate row-blocking for cache-efficient multiplication  
- Add full LU decomposition with pivoting  
- Create automated test suite with expected outputs  
- Migrate timing collection into Python benchmarking scripts  

---

## License

This project is licensed under the **MIT License**. See `LICENSE` for details.

