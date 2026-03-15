# Matrix — Linear Algebra in C with OpenMPI Benchmarks

This project implements core matrix operations in pure C with OpenMPI-aware test
drivers for matrix multiplication, REF, RREF, and LU decomposition. The matrix
storage is now heap-backed, so the old fixed `2100 x 2100` cap is gone and the
practical limit is determined by available RAM, MPI count limits, and runtime.

The repository includes the core implementation, interactive MPI test binaries,
and a repeatable benchmark suite that generates CSV timing tables, PNG figures,
and a short LaTeX report.

---

## Project Structure

```
Matrix/
│
├── benchmarks/       # New benchmark data, figures, scripts, and LaTeX report
├── src/              # Core library: MatrixC.c, MatrixC.h
├── tests/            # Standalone test drivers (luTest.c, refTest.c, etc.)
├── results/          # Older timing spreadsheets (.ods)
├── figs/             # Older performance plots (.png)
├── docs/             # Original project report and slides
├── makefile          # Build rules for all test programs
└── README.md         # This file
```

---

## Features

### **Matrix Operations**
- Matrix multiplication  
- Row-Echelon Form (REF)  
- Reduced Row-Echelon Form (RREF)  
- LU decomposition  
- Heap-backed matrix storage with allocation checks

### **Design Philosophy**
- Pure C with explicit memory management  
- Contiguous row-major storage for cache-friendly access  
- Clean separation between library and tests  
- Easy to benchmark and extend into stronger MPI implementations

---

## Compilation & Usage

### **Build All Tests**
Run:

```bash
make
```

### How to Run Tests

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

### Benchmark Suite

To regenerate the new benchmark data, plots, and summary CSV files, run:

```bash
python3 benchmarks/scripts/run_benchmarks.py
```

This benchmark sweep currently covers:

- matrix sizes `250`, `500`, `1000`, and `1500`
- MPI process counts `1`, `2`, `4`, `8`, and `16`
- all four shipped executables

Generated artifacts are written to:

- `benchmarks/data/` for timing tables and raw logs
- `benchmarks/figures/` for timing and speedup plots
- `benchmarks/report/` for the LaTeX benchmark report

### **Test Executables**

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

This repository now contains two sets of performance artifacts:

- the original spreadsheets and figures in `results/` and `figs/`
- the new reproducible benchmark suite in `benchmarks/`

### New Benchmark Outputs

Located in `benchmarks/`:

- `data/timings.csv` — all benchmark rows with kernel time, wall time, and log path
- `data/summary.csv` — best MPI result at each tested matrix size
- `data/logs/` — raw stdout/stderr from each benchmark run
- `figures/multiply_timing.png` — multiplication timing curves
- `figures/ref_timing.png` — REF timing curves
- `figures/rref_timing.png` — RREF timing curves
- `figures/lu_timing.png` — LU timing curves
- `figures/speedup_summary.png` — combined speedup overview
- `report/Benchmark_Timing_Report.tex` — short LaTeX benchmark writeup
- `report/Benchmark_Timing_Report.pdf` — compiled benchmark report

### Original Results

Located in `results/` and `figs/`:

- `multiply.ods`
- `ref.ods`
- `rref.ods`
- `lu.ods`
- `matrix.png`
- `ref.png`
- `rref.png`
- `lu.png`

The new benchmark suite uses the timings printed by the test executables
themselves, which makes it easy to rerun after code changes.

---

## Documentation

Documentation is split between the original project writeup in `docs/` and the
new benchmark report in `benchmarks/report/`.

### Original Documentation in `docs/`

- **Matrix_Project_Report.pdf** — Complete writeup  
- **Matrix_Project_Slides.pdf** — Presentation / overview  
- **Matrix_Project_Report.tex** — LaTeX source  

### New Benchmark Documentation in `benchmarks/report/`

- **Benchmark_Timing_Report.pdf** — Short benchmark summary
- **Benchmark_Timing_Report.tex** — LaTeX source

---

## Source Code Details

The core of the project lives in `src/MatrixC.c` and `src/MatrixC.h`.

### Highlights:

- `Matrix` uses heap-backed contiguous storage.
- Each operation is isolated and testable.
- The MPI paths allocate local work buffers based on the real matrix size.
- Large temporary buffers have been moved off the stack.
- The benchmark suite makes it easy to rerun timing studies after changes.

---

## Test Programs

Each test program in `tests/` focuses on one operation:

- `luTest.c`  
- `refTest.c`  
- `rrefTest.c`  
- `multiplyTest.c`

These programs validate the core matrix operations (LU, REF, RREF, and multiplication).

---

## Future Work

- Add a truly distributed MPI LU decomposition
- Add partial pivoting to LU, REF, and RREF
- Add assertion-based correctness tests
- Expand the benchmark suite to larger problem sizes and repeated trials
- Add CSV-to-LaTeX table generation for the benchmark report

---

## License

This project is licensed under the **MIT License**. See `LICENSE` for details.
