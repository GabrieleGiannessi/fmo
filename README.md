# Fluence Map Optimization (FMO) - NSGA-II

A multi-objective Fluence Map Optimization framework for Intensity-Modulated Radiation Therapy (IMRT) using NSGA-II, implemented in C++ across Sequential, OpenMP, and FastFlow execution models.

---

## 1. Prerequisites & Dependencies

### System Packages (Ubuntu/Debian)
```bash
sudo apt update && sudo apt install -y build-essential cmake zlib1g-dev libomp-dev git curl unzip
```

### Git Dependencies
Clone the required third-party repositories into the **parent directory** (`..`) of `fmo`:

```bash
cd ..
git clone https://gitlab.com/libeigen/eigen.git
git clone https://github.com/tbeu/matio.git
git clone https://github.com/fastflow/fastflow.git # if not already present
cd fmo
```

### Dataset
Download the required influence matrix datasets into `../data`:

```bash
bash scripts/download_data.sh
```

---

## 2. Build

Compile all targets using the build script:

```bash
bash scripts/build.sh
```

*(Alternatively: `mkdir -p build && cd build && cmake .. && make -j$(nproc) && cd ..`)*

---

## 3. Testing & Execution

### Automated Tests
Run unit tests and CLI verification via CTest:

```bash
ctest --test-dir build --output-on-failure
```

### Running NSGA-II Implementations

Execute commands from the `fmo` root directory to ensure dataset resolution (`../data/prostate`).

#### 1. Sequential (`nsga2-seq`)
```bash
./build/nsga2-seq [pop_size] [generations]

# Example (pop=256, gen=50):
./build/nsga2-seq 256 50
```

#### 2. OpenMP (`nsga2-omp`)
```bash
./build/nsga2-omp <workers> [pop_size] [generations]

# Example (8 threads, pop=256, gen=50):
./build/nsga2-omp 8 256 50
```

#### 3. FastFlow (`nsga2-ff`)
Modes: `0` = ParallelFor (`parfor`), `1` = Farm (`farm`).
```bash
./build/nsga2-ff <workers> <mode> [pop_size] [generations]

# FastFlow ParallelFor (8 workers):
./build/nsga2-ff 8 0 256 50

# FastFlow Farm (8 workers):
./build/nsga2-ff 8 1 256 50
```

---

## 4. Automated Benchmarking

Run the automated scalability and stability benchmark across all variants:

```bash
# Quick smoke test (OpenMP, FastFlow parfor/farm, Sequential):
./scripts/benchmark.sh -v all -w "1 2 4 8" -r 1 -p 64 -g 10

# Full benchmark run:
./scripts/benchmark.sh
```
