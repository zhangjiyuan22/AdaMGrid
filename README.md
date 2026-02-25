# AdaMGrid

### Adaptive Magnification Map Grid Search code for binary/triple lensing running on CPU/GPU 

### (binary lensing grid search on CPU is released; triple grid on GPU is in development)

Multiple-lens microlensing events (binary and triple lenses) are key probes of cold planets but are difficult to model because their light curves depend on a complex, highly degenerate parameter space. Robust interpretation therefore requires near-full parameter-space searches that are already computationally demanding for binary lenses and, given the larger number of parameters and slower magnification evaluation, not yet feasible for triple lenses. Precomputing reusable magnification maps on a source-plane grid and then evaluating magnifications via interpolation can in principle accelerate searches, but as magnification varies non-uniformly across the source plane, uniform grids are either too coarse near caustics or waste most resolution in smooth regions. 

We present AdaMGrid, a grid-search framework that evaluates light curves by interpolating on precomputed, adaptively refined magnification maps. Starting from a coarse source-plane grid, AdaMGrid iteratively subdivides only those cells where the discrepancy between interpolated and numerically computed magnifications exceeds a threshold, yielding an adaptive grid that is dense around caustics and sparse in smooth regions. By tying this threshold to a fixed fraction of the expected photometric noise, AdaMGrid keeps interpolation uncertainties well below the photometric uncertainties across the source plane and achieves ≲0.1% relative precision in high-magnification regions, while accelerating magnification evaluation by three orders of magnitude relative to contour integration. 

Applied to binary-lens events, AdaMGrid reduces a grid search to ≲30 minutes on a 100–core node using CPU-based interpolation and MCMC, and a forthcoming GPU implementation should reduce this to a few minutes. We further outline a path toward triple-lens parameter space exploration: generating of order 10⁸ reusable triple-lens adaptive magnification maps with VBMicrolensing on a CPU cluster, then exploiting GPU-accelerated interpolation and MCMC to bring triple-lens grid search to timescales of several days per event. These capabilities would enable systematic searches for planets in binary systems in archival data, robust re-analyses of known and candidate triple-lens events, and scalable modeling pipelines for upcoming microlensing surveys with Roman, ET, and other facilities. AdaMGrid is released publicly. 

## How to generate a binary-lens magnification map set and run a binary grid search

### Prerequisites

* Python **≥ 3.7**
* Python packages: `numpy`, `matplotlib`, `ctypes`, `emcee`
* `VBMicrolensing` (installed via pip; see below)

---

### 1) Clone the repository

```bash
git clone https://github.com/zhangjiyuan22/AdaMGrid.git
cd AdaMGrid
```

Make sure your Python environment has the required packages (`numpy/matplotlib/ctypes/emcee`) installed.

---

### 2) Install VBMicrolensing

```bash
pip install VBMicrolensing
```

---

### 3) Compile the C interpolation code

```bash
chmod a+x compile_all
./compile_all
```

This builds the interpolation code.

---

### 4) Generate the reusable magnification map set

```bash
python adaptive_map_generator_VBMicrolensingPython_BinaryMag2.py
```

#### Default map-set configuration

With the default settings, the script generates a grid of maps over:

* `logs ∈ [-1.5, 1.5]` with `d_logs = 0.05`
* `logq ∈ [-6.0, 4.0]` with `d_logq = 0.1`
* `logrho ∈ [-4.0, -1.6]` with `d_logrho = 0.3`

Each map is:

* a **7 × 7 θE** square,
* centered on the **magnification center**.

Later for points **outside** the map boundary, the interpolation code falls back to a **single-lens approximation**.

**Resource footprint (default):**

* Total map-set size: **~63 GB**
* Typical runtime: **~1 hour on a 64-core node**

---

## Common configuration edits (recommended)

### Change the map-set name

* Edit **line 14**

### Change parameter ranges / resolutions

* Edit **lines 224–231** to modify:

  * `logs` range and step
  * `logq` range and step
  * `logrho` range and step

### Change the number of CPU cores

* Edit **line 301**

---

## Advanced controls (usually not needed)

### Change map size

* Edit **line 12** (0.5 * map side length in units of θE)

### Change interpolation accuracy target

* Edit **line 104** via `threshold_coefficient`

The adaptive refinement aims to satisfy:

* `|A_interp - A_VBML| < threshold_coefficient * sqrt(A)`

### Tune VBMicrolensing numerical settings

* Edit **lines 141–143** to adjust:

  * absolute tolerance
  * relative tolerance
  * limb-darkening coefficient
    
### Change the maximum refinement depth (map layers)

* Edit **line 173** (`max_layer`)

Notes:

* Up to **16 layers** (with a 7 × 7 θE map) yields a finest spatial resolution of: `dx ≈ 7 / 2^16 ≈ 1e-4 θE`

