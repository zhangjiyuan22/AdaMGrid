# AdaMGrid

### Adaptive Magnification Map Grid Search code for binary/triple lensing running on CPU/GPU 

### (binary lensing grid search on CPU is released; triple grid on GPU is in development)

Multiple-lens microlensing events (binary and triple lenses) are key probes of cold planets but are difficult to model because their light curves depend on a complex, highly degenerate parameter space. Robust interpretation therefore requires near-full parameter-space searches that are already computationally demanding for binary lenses and, given the larger number of parameters and slower magnification evaluation, not yet feasible for triple lenses. Precomputing reusable magnification maps on a source-plane grid and then evaluating magnifications via interpolation can in principle accelerate searches, but as magnification varies non-uniformly across the source plane, uniform grids are either too coarse near caustics or waste most resolution in smooth regions. 

We present AdaMGrid, a grid-search framework that evaluates light curves by interpolating on precomputed, adaptively refined magnification maps. Starting from a coarse source-plane grid, AdaMGrid iteratively subdivides only those cells where the discrepancy between interpolated and numerically computed magnifications exceeds a threshold, yielding an adaptive grid that is dense around caustics and sparse in smooth regions. By tying this threshold to a fixed fraction of the expected photometric noise, AdaMGrid keeps interpolation uncertainties well below the photometric uncertainties across the source plane and achieves ≲0.1% relative precision in high-magnification regions, while accelerating magnification evaluation by three orders of magnitude relative to contour integration. 

Applied to binary-lens events, AdaMGrid reduces a grid search to ≲30 minutes on a 100–core node using CPU-based interpolation and MCMC, and a forthcoming GPU implementation should reduce this to a few minutes. We further outline a path toward triple-lens parameter space exploration: generating of order 10⁸ reusable triple-lens adaptive magnification maps with VBMicrolensing on a CPU cluster, then exploiting GPU-accelerated interpolation and MCMC to bring triple-lens grid search to timescales of several days per event. These capabilities would enable systematic searches for planets in binary systems in archival data, robust re-analyses of known and candidate triple-lens events, and scalable modeling pipelines for upcoming microlensing surveys with Roman, ET, and other facilities. AdaMGrid is released publicly. 

## How to generate binary-lens magnification map set and run a binary grid search?
0. git clone https://github.com/zhangjiyuan22/AdaMGrid.git, and make sure you have numpy/matplotlib/ctypes/emcee/ installed
1. pip install VBMicrolensing (require python >= 3.7)
2. chmod a+x compile_all
3. ./compile_all (to compile the interpolation code written in C)
4. python adaptive_map_generator_VBMicrolensingPython_BinaryMag2.py to generate resuable binary-lens magnification map set;
  
the default setting leads to a map set with logs in [-1.5, 1.5] with d_logs = 0.05, logq in [-6.0, 4.0] with d_logq = 0.1, logrho in [-4.0, -1.6] with d_logrho = 0.3; each map is a square of 7*7 thetaE, centered at magnification center; for points outside the map, we use single lens approximation; the default map set has size of 63 GB, and need 1 hour to generate on a 64-core node

change map set name at line 14; change logs/logq/logrho range and resolution from line 224 to 231; change number of CPU core used at line 301

detail control (generally no need to change): change map size at line 12; change expected interpolation accuracy of map at line 104 by changing threshold_coefficient, where map will satisfy |A_interpolation - A_VBML| < threshold_coefficient * sqrt(A); change absolute tolerance, relative tolerance, and limb-darkening coefficients of VBML from line 141 to 143; change max layer of map at line 173 (at most 16 layer together with map size of 7*7 thetaE, leads to map finest x/y resolution of 7/(2^16) = 1e-4 thetaE) 

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

* Edit **line 12** (map side length in units of θE)

### Change interpolation accuracy target

* Edit **line 104** via `threshold_coefficient`

The adaptive refinement aims to satisfy:
[
|A_{\rm interp} - A_{\rm VBML}| < \texttt{threshold_coefficient},\sqrt{A}.
]

### Tune VBMicrolensing numerical settings

* Edit **lines 141–143** to adjust:

  * absolute tolerance
  * relative tolerance
  * limb-darkening coefficients

### Change the maximum refinement depth (map layers)

* Edit **line 173** (`max_layer`)

Notes:

* Up to **16 layers** (with a 7 × 7 θE map) yields a finest spatial resolution of:
  [
  \Delta x \approx \frac{7}{2^{16}} \approx 1\times10^{-4},\theta_{\rm E}.
  ]

---

If you want, I can also rewrite this into a more “README-like” style with a short overview, directory structure, and a minimal “quick start” block at the top.



