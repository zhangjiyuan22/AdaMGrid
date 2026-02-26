
# AdaMGrid

#### Adaptive magnification-map grid search for binary/triple microlensing on CPU/GPU  
#### *(Binary-lens grid search on CPU is released; triple-lens grid search on GPU is in development)*

Multiple-lens microlensing events (binary and triple lenses) are key probes of cold planets but are difficult to model because their light curves depend on a complex, highly degenerate parameter space. Robust interpretation therefore requires near-full parameter-space searches that are already computationally demanding for binary lenses and, given the larger number of parameters and slower magnification evaluation, not yet feasible for triple lenses. Precomputing reusable magnification maps on a source-plane grid and then evaluating magnifications via interpolation can in principle accelerate searches, but as magnification varies non-uniformly across the source plane, uniform grids are either too coarse near caustics or waste most resolution in smooth regions.

We present AdaMGrid, a grid-search framework that evaluates light curves by interpolating on precomputed, adaptively refined magnification maps. Starting from a coarse source-plane grid, AdaMGrid iteratively subdivides only those cells where the discrepancy between interpolated and numerically computed magnifications exceeds a threshold, yielding an adaptive grid that is dense around caustics and sparse in smooth regions. By tying this threshold to a fixed fraction of the expected photometric noise, AdaMGrid keeps interpolation uncertainties well below the photometric uncertainties across the source plane and achieves ≲0.1% relative precision in high-magnification regions, while accelerating magnification evaluation by three orders of magnitude relative to contour integration.

Applied to binary-lens events, AdaMGrid reduces a grid search to ≲30 minutes on a 100–core node using CPU-based interpolation and MCMC, and a forthcoming GPU implementation should reduce this to a few minutes. We further outline a path toward triple-lens parameter space exploration: generating of order 10⁸ reusable triple-lens adaptive magnification maps with VBMicrolensing on a CPU cluster, then exploiting GPU-accelerated interpolation and MCMC to bring triple-lens grid search to timescales of several days per event. These capabilities would enable systematic searches for planets in binary systems in archival data, robust re-analyses of known and candidate triple-lens events, and scalable modeling pipelines for upcoming microlensing surveys with Roman, ET, and other facilities. AdaMGrid is released publicly.

## Generate a resuable binary-lens magnification map set and run a binary grid search
(point out we only need to generate map set for one time, and reuse it for each binary lens grid search)

### Prerequisites
- Python **≥ 3.7**
- Python packages: `numpy`, `scipy`, `matplotlib`, `ctypes`, `emcee`, `PyAstronomy`
- `VBMicrolensing` (installed via pip; see below)

---

### 1) Clone the repository
```bash
git clone https://github.com/zhangjiyuan22/AdaMGrid.git
cd AdaMGrid
````

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

---

### 4) Generate the reusable magnification map set

```bash
python adaptive_map_generator_VBMicrolensingPython_BinaryMag2.py
```

#### 4.1 Default map-set configuration

With the default settings, the script generates a grid of maps over:

* `logs ∈ [-1.5, 1.5]` with `d_logs = 0.05`
* `logq ∈ [-6.0, 4.0]` with `d_logq = 0.1`
* `logrho ∈ [-4.0, -1.6]` with `d_logrho = 0.3`

Each map is:

* a **7 × 7 θE** square,
* centered on the **magnification center**.

For points **outside** the map boundary, the interpolation code falls back to a **single-lens approximation**.

**Resource footprint (default):**

* Total map-set size: **~63 GB**
* Typical runtime: **~4 hour on a 96-core node**(might be a tail that a tiny portion of map need much longer time than others, so you can terminate the program earlier like if 99.9% maps have been generated; by checking the real-time file number in the ./map_set_VBMicrolensing5p0Python_logs_minus1p5_to_1p5_dlogs_0p05_logq_minus6_to_4_dlogq_0p1_logrho_minus4p0_to_minus1p6_dlogrho_0p3_layer_16_boxsize_3p5/ through ls -l | grep "^-" | wc -l command, the total number of map for the default map set is 33993)

#### 4.2 Configuration edits (apply to Step 4)

##### Common edits (recommended)

* **Change the map-set name:** edit **line 14**
* **Change parameter ranges / resolutions:** edit **lines 224–231**

  * `logs` range and step
  * `logq` range and step
  * `logrho` range and step
* **Change the number of CPU cores:** edit **line 301**

##### Advanced controls (usually not needed)

* **Change map size:** edit **line 12** *(0.5 × map side length in units of θE)*
* **Change interpolation accuracy target:** edit **line 104** via `threshold_coefficient`

  * Target condition: `|A_interp - A_VBML| < threshold_coefficient * sqrt(A)`
* **Tune VBMicrolensing numerical settings:** edit **lines 141–143**

  * absolute tolerance
  * relative tolerance
  * limb-darkening coefficient
* **Change the maximum refinement depth (map layers):** edit **line 173** (`max_layer`)
  * Note: up to **16 layers** (with a 7 × 7 θE map) yields a finest spatial resolution of `dx ≈ 7 / 2^16 ≈ 1e-4 θE`

 To run a grid search for a binary lens event: 
 first prepare the light curve data, put into ./data/xxxxxx(name_of_your_event)/; you can bin the data for non-anomaly region in advance to speed up the grid search; 
 then for each event copy a new binary_grid_xxxxxx(name_of_your_event)_fast_struct.py from the current binary_grid_kb240697_fast_struct.py; 
 in that file, change which map set to use at line 14; 
 change alpha initial guess at line 205 and 206, currently 16 alpha initial guesses equally spaced in [0,2pi); 
 change MCMC nburn_in and nsample at line 223 and 224, default to 300+500 step; 
 chenge the event name at line 294, change event's alpha/dec at line 295/296 (in degree), change light curve files at line 298, the format should be in column 0/1/2 is hjd/magnitude/magnitude_error, while hjd-2450000 or hjd-2400000 is automatically converted, while hjd/flux/flux_error is also allowed by repeting the light curve file name at line 299, and jd is allowed by specifying at line 304; set the k and emin (Yee et al. 2012) for each light curve file at line 301 and 302; 

 change the initial guess of t0/u0/tE at line 308-310, which typically from single fit with the anomaly excluded; for some stellar binary events that u0 and tE can not be robustly obtained from single fit, multiple grid search with different u0 intial guess might be needed; 

 change the time window for the light curve data used for grid search at line 314-315, which is in format of HJD-2450000; 

 change whether use fmin+MCMC or just fmin at line 318; using fmin+MCMC is dozens times slower than just fmin (but fmin+MCMC is still already fast, like 30 min per event), but is strongly recommended, as fmin may fail catastrophyilly for bad initial guess, while MCMC is much more robust; the nburn_in and nsample of MCMC might be adjusted (line 223 and 224); we recommend hot MCMC, change the hot MCMC factor at line 321, where factor=0.1 means the chi2 is multiplied by 0.1;

 change the logs/logq/logrho range and resolution used from line 336 to 342

 change the number of CPU core used and file's save name at line 481 and 488
 
 

