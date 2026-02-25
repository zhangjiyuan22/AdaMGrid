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
4. python adaptive_map_generator_VBMicrolensingPython_BinaryMag2.py to generate resuable binary-lens magnification map set; the default setting leads to a map set with logs in [-1.5, 1.5] with d_logs = 0.05, logq in [-6.0, 4.0] with d_logq = 0.1, logrho in [-4.0, -1.6] with d_logrho = 0.3; each map is a square of 7*7 thetaE; the default map set has size of 63 GB



