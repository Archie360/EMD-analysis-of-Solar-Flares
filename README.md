
## Why use Empirical Mode Decomposition for flares?

- Solar and stellar flare light curves are typically non-stationary and may contain multiple overlapping processes: impulsive rise, quasi-periodic pulsations (QPPs), oscillatory relaxation, and slow trends.
- EMD is a data-driven, adaptive method that decomposes a time series into intrinsic mode functions (IMFs) without imposing a fixed basis (like Fourier or wavelets). This makes it well suited to capture transient oscillations and evolving frequencies intrinsic to flare signals.

## scientific workflow

1. Select event and data: choose a flare with sufficient signal-to-noise and cadence to resolve the timescales of interest.
2. Preprocess: correct for instrumental effects, handle gaps, optionally detrend or normalize (but keep originals for comparison).
3. Decompose: run EMD (or EEMD) to obtain IMFs + residual trend.
4. Analyze IMFs: compute instantaneous amplitude/frequency (via Hilbert transform), power spectra, and map IMFs to physical timescales.
5. Validate: use surrogate testing (red/white noise surrogates), cross-validate with wavelet analysis, and check robustness against EMD parameters.
6. Report: present stacked IMFs, Hilbert spectrum, marginal spectrum, and statistical significance for any claimed periodicities.


## EMD-specific choices that matter scientifically

- Stopping criterion for sifting: the convergence threshold (e.g., standard deviation criterion) affects IMF shapes — report it.
- Ensemble EMD (EEMD) or CEEMDAN: helps mitigate mode mixing; report ensemble size and added-noise amplitude.
- End effects and padding: use mirror-padding or tapering and test sensitivity to padding choices.



## references

- Huang et al. 1998 — "The empirical mode decomposition and the Hilbert spectrum for nonlinear and non-stationary time series analysis".
- Wu & Huang 2009 — "Ensemble empirical mode decomposition: a noise-assisted data analysis method".
- Torrence & Compo 1998 — "A practical guide to wavelet analysis" (for cross-validation and significance testing).


