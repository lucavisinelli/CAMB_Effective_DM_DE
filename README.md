# CAMB Effective DM→DE Model (CAMB 1.6.7)

This repository implements the effective Dark Matter → Dark Energy (DM→DE) transition model described in:

**L. Visinelli (2019)**  
https://arxiv.org/abs/1906.11255

---

## Overview

This code provides a modern implementation of the DM→DE transition model within CAMB 1.6.x.

### Key features

- Analytic implementation of the equation of state `w(a)`
- Analytic dark energy density evolution `rho_DE(a)`
- Consistent perturbation evolution including `dw/d ln a`
- Compatible with CAMB 1.6.x
- Validated against CAMB's tabulated `w(a)` implementation

---

## Model

The equation of state is

```text
w(a) = w0 / (1 + (a/a_transition)^(-2/delta))
```

where:

- `w0` is the late-time equation of state, typically close to `-1`
- `a_transition` is the transition scale factor
- `delta` controls the transition width

The model interpolates between:

- dark matter-like behaviour at early times, `w(a) ≈ 0`
- dark energy-like behaviour at late times, `w(a) → w0`

For `w0 = -1`, the late-time limit is cosmological-constant-like.

---

## Files modified

The implementation modifies two CAMB files:

```text
fortran/DarkEnergyFluid.f90
camb/dark_energy.py
```

The Fortran side adds the analytic background and perturbation evolution, while the Python side exposes the new class:

```python
from camb.dark_energy import DarkEnergyDMDE
```

---

## Basic usage in Python

```python
import camb
from camb.dark_energy import DarkEnergyDMDE

pars = camb.CAMBparams()
pars.set_cosmology(H0=67.5, ombh2=0.0224, omch2=0.12)
pars.InitPower.set_params(As=2.1e-9, ns=0.965)

pars.DarkEnergy = DarkEnergyDMDE()
pars.DarkEnergy.set_params(
    w0=-1.0,
    a_transition=0.3,
    delta=0.3,
)

results = camb.get_results(pars)
print(results.get_derived_params())
```

---

## Validation

The analytic implementation was validated against CAMB's tabulated `w(a)` fluid implementation.

For the test case:

```text
w0 = -1.0
a_transition = 0.3
delta = 0.3
```

we find:

- background quantities agree at numerical precision
- CMB TT spectra agree at approximately the `10^-6` level

This confirms that the analytic Fortran implementation reproduces the tabulated CAMB fluid result while avoiding spline overhead.

---

## Example: DMDE vs ΛCDM

```python
import camb
from camb.dark_energy import DarkEnergyDMDE


def base_pars():
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=67.5, ombh2=0.0224, omch2=0.12)
    pars.InitPower.set_params(As=2.1e-9, ns=0.965)
    pars.set_for_lmax(2500, lens_potential_accuracy=1)
    return pars


pars_lcdm = base_pars()

pars_dmde = base_pars()
pars_dmde.DarkEnergy = DarkEnergyDMDE()
pars_dmde.DarkEnergy.set_params(
    w0=-1.0,
    a_transition=0.3,
    delta=0.3,
)

res_lcdm = camb.get_results(pars_lcdm)
res_dmde = camb.get_results(pars_dmde)

cl_lcdm = res_lcdm.get_cmb_power_spectra(pars_lcdm, CMB_unit="muK")["total"]
cl_dmde = res_dmde.get_cmb_power_spectra(pars_dmde, CMB_unit="muK")["total"]

for ell in [2, 10, 100, 500, 1000, 2000]:
    rel = (cl_dmde[ell, 0] / cl_lcdm[ell, 0] - 1) * 100
    print(f"ell={ell:4d}: TT difference = {rel:+.6f}%")
```

---

## Example Cobaya configuration

Below is a minimal Cobaya-style example showing how one could expose the DMDE parameters in an MCMC run. The exact likelihood block should be adapted to the data combination used in the analysis.

```yaml
theory:
  camb:
    extra_args:
      dark_energy_model: DMDE

params:
  H0:
    prior:
      min: 50
      max: 90
    ref:
      dist: norm
      loc: 67.5
      scale: 2.0
    proposal: 0.5

  ombh2:
    prior:
      min: 0.005
      max: 0.1
    ref:
      dist: norm
      loc: 0.0224
      scale: 0.0002
    proposal: 0.0001

  omch2:
    prior:
      min: 0.001
      max: 0.99
    ref:
      dist: norm
      loc: 0.12
      scale: 0.002
    proposal: 0.001

  logA:
    prior:
      min: 1.61
      max: 3.91
    ref:
      dist: norm
      loc: 3.05
      scale: 0.001
    proposal: 0.001

  ns:
    prior:
      min: 0.8
      max: 1.2
    ref:
      dist: norm
      loc: 0.965
      scale: 0.004
    proposal: 0.002

  tau:
    prior:
      min: 0.01
      max: 0.8
    ref:
      dist: norm
      loc: 0.055
      scale: 0.006
    proposal: 0.003

  DMDE_w0:
    prior:
      min: -1.0
      max: -0.3
    ref: -1.0
    proposal: 0.02

  DMDE_a_transition:
    prior:
      min: 0.001
      max: 0.8
    ref: 0.3
    proposal: 0.02

  DMDE_delta:
    prior:
      min: 0.02
      max: 2.0
    ref: 0.3
    proposal: 0.02

sampler:
  mcmc:
    Rminus1_stop: 0.01
    max_tries: 1000
```

### Important Cobaya note

The YAML block above assumes that the CAMB-Cobaya interface has been extended to pass the DMDE parameters to CAMB, e.g. by constructing

```python
from camb.dark_energy import DarkEnergyDMDE

pars.DarkEnergy = DarkEnergyDMDE()
pars.DarkEnergy.set_params(
    w0=DMDE_w0,
    a_transition=DMDE_a_transition,
    delta=DMDE_delta,
)
```

If this mapping has not yet been added, the model will work directly in Python CAMB, but Cobaya will not automatically know how to pass the new parameters.

---

## Installation notes

A development installation of CAMB is recommended:

```bash
git clone https://github.com/lucavisinelli/CAMB_Effective_DM_DE.git
cd CAMB_Effective_DM_DE
git checkout camb-1.6.7-dmde
python3 -m pip install -e . --no-build-isolation
```

A Fortran compiler such as `gfortran` is required.

---

## Citation

If you use this implementation, please cite:

```bibtex
@article{Visinelli:2019,
  author = {Visinelli, Luca},
  title = {Effective Dark Matter to Dark Energy Transition},
  eprint = {1906.11255},
  archivePrefix = {arXiv},
  primaryClass = {astro-ph.CO}
}
```

Paper link:

https://arxiv.org/abs/1906.11255

---

## Disclaimer

This is a research implementation intended for cosmological model comparison and parameter studies. Users should validate the model against their own numerical settings, priors, and likelihood choices before using it in production analyses.
