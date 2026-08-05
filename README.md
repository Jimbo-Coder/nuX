<!-- <img align="top" src="Docs/figures/nux.png" width="140"> -->

[![CI](https://github.com/jaykalinani/nuX/actions/workflows/ci.yml/badge.svg?branch=odesolvers)](https://github.com/jaykalinani/nuX/actions/workflows/ci.yml)
[![AMReX](https://amrex-codes.github.io/badges/powered%20by-AMReX-red.svg)](https://amrex-codes.github.io)
[![Ask DeepWiki](https://deepwiki.com/badge.svg)](https://deepwiki.com/jaykalinani/nuX)

**nuX** is a GPU-accelerated neutrino-transport module for dynamical spacetimes, written in C++ and designed to run within the [Einstein Toolkit](https://einsteintoolkit.org/) using the [CarpetX](https://github.com/EinsteinToolkit/CarpetX) driver. **CarpetX** itself is built atop [AMReX](https://amrex-codes.github.io), a framework for block-structured adaptive mesh refinement (AMR).  
**nuX** provides a two-moment (M1) evolution of radiation variables with analytic closures and multi-species support, interfacing with our GRMHD code [AsterX](https://github.com/EinsteinToolkit/AsterX).

---

## Overview

- Two-moment (**M1**) neutrino transport with analytic closures (e.g., Minerbo/Levermore)  
- Port of the `THC_M1` code of Radice et al. 2022
- Multi-species support (typically $ν_e$, $\barν_e$, $ν_x$)  
- Source terms for neutrino–matter coupling and stress–energy feedback ($T^{μν}_{\rm rad}$) coupled with **AsterX** 
- Robust floors, masking, and diagnostics for production-quality BNS/CCSN simulations

---

## Available Modules

- `nuX_M1` — Core M1 evolution: fluxes, closures, source terms, analysis/diagnostics, stress–energy output
- `nuX_NuRates` — Physical interaction rates & opacities via `bns_nurates`
- `nuX_WeakRates` — Independent GPU-capable WeakRates implementation for transport and THC_M1 cross-validation
- `nuX_FakeRates` — Lightweight fake/constant rates for testing
- `nuX_NuRatesToy` — Toy module for testing
- `nuX_Seeds` — Initial data for radiation fields as well as MHD variables
- `nuX_Utils` — Tensor utilities, metric helpers, math wrappers

---

## Getting Started

* Instructions for downloading and building nuX with the Einstein Toolkit are available [here](https://github.com/EinsteinToolkit/CarpetX/wiki/Getting-Started).  
* Thornlist is available [here](https://github.com/jaykalinani/nuX/blob/development/Docs/thornlist/nuX.th)
* Simfactory files for various clusters and setup instructions can be found [here](https://github.com/lwJi/ETK-Compile-Guides).    

---

## Useful Repositories

* [AsterX](https://github.com/EinsteinToolkit/AsterX) - GRMHD code
* [CarpetX](https://github.com/EinsteinToolkit/CarpetX) – Next-generation driver for the Einstein Toolkit  
* [SpacetimeX](https://github.com/EinsteinToolkit/SpacetimeX) – Modules for spacetime evolution 
* [BNSTools](https://github.com/jaykalinani/BNSTools) – Utilities supporting BNS merger simulations  

---

## Additional Resources

* Radice et al 2022: [MNRAS](https://academic.oup.com/mnras/article/512/1/1499/6542449) [arXiv](https://arxiv.org/abs/2111.14858)
* Shibata et al 2021: [JPS](https://academic.oup.com/ptp/article/125/6/1255/1861577) [arXiv](https://arxiv.org/abs/1104.3937)
* Introductory notebooks on nuX_M1, truncated moments formalism [here](https://github.com/jaykalinani/nuX_Docs) 
