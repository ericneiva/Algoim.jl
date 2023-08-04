# Algoim

[![DOI](https://zenodo.org/badge/621924496.svg)](https://zenodo.org/badge/latestdoi/621924496)
[![Build Status](https://github.com/ericneiva/Algoim.jl/workflows/CI/badge.svg?branch=main)](https://github.com/ericneiva/Algoim.jl/actions?query=workflow%3ACI)

A Julia wrapper for [algoim](https://github.com/algoim/algoim), providing algorithms for implicitly defined geometry, level set methods, and Voronoi implicit interface methods

**If you use this package, please cite algoim, according to the guidelines on the [algoim Github page](https://algoim.github.io).**

This code wraps algoim's 

 - high-order quadrature algorithms for domains implicitly-defined by multivariate polynomials and 
 - high-order accurate algorithms for computing closest points on implicitly-defined surfaces.

When using the **quadrature algorithms**, please cite:

> [R. I. Saye, High-order quadrature on multi-component domains implicitly defined by multivariate polynomials, Journal of Computational Physics, 448, 110720 (2022).](https://doi.org/10.1016/j.jcp.2021.110720)

When using the **closest-point algorithms**, please cite:

> [R. I. Saye, High-order methods for computing distances to implicitly defined surfaces, Communications in Applied Mathematics and Computational Science, 9(1), 107-141 (2014).](http://dx.doi.org/10.2140/camcos.2014.9.107)
