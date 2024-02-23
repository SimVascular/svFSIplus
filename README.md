# svFSIplus 

<div align="center">

[![Build Status](https://github.com/SimVascular/svFSIplus/actions/workflows/test.yml/badge.svg)](https://github.com/SimVascular/svFSIplus/actions)
[![codecov](https://codecov.io/github/SimVascular/svFSIplus/graph/badge.svg?token=I848DNIHSP)](https://codecov.io/github/SimVascular/svFSIplus)
![Latest Release](https://img.shields.io/github/v/release/SimVascular/svFSIplus?label=latest)
![Platform](https://img.shields.io/badge/platform-macOS%20|%20linux-blue)
[![DOI](https://img.shields.io/badge/DOI-10.21105%2Fjoss.04118-green)](https://doi.org/10.21105/joss.04118)

</div>

## Introduction

svFSIplus is a C++ implementation of the Fortran [svFSI](https://github.com/SimVascular/svFSI) multi-physics finite element solver designed for computational modeling of the cardiovascular system. It represents the <i>First Stage</i> in the development of a C++ multi-physics finite element solver and is essentially a direct line-by-line translation of the [svFSI](https://github.com/SimVascular/svFSI) Fortran code. This was done to maintain a simple mapping between the code of the two versions and to facilitate debugging by comparing intermediate data (i.e. element stiffness matrices) and simulation results.

The *Second Stage* of the solver development will be an entirely new implementation encapsulating portions of the *First Stage* code into C++ classes forming a clean and simple abstraction that can be enhanced by any developer.

## Getting started

Please see the documentation for [getting started](https://simvascular.github.io/svFSIplus/index.html) and [implementation details](https://simvascular.github.io/svFSIplus/implementation.html). To run examples, have a look at our [testing guide](https://simvascular.github.io/svFSIplus/testing.html).
