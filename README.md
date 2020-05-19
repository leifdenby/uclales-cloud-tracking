# Cloud-tracking for UCLALES

[![Build
Status](https://travis-ci.com/leifdenby/uclales-cloud-tracking.svg?branch=master)](https://travis-ci.com/leifdenby/uclales-cloud-tracking)

This repo contains tools for tracking convective clouds in simulations run with
[UCLALES](https://github.com/uclales/uclales/). It is a separate fork and
rewrite of large part of the original cloud tracking code developed by Thijs
Heus.

## Installation

To compile `netcdf4` is required (e.g. load with `module load netcdf4` or
install dependencies with conda `python .conda/environment.yml.meta`)

Build with cmake, e.g.

```bash
> mkdir build
> cd build/
> cmake ..
> make
```

To make the documentation `Doxygen` will be required, build with

```bash
> mkdir build
> cd build/
> cmake -DBUILD_DOC=ON ..
> make
```

## Testing

Tests a written in python including util to generate synthetic input files. All
tests are in `tests/`. Run the tests with

```bash
> mkdir build
> cd build/
> cmake ..
> make check
