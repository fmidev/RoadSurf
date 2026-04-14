# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

RoadSurf is a Fortran library for predicting road surface conditions (temperature, water/ice/snow/deposit storage). Developed by the Finnish Meteorological Institute (FMI). The Fortran library is called from C++ via ISO_C_BINDING. Physics documented in: Karsisto 2024, Geosci. Model Dev., 17, 4837-4853.

## Build commands

```bash
make                    # Build libroadsurf.so (release, gfortran with -O2 -Ofast)
make debug              # Build with -O0 -Wall
make clean              # Remove libroadsurf.so and obj/
make install            # Install to PREFIX=/usr (lib64 + include/roadsurf/)
make rpm                # Build RPM package
make examples           # Build both examples (requires libroadsurf.so built first)
make example1           # Build example1 only (standalone, needs jsoncpp)
make example2           # Build example2 only (needs SmartMet libraries installed)
```

Sanitizers: `make ASAN=yes`, `make TSAN=yes`, `make BSAN=yes` (address, thread, bounds).

There is no test suite. The examples in `examples/` serve as integration tests.

## Source layout

```
src/                        # Fortran 90 library source
  RoadSurfVariables.f90     # Root module: includes all type definitions from .inc files
  RoadSurf.f90              # Main module: declares public interface (all subroutines)
  *.f90                     # Submodules implementing each subroutine
  *.f90.inc                 # Type definitions (derived types for inputs, outputs, settings, physics)
  Constants.h               # Preprocessor constants (precipitation/snow types)
examples/
  example1/                 # Standalone C++ wrapper (JSON I/O, jsoncpp, no SmartMet deps)
  example2/                 # SmartMet-integrated C++ wrapper (QueryData, PostgreSQL, parallel)
```

Both examples produce a binary called `roadrunner`.

## Architecture

The library is a single shared object (`libroadsurf.so`) with a Fortran module hierarchy:

- **RoadSurfVariables** (root) - All derived types (`InputArrays`, `OutputArrays`, `ModelSettings`, `PhysicalParameters`, etc.) defined via `.inc` files
- **RoadSurf** (main interface) - Declares all public subroutines as module interfaces
- **Submodules** - Each `.f90` file (except the two above) is a submodule implementing one or more subroutines

### Simulation flow (called from C++)

1. `ConnectFortran2Carrays()` - Maps C pointers to Fortran arrays via `C_F_POINTER`
2. `Initialization()` - Sets up model parameters, ground layers, initial temperatures
3. Per-timestep loop:
   - `CouplingOperations1()` - Handle observation coupling (radiation coefficient adjustment)
   - `RelaxationOperations()` - Smooth transition from observations to forecast
   - `SetCurrentValues()` - Load atmospheric data for current timestep
   - `CheckValues()` - Validate input data
   - `ModRadiationBySurroundings()` - Adjust radiation using sky view factor / horizon angles
   - `CalcAlbedo()` - Surface albedo from snow/ice state
   - `BalanceModelOneStep()` - Core heat balance computation
   - `WearFactors()` + `RoadCond()` - Traffic wear and road condition determination
   - `PrecipitationToStorage()` - Phase determination and storage update
   - `CheckEndCoupling()` - End-of-coupling-period checks
   - `SaveOutput()` - Write results to output arrays
4. C++ side reads output arrays (surface temperature, snow, water, ice, deposit, ice2)

### Module dependency chain

`RoadSurfVariables` -> `RoadSurf` -> all submodules (BalanceModel, BoundaryLayer, Cond, Coupling, etc.)

This is explicitly declared in the Makefile dependency rules (lines 76-87).

## Language and compiler notes

- **Fortran 90/95** with Fortran 2008 submodules, compiled with `gfortran`
- Uses `-cpp` flag for C preprocessor (`#include`, `#define` in Constants.h)
- Aggressive math optimizations enabled by default (unsafe-math, no-errno, reciprocal-math, etc.) - NaN handling is intentionally preserved (`-fno-finite-math-only` is NOT used)
- Module files (`.mod`) are compiled to `obj/` directory (`-Jobj`)
- C interoperability through `ISO_C_BINDING` and `BIND(C)` (see ConnectFortran2Carrays.f90)

## RPM packaging

- Library package: `roadsurf` (installs `libroadsurf.so`)
- Devel package: `roadsurf-devel` (installs `.mod` files and `Constants.h` to `/usr/include/roadsurf/`)
- Current version tracked in `roadsurf.spec`
