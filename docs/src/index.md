# XraySpectra.jl

XraySpectra loads OGIP X-ray spectra and instrument responses into Julia. PHA
data is stored using types from
[SpectrumBase.jl](https://github.com/JuliaAstro/SpectrumBase.jl), while
XraySpectra provides the X-ray-specific response and ancillary response types.

!!! warning
    XraySpectra is a work in progress and its API may still change.

## Installation

The package and its current dependencies are not all registered yet. Install
them from GitHub:

```julia
import Pkg
Pkg.add(url="https://github.com/JuliaAstro/SpectrumBase.jl")
Pkg.add(url="https://github.com/JuliaAstro/FITSFiles.jl")
Pkg.add(url="https://github.com/JuliaAstro/XraySpectra.jl")
```

## Quick start

```julia
using XraySpectra
using SpectrumBase: flux_axis, spectral_axis

data = read_dataset(
    "source.pha";
    read_background = false,
)

length(data.spectrum)
size(data.response.matrix)
length(effective_area(data.ancillary))
```

The PHA spectrum is stored in detector channels. Its RMF describes how photons
from model-energy bins are redistributed into those detector channels, while
the ARF gives the effective area on the model-energy grid.

See [Loading OGIP data](@ref) for partial datasets and linked files, then the
[folding guide](folding.md) for applying responses to a model vector.
