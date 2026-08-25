# Loading OGIP data

An OGIP observation may be split across a source PHA file, an RMF or RSP
response, an ARF ancillary response, and an optional background PHA file.

## Read one product

Each product can be read separately:

```julia
using XraySpectra

spec = read_pha("source.pha")
resp = read_rmf("source.rmf")
arf = read_ancillary_response("source.arf")
background = read_background("background.pha")
```

`read_pha` returns a SpectrumBase spectrum whose spectral axis contains detector
channel numbers. `read_background` reads the background through the same PHA
reader. `read_rmf` returns a [`ResponseMatrix`](@ref), and
`read_ancillary_response` returns an [`AncillaryResponse`](@ref).

All readers use `Float64` by default. Pass `T = Float32` to request another
floating-point type:

```julia
spec = read_pha("source.pha"; T = Float32)
```

## Read a complete dataset

[`read_dataset`](@ref) starts with a PHA file and follows the OGIP `RESPFILE`,
`ANCRFILE`, and `BACKFILE` headers:

```julia
data = read_dataset(
    "source.pha";
    read_response = true,
    read_ancillary = true,
    read_background = false,
)
```

The result is a named tuple:

```julia
data.spectrum
data.response
data.ancillary
data.background
data.paths
```

Response and ancillary loading default to `true`. Background loading defaults
to `false`, because a background file is optional in many workflows. When a
requested companion file is missing, the loader raises an error rather than
silently returning incomplete data.

Set a flag to `false` when that product is absent or not needed:

```julia
pha_only = read_dataset(
    "source.pha";
    read_response = false,
    read_ancillary = false,
    read_background = false,
)
```

For an RSP file that already includes the effective area, do not request a
separate ARF:

```julia
data = read_dataset("source.pha"; read_ancillary = false)
```

## Inspect linked paths

Use [`read_paths_from_spectrum`](@ref) to resolve the companion paths without
loading their arrays:

```julia
paths = read_paths_from_spectrum("source.pha")

paths.response
paths.ancillary
paths.background
```

Relative filenames are resolved from the PHA directory. The OGIP values
`NONE`, `none`, and `%match%` are handled by the same path lookup.

## PHA metadata

The loaded PHA keeps common OGIP information in its SpectrumBase metadata:

```julia
spec.units
spec.errors
spec.error_statistics
spec.poisson_errors
spec.exposure_time
spec.quality
spec.grouping
spec.telescope
spec.instrument
```

`spec.units` is currently `:counts` or `:rate`; XraySpectra does not attach
Unitful units or perform unit conversion during OGIP loading.
