# Responses

An X-ray response connects two different grids:

- **Input/model-energy bins:** where a theoretical photon model is evaluated.
- **Output/channel bins:** where the detector records events.

For a response matrix with `m` detector channels and `n` model-energy bins,
`size(resp.matrix) == (m, n)`.

```julia
resp = read_rmf("source.rmf")

size(resp.matrix)
response_bins(resp)[1:5, :]
channel_bins(resp)[1:5, :]
```

Each row of `response_bins(resp)` contains the low and high edge of one input
energy bin. Each row of `channel_bins(resp)` contains the energy bounds for one
detector channel from EBOUNDS.

The edge-vector forms are also available:

```julia
response_energy(resp)
folded_energy(resp)
```

Each contains one more element than the corresponding number of bins.

## RMF and RSP responses

A redistribution-only RMF is represented by [`RedistributionResponse`](@ref).
It can be combined with a separate ARF. An RSP, or a response whose OGIP
`HDUCLAS3` is `FULL`, is represented by [`FullResponse`](@ref) because its
effective area is already included.

```julia
response_kind(resp)
arf_folded(resp)
```

XraySpectra rejects attempts to combine another ARF with a full response. This
prevents applying the effective area twice.

## Ancillary responses

An ARF stores one effective-area value for every input energy bin:

```julia
arf = read_ancillary_response("source.arf")

ancillary_bins(arf)[1:5, :]
effective_area(arf)[1:5]
```

The ARF and RMF input grids must match before they can be combined. XraySpectra
currently reports mismatched grids as an error; it does not automatically
rebin them.

## Attach channel-energy bins to a PHA

PHA files normally store channel numbers rather than physical energy edges.
After loading the response, [`energy_binned_spectrum`](@ref) can produce a
binned SpectrumBase spectrum using the RMF EBOUNDS:

```julia
energy_spec = energy_binned_spectrum(data.spectrum, data.response)
```

The original detector channel numbers remain in `energy_spec.channels`.
