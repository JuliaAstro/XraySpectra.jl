# Folding

Response folding predicts what the detector would record for a model evaluated
on the response input-energy grid.

If `R` is an `m × n` response matrix and `f` is a model vector of length `n`,
folding computes

```math
y = Rf,
```

where `y` has length `m`, one value per detector channel.

## Fold through an RMF and ARF

```julia
data = read_dataset("source.pha"; read_background = false)

model_flux = ones(size(data.response.matrix, 2))
predicted = fold(
    data.response,
    model_flux;
    ancillary = data.ancillary,
)
```

Supplying `ancillary` first scales every input-energy column of the RMF by the
matching effective-area value. [`combine`](@ref) exposes that intermediate
matrix directly:

```julia
combined_matrix = combine(data.response, data.ancillary)
size(combined_matrix) == size(data.response.matrix)
```

For a full RSP response, omit the ancillary argument:

```julia
predicted = fold(rsp, model_flux)
```

## Preallocated forms

The functions ending in `!` write into an existing output array:

```julia
combined = copy(data.response.matrix)
combine!(combined, data.response, data.ancillary)

predicted = zeros(size(data.response.matrix, 1))
fold!(
    predicted,
    data.response,
    model_flux;
    ancillary = data.ancillary,
)
```

The mutating forms can reuse storage in repeated calculations. The caller must
provide an output with the expected shape.

## Current physical assumptions

`fold` and `combine` perform numerical matrix operations. They validate vector
lengths, matrix shapes, and RMF/ARF bin compatibility, but they do not perform:

- unit conversion,
- exposure-time scaling,
- model integration over energy bins, or
- automatic energy-grid rebinning.

The model vector must already represent the quantity expected by the response
and must already be evaluated on `response_bins(data.response)`.
