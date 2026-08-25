# Channel rebinning

Channel rebinning combines adjacent detector channels. It is explicit: loading
a dataset never rebins it automatically.

## Rebin by a factor

The simplest interface combines every `factor` adjacent channels:

```julia
rebinned = rebin_channels(data; factor = 4)
```

This rebins `data.spectrum`, `data.response`, and a loaded `data.background` in
the same way. The ARF is unchanged because it belongs to the response input
energy grid rather than the detector-channel grid.

The last group may contain fewer channels when the original length is not
divisible by the factor.

## Use an explicit grouping

[`channel_grouping`](@ref) creates an OGIP-style grouping vector. A value of `1`
starts a new output bin, while `0` or `-1` continues the current bin:

```julia
grouping = channel_grouping(length(data.spectrum), 4)
rebinned = rebin_channels(data, grouping)
```

For example:

```julia
grouping = [1, 0, 0, 1, -1]
```

creates groups covering channels `1:3` and `4:5`.

## Choose which products to rebin

The dataset method accepts flags for companion products:

```julia
rebinned = rebin_channels(
    data;
    factor = 4,
    rebin_response = false,
    rebin_background = false,
)
```

The source spectrum is always rebinned. A response or background is rebinned
only when its flag is `true` and the object is present.

## What is combined

- PHA counts or rates are summed.
- A rebinned response row is the sum of its original detector-channel rows.
- Channel-energy bounds use the first low edge and last high edge in the group.
- If any original quality value is nonzero, the grouped channel is marked bad.
- Explicit numeric errors are combined in quadrature.
- Poisson errors are recalculated from grouped counts. For rates, the exposure
  time is used to recover counts before calculating the grouped error.

Channel rebinning does not change the response input-energy bins or ARF bins.
