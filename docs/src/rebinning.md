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

## Rebin by minimum counts

To group channels until each bin accumulates at least a target number of counts (useful for ensuring Poisson-to-Gaussian approximation for $\chi^2$ fitting):

```julia
rebinned = rebin_channels(data; min_counts = 20)
```

By default, any leftover channels at the end of the spectrum that do not reach `min_counts` are merged into the preceding group (`merge_trailing = true`). Pass `merge_trailing = false` to keep them as a distinct trailing group.

When channel quality flags are present in the spectrum metadata, channels marked as bad (`quality != 0`) are excluded from count accumulation and isolated into separate bins (`ignore_bad = true`), preventing bad channels from contaminating good channels.

## Rebin by minimum signal-to-noise ratio (S/N)

Channels can also be grouped until each bin achieves a target signal-to-noise ratio ($S/N$), following the XMM-Newton SAS `specgroup` definition:

```math
S/N = \frac{\mathrm{source} - \mathrm{background} \times \mathrm{areanorm}}{\sqrt{\mathrm{source} + \mathrm{background} \times \mathrm{areanorm}^2}}
```

where $\mathrm{source}$ and $\mathrm{background}$ are the accumulated counts in the group, and $\mathrm{areanorm}$ is the scaling factor accounting for the ratios of `BACKSCAL`, `AREASCAL`, and exposure times between the source and background extractions:

```julia
rebinned = rebin_channels(data; min_sn = 6.0)
```

Alternatively, `channel_grouping` can be called directly with source and background spectra:

```julia
grouping = channel_grouping(src_spec, bg_spec; min_sn = 6.0)
rebinned = rebin_channels(src_spec, grouping)
```

If an explicit area normalization factor is not provided, [`calculate_areanorm`](@ref) automatically derives it from the metadata of both spectra. As with count rebinning, bad channels (`quality != 0`) are isolated, and leftover trailing channels below `min_sn` are merged into the preceding group by default (`merge_trailing = true`).

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
