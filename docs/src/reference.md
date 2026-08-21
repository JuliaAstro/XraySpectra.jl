# API reference

## Loading

```@docs
read_pha
read_rmf
read_ancillary_response
read_background
read_dataset
read_paths_from_spectrum
```

## Responses

```@docs
ResponseMatrix
AncillaryResponse
RedistributionResponse
FullResponse
response_kind
arf_folded
response_bins
response_bins_low
response_bins_high
channel_bins
channel_bins_low
channel_bins_high
ancillary_bins
ancillary_bins_low
ancillary_bins_high
effective_area
response_energy
folded_energy
energy_binned_spectrum
```

## Folding functions

```@docs
combine
combine!
fold
fold!
```

## Channel rebinning

```@docs
channel_grouping
rebin_channels
```
