module XraySpectra

import FITSFiles
import SpectrumBase
import Distributions

using SpecialFunctions: gamma_inc_inv
using SparseArrays

export PHA, RMF, ARF, read_pha, read_rmf
export read_ancillary_response, read_background, read_dataset
export read_paths_from_spectrum
export ResponseMatrix, response_bins, channel_bins
export response_bins_low, response_bins_high, channel_bins_low, channel_bins_high
export response_energy, folded_energy, energy_binned_spectrum
export AncillaryResponse, ancillary_bins, ancillary_bins_low, ancillary_bins_high, effective_area

include("response.jl")
include("ogip.jl")

end # module XraySpectra
