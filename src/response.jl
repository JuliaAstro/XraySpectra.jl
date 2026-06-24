mutable struct ResponseMatrix{T}
    matrix::SparseMatrixCSC{T,Int}
    channels::Vector{Int}
    channel_bins::Matrix{T}
    bins::Matrix{T}
end

struct AncillaryResponse{T}
    bins::Matrix{T}
    effective_area::Vector{T}
end

"""
    response_bins(response::ResponseMatrix)

Return the low/high energy bins for the input domain of the response matrix.
"""
response_bins(resp::ResponseMatrix) = resp.bins

response_bins_low(resp::ResponseMatrix) = @view resp.bins[:, 1]
response_bins_high(resp::ResponseMatrix) = @view resp.bins[:, 2]

"""
    channel_bins(response::ResponseMatrix)

Return the low/high channel-energy bins for the output domain of the response matrix.
"""
channel_bins(resp::ResponseMatrix) = resp.channel_bins

channel_bins_low(resp::ResponseMatrix) = @view resp.channel_bins[:, 1]
channel_bins_high(resp::ResponseMatrix) = @view resp.channel_bins[:, 2]

ancillary_bins(arf::AncillaryResponse) = arf.bins
ancillary_bins_low(arf::AncillaryResponse) = @view arf.bins[:, 1]
ancillary_bins_high(arf::AncillaryResponse) = @view arf.bins[:, 2]
effective_area(arf::AncillaryResponse) = arf.effective_area

"""
    response_energy(response::ResponseMatrix)

Get the contiguously binned energy corresponding to the *input domain* of the
response matrix. This is equivalent to domain on which a theoretical model would be evaluated on.
"""
response_energy(resp::ResponseMatrix) = [response_bins_low(resp); response_bins_high(resp)[end]]

"""
    folded_energy(response::ResponseMatrix)

Get the contiguously binned energy corresponding to the *output (folded) domain*
of the response matrix. That is, the channel energies as used by the spectrum.
"""
folded_energy(resp::ResponseMatrix) = [channel_bins_low(resp); channel_bins_high(resp)[end]]

function energy_binned_spectrum(spec::SpectrumBase.AbstractSpectrum, resp::ResponseMatrix)
    channels = Int.(SpectrumBase.spectral_axis(spec))
    @assert issorted(resp.channels)
    bins = Matrix{eltype(resp.channel_bins)}(undef, length(channels), 2)
    ebounds = channel_bins(resp)

    for (i, channel) in pairs(channels)
        index = searchsortedfirst(resp.channels, channel)
        if index > length(resp.channels) || resp.channels[index] != channel
            throw(ArgumentError("Channel $channel is not present in the response EBOUNDS."))
        end
        bins[i, :] .= ebounds[index, :]
    end

    metadata = copy(SpectrumBase.meta(spec))
    metadata[:channels] = channels
    SpectrumBase.Spectrum(bins, SpectrumBase.flux_axis(spec), metadata)
end
