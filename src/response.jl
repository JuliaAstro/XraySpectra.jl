mutable struct ResponseMatrix{T}
    matrix::SparseMatrixCSC{T,Int}
    channels::Vector{Int}
    channel_bins::Matrix{T}
    bins::Matrix{T}
    arf_folded::Bool
end

ResponseMatrix(
    matrix::SparseMatrixCSC{T,Int},
    channels::Vector{Int},
    channel_bins::Matrix{T},
    bins::Matrix{T},
) where {T} = ResponseMatrix(matrix, channels, channel_bins, bins, false)

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
arf_folded(resp::ResponseMatrix) = resp.arf_folded

_bins_match(left, right) = isapprox(left, right; rtol = 1e-5, atol = 1e-6)

function _check_ancillary_compatible(resp::ResponseMatrix, arf::AncillaryResponse)
    if arf_folded(resp)
        throw(ArgumentError("Response matrix already has an ancillary response folded in."))
    elseif length(effective_area(arf)) != size(resp.matrix, 2)
        throw(ArgumentError(
            "Ancillary response length must match the response matrix input bins: " *
            "length(effective_area(arf)) != size(resp.matrix, 2) " *
            "($(length(effective_area(arf))) != $(size(resp.matrix, 2)))",
        ))
    elseif size(ancillary_bins(arf)) != size(response_bins(resp))
        throw(ArgumentError(
            "Ancillary response bins must match the response matrix input bins: " *
            "size(ancillary_bins(arf)) != size(response_bins(resp)) " *
            "($(size(ancillary_bins(arf))) != $(size(response_bins(resp))))",
        ))
    elseif !_bins_match(ancillary_bins(arf), response_bins(resp))
        throw(ArgumentError("Ancillary response bins do not match the response matrix input bins."))
    end
    nothing
end

function _check_combine_output_compatible(resp::ResponseMatrix, output)
    if size(output) != size(resp.matrix)
        throw(ArgumentError(
            "Output size must match the response matrix size: " *
            "size(output) != size(resp.matrix) ($(size(output)) != $(size(resp.matrix)))",
        ))
    end
    nothing
end

function _check_fold_compatible(resp::ResponseMatrix, flux)
    if length(flux) != size(resp.matrix, 2)
        throw(ArgumentError(
            "Flux length must match the response matrix input bins: " *
            "length(flux) != size(resp.matrix, 2) ($(length(flux)) != $(size(resp.matrix, 2)))",
        ))
    end
    nothing
end

function _check_fold_output_compatible(resp::ResponseMatrix, output)
    if length(output) != size(resp.matrix, 1)
        throw(ArgumentError(
            "Output length must match the response matrix output bins: " *
            "length(output) != size(resp.matrix, 1) ($(length(output)) != $(size(resp.matrix, 1)))",
        ))
    end
    nothing
end

"""
    combine(response::ResponseMatrix, ancillary::AncillaryResponse)

Fold the ancillary effective area into the response matrix.

The ancillary bins must match the response input bins.
"""
function combine(resp::ResponseMatrix, arf::AncillaryResponse)
    output = copy(resp.matrix)
    combine!(output, resp, arf)
    return output
end

"""
    combine!(output, response::ResponseMatrix, ancillary::AncillaryResponse)

Write `combine(response, ancillary)` into a preallocated output matrix.
"""
function combine!(output, resp::ResponseMatrix, arf::AncillaryResponse)
    _check_ancillary_compatible(resp, arf)
    _check_combine_output_compatible(resp, output)
    output .= effective_area(arf)' .* resp.matrix
end

function combine!(output::SparseMatrixCSC, resp::ResponseMatrix, arf::AncillaryResponse)
    _check_ancillary_compatible(resp, arf)
    _check_combine_output_compatible(resp, output)

    copyto!(output, resp.matrix)
    area = effective_area(arf)

    for col in axes(output, 2)
        for i in nzrange(output, col)
            output.nzval[i] *= area[col]
        end
    end

    return output
end

"""
    fold(response::ResponseMatrix, flux; ancillary=nothing)

Fold a vector through the response matrix.

If `ancillary` is provided, the ancillary response is first combined with the
response matrix. The flux vector must already be on the response input bins.
"""
function fold(resp::ResponseMatrix, flux; ancillary = nothing)
    _check_fold_compatible(resp, flux)
    matrix = isnothing(ancillary) ? resp.matrix : combine(resp, ancillary)
    matrix * flux
end

"""
    fold!(output, response::ResponseMatrix, flux; ancillary=nothing)

Mutating form of [`fold`](@ref), writing the folded vector into `output`.
The output vector must already match the response output bins.
"""
function fold!(output, resp::ResponseMatrix, flux; ancillary = nothing)
    _check_fold_compatible(resp, flux)
    _check_fold_output_compatible(resp, output)
    matrix = isnothing(ancillary) ? resp.matrix : combine(resp, ancillary)
    mul!(output, matrix, flux)
end

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
