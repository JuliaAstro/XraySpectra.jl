function channel_grouping(n::Integer, factor::Integer)
    n > 0 || throw(ArgumentError("Number of channels must be positive."))
    factor > 0 || throw(ArgumentError("Channel grouping factor must be positive."))

    grouping = zeros(Int, n)
    grouping[1:factor:end] .= 1
    grouping
end

function rebin_channels(obj; factor::Integer)
    rebin_channels(obj, channel_grouping(_channel_count(obj), factor))
end

_channel_count(spec::SpectrumBase.AbstractSpectrum) = length(spec)
_channel_count(resp::ResponseMatrix) = size(resp.matrix, 1)
_channel_count(data::NamedTuple) = _channel_count(data.spectrum)

function _group_spans(grouping)
    isempty(grouping) && throw(ArgumentError("Channel grouping cannot be empty."))
    first(grouping) == 1 || throw(ArgumentError("Channel grouping must start with 1."))

    starts = findall(==(1), grouping)
    map(enumerate(starts)) do (new_index, first_index)
        last_index = new_index == length(starts) ? length(grouping) : starts[new_index + 1] - 1
        (new_index, first_index, last_index)
    end
end

function _check_grouping_length(name, n, grouping)
    if length(grouping) != n
        throw(ArgumentError(
            "Channel grouping length must match $name channels: " *
            "length(grouping) != $n ($(length(grouping)) != $n)",
        ))
    end
    nothing
end

function _combine_quality(quality, first_index, last_index)
    any(!=(0), @view quality[first_index:last_index]) ? 1 : 0
end

function _combine_errors(value, errors, stats, units, exposure_time, first_index, last_index)
    if stats == :poisson
        if units == :rate
            counts = value * exposure_time
            count_error(counts, 1.0) / exposure_time
        else
            count_error(value, 1.0)
        end
    elseif stats == :numeric
        sqrt(sum(abs2, @view errors[first_index:last_index]))
    else
        zero(eltype(errors))
    end
end

function _rebin_spectral_axis(axis::AbstractVector, spans)
    [axis[first_index] for (_, first_index, _) in spans]
end

function _rebin_spectral_axis(axis::AbstractMatrix, spans)
    rebinned = Matrix{eltype(axis)}(undef, length(spans), 2)
    for (new_index, first_index, last_index) in spans
        rebinned[new_index, 1] = axis[first_index, 1]
        rebinned[new_index, 2] = axis[last_index, 2]
    end
    rebinned
end

function rebin_channels(spec::SpectrumBase.AbstractSpectrum, grouping)
    _check_grouping_length("spectrum", length(spec), grouping)
    spans = _group_spans(grouping)

    SpectrumBase.flux_axis(spec) isa AbstractVector ||
        throw(ArgumentError("Channel rebinning currently expects a vector flux axis."))

    values = collect(SpectrumBase.flux_axis(spec))
    rebinned_values = Vector{eltype(values)}(undef, length(spans))
    for (new_index, first_index, last_index) in spans
        rebinned_values[new_index] = sum(@view values[first_index:last_index])
    end

    metadata = copy(SpectrumBase.meta(spec))
    metadata[:channel_groups] = [(first_index, last_index) for (_, first_index, last_index) in spans]

    if haskey(metadata, :quality)
        quality = metadata[:quality]
        _check_grouping_length("quality", length(quality), grouping)
        metadata[:quality] = [_combine_quality(quality, first_index, last_index) for (_, first_index, last_index) in spans]
    end

    metadata[:grouping] = ones(Int, length(spans))

    if haskey(metadata, :errors)
        errors = metadata[:errors]
        _check_grouping_length("errors", length(errors), grouping)
        stats = get(metadata, :error_statistics, :unknown)
        units = get(metadata, :units, :unknown)
        exposure_time = get(metadata, :exposure_time, one(eltype(rebinned_values)))
        metadata[:errors] = [
            _combine_errors(
                rebinned_values[new_index],
                errors,
                stats,
                units,
                exposure_time,
                first_index,
                last_index,
            )
            for (new_index, first_index, last_index) in spans
        ]
    end

    axis = _rebin_spectral_axis(SpectrumBase.spectral_axis(spec), spans)
    SpectrumBase.Spectrum(axis, rebinned_values, metadata)
end

function rebin_channels(resp::ResponseMatrix, grouping)
    _check_grouping_length("response", size(resp.matrix, 1), grouping)
    spans = _group_spans(grouping)

    row_map = zeros(Int, size(resp.matrix, 1))
    for (new_index, first_index, last_index) in spans
        row_map[first_index:last_index] .= new_index
    end

    rows = rowvals(resp.matrix)
    values = nonzeros(resp.matrix)
    new_rows = Int[]
    new_cols = Int[]
    new_values = eltype(values)[]

    for col in axes(resp.matrix, 2)
        for i in nzrange(resp.matrix, col)
            push!(new_rows, row_map[rows[i]])
            push!(new_cols, col)
            push!(new_values, values[i])
        end
    end

    rebinned_matrix = sparse(new_rows, new_cols, new_values, length(spans), size(resp.matrix, 2))
    rebinned_channels = [resp.channels[first_index] for (_, first_index, _) in spans]

    rebinned_channel_bins = Matrix{eltype(resp.channel_bins)}(undef, length(spans), 2)
    for (new_index, first_index, last_index) in spans
        rebinned_channel_bins[new_index, 1] = resp.channel_bins[first_index, 1]
        rebinned_channel_bins[new_index, 2] = resp.channel_bins[last_index, 2]
    end

    ResponseMatrix(
        rebinned_matrix,
        rebinned_channels,
        rebinned_channel_bins,
        copy(resp.bins),
        arf_folded(resp),
    )
end

function rebin_channels(data::NamedTuple, grouping)
    spectrum = rebin_channels(data.spectrum, grouping)
    response = isnothing(data.response) ? nothing : rebin_channels(data.response, grouping)
    background = isnothing(data.background) ? nothing : rebin_channels(data.background, grouping)

    merge(data, (;
        spectrum = spectrum,
        response = response,
        background = background,
    ))
end
