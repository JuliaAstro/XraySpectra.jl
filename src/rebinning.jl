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

function rebin_channels(data::NamedTuple; factor::Integer, kwargs...)
    rebin_channels(data, channel_grouping(_channel_count(data), factor); kwargs...)
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

_error_statistic(stats) = stats in (:poisson, :numeric) ? Val(stats) : Val(:unknown)
_error_units(units) = units in (:counts, :rate) ? Val(units) : Val(:unknown)

function _combine_errors(value, errors, stats, units, exposure_time, first_index, last_index)
    _combine_errors(
        _error_statistic(stats),
        value,
        errors,
        _error_units(units),
        exposure_time,
        first_index,
        last_index,
    )
end

function _combine_errors(::Val{:poisson}, value, errors, ::Val{:rate}, exposure_time, first_index, last_index)
    counts = value * exposure_time
    count_error(counts, 1.0) / exposure_time
end

function _combine_errors(::Val{:poisson}, value, errors, ::Val{:counts}, exposure_time, first_index, last_index)
    count_error(value, 1.0)
end

function _combine_errors(::Val{:poisson}, value, errors, ::Val{:unknown}, exposure_time, first_index, last_index)
    throw(ArgumentError("Cannot group Poisson errors without count or rate units."))
end

function _combine_errors(::Val{:numeric}, value, errors, units, exposure_time, first_index, last_index)
    sqrt(sum(abs2, @view errors[first_index:last_index]))
end

function _combine_errors(::Val{:unknown}, value, errors, units, exposure_time, first_index, last_index)
    zero(eltype(errors))
end

function _rebin_spectral_axis(axis::AbstractVector, spans)
    collect(1:length(spans))
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
    rebinned_channels = collect(1:length(spans))

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
        response_kind(resp),
    )
end

function rebin_channels(
    data::NamedTuple,
    grouping;
    rebin_response::Bool = true,
    rebin_background::Bool = true,
)
    spectrum = rebin_channels(data.spectrum, grouping)
    response = if rebin_response && !isnothing(data.response)
        rebin_channels(data.response, grouping)
    else
        data.response
    end
    background = if rebin_background && !isnothing(data.background)
        rebin_channels(data.background, grouping)
    else
        data.background
    end

    merge(data, (;
        spectrum = spectrum,
        response = response,
        background = background,
    ))
end
