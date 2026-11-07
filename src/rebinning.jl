"""
    channel_grouping(n, factor)

Create a grouping vector for `n` channels in adjacent groups of at most
`factor` channels. A value of `1` starts each group.
"""
function channel_grouping(n::Integer, factor::Integer)
    n > 0 || throw(ArgumentError("Number of channels must be positive."))
    factor > 0 || throw(ArgumentError("Channel grouping factor must be positive."))

    grouping = zeros(Int, n)
    grouping[1:factor:end] .= 1
    grouping
end

"""
    channel_grouping(counts::AbstractVector; min_counts::Real, quality=nothing, merge_trailing::Bool=true)

Create an OGIP-style grouping vector where adjacent channels are combined until each group
accumulates at least `min_counts` counts.

If `quality` is provided, channels flagged as bad (`quality != 0`) are excluded from count
accumulation and isolated into separate groups so they do not contaminate good channels.

A value of `1` starts each group. If `merge_trailing` is `true` (default), leftover channels
at the end of the spectrum that do not reach `min_counts` are merged into the preceding group.
"""
function channel_grouping(
    counts::AbstractVector{<:Real};
    min_counts::Real,
    quality::Union{AbstractVector{<:Integer},Nothing} = nothing,
    merge_trailing::Bool = true,
)
    min_counts > 0 || throw(ArgumentError("min_counts must be positive."))
    isempty(counts) && throw(ArgumentError("Channel counts vector cannot be empty."))
    if !isnothing(quality) && length(quality) != length(counts)
        throw(ArgumentError("quality length must match counts length: $(length(quality)) != $(length(counts))."))
    end

    n = length(counts)
    grouping = zeros(Int, n)
    grouping[1] = 1

    accum = zero(eltype(counts))
    last_group_start = 1
    in_bad_group = !isnothing(quality) && quality[1] != 0

    for i in 1:n
        is_bad = !isnothing(quality) && quality[i] != 0

        # Transition between good and bad channels: always start a new group boundary
        if !isnothing(quality) && i > 1
            was_bad = quality[i - 1] != 0
            if is_bad != was_bad
                grouping[i] = 1
                last_group_start = i
                accum = zero(eltype(counts))
                in_bad_group = is_bad
            end
        end

        if !is_bad
            accum += counts[i]
            if accum >= min_counts && i < n
                next_is_bad = !isnothing(quality) && quality[i + 1] != 0
                if !next_is_bad
                    grouping[i + 1] = 1
                    last_group_start = i + 1
                    accum = zero(eltype(counts))
                end
            end
        end
    end

    # If trailing good channels don't reach min_counts, merge them into previous good group
    if merge_trailing && last_group_start > 1 && accum > 0 && accum < min_counts && !in_bad_group
        prev_is_good = isnothing(quality) || quality[last_group_start - 1] == 0
        if prev_is_good
            grouping[last_group_start] = 0
        end
    end

    grouping
end

"""
    channel_grouping(spec::SpectrumBase.AbstractSpectrum; min_counts::Real, ignore_bad::Bool=true, merge_trailing::Bool=true)

Create an OGIP grouping vector for a spectrum based on a minimum number of counts per bin.
If `ignore_bad` is `true` (default) and quality metadata is present, bad channels (`quality != 0`)
are excluded from count accumulation and isolated so good channels are not contaminated.
If the spectrum is stored in count rates, the exposure time from metadata is used to recover counts.
"""
function channel_grouping(
    spec::SpectrumBase.AbstractSpectrum;
    min_counts::Real,
    ignore_bad::Bool = true,
    merge_trailing::Bool = true,
)
    flux = SpectrumBase.flux_axis(spec)
    meta = SpectrumBase.meta(spec)
    units = get(meta, :units, :counts)
    quality = if ignore_bad && haskey(meta, :quality)
        meta[:quality]
    else
        nothing
    end

    counts = if units == :rate
        exposure = get(meta, :exposure_time, one(eltype(flux)))
        flux .* exposure
    else
        flux
    end

    channel_grouping(
        counts;
        min_counts = min_counts,
        quality = quality,
        merge_trailing = merge_trailing,
    )
end

channel_grouping(data::NamedTuple; kwargs...) = channel_grouping(data.spectrum; kwargs...)

"""
    calculate_areanorm(src_spec::SpectrumBase.AbstractSpectrum, bg_spec::SpectrumBase.AbstractSpectrum)

Calculate the background-to-source area normalization factor (`areanorm`):

```math
\\mathrm{areanorm} = \\frac{\\mathrm{AREASCAL}_{\\mathrm{src}} \\times \\mathrm{BACKSCAL}_{\\mathrm{src}} \\times \\mathrm{EXPOSURE}_{\\mathrm{src}}}{\\mathrm{AREASCAL}_{\\mathrm{bg}} \\times \\mathrm{BACKSCAL}_{\\mathrm{bg}} \\times \\mathrm{EXPOSURE}_{\\mathrm{bg}}}
```

Scales background counts to the source extraction region for background-subtracted signal-to-noise ratio calculations.
"""
function calculate_areanorm(src_spec::SpectrumBase.AbstractSpectrum, bg_spec::SpectrumBase.AbstractSpectrum)
    src_meta = SpectrumBase.meta(src_spec)
    bg_meta = SpectrumBase.meta(bg_spec)

    src_back = get(src_meta, :background_scale, 1.0)
    bg_back = get(bg_meta, :background_scale, 1.0)

    src_area = get(src_meta, :area_scale, 1.0)
    bg_area = get(bg_meta, :area_scale, 1.0)

    src_exp = get(src_meta, :exposure_time, 1.0)
    bg_exp = get(bg_meta, :exposure_time, 1.0)

    (src_back .* src_area .* src_exp) ./ (bg_back .* bg_area .* bg_exp)
end

"""
    channel_grouping(source::AbstractVector{<:Real}, background::AbstractVector{<:Real};
                     min_sn::Real, areanorm::Union{Real,AbstractVector{<:Real}}=1.0,
                     quality::Union{AbstractVector{<:Integer},Nothing}=nothing,
                     merge_trailing::Bool=true)

Create an OGIP-style grouping vector where adjacent channels are combined until each group
achieves a minimum signal-to-noise ratio (`min_sn`), following the XMM-Newton SAS `specgroup`
formulation:

```math
S/N = \\frac{\\mathrm{source} - \\mathrm{background} \\times \\mathrm{areanorm}}{\\sqrt{\\mathrm{source} + \\mathrm{background} \\times \\mathrm{areanorm}^2}}
```

Channels with non-zero quality flags (`quality != 0`) are excluded from S/N accumulation and isolated
into separate groups so they do not contaminate adjacent good channels.

If `merge_trailing` is `true` (default), leftover channels at the end that do not reach `min_sn`
are merged into the preceding good group.
"""
function channel_grouping(
    source::AbstractVector{<:Real},
    background::AbstractVector{<:Real};
    min_sn::Real,
    areanorm::Union{Real,AbstractVector{<:Real}} = 1.0,
    quality::Union{AbstractVector{<:Integer},Nothing} = nothing,
    merge_trailing::Bool = true,
)
    min_sn > 0 || throw(ArgumentError("min_sn must be positive."))
    n = length(source)
    n > 0 || throw(ArgumentError("Source counts vector cannot be empty."))
    length(background) == n || throw(ArgumentError("Source and background vectors must have the same length: $(length(source)) != $(length(background))."))

    if areanorm isa AbstractVector
        length(areanorm) == n || throw(ArgumentError("areanorm length must match channels: $(length(areanorm)) != $n."))
        all(>=(0), areanorm) || throw(ArgumentError("areanorm values must be non-negative."))
    else
        areanorm >= 0 || throw(ArgumentError("areanorm must be non-negative."))
    end

    if !isnothing(quality) && length(quality) != n
        throw(ArgumentError("quality length must match channels: $(length(quality)) != $n."))
    end

    grouping = zeros(Int, n)
    grouping[1] = 1

    accum_s = zero(Float64)
    accum_b_scale = zero(Float64)
    accum_b_var = zero(Float64)
    last_group_start = 1
    in_bad_group = !isnothing(quality) && quality[1] != 0

    for i in 1:n
        is_bad = !isnothing(quality) && quality[i] != 0

        # Boundary between good and bad channels
        if !isnothing(quality) && i > 1
            was_bad = quality[i - 1] != 0
            if is_bad != was_bad
                grouping[i] = 1
                last_group_start = i
                accum_s = zero(Float64)
                accum_b_scale = zero(Float64)
                accum_b_var = zero(Float64)
                in_bad_group = is_bad
            end
        end

        if !is_bad
            a_i = areanorm isa AbstractVector ? areanorm[i] : areanorm
            s_i = Float64(source[i])
            b_i = Float64(background[i])

            accum_s += s_i
            accum_b_scale += a_i * b_i
            accum_b_var += (a_i^2) * b_i

            net_signal = accum_s - accum_b_scale
            variance = accum_s + accum_b_var

            current_sn = if variance > 0 && net_signal > 0
                net_signal / sqrt(variance)
            else
                zero(Float64)
            end

            if current_sn >= min_sn && i < n
                next_is_bad = !isnothing(quality) && quality[i + 1] != 0
                if !next_is_bad
                    grouping[i + 1] = 1
                    last_group_start = i + 1
                    accum_s = zero(Float64)
                    accum_b_scale = zero(Float64)
                    accum_b_var = zero(Float64)
                end
            end
        end
    end

    # If trailing group does not reach min_sn, merge it into the previous group
    if merge_trailing && last_group_start > 1 && !in_bad_group
        net_signal = accum_s - accum_b_scale
        variance = accum_s + accum_b_var
        final_sn = (variance > 0 && net_signal > 0) ? (net_signal / sqrt(variance)) : zero(Float64)
        if final_sn < min_sn
            prev_is_good = isnothing(quality) || quality[last_group_start - 1] == 0
            if prev_is_good
                grouping[last_group_start] = 0
            end
        end
    end

    grouping
end

"""
    channel_grouping(src_spec::SpectrumBase.AbstractSpectrum, bg_spec::SpectrumBase.AbstractSpectrum;
                     min_sn::Real, areanorm=nothing, ignore_bad::Bool=true, merge_trailing::Bool=true)

Create an OGIP grouping vector for a source spectrum and companion background spectrum based on
a minimum signal-to-noise ratio (`min_sn`). If `areanorm` is omitted or `nothing`, it is automatically
derived from the `BACKSCAL`, `AREASCAL`, and exposure times of both spectra.
"""
function channel_grouping(
    src_spec::SpectrumBase.AbstractSpectrum,
    bg_spec::SpectrumBase.AbstractSpectrum;
    min_sn::Real,
    areanorm::Union{Real,AbstractVector{<:Real},Nothing} = nothing,
    ignore_bad::Bool = true,
    merge_trailing::Bool = true,
)
    length(src_spec) == length(bg_spec) ||
        throw(ArgumentError("Source and background spectra must have the same length: $(length(src_spec)) != $(length(bg_spec))."))

    src_flux = SpectrumBase.flux_axis(src_spec)
    src_meta = SpectrumBase.meta(src_spec)
    src_units = get(src_meta, :units, :counts)
    src_counts = if src_units == :rate
        exposure = get(src_meta, :exposure_time, one(eltype(src_flux)))
        src_flux .* exposure
    else
        src_flux
    end

    bg_flux = SpectrumBase.flux_axis(bg_spec)
    bg_meta = SpectrumBase.meta(bg_spec)
    bg_units = get(bg_meta, :units, :counts)
    bg_counts = if bg_units == :rate
        exposure = get(bg_meta, :exposure_time, one(eltype(bg_flux)))
        bg_flux .* exposure
    else
        bg_flux
    end

    actual_areanorm = if isnothing(areanorm)
        calculate_areanorm(src_spec, bg_spec)
    else
        areanorm
    end

    quality = if ignore_bad
        src_qual = get(src_meta, :quality, nothing)
        bg_qual = get(bg_meta, :quality, nothing)
        if !isnothing(src_qual) && !isnothing(bg_qual)
            [ (sq != 0 || bq != 0) ? 1 : 0 for (sq, bq) in zip(src_qual, bg_qual) ]
        elseif !isnothing(src_qual)
            src_qual
        elseif !isnothing(bg_qual)
            bg_qual
        else
            nothing
        end
    else
        nothing
    end

    channel_grouping(
        src_counts,
        bg_counts;
        min_sn = min_sn,
        areanorm = actual_areanorm,
        quality = quality,
        merge_trailing = merge_trailing,
    )
end

"""
    rebin_channels(object; factor=nothing, min_counts=nothing, min_sn=nothing, background=nothing, areanorm=nothing, merge_trailing=true)
    rebin_channels(object, grouping)

Combine adjacent detector channels in a spectrum, response matrix, or dataset.
Channels can be grouped by a fixed integer `factor`, a minimum count threshold `min_counts`,
or a minimum signal-to-noise ratio `min_sn` (requiring companion `background`).
"""
function rebin_channels(
    obj;
    factor::Union{Integer,Nothing} = nothing,
    min_counts::Union{Real,Nothing} = nothing,
    min_sn::Union{Real,Nothing} = nothing,
    background = nothing,
    areanorm = nothing,
    merge_trailing::Bool = true,
)
    specified = count(!isnothing, (factor, min_counts, min_sn))
    if specified > 1
        throw(ArgumentError("Cannot specify more than one of `factor`, `min_counts`, and `min_sn`."))
    elseif specified == 0
        throw(ArgumentError("Must specify one of `factor`, `min_counts`, or `min_sn`."))
    end

    if !isnothing(factor)
        rebin_channels(obj, channel_grouping(_channel_count(obj), factor))
    elseif !isnothing(min_counts)
        rebin_channels(obj, channel_grouping(obj; min_counts = min_counts, merge_trailing = merge_trailing))
    else # !isnothing(min_sn)
        if isnothing(background)
            throw(ArgumentError("Rebinning by `min_sn` requires companion `background`."))
        end
        grouping = channel_grouping(obj, background; min_sn = min_sn, areanorm = areanorm, merge_trailing = merge_trailing)
        rebin_channels(obj, grouping)
    end
end

function rebin_channels(
    data::NamedTuple;
    factor::Union{Integer,Nothing} = nothing,
    min_counts::Union{Real,Nothing} = nothing,
    min_sn::Union{Real,Nothing} = nothing,
    areanorm = nothing,
    merge_trailing::Bool = true,
    kwargs...,
)
    specified = count(!isnothing, (factor, min_counts, min_sn))
    if specified > 1
        throw(ArgumentError("Cannot specify more than one of `factor`, `min_counts`, and `min_sn`."))
    elseif specified == 0
        throw(ArgumentError("Must specify one of `factor`, `min_counts`, or `min_sn`."))
    end

    if !isnothing(factor)
        rebin_channels(data, channel_grouping(_channel_count(data), factor); kwargs...)
    elseif !isnothing(min_counts)
        rebin_channels(data, channel_grouping(data.spectrum; min_counts = min_counts, merge_trailing = merge_trailing); kwargs...)
    else # !isnothing(min_sn)
        if isnothing(data.background)
            throw(ArgumentError("Dataset must have a non-nothing `background` spectrum to rebin by `min_sn`."))
        end
        grouping = channel_grouping(data.spectrum, data.background; min_sn = min_sn, areanorm = areanorm, merge_trailing = merge_trailing)
        rebin_channels(data, grouping; kwargs...)
    end
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

function _rebin_bins(bins, spans)
    rebinned = Matrix{eltype(bins)}(undef, length(spans), 2)
    for (new_index, first_index, last_index) in spans
        rebinned[new_index, 1] = bins[first_index, 1]
        rebinned[new_index, 2] = bins[last_index, 2]
    end
    rebinned
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
    _rebin_bins(axis, spans)
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
