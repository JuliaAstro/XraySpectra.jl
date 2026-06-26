"""
    AbstractXrayProduct

The abstract type of all X-ray observation products (PHA, RMF, etc.) that are handled by this package.
"""
abstract type AbstractXrayProduct end

struct PHA <: AbstractXrayProduct end
struct RMF <: AbstractXrayProduct end
struct ARF <: AbstractXrayProduct end

struct MissingHeader <: Exception
    header::String
end
Base.showerror(io::IO, e::MissingHeader) = print(io, "Header: '$(e.header)' is not defined")

struct RMFHeader
    first_channel::Int
    num_channels::Int
end

struct RMFMatrix{V,T,M}
    f_chan::V
    n_chan::V
    bins_low::Vector{T}
    bins_high::Vector{T}
    matrix::M
    header::RMFHeader
end

struct RMFChannels{T}
    channels::Vector{Int}
    bins_low::Vector{T}
    bins_high::Vector{T}
end

function read_pha(path; T::Type = Float64)
    fits = FITSFiles.fits(path)
    parse_hdu(PHA, fits[2]; T = T)
end

function parse_hdu(::Type{PHA}, hdu::FITSFiles.HDU; T::Type = Float64)
    column_names = keys(hdu.data)
    _require_column(column_names, :CHANNEL)

    channels = _read_vector(Int, hdu.data.CHANNEL, :CHANNEL)
    values, units = if :RATE in column_names
        _read_vector(T, hdu.data.RATE, :RATE), :rate
    else
        _require_column(column_names, :COUNTS)
        _read_vector(T, hdu.data.COUNTS, :COUNTS), :counts
    end

    quality = if :QUALITY in column_names
        _read_vector(Int, hdu.data.QUALITY, :QUALITY)
    else
        zeros(Int, length(channels))
    end

    grouping = if :GROUPING in column_names
        _read_vector(Int, hdu.data.GROUPING, :GROUPING)
    else
        ones(Int, length(channels))
    end

    is_poisson = _string_boolean(get(hdu.cards, "POISSERR", false))
    flux, error_statistics, errors = if :STAT_ERR in column_names
        stat_err = _read_vector(T, hdu.data.STAT_ERR, :STAT_ERR)
        values, :numeric, stat_err
    elseif is_poisson
        values, :poisson, @. T(count_error(values, 1.0))
    else
        values, :unknown, zeros(T, length(values))
    end

    meta = Dict{Symbol,Any}(
        :quality => quality,
        :grouping => grouping,
        :units => units,
        :errors => errors,
        :error_statistics => error_statistics,
        :poisson_errors => is_poisson,
        :exposure_time => T(_get_exposure_time(hdu.cards)),
        :background_scale => _get_stable(T, hdu.cards, "BACKSCAL", one(T)),
        :area_scale => _get_stable(T, hdu.cards, "AREASCAL", one(T)),
        :systematic_error => _get_stable(T, hdu.cards, "SYS_ERR", zero(T)),
        :telescope => strip(String(get(hdu.cards, "TELESCOP", ""))),
        :instrument => strip(String(get(hdu.cards, "INSTRUME", ""))),
        :cards => hdu.cards,
    )

    SpectrumBase.Spectrum(channels, flux, meta)
end

function read_rmf(path::AbstractString; T::Type = Float64)
    (header, rmf, channels::RMFChannels{T}) = _read_fits_and_close(path) do fits
        rmf_index = findfirst(fits) do hdu
            extname = get(hdu.cards, "EXTNAME", "")
            occursin("RESP", extname) || occursin("MATRIX", extname)
        end
        if isnothing(rmf_index)
            throw(MissingHeader("MATRIX"))
        end
        rmf_hdu = fits[rmf_index]
        hdr = parse_rmf_header(rmf_hdu)
        _rmf = read_rmf_matrix(rmf_hdu, hdr, T)
        _channels = read_rmf_channels(fits["EBOUNDS"], T)
        (hdr, _rmf, _channels)
    end

    _build_response_matrix(header, rmf, channels, T; arf_folded = _is_rsp(path))
end

_is_rsp(path::AbstractString) = lowercase(splitext(String(path))[2]) == ".rsp"

function read_ancillary_response(path::AbstractString; T::Type = Float64)
    _read_fits_and_close(path) do fits
        parse_hdu(ARF, fits["SPECRESP"]; T = T)
    end
end

function read_background(path::AbstractString; T::Type = Float64)
    read_pha(path; T = T)
end

function read_dataset(
    path::AbstractString;
    read_response::Bool = true,
    read_ancillary::Bool = true,
    read_background::Bool = false,
    T::Type = Float64,
)
    paths = read_paths_from_spectrum(path)
    spectrum = read_pha(path; T = T)

    response = if read_response
        response_path = _required_ogip_file(paths.response, "response")
        read_rmf(response_path; T = T)
    else
        nothing
    end

    ancillary = if read_ancillary
        ancillary_path = _required_ogip_file(paths.ancillary, "ancillary")
        read_ancillary_response(ancillary_path; T = T)
    else
        nothing
    end

    background = if read_background
        background_path = _required_ogip_file(paths.background, "background")
        read_pha(background_path; T = T)
    else
        nothing
    end

    (
        spectrum = spectrum,
        response = response,
        ancillary = ancillary,
        background = background,
        paths = paths,
    )
end

function parse_hdu(
    ::Type{RMF},
    matrix_hdu::FITSFiles.HDU,
    ebounds_hdu::FITSFiles.HDU;
    T::Type = Float64,
)
    header = parse_rmf_header(matrix_hdu)
    rmf = read_rmf_matrix(matrix_hdu, header, T)
    channels = read_rmf_channels(ebounds_hdu, T)
    _build_response_matrix(header, rmf, channels, T)
end

function parse_hdu(::Type{ARF}, hdu::FITSFiles.HDU; T::Type = Float64)
    column_names = keys(hdu.data)
    _require_column(column_names, :ENERG_LO; table = "ARF")
    _require_column(column_names, :ENERG_HI; table = "ARF")
    _require_column(column_names, :SPECRESP; table = "ARF")

    energy_low = _read_vector(T, hdu.data.ENERG_LO, :ENERG_LO; table = "ARF")
    energy_high = _read_vector(T, hdu.data.ENERG_HI, :ENERG_HI; table = "ARF")
    area = _read_vector(T, hdu.data.SPECRESP, :SPECRESP; table = "ARF")

    AncillaryResponse(hcat(energy_low, energy_high), area)
end

function _read_fits_and_close(f, path)
    fits_file = FITSFiles.fits(path)
    f(fits_file)
end

function read_paths_from_spectrum(path::AbstractString)
    header = _read_fits_and_close(path) do fits
        fits[2].cards
    end

    possible_ext = splitext(path)[2]
    response_path = read_filename(header, "RESPFILE", path, ".rmf", ".rsp")
    ancillary_path = read_filename(header, "ANCRFILE", path, possible_ext)
    background_path = read_filename(header, "BACKFILE", path, possible_ext)

    (
        spectrum = path,
        response = response_path,
        ancillary = ancillary_path,
        background = background_path,
    )
end

function read_filename(header, entry, parent, exts...)
    data_directory = dirname(parent)
    parent_name = basename(parent)
    if haskey(header, entry)
        path = strip(String(get(header, entry)))
        if path == "NONE"
            return nothing
        end
        name = find_file(data_directory, path, parent_name, exts)
        if !isnothing(name)
            return name
        end
    end
    nothing
end

function find_file(dir, name, parent, extensions)
    if length(name) == 0
        return nothing
    elseif match(r"%match%", name) !== nothing
        base = splitext(parent)[1]
        for ext in extensions
            testfile = joinpath(dir, base * ext)
            if isfile(testfile)
                return testfile
            end
        end
        @warn "Missing! Could not find file '%match%': tried $extensions"
        return nothing
    elseif match(r"^none\b", name) !== nothing
        return nothing
    end
    joinpath(dir, name)
end

function _required_ogip_file(path, label)
    if isnothing(path)
        throw(ArgumentError("No $label file found in the PHA header."))
    elseif !isfile(path)
        throw(ArgumentError("The $label file '$path' does not exist."))
    end
    path
end

function parse_rmf_header(table::FITSFiles.HDU)
    findex = findfirst(==(:F_CHAN), keys(table.data))
    if isnothing(findex)
        throw(MissingHeader("F_CHAN"))
    end

    tlindex = "TLMIN$findex"
    first_channel = haskey(table.cards, tlindex) ? _parse_any(Int, get(table.cards, tlindex)) : begin
        @warn "No TLMIN key set in RMF header. Assuming channels start at 1."
        1
    end

    num_channels = haskey(table.cards, "DETCHANS") ? _parse_any(Int, get(table.cards, "DETCHANS")) : begin
        @warn "DETCHANS is not set in RMF header. Infering channel count from table length."
        -1
    end
    RMFHeader(first_channel, num_channels)
end

function read_rmf_channels(table::FITSFiles.HDU, T::Type)
    channels = _parse_any.(Int, table.data.CHANNEL)
    energy_low = _parse_any.(T, table.data.E_MIN)
    energy_high = _parse_any.(T, table.data.E_MAX)
    RMFChannels(channels, energy_low, energy_high)
end

function read_rmf_matrix(table::FITSFiles.HDU, header::RMFHeader, T::Type)
    energy_low = convert.(T, table.data.ENERG_LO)
    energy_high = convert.(T, table.data.ENERG_HI)
    f_chan_raw = table.data.F_CHAN
    n_chan_raw = table.data.N_CHAN
    matrix_raw = table.data.MATRIX

    # type stable: convert to common vector of vector format
    f_chan::Vector{Vector{Int}} = _translate_channel_array(f_chan_raw)
    n_chan::Vector{Vector{Int}} = _translate_channel_array(n_chan_raw)

    RMFMatrix(
        f_chan,
        n_chan,
        energy_low,
        energy_high,
        _adapt_matrix_type(T, matrix_raw),
        header,
    )
end

function _build_response_matrix(
    header::RMFHeader,
    rmf::RMFMatrix,
    channels::RMFChannels,
    T::Type;
    arf_folded::Bool = false,
)
    R = build_response_matrix(
        rmf.f_chan,
        rmf.n_chan,
        rmf.matrix,
        header.num_channels,
        header.first_channel,
        T,
    )
    ResponseMatrix(
        R,
        channels.channels,
        hcat(channels.bins_low, channels.bins_high),
        hcat(rmf.bins_low, rmf.bins_high),
        arf_folded,
    )
end

function build_response_matrix(
    f_chan::Vector,
    n_chan::Vector,
    matrix_rows::Vector,
    num_cols::Int,
    first_channel,
    T::Type,
)
    ptrs = Int[1]
    indices = Int[]
    matrix = Float64[]

    prev = first(ptrs)

    for i in eachindex(f_chan)
        row_values = matrix_rows[i]

        row_offset = 0
        for (first_chan, n_chan_group) in zip(f_chan[i], n_chan[i])
            if n_chan_group == 0
                # advance row
                break
            end
            first_index = (first_chan - first_channel) + 1
            # append all of the indices
            for j = 0:(n_chan_group-1)
                push!(indices, j + first_index)
            end
            append!(matrix, row_values[(row_offset+1):(row_offset+n_chan_group)])
            row_offset += n_chan_group
        end

        next = row_offset + prev
        push!(ptrs, next)
        prev = next
    end

    SparseMatrixCSC{T,Int}(num_cols, length(f_chan), ptrs, indices, matrix)
end

function _chan_to_vectors(chan::AbstractMatrix)
    map(eachcol(chan)) do column
        zero_index = findfirst(==(0), column)
        # if no zeroes, return full column
        if isnothing(zero_index)
            column
        else
            # exclude the zero from our selection
            last_nonzero = zero_index > 1 ? zero_index - 1 : zero_index
            column[1:last_nonzero]
        end
    end
end

function _translate_channel_array(channel)
    if channel isa AbstractMatrix
        # transpose the matrix, since FITSFiles reads things in Julia-style
        _chan_to_vectors(transpose(channel))
    elseif eltype(channel) <: AbstractVector
        # reorder for the same reason as above transpose
        channel
    else
        map(i -> [i], channel)
    end
end

function _adapt_matrix_type(T::Type, mat::M) where {M}
    if eltype(M) <: AbstractVector
        map(row -> convert.(T, row), mat)
    elseif M <: AbstractMatrix
        map(row -> convert.(T, row), eachcol(mat))
    else
        # handle simple vector format where each energy bin maps to a single
        # value (e.g., for responses produced using ftflx2xsp)
        map(val -> [convert(T, val)], mat)
    end
end

function _require_column(column_names, column::Symbol; table = "PHA")
    column in column_names && return nothing
    throw(ArgumentError("$table table is missing required column $column"))
end

function _read_vector(::Type{T}, data, column::Symbol; table = "PHA") where {T}
    if data isa AbstractMatrix
        throw(ArgumentError("PHA II is not supported yet"))
    elseif !(data isa AbstractVector)
        throw(ArgumentError("$table column $column must be a vector"))
    end
    convert.(T, data)
end

function _parse_any(::Type{T}, value::AbstractString)::T where {T}
    parse(T, value)
end

function _parse_any(::Type{T}, value)::T where {T}
    convert(T, value)
end

function _string_boolean(@nospecialize(value::V))::Bool where {V}
    if V <: AbstractString
        if value == "F"
            false
        elseif value == "T"
            true
        else
            @warn("Unknown boolean string: $(value)")
            false
        end
    else
        value
    end
end

function _get_exposure_time(header)
    if haskey(header, "EXPOSURE")
        return get(header, "EXPOSURE")
    end
    if haskey(header, "TELAPSE")
        return get(header, "TELAPSE")
    end
    if (haskey(header, "TSTART")) && (haskey(header, "TSTOP"))
        return get(header, "TSTOP") - get(header, "TSTART")
    end
    @warn "Cannot find or infer exposure time."
    0.0
end

function _get_stable(::Type{T}, header, key, default)::T where {T}
    get(header, key, T(default))
end

function count_error(k, σ)
    p = Distributions.cdf(Distributions.Normal(), σ)
    kₑ = gamma_inc_inv(k + 1, p, 1 - p)
    abs(k - kₑ)
end
