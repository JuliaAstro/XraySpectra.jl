using Test
using Artifacts
using FITSFiles
using SparseArrays
using SpectrumBase
using XraySpectra

const TEST_DATA = artifact"test_data"
const NUSTAR_PHA = joinpath(TEST_DATA, "nustar", "nu60001047002A01_sr_grp_simple.pha")
const NUSTAR_RMF = joinpath(TEST_DATA, "nustar", "nu60001047002A01_sr.rmf")
const NUSTAR_ARF = joinpath(TEST_DATA, "nustar", "nu60001047002A01_sr.arf")
const CHANDRA_PHA2 = joinpath(TEST_DATA, "chandra", "acisf13839N002_pha2.fits")
const EXOSAT_RSP = joinpath(TEST_DATA, "rsp", "s54405.rsp")
const GINGA_RSP = joinpath(TEST_DATA, "rsp", "ginga_lac.rsp")

@testset "PHA I reader" begin
    spec = read_pha(NUSTAR_PHA)

    @test spec isa SingleSpectrum
    @test length(spec) == 4096
    @test first(spectral_axis(spec)) == 0
    @test last(spectral_axis(spec)) == 4095
    @test length(flux_axis(spec)) == length(spectral_axis(spec))
    @test first(flux_axis(spec)) == 6.0

    @test spec.units == :counts
    @test spec.error_statistics == :poisson
    @test spec.telescope == "NuSTAR"
    @test spec.instrument == "FPMA"
    @test length(spec.quality) == length(spec)
    @test length(spec.grouping) == length(spec)
    @test length(spec.errors) == length(spec)
end

@testset "explicit STAT_ERR stays in metadata" begin
    fake_hdu = FITSFiles.HDU(
        (;
            CHANNEL = [1, 2, 3],
            COUNTS = [10.0, 11.0, 12.0],
            STAT_ERR = [1.0, 1.1, 1.2],
        ),
        FITSFiles.Card[FITSFiles.Card("EXPOSURE", 1.0)],
    )

    spec = XraySpectra.parse_hdu(PHA, fake_hdu)

    @test spec.error_statistics == :numeric
    @test eltype(flux_axis(spec)) == Float64
    @test flux_axis(spec) == fake_hdu.data.COUNTS
    @test spec.errors == fake_hdu.data.STAT_ERR
end

@testset "PHA II is out of scope for now" begin
    @test_throws ArgumentError read_pha(CHANDRA_PHA2)
end

@testset "RMF reader" begin
    response = read_rmf(NUSTAR_RMF)

    @test response isa ResponseMatrix
    @test response_kind(response) == RedistributionResponse
    @test !arf_folded(response)
    @test length(response.channels) == 4096
    @test first(response.channels) == 0
    @test last(response.channels) == 4095
    @test size(response.channel_bins) == (4096, 2)
    @test size(response.bins, 2) == 2
    @test channel_bins_low(response)[1] ≈ 1.60
    @test channel_bins_high(response)[1] ≈ 1.64
    @test size(response.matrix) == (4096, 4096)
    @test nnz(response.matrix) > 0
    @test response.matrix[1, 1] ≈ 3.7741076084785163e-5
    @test response.matrix[40, 116] ≈ 0.00010734131501521915
    @test response.matrix[4096, 4096] ≈ 0.004404161591082811
    @test response_bins(response)[1, 1] ≈ response_bins_low(response)[1]
    @test response_bins(response)[1, 2] ≈ response_bins_high(response)[1]
    @test channel_bins(response)[1, 1] ≈ channel_bins_low(response)[1]
    @test channel_bins(response)[1, 2] ≈ channel_bins_high(response)[1]
    @test response_energy(response) == [response_bins_low(response); response_bins_high(response)[end]]
    @test folded_energy(response) == [channel_bins_low(response); channel_bins_high(response)[end]]
end

@testset "RSP reader" begin
    exosat = read_rmf(EXOSAT_RSP)
    ginga = read_rmf(GINGA_RSP)

    @test exosat isa ResponseMatrix
    @test response_kind(exosat) == FullResponse
    @test arf_folded(exosat)
    @test size(exosat.matrix) == (128, 128)
    @test length(exosat.channels) == 128
    @test size(response_bins(exosat)) == (128, 2)
    @test size(channel_bins(exosat)) == (128, 2)
    @test nnz(exosat.matrix) > 0

    @test ginga isa ResponseMatrix
    @test response_kind(ginga) == FullResponse
    @test arf_folded(ginga)
    @test size(ginga.matrix) == (48, 700)
    @test length(ginga.channels) == 48
    @test size(response_bins(ginga)) == (700, 2)
    @test size(channel_bins(ginga)) == (48, 2)
    @test nnz(ginga.matrix) > 0

    rsp_arf = AncillaryResponse(response_bins(exosat), ones(size(exosat.matrix, 2)))
    @test_throws ArgumentError combine(exosat, rsp_arf)
end

@testset "ARF reader" begin
    fake_hdu = FITSFiles.HDU(
        (;
            ENERG_LO = [1.0, 2.0, 3.0],
            ENERG_HI = [2.0, 3.0, 4.0],
            SPECRESP = [10.0, 20.0, 30.0],
        ),
        FITSFiles.Card[FITSFiles.Card("EXTNAME", "SPECRESP")],
    )

    arf = XraySpectra.parse_hdu(ARF, fake_hdu)

    @test arf isa AncillaryResponse
    @test ancillary_bins(arf) == [1.0 2.0; 2.0 3.0; 3.0 4.0]
    @test ancillary_bins_low(arf) == [1.0, 2.0, 3.0]
    @test ancillary_bins_high(arf) == [2.0, 3.0, 4.0]
    @test effective_area(arf) == [10.0, 20.0, 30.0]
end

@testset "response folding" begin
    response = ResponseMatrix(
        sparse([1.0 2.0 0.0; 0.0 3.0 4.0]),
        [1, 2],
        [0.1 0.2; 0.2 0.3],
        [1.0 2.0; 2.0 3.0; 3.0 4.0],
    )
    arf = AncillaryResponse(
        [1.0 2.0; 2.0 3.0; 3.0 4.0],
        [10.0, 20.0, 30.0],
    )
    flux = [1.0, 2.0, 3.0]

    expected_matrix = [10.0 40.0 0.0; 0.0 60.0 120.0]

    @test combine(response, arf) isa SparseMatrixCSC
    @test Matrix(combine(response, arf)) == expected_matrix

    combined = zeros(2, 3)
    @test combine!(combined, response, arf) == expected_matrix
    @test combined == expected_matrix

    sparse_combined = similar(response.matrix)
    @test combine!(sparse_combined, response, arf) === sparse_combined
    @test sparse_combined isa SparseMatrixCSC
    @test Matrix(sparse_combined) == expected_matrix

    @test_throws ArgumentError combine!(zeros(1, 3), response, arf)

    @test fold(response, flux) == response.matrix * flux
    @test fold(response, flux; ancillary = arf) == expected_matrix * flux

    output = zeros(2)
    @test fold!(output, response, flux) == response.matrix * flux
    @test output == response.matrix * flux

    folded = zeros(2)
    @test fold!(folded, response, flux; ancillary = arf) == expected_matrix * flux
    @test folded == expected_matrix * flux

    @test_throws ArgumentError fold(response, [1.0, 2.0])
    @test_throws ArgumentError fold!(zeros(1), response, flux)

    short_arf = AncillaryResponse([1.0 2.0; 2.0 3.0], [10.0, 20.0])
    shifted_arf = AncillaryResponse([1.0 2.0; 2.0 3.0; 3.1 4.1], [10.0, 20.0, 30.0])

    @test_throws ArgumentError combine(response, short_arf)
    @test_throws ArgumentError combine(response, shifted_arf)

    try
        combine(response, short_arf)
    catch err
        @test err isa ArgumentError
        @test occursin("length(effective_area(arf)) != size(resp.matrix, 2)", err.msg)
        @test occursin("(2 != 3)", err.msg)
    end
end

@testset "channel rebinning" begin
    grouping = [1, 0, 1, 0]

    @test channel_grouping(5, 2) == [1, 0, 1, 0, 1]
    @test channel_grouping(4, 10) == [1, 0, 0, 0]
    @test_throws ArgumentError channel_grouping(0, 2)
    @test_throws ArgumentError channel_grouping(4, 0)

    spec = SpectrumBase.Spectrum(
        [0, 1, 2, 3],
        [10.0, 20.0, 30.0, 40.0],
        Dict{Symbol,Any}(
            :quality => [0, 1, 0, 0],
            :grouping => grouping,
            :errors => [1.0, 2.0, 3.0, 4.0],
            :error_statistics => :numeric,
            :units => :counts,
        ),
    )
    rebinned_spec = rebin_channels(spec, grouping)

    @test spectral_axis(rebinned_spec) == [1, 2]
    @test flux_axis(rebinned_spec) == [30.0, 70.0]
    @test rebinned_spec.quality == [1, 0]
    @test rebinned_spec.grouping == [1, 1]
    @test rebinned_spec.errors ≈ [sqrt(5.0), 5.0]
    @test rebinned_spec.channel_groups == [(1, 2), (3, 4)]
    @test flux_axis(rebin_channels(spec, [1, -1, 1, -1])) == [30.0, 70.0]

    poisson_spec = SpectrumBase.Spectrum(
        [0, 1, 2, 3],
        [1.0, 2.0, 3.0, 4.0],
        Dict{Symbol,Any}(
            :errors => zeros(4),
            :error_statistics => :poisson,
            :units => :counts,
        ),
    )
    rebinned_poisson = rebin_channels(poisson_spec; factor = 2)

    @test flux_axis(rebinned_poisson) == [3.0, 7.0]
    @test rebinned_poisson.errors ≈ XraySpectra.count_error.(flux_axis(rebinned_poisson), 1.0)

    poisson_rate_spec = SpectrumBase.Spectrum(
        [0, 1],
        [0.5, 1.5],
        Dict{Symbol,Any}(
            :errors => zeros(2),
            :error_statistics => :poisson,
            :units => :rate,
            :exposure_time => 10.0,
        ),
    )
    rebinned_poisson_rate = rebin_channels(poisson_rate_spec, [1, 0])

    @test flux_axis(rebinned_poisson_rate) == [2.0]
    @test only(rebinned_poisson_rate.errors) ≈ XraySpectra.count_error(20.0, 1.0) / 10.0

    poisson_unknown_units = SpectrumBase.Spectrum(
        [0, 1],
        [1.0, 2.0],
        Dict{Symbol,Any}(
            :errors => zeros(2),
            :error_statistics => :poisson,
        ),
    )

    @test_throws ArgumentError rebin_channels(poisson_unknown_units, [1, 0])

    binned_spec = SpectrumBase.Spectrum(
        [0.1 0.2; 0.2 0.3; 0.3 0.4; 0.4 0.5],
        [10.0, 20.0, 30.0, 40.0],
        Dict{Symbol,Any}(),
    )
    rebinned_binned_spec = rebin_channels(binned_spec, grouping)

    @test spectral_axis(rebinned_binned_spec) == [0.1 0.3; 0.3 0.5]
    @test flux_axis(rebinned_binned_spec) == [30.0, 70.0]

    response = ResponseMatrix(
        sparse([1.0 0.0 2.0; 0.0 3.0 0.0; 4.0 0.0 5.0; 0.0 6.0 0.0]),
        [0, 1, 2, 3],
        [0.1 0.2; 0.2 0.3; 0.3 0.4; 0.4 0.5],
        [1.0 2.0; 2.0 3.0; 3.0 4.0],
        FullResponse,
    )
    rebinned_response = rebin_channels(response, grouping)

    @test rebinned_response.matrix isa SparseMatrixCSC
    @test Matrix(rebinned_response.matrix) == [1.0 3.0 2.0; 4.0 6.0 5.0]
    @test rebinned_response.channels == [1, 2]
    @test channel_bins(rebinned_response) == [0.1 0.3; 0.3 0.5]
    @test response_bins(rebinned_response) == response_bins(response)
    @test arf_folded(rebinned_response)
    @test response_kind(rebinned_response) == FullResponse

    ancillary = AncillaryResponse([1.0 2.0; 2.0 3.0; 3.0 4.0], [10.0, 20.0, 30.0])
    data = (;
        spectrum = spec,
        response = response,
        ancillary = ancillary,
        background = spec,
        paths = (; spectrum = "source.pha"),
    )
    rebinned_data = rebin_channels(data, grouping)

    @test flux_axis(rebinned_data.spectrum) == [30.0, 70.0]
    @test Matrix(rebinned_data.response.matrix) == [1.0 3.0 2.0; 4.0 6.0 5.0]
    @test flux_axis(rebinned_data.background) == [30.0, 70.0]
    @test rebinned_data.ancillary === ancillary
    @test rebinned_data.paths === data.paths

    spectrum_only_rebinned = rebin_channels(
        data;
        factor = 2,
        rebin_response = false,
        rebin_background = false,
    )

    @test flux_axis(spectrum_only_rebinned.spectrum) == [30.0, 70.0]
    @test spectrum_only_rebinned.response === response
    @test spectrum_only_rebinned.background === spec

    @test_throws ArgumentError rebin_channels(response, [1, 0])
    @test_throws ArgumentError rebin_channels(spec, [0, 1, 0, 1])
end

@testset "NuSTAR response and ancillary" begin
    response = read_rmf(NUSTAR_RMF)
    ancillary = read_ancillary_response(NUSTAR_ARF)

    combined = combine(response, ancillary)
    flux = ones(size(response.matrix, 2))

    @test size(ancillary_bins(ancillary)) == size(response_bins(response))
    @test length(effective_area(ancillary)) == size(response.matrix, 2)
    @test size(combined) == size(response.matrix)
    @test combined isa SparseMatrixCSC
    @test fold(response, flux; ancillary = ancillary) == combined * flux
end

@testset "OGIP spectrum paths" begin
    paths = read_paths_from_spectrum(NUSTAR_PHA)

    @test paths.spectrum == NUSTAR_PHA
    @test basename(paths.response) == "nu60001047002A01_sr.rmf"
    @test basename(paths.ancillary) == "nu60001047002A01_sr.arf"
    @test basename(paths.background) == "nu60001047002A01_bk.pha"
    @test isfile(paths.response)
    @test isfile(paths.ancillary)
    @test !isfile(paths.background)
end

@testset "read_dataset flags" begin
    data = read_dataset(NUSTAR_PHA; read_background = false)

    @test data.spectrum isa SingleSpectrum
    @test data.response isa ResponseMatrix
    @test data.ancillary isa AncillaryResponse
    @test isnothing(data.background)
    @test data.paths.response == NUSTAR_RMF
    @test data.paths.ancillary == NUSTAR_ARF

    pha_only = read_dataset(
        NUSTAR_PHA;
        read_response = false,
        read_ancillary = false,
        read_background = false,
    )

    @test pha_only.spectrum isa SingleSpectrum
    @test isnothing(pha_only.response)
    @test isnothing(pha_only.ancillary)
    @test isnothing(pha_only.background)

    @test spectral_axis(read_background(NUSTAR_PHA)) == spectral_axis(read_pha(NUSTAR_PHA))
    @test_throws ArgumentError read_dataset(NUSTAR_PHA; read_background = true)
end

@testset "PHA plus RMF energy bins" begin
    spec = read_pha(NUSTAR_PHA)
    response = read_rmf(NUSTAR_RMF)
    binned = energy_binned_spectrum(spec, response)

    @test binned isa SpectrumBase.Spectrum{S,F,2,1} where {S,F}
    @test size(spectral_axis(binned)) == (length(spec), 2)
    @test spectral_axis(binned)[1, 1] ≈ 1.60
    @test spectral_axis(binned)[1, 2] ≈ 1.64
    @test flux_axis(binned) == flux_axis(spec)
    @test binned.channels == spectral_axis(spec)
    @test binned.units == spec.units
    @test binned.error_statistics == spec.error_statistics

    bad_response = ResponseMatrix(
        response.matrix[2:end, :],
        response.channels[2:end],
        channel_bins(response)[2:end, :],
        response_bins(response),
    )
    @test_throws ArgumentError energy_binned_spectrum(spec, bad_response)
end
