using Test
using Artifacts
using FITSFiles
using SparseArrays
using SpectrumBase
using XraySpectra

const TEST_DATA = artifact"test_data"
const NUSTAR_PHA = joinpath(TEST_DATA, "nustar", "nu60001047002A01_sr_grp_simple.pha")
const NUSTAR_RMF = joinpath(TEST_DATA, "nustar", "nu60001047002A01_sr.rmf")
const CHANDRA_PHA2 = joinpath(TEST_DATA, "chandra", "acisf13839N002_pha2.fits")

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

    spec = parse_hdu(PHA, fake_hdu)

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

@testset "ARF reader" begin
    fake_hdu = FITSFiles.HDU(
        (;
            ENERG_LO = [1.0, 2.0, 3.0],
            ENERG_HI = [2.0, 3.0, 4.0],
            SPECRESP = [10.0, 20.0, 30.0],
        ),
        FITSFiles.Card[FITSFiles.Card("EXTNAME", "SPECRESP")],
    )

    arf = parse_hdu(ARF, fake_hdu)

    @test arf isa AncillaryResponse
    @test ancillary_bins(arf) == [1.0 2.0; 2.0 3.0; 3.0 4.0]
    @test ancillary_bins_low(arf) == [1.0, 2.0, 3.0]
    @test ancillary_bins_high(arf) == [2.0, 3.0, 4.0]
    @test effective_area(arf) == [10.0, 20.0, 30.0]
end

@testset "OGIP spectrum paths" begin
    paths = read_paths_from_spectrum(NUSTAR_PHA)

    @test paths.spectrum == NUSTAR_PHA
    @test basename(paths.response) == "nu60001047002A01_sr.rmf"
    @test basename(paths.ancillary) == "nu60001047002A01_sr.arf"
    @test basename(paths.background) == "nu60001047002A01_bk.pha"
    @test isfile(paths.response)
    @test !isfile(paths.ancillary)
    @test !isfile(paths.background)
end

@testset "read_spectrum flags" begin
    data = read_spectrum(NUSTAR_PHA; read_ancillary = false, read_background = false)

    @test data.spectrum isa SingleSpectrum
    @test data.response isa ResponseMatrix
    @test isnothing(data.ancillary)
    @test isnothing(data.background)
    @test data.paths.response == NUSTAR_RMF

    pha_only = read_spectrum(
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
    @test_throws ArgumentError read_spectrum(NUSTAR_PHA)
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
