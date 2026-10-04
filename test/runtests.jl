using SFHWorkflows
using SFHWorkflows.SFHFitting.ASTs: snr_model, fill_nan
using SFHWorkflows.SFHFitting.Parsing: parse_filter_models, check_filter_models, parse_binaries, parse_imf
import StarFormationHistories as SFH
using Test

@testset "snr_model" begin
    # log10(SNR) falls linearly from 2 at mag 20 to 1 at mag 25
    c, b, e = snr_model([20.0, 25.0], [100.0, 10.0])
    @test e(22.5) ≈ SFH.magerr_snr(exp10(1.5))
    @test e(30.0) ≈ SFH.magerr_snr(exp10(0.0)) # Linear extrapolation in log10(SNR) past the faint end
    @test b(27.0) == 0
    m50 = 20 + (2 - log10(5)) / 0.2 # SNR = snr50 = 5
    @test c(m50) ≈ 0.5
    @test c(20.0) ≈ 1
    @test c(35.0) < 1e-6 # Completeness keeps falling past the table
    @test snr_model([20.0, 25.0], [100.0, 10.0]; minerr=0.05).err(22.5) == 0.05
    @test snr_model([25.0, 20.0], [10.0, 100.0]).err(22.5) ≈ e(22.5) # Input order does not matter
    @test_throws "must not increase toward fainter" snr_model([20.0, 25.0, 26.0], [100.0, 10.0, 20.0])
    @test_throws "between the two faintest" snr_model([20.0, 25.0, 26.0], [100.0, 10.0, 10.0])
    @test_throws "must be positive" snr_model([20.0, 25.0], [100.0, 0.0])
    bias = snr_model([25.0, 20.0], [10.0, 100.0]; bias=[-0.02, 0.0]).bias
    @test (bias(22.5), bias(19.0), bias(30.0)) == (-0.01, 0.0, -0.02) # Linear inside the table, constant outside
    @test_throws "one value per entry of `mag`" snr_model([20.0, 25.0], [100.0, 10.0]; bias=[0.0])
end

@testset "fill_nan" begin
    # A row without stars (NaN) takes the nearest row; a column without stars takes the nearest column
    @test fill_nan([1.0 NaN 3.0; NaN NaN NaN; 7.0 NaN 9.0]) == [1.0 1.0 3.0; 1.0 1.0 3.0; 7.0 7.0 9.0]
    @test isequal(fill_nan(fill(NaN, 2, 2)), fill(NaN, 2, 2))
end

@testset "parse_filter_models" begin
    @test isempty(parse_filter_models(Dict()))
    d = Dict("filters" => Dict("F606W" => Dict("mag" => [20, 25], "snr" => [100, 10])))
    @test parse_filter_models(d)["F606W"] == (mag=[20.0, 25.0], snr=[100.0, 10.0], bias=nothing, minerr=0.0, snr50=5.0, width=1.0)
    d = Dict("filters" => Dict("F606W" => Dict("mag" => [20, 25], "mag_err" => [0.01, 0.1], "bias" => [0, -0.02], "minerr" => 0.02, "completeness" => Dict("snr50" => 3, "width" => 0.5))))
    m = parse_filter_models(d)["F606W"]
    @test m.snr ≈ SFH.snr_magerr.([0.01, 0.1])
    @test (m.bias, m.minerr, m.snr50, m.width) == ([0.0, -0.02], 0.02, 3.0, 0.5)
    @test_throws "exactly one of `snr` or `mag_err`" parse_filter_models(Dict("filters" => Dict("F606W" => Dict("mag" => [20, 25]))))
    @test_throws "exactly one of `snr` or `mag_err`" parse_filter_models(Dict("filters" => Dict("F606W" => Dict("mag" => [20, 25], "snr" => [100, 10], "mag_err" => [0.01, 0.1]))))
end

@testset "check_filter_models" begin
    @test isnothing(check_filter_models(["F475W", "F606W", "F814W"], ["F475W", "F814W"], ["F606W"]))
    @test isnothing(check_filter_models(["F606W"], String[], ["F606W"]))
    @test_throws "has both an AST model" check_filter_models(["F475W", "F814W"], ["F475W", "F814W"], ["F814W"])
    @test_throws "AST filter F475W is not used in the binning" check_filter_models(["F606W", "F814W"], ["F475W", "F814W"], ["F606W"])
    @test_throws "no observational model for filter F606W" check_filter_models(["F475W", "F606W", "F814W"], ["F475W", "F814W"], String[])
end

@testset "parse_binaries" begin
    @test parse_binaries(Dict("binaries" => Dict("model" => "NoBinaries"))) isa SFH.NoBinaries
    @test parse_binaries(Dict("binaries" => Dict("model" => "RandomBinaryPairs", "binary_fraction" => 0.4))) == SFH.RandomBinaryPairs(0.4)
    bm = parse_binaries(Dict("binaries" => Dict("model" => "BinaryMassRatio", "binary_fraction" => 0.3)))
    @test bm isa SFH.BinaryMassRatio
    @test bm.fraction == 0.3
    @test_throws "valid options are NoBinaries, RandomBinaryPairs, and BinaryMassRatio" parse_binaries(Dict("binaries" => Dict("model" => "Pairs")))
end

@testset "parse_imf" begin
    imf = parse_imf(Dict("imf" => Dict("model" => "Kroupa2001", "mmin" => 0.1, "mmax" => 50.0)))
    @test extrema(imf) == (0.1, 50.0)
    @test_throws "IMF model Kroupa invalid" parse_imf(Dict("imf" => Dict("model" => "Kroupa")))
end
