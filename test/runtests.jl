using SFHWorkflows
using SFHWorkflows.SFHFitting.ASTs: snr_model, fill_nan
using SFHWorkflows.SFHFitting.Parsing: parse_gates, parse_filter_models, check_filter_models, parse_binaries, parse_imf, parse_metallicity,
    parse_prior, age_prior, distance_prior, restrict, parse_parameters, TransformedPrior, parse_background
import Distributions
using SFHWorkflows.Simulate: sfh_mass_fractions
using SFHWorkflows.SFHFitting.Plotting: sparse_mask
using SFHWorkflows.SFHFitting.Systematics: optim_status
using SFHWorkflows.SSPFitting: column_format, write_table, read_best
using TypedTables: Table
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

@testset "parse_gates" begin
    @test isempty(parse_gates(Dict())) && isempty(parse_gates(Dict("gates" => nothing)))
    @test parse_gates(Dict("gates" => [[[1, 20], [2.5, 20], [2, 24]]])) == [[(1.0, 20.0), (2.5, 20.0), (2.0, 24.0)]]
    @test_throws "entry 2 must be a list of at least 3" parse_gates(Dict("gates" => [[[1, 20], [2, 20], [2, 24]], [[1, 20], [2, 20]]]))
    @test_throws "entry 1 must be a list of at least 3" parse_gates(Dict("gates" => [[[1, 20, 3], [2, 20], [2, 24]]]))
end

@testset "parse_background" begin
    cfg(b) = Dict("data" => Dict("path" => "dir", "background" => b))
    @test parse_background(Dict("data" => Dict("path" => "dir"))) == (photometry_file=nothing, hess_file=nothing, floor=0.05, enabled=true)
    @test !parse_background(cfg("none")).enabled
    @test parse_background(cfg(Dict("hess_file" => "f.txt", "floor" => 0.1))) == (photometry_file=nothing, hess_file=joinpath("dir", "f.txt"), floor=0.1, enabled=true)
    @test_throws "exactly one of photometry_file or hess_file" parse_background(cfg(Dict("floor" => 0.1)))
    @test_throws "must be `none` or a section" parse_background(cfg("flat"))
    @test_throws "floor must be between 0 and 1" parse_background(cfg(Dict("hess_file" => "f.txt", "floor" => 2)))
end

@testset "check_filter_models" begin
    @test isnothing(check_filter_models(["F475W", "F606W", "F814W"], ["F475W", "F814W"], ["F606W"]))
    @test isnothing(check_filter_models(["F606W"], String[], ["F606W"]))
    @test_throws "has both an AST model" check_filter_models(["F475W", "F814W"], ["F475W", "F814W"], ["F814W"])
    @test_throws "AST filter F475W is not used in the binning" check_filter_models(["F606W", "F814W"], ["F475W", "F814W"], ["F606W"])
    @test_throws "no observational model for filter F606W" check_filter_models(["F475W", "F606W", "F814W"], ["F475W", "F814W"], String[])
end

@testset "parse_prior" begin
    @test parse_prior(24.95) === 24.95
    @test parse_prior("Normal(24.95, 0.1)") == Distributions.Normal(24.95, 0.1)
    @test parse_prior("Uniform(-2, -0.5)") == Distributions.Uniform(-2.0, -0.5)
    d = parse_prior("truncated(Normal(0.1, 0.05), 0, Inf)")
    @test (minimum(d), maximum(d)) == (0.0, Inf)
    @test d.untruncated == Distributions.Normal(0.1, 0.05)
    @test_throws "Invalid prior \"Poisson(3)\"" parse_prior("Poisson(3)")
    @test_throws "Invalid prior \"run(`ls`)\"" parse_prior("run(`ls`)")
    @test_throws "Invalid prior \"Normal(1, \"" parse_prior("Normal(1, ")
    @test_throws "Invalid prior \"Normal(0, -1)\"" parse_prior("Normal(0, -1)")
end

@testset "TransformedPrior" begin
    # A Gaussian prior on age in Gyr, as a prior on logAge
    base = Distributions.Normal(5.0, 0.5)
    d = age_prior(base)
    @test d isa TransformedPrior
    la = 9.7
    @test Distributions.logpdf(d, la) ≈ Distributions.logpdf(base, exp10(la - 9)) + log(exp10(la - 9) * log(10))
    grid = range(9.0, 10.2; length=20_001)
    @test sum(Distributions.pdf.(d, grid)) * step(grid) ≈ 1 rtol=1e-4 # Normalized in logAge
    @test Distributions.cdf(d, Distributions.quantile(d, 0.3)) ≈ 0.3
    @test Distributions.quantile(d, 0.5) ≈ log10(5.0) + 9
    @test age_prior(5.0) ≈ log10(5.0) + 9
    @test distance_prior(1e6) ≈ 25.0 # 1 Mpc
    dd = distance_prior(Distributions.Normal(1e6, 5e4))
    @test Distributions.quantile(dd, 0.5) ≈ 25.0
    dgrid = range(24.0, 26.0; length=20001)
    @test sum(Distributions.pdf.(dd, dgrid)) * step(dgrid) ≈ 1 rtol=1e-4 # Normalized in distance modulus
    # Truncation of a transformed prior is applied to the underlying age distribution
    t = restrict(d, 9.6, 10.13, "logAge")
    @test all((minimum(t), maximum(t)) .≈ (9.6, 10.13))
    @test Distributions.logpdf(t, 9.5) == -Inf
    @test restrict(Distributions.Uniform(9.0, 10.0), 9.0, 10.13, "logAge") == Distributions.Uniform(9.0, 10.0)
    @test (minimum(restrict(Distributions.Normal(10.0, 0.3), 9.0, 10.13, "logAge")), maximum(restrict(Distributions.Normal(10.0, 0.3), 9.0, 10.13, "logAge"))) == (9.0, 10.13)
    @test_throws "outside the valid range" restrict(10.5, 9.0, 10.13, "logAge")
end

@testset "parse_parameters" begin
    p = Dict("logAge" => "Uniform(9.0, 10.1)", "MH" => -1.2, "distance_modulus" => "Normal(24.95, 0.1)", "Av" => 0.1, "binary_fraction" => "Uniform(0, 1)")
    r = parse_parameters(Dict("parameters" => p), SFH.BinaryMassRatio(0.0))
    @test keys(r) == (:logAge, :MH, :dmod, :Av, :binary_fraction, :background_fraction)
    @test r.background_fraction == Distributions.Uniform(0.0, 1.0)
    @test r.MH === -1.2
    r2 = parse_parameters(Dict("parameters" => merge(Dict(k => v for (k, v) in p if k ∉ ("logAge", "distance_modulus")), Dict("age" => 5.0, "distance" => "Normal(1e6, 5e4)"))), SFH.BinaryMassRatio(0.0))
    @test r2.logAge ≈ log10(5.0) + 9
    @test r2.dmod isa TransformedPrior
    @test_throws "exactly one of parameters.logAge or parameters.age" parse_parameters(Dict("parameters" => merge(p, Dict("age" => 5.0))), SFH.NoBinaries())
    @test_throws "parameters.MH is required" parse_parameters(Dict("parameters" => Dict(k => v for (k, v) in p if k != "MH")), SFH.BinaryMassRatio(0.0))
    @test_throws "requires binaries.model" parse_parameters(Dict("parameters" => p), SFH.NoBinaries())
    @test parse_parameters(Dict("parameters" => Dict(k => v for (k, v) in p if k != "binary_fraction")), SFH.NoBinaries()).binary_fraction == 0
    @test_throws "parameters.binary_fraction is required" parse_parameters(Dict("parameters" => Dict(k => v for (k, v) in p if k != "binary_fraction")), SFH.RandomBinaryPairs(0.0))
    # data.background: none fixes the background fraction to 0
    @test parse_parameters(Dict("parameters" => p), SFH.BinaryMassRatio(0.0); background=false).background_fraction === 0.0
    @test parse_parameters(Dict("parameters" => merge(p, Dict("background_fraction" => 0))), SFH.BinaryMassRatio(0.0); background=false).background_fraction == 0
    for bf in ("Uniform(0, 0.5)", 0.1)
        @test_throws "background_fraction must be omitted or 0" parse_parameters(Dict("parameters" => merge(p, Dict("background_fraction" => bf))), SFH.BinaryMassRatio(0.0); background=false)
    end
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

@testset "parse_metallicity" begin
    d = Dict("name" => "LinearAMR", "alpha" => Dict("x0" => 0.1, "free" => false), "beta" => -2.5, "std" => 0.1, "T_max" => 13.7)
    m, disp = parse_metallicity(Dict("metallicity" => d))
    @test m == SFH.LinearAMR(0.1, -2.5, 13.7, (false, true))
    @test disp == SFH.GaussianDispersion(0.1, (false,))
    # Constraints, with T_max from the keyword when the section has none
    c = Dict("name" => "LogarithmicAMR", "constraints" => [[-2.5, 13.7], [-1.0, 0.0]], "std" => 0.1)
    @test first(parse_metallicity(c; T_max=13.7)) == SFH.LogarithmicAMR((-2.5, 13.7), (-1.0, 0.0), 13.7)
    @test first(parse_metallicity(merge(c, Dict("name" => "LinearAMR")); T_max=13.7)) == SFH.LinearAMR((-2.5, 13.7), (-1.0, 0.0), 13.7)
    @test first(parse_metallicity(Dict("name" => "PowerLawMZR", "alpha" => 0.3, "beta" => -1.5, "mstar0" => 1e6, "std" => 0.1))) == SFH.PowerLawMZR(0.3, -1.5, 6.0)
    @test_throws "requires `T_max`" parse_metallicity(c)
    @test_throws "supported only for LinearAMR and LogarithmicAMR" parse_metallicity(merge(c, Dict("alpha" => 0.1)); T_max=13.7)
    @test_throws "must be two" parse_metallicity(merge(c, Dict("constraints" => [[-2.5, 13.7]])); T_max=13.7)
end

@testset "sfh_mass_fractions" begin
    logAge = [9.0, 9.5, 10.0]
    t = [1.0, exp10(0.5), 10.0, 13.0] # Bin edges in lookback time [Gyr]
    # Constant SFR: mass in each bin is proportional to its duration
    f = sfh_mass_fractions(Dict("T_max" => 13.0, "model" => "constant"), logAge)
    @test f ≈ diff(t) ./ 12
    # Cumulative SFH: 40% formed between 13 and 10 Gyr, the rest between 10 Gyr and exp10(0.5) Gyr
    sfh = Dict("T_max" => 13.0, "model" => "cumulative", "logAge" => [9.5, 10.0], "cum_sfh" => [1.0, 0.4])
    @test sfh_mass_fractions(sfh, logAge) ≈ [0, 0.6, 0.4]
    # Linear in lookback time (constant SFR) between 13 Gyr and 1 Gyr, so half of the mass forms before 7 Gyr
    @test sfh_mass_fractions(merge(sfh, Dict("logAge" => [9.0], "cum_sfh" => [1.0])), [9.0, log10(7e9)]) ≈ [0.5, 0.5]
    @test_throws "must increase toward the present" sfh_mass_fractions(merge(sfh, Dict("cum_sfh" => [0.4, 1.0])), logAge)
    @test_throws "younger than the youngest" sfh_mass_fractions(merge(sfh, Dict("logAge" => [8.5, 10.0])), logAge)
    @test_throws "must be older than the oldest" sfh_mass_fractions(Dict("T_max" => 5.0, "model" => "constant"), logAge)
end

@testset "output helpers" begin
    # 30 points at one spot and 2 isolated points; only the isolated points are in sparse bins
    x, y = vcat(fill(0.5, 30), 0.0, 1.0), vcat(fill(0.5, 30), 0.0, 1.0)
    @test sparse_mask(x, y; bins=4, threshold=10) == vcat(falses(30), true, true)
    @test column_format(:background_fraction) == column_format(:mass) == "%.4e"
    @test column_format(:logAge) == "%.5f"
    @test column_format(:walker) == "%d"
    # The best fit is read back from a summary table as written by fit_ssp
    f = tempname()
    write_table(f, Table([(name="A", lp_best=-1.0, logAge_best=9.7, logAge_lower=9.6, mass_best=5e5)]), ["comment"])
    @test read_best(f, "A") == (logAge=9.7, mass=5e5)
    # Convergence of an optimization that converges and of one stopped by its iteration limit
    quad(x) = sum(abs2, x .- (1, 2))
    ok = optim_status(SFH.Optim.optimize(quad, [0.0, 0.0], SFH.Optim.BFGS()))
    @test ok.converged && ok.iterations > 0 && ok.g_residual < 1e-8
    stopped = optim_status(SFH.Optim.optimize(x -> sum(abs2, x .- (1, 2)) + x[1]^4, [10.0, 10.0], SFH.Optim.BFGS(), SFH.Optim.Options(iterations=1)))
    @test !stopped.converged && stopped.iterations == 1 && stopped.termination == "Iterations"
end
