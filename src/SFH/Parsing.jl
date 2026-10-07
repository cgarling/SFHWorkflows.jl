"""Module containing code to parse input configuration file and construct input for Systematics module."""
module Parsing

export parse_config, parse_ssp_config

import Distributions
import YAML
using OrderedCollections: OrderedDict
using InitialMassFunctions
using StarFormationHistories: NoBinaries, RandomBinaryPairs, BinaryMassRatio, MH_from_Z, dMH_dZ, PowerLawMZR, LinearAMR, LogarithmicAMR, GaussianDispersion, snr_magerr, match_extinction
using StellarTracks: PARSECLibrary, MISTv1Library, MISTv2Library, BaSTIv1Library, BaSTIv2Library
using BolometricCorrections: YBCGrid, MISTv1BCGrid, MISTv2BCGrid

strip_whitespace(s::AbstractString) = replace(s, r"\s+" => "")

"""
    parse_to_vector(input::AbstractString)
Parse a string like `"F475W, F606W, F814W" to a `Vector{String}`.
"""
function parse_to_vector(input::AbstractString)
    return String.(split(strip_whitespace(input), ","))
end

"""
    parse_range(input::AbstractString)::Vector{Float64}
Parse a range specified as a string (e.g., input = "range(0, 10; step=0.1)").
"""
parse_range(input::AbstractString) = eval(Meta.parse(input))::StepRangeLen

"""
    parse_imf(dict)
Parses the IMF portion of the configuration dictionary and returns an instantiated IMF object.
"""
function parse_imf(dict)
    imf = dict["imf"]
    model = strip_whitespace(imf["model"])
    valid_models = ("Kroupa2001", "Chabrier2001BPL", "Chabrier2001LogNormal", "Chabrier2003", "Chabrier2003System", "Salpeter1955")
    if !(model ∈ valid_models)
        error("IMF model $model invalid; valid IMF models are $(join(valid_models, ", ")).")
    end
    imf_type = eval(Meta.parse(model))
    default = imf_type() # Create instance with default limits
    mmin = get(imf, "mmin", minimum(default))::Float64
    mmax = get(imf, "mmax", maximum(default))::Float64
    return imf_type(mmin, mmax)
end

function parse_binaries(dict)
    # binary_model = eval(Meta.parse(dict["binaries"]["model"]))
    binary_model = strip_whitespace(dict["binaries"]["model"])
    binary_model = if binary_model == "NoBinaries"
        NoBinaries()
    elseif binary_model == "RandomBinaryPairs"
        RandomBinaryPairs(dict["binaries"]["binary_fraction"])
    elseif binary_model == "BinaryMassRatio" # ponytail: default Uniform(0.1, 1) mass-ratio distribution only; parse `qdist` when needed
        BinaryMassRatio(dict["binaries"]["binary_fraction"])
    else
        error("Binary model $binary_model unrecognized; valid options are NoBinaries, RandomBinaryPairs, and BinaryMassRatio.")
    end
    return binary_model
end

"""
    parse_extinction(dict)
Returns a function of `logAge` giving the distribution of the V-band extinction of stars of that age under the model of
MATCH (`StarFormationHistories.match_extinction`): uniform over `[Av, Av + dAv]` plus an independent uniform
spread up to `dAvy` for young stars, which tapers from full at age `dAvy_t1` to zero at `dAvy_t2` [Gyr]. The function
returns `nothing` for ages at which every star has extinction `Av`.
"""
function parse_extinction(dict)
    p = dict["properties"]
    Av, dAv, dAvy = Float64(p["Av"]), Float64(get(p, "dAv", 0.0)), Float64(get(p, "dAvy", 0.0))
    t1, t2 = 1e9 * Float64(get(p, "dAvy_t1", 0.04)), 1e9 * Float64(get(p, "dAvy_t2", 0.1))
    (dAv >= 0 && dAvy >= 0) || error("Invalid configuration: properties.dAv and properties.dAvy must be non-negative.")
    t1 < t2 || error("Invalid configuration: properties.dAvy_t1 must be less than properties.dAvy_t2.")
    return logAge -> (dAv > 0 || (dAvy > 0 && exp10(logAge) < t2)) ? match_extinction(Av, dAv, dAvy, logAge; t1, t2) : nothing
end

function parse_tracks(dict)
    # Recursion: If this is top-level dict, call again with sub-dictionaries as argument 
    if "stellartracks" ∈ keys(dict)
        # return [parse_tracks(track) for track in dict["stellartracks"]]
        st = dict["stellartracks"]
        # Filter out only keys that contain "track" and sort so that order is track1, track2, track3, etc.
        goodkeys = sort([key for key in keys(st) if occursin("track", key)])
        return [parse_tracks(st[key]) for key in goodkeys]
    end
    valid_models = ("parsec", "mistv1", "mistv2", "bastiv1", "bastiv2")
    name = lowercase(dict["name"])
    if name == "mist"
        @warn "Stellar track name \"mist\" is deprecated; please use \"mistv1\" for MIST v1.2 or \"mistv2\" for MIST v2.5. Assuming \"mistv1\"."
        name = "mistv1"
    end
    if !(name ∈ valid_models)
        error("Stellar track model $name invalid; valid stellar track models are $(join(valid_models, ", ")).")
    end
    if name == "parsec"
        return PARSECLibrary()
    elseif name == "mistv1"
        vvcrit = get(dict, "vvcrit", 0.0)
        return MISTv1Library(vvcrit)
    elseif name == "mistv2"
        vvcrit = get(dict, "vvcrit", 0.0)
        afe = get(dict, "alpha_fe", 0.0)::Float64
        return MISTv2Library(vvcrit, afe)
    elseif name == "bastiv1"
        α_fe = get(dict, "alpha_fe", 0.0)::Float64
        canonical = get(dict, "canonical", false)::Bool
        agb = get(dict, "agb", false)::Bool
        eta = get(dict, "eta", 0.4)::Float64
        return BaSTIv1Library(α_fe, canonical, agb, eta)
    elseif name == "bastiv2"
        α_fe = get(dict, "alpha_fe", 0.0)::Float64
        canonical = get(dict, "canonical", false)::Bool
        diffusion = get(dict, "diffusion", true)::Bool
        yp = get(dict, "yp", 0.247)::Float64
        eta = get(dict, "eta", 0.3)::Float64
        return BaSTIv2Library(α_fe, canonical, diffusion, yp, eta)
    end
end

function parse_bcs(dict)
    # Recursion: If this is top-level dict, call again with sub-dictionaries as argument 
    if "bolometriccorrections" ∈ keys(dict)
        bc = dict["bolometriccorrections"]
        return [parse_bcs(bc[key]) for key in keys(bc)]
    end
    valid_models = ("ybc", "mistv1", "mistv2")
    name = lowercase(dict["name"])
    if name == "mist"
        @warn "Bolometric correction name \"mist\" is deprecated; please use \"mistv1\" for MIST v1.2 or \"mistv2\" for MIST v2.5. Assuming \"mistv1\"."
        name = "mistv1"
    end
    if !(name ∈ valid_models)
        error("Bolometric correction grid $name invalid; valid BC grids are $(join(valid_models, ", ")).")
    end
    filterset = dict["filterset"]
    return if name == "ybc"
        try
            YBCGrid(filterset)
        catch e
            println("Failed to initialize YBC bolometric correction grid with filterset $filterset; error shown below")
            rethrow(e)
        end
    elseif name == "mistv1"
        try
            MISTv1BCGrid(filterset)
        catch e
            println("Failed to initialize MIST v1.2 bolometric correction grid with filterset $filterset; error shown below")
            rethrow(e)
        end
    elseif name == "mistv2"
        try
            MISTv2BCGrid(filterset)
        catch e
            println("Failed to initialize MIST v2.5 bolometric correction grid with filterset $filterset; error shown below")
            rethrow(e)
        end
    end
end

# metallicity: # Options are LinearAMR, LogarithmicAMR, PowerLawMZR, examples below
#   name: PowerLawMZR # [M/H] = alpha * (log10(M_*(t)) - log10(mstar0)) + beta
#   alpha: # power-law slope parameter
#     x0: 1.0 # Initial guess
#     free: true # Whether to parameter should be free to vary (true) or fixed (false)
#   beta: # power-law intercept parameter
#     x0: -1.5 # Initial guess
#     free: true
#   mstar0: 1e6 # Stellar mass normalization; by definition, metallicity is beta at mstar0. Generally leave this be.
#   std: 0.1 # Gaussian σ for the spread in metallicity at fixed time

# A model parameter is either `{x0, free}` (initial guess for fit_sfh) or a plain value, which is free to vary in fit_sfh
parse_param(d::AbstractDict) = (Float64(d["x0"]), Bool(get(d, "free", true)))
parse_param(d) = (Float64(d), true)

"""
    parse_metallicity(dict; T_max=nothing)
Parses a metallicity model section into `(MH_model, disp_model)`. `alpha` and `beta` may each be `{x0, free}` or a plain
value; the AMRs alternatively accept `constraints: [[MH1, t1], [MH2, t2]]` ([M/H] at two lookback times in Gyr), which
are free in fit_sfh. AMRs use `dict["T_max"]` if present and the keyword `T_max` otherwise.
"""
function parse_metallicity(dict; T_max=nothing)
    # Recursion: If this is top-level dict, call again with sub-dictionaries as argument
    if "metallicity" ∈ keys(dict)
        return parse_metallicity(dict["metallicity"]; T_max)
    end
    valid_models = ("PowerLawMZR", "LinearAMR", "LogarithmicAMR")
    name = lowercase(dict["name"])
    if !(name ∈ lowercase.(valid_models))
        error("Metallicity model $name invalid; valid metallicity models are $(join(valid_models, ", ")).")
    end
    T_max = get(dict, "T_max", T_max)
    if name != "powerlawmzr" && isnothing(T_max)
        error("Invalid configuration: metallicity model $(dict["name"]) requires `T_max`.")
    end
    MH_model0 = if haskey(dict, "constraints")
        if name == "powerlawmzr" || haskey(dict, "alpha") || haskey(dict, "beta")
            error("Invalid configuration: metallicity `constraints` are supported only for LinearAMR and LogarithmicAMR, and replace `alpha` and `beta`.")
        end
        c = dict["constraints"]
        if length(c) != 2 || any(x -> length(x) != 2, c)
            error("Invalid configuration: metallicity `constraints` must be two [[M/H], lookback time [Gyr]] pairs.")
        end
        c1, c2 = (Tuple(Float64.(x)) for x in c)
        name == "linearamr" ? LinearAMR(c1, c2, T_max) : LogarithmicAMR(c1, c2, T_max)
    else
        (α, αfree), (β, βfree) = parse_param(dict["alpha"]), parse_param(dict["beta"])
        if name == "powerlawmzr"
            PowerLawMZR(α, β, log10(dict["mstar0"]), (αfree, βfree))
        elseif name == "linearamr"
            LinearAMR(α, β, T_max, (αfree, βfree))
        else
            LogarithmicAMR(α, β, T_max, MH_from_Z, dMH_dZ, (αfree, βfree))
        end
    end
    disp_model0 = GaussianDispersion(dict["std"], (false,))
    return MH_model0, disp_model0
end

"""
    parse_filter_models(dict)
Parses the per-filter SNR(mag) or mag_err(mag) models (with optional bias(mag)) in `dict["filters"]` (e.g., the `data` section of a fit_sfh
configuration) into an `OrderedDict` mapping filter name to a `NamedTuple` of arguments for `ASTs.snr_model`.
Returns an empty `OrderedDict` if `dict` has no `filters` entry.
"""
function parse_filter_models(dict)
    models = OrderedDict{String, NamedTuple}()
    for (name, d) in get(dict, "filters", OrderedDict())
        name = string(name)
        if haskey(d, "snr") == haskey(d, "mag_err")
            error("Invalid configuration: the model for filter $name must give exactly one of `snr` or `mag_err`.")
        end
        mag = Float64.(d["mag"])
        snr = haskey(d, "snr") ? Float64.(d["snr"]) : snr_magerr.(Float64.(d["mag_err"]))
        c = get(d, "completeness", OrderedDict())
        bias = haskey(d, "bias") ? Float64.(d["bias"]) : nothing
        models[name] = (mag=mag, snr=snr, bias=bias, minerr=Float64(get(d, "minerr", 0.0)), snr50=Float64(get(c, "snr50", 5.0)), width=Float64(get(c, "width", 1.0)))
    end
    return models
end

"""
    parse_gates(binning)
Parses the optional `gates` entry of `binning` (the `data.binning` section of a fit_sfh configuration), a list of polygons
each given as a list of at least 3 `[color, magnitude]` vertices, into a `Vector{Vector{NTuple{2, Float64}}}`.
Hess diagram bins whose centers lie inside a gate are excluded from the fit. Returns an empty vector if there are no gates.
"""
function parse_gates(binning)
    gates = Vector{NTuple{2, Float64}}[]
    for (i, g) in enumerate(something(get(binning, "gates", nothing), []))
        if !(g isa AbstractVector && length(g) >= 3 && all(v -> v isa AbstractVector && length(v) == 2 && all(x -> x isa Real, v), g))
            error("Invalid configuration: data.binning.gates entry $i must be a list of at least 3 [color, magnitude] vertices.")
        end
        push!(gates, [(Float64(v[1]), Float64(v[2])) for v in g])
    end
    return gates
end

"""
    check_filter_models(needed, ast_filters, model_filters)
Throws an error if a filter has both an AST model and a per-filter model, if any filter in `needed` has neither, or if
an AST filter is not in `needed` (the AST models are joint functions of both AST filters).
"""
function check_filter_models(needed, ast_filters, model_filters)
    for f in ast_filters
        if f ∉ needed
            error("Invalid configuration: AST filter $f is not used in the binning; the AST models are joint functions of both AST filters, so both must be used in `yfilter` or `xcolor`.")
        end
    end
    for f in intersect(ast_filters, model_filters)
        error("Invalid configuration: filter $f has both an AST model and a per-filter `filters` model; give only one.")
    end
    for f in needed
        if f ∉ ast_filters && f ∉ model_filters
            error("Invalid configuration: no observational model for filter $f; provide one through `ASTs` or `filters`.")
        end
    end
end

"""
    parse_data(config)
Parses the `data` section of a fit_sfh or fit_ssp configuration: the photometry and artificial star test files, the
per-filter observational models, and the Hess diagram binning and gates.
"""
function parse_data(config)
    data_path = config["data"]["path"]
    phot_file = joinpath(data_path, config["data"]["photometry"]["photometry_file"])
    filters = parse_to_vector(config["data"]["photometry"]["filters"])
    asts = get(config["data"], "ASTs", nothing)
    ast_file = isnothing(asts) ? nothing : joinpath(data_path, asts["ast_file"])
    # The two AST columns default to the first two photometry filters
    ast_filters = isnothing(asts) ? String[] : (haskey(asts, "filters") ? parse_to_vector(asts["filters"]) : first(filters, 2))
    if !isnothing(asts) && length(ast_filters) != 2
        error("Invalid configuration: config[\"data\"][\"ASTs\"][\"filters\"] must list exactly two filters.")
    end
    filter_models = parse_filter_models(config["data"])
    for f in Iterators.flatten((ast_filters, keys(filter_models)))
        if f ∉ filters
            error("Invalid configuration: filter $f with an observational model not found in provided config[\"data\"][\"photometry\"][\"filters\"] $(join(filters, ", ")).")
        end
    end
    ystring = config["data"]["binning"]["yfilter"]
    if ystring ∉ filters
        error("Invalid configuration: config[\"data\"][\"binning\"][\"yfilter\"] $ystring not found in provided config[\"data\"][\"photometry\"][\"filters\"] $(join(filters, ", ")).")
    end
    xstrings = string.(split(strip_whitespace(config["data"]["binning"]["xcolor"]), ","))
    for xstring in xstrings
        if xstring ∉ filters
            error("Invalid configuration: config[\"data\"][\"binning\"][\"xcolor\"] entry $xstring not found in provided config[\"data\"][\"photometry\"][\"filters\"] $(join(filters, ", ")).")
        end
    end
    check_filter_models(unique(vcat(ystring, xstrings)), ast_filters, keys(filter_models))
    ybins = parse_range(config["data"]["binning"]["ybins"])
    xbins = parse_range(config["data"]["binning"]["xbins"])
    gates = parse_gates(config["data"]["binning"])
    badval = isnothing(asts) ? 99.999 : asts["badval"]
    maxerr = isnothing(asts) ? Inf : get(asts, "maxerr", Inf)::Float64 # If maxerr not provided, use Inf
    minerr = isnothing(asts) ? 0.0 : get(asts, "minerr", 0.0)::Float64
    return (; phot_file, ast_file, ast_filters, filter_models, filters, badval, maxerr, minerr, ystring, xstrings, xbins, ybins, gates)
end

# This function will parse YAML file to dictionary, then call below function
# that takes input dictionary. This way, if you want, you can load the dict from 
# file, programmatically update the dict, and pass the altered dict to parse_config
# to easily run different variations on the same YAML without manually editing it.
function parse_config(file::AbstractString)
    @info "Parsing config"
    if !isfile(file)
        throw(ArgumentError("Config file $file not found."))
    end

    config = try
        YAML.load_file(file; dicttype=OrderedDict{String, Any})
    catch e
        println("Failed to parse configuration YAML file $file with error: ")
        rethrow(e)
    end
    return parse_config(config)
end

function parse_config(config::AbstractDict)
    output_path = get(config["output"], "path", ".")
    if !isdir(output_path)
        try
            mkdir(output_path)
        catch e
            "Requested output_path $output_path does not exist and attempt to create directory failed. Ensure you have write permissions to this path. Full error: "
            rethrow(e)
        end
    end
    # Save copy of input config to output path
    YAML.write_file(joinpath(output_path, "input.yml"), config)

    data = parse_data(config)
    background = parse_background(config)
    imf = parse_imf(config)
    binary_model = parse_binaries(config)
    @info "Loading stellar tracks"
    stellar_tracks = parse_tracks(config)
    @info "Loading bolometric corrections"
    bcs = parse_bcs(config)
    # Lookback time [Gyr] when star formation begins (right edge of the oldest logAge bin); the AMRs use the same T_max
    T_max = Float64(get(config["properties"], "T_max", get(config["metallicity"], "T_max", 13.7)))
    if haskey(config["metallicity"], "T_max") && config["metallicity"]["T_max"] != T_max
        error("Invalid configuration: metallicity.T_max $(config["metallicity"]["T_max"]) differs from properties.T_max $T_max; give T_max once, in properties.")
    end
    MH_model0, disp_model0 = parse_metallicity(config; T_max)
    logAge = eval(Meta.parse(config["stellartracks"]["logAge"]))
    MH = eval(Meta.parse(config["stellartracks"]["MH"]))

    return (; data..., background, extinction=parse_extinction(config), plot_diagnostics=config["plotting"]["diagnostics"], imf=imf, binary_model=binary_model, Av=config["properties"]["Av"], dmod=config["properties"]["distance_modulus"], Mstar=config["properties"]["Mstar"], stellar_tracks=stellar_tracks, bcs=bcs, MH_model0=MH_model0, disp_model0=disp_model0, output_path=output_path, output_filename=config["output"]["filename"], logAge=logAge, MH=MH, T_max=T_max)
end


#################################
# fit_ssp configuration

# Distributions allowed in prior strings, e.g., "Normal(24.95, 0.1)" or "truncated(Normal(0.1, 0.05), 0, Inf)"
const PRIOR_DISTRIBUTIONS = Dict(:Normal => Distributions.Normal, :LogNormal => Distributions.LogNormal,
                                 :Uniform => Distributions.Uniform, :Beta => Distributions.Beta, :Gamma => Distributions.Gamma,
                                 :Exponential => Distributions.Exponential, :truncated => Distributions.truncated)

"""
    parse_prior(x)
Parses a fit_ssp parameter: a number is a fixed value, returned as a `Float64`, and a string is a prior distribution from
`PRIOR_DISTRIBUTIONS` with numeric arguments (e.g., `"Normal(24.95, 0.1)"`), built without evaluating arbitrary code.
"""
parse_prior(x::Real) = Float64(x)
function parse_prior(s::AbstractString)
    invalid() = error("Invalid prior \"$s\": a prior is a number or one of the distributions $(join(sort(string.(keys(PRIOR_DISTRIBUTIONS))), ", ")) with numeric arguments, e.g., \"Normal(24.95, 0.1)\".")
    function build(ex)
        ex isa Real && return Float64(ex)
        ex === :Inf && return Inf
        if ex isa Expr && ex.head === :call
            f, args = ex.args[1], ex.args[2:end]
            f === :- && length(args) == 1 && return -build(args[1])
            f isa Symbol && haskey(PRIOR_DISTRIBUTIONS, f) && return PRIOR_DISTRIBUTIONS[f](build.(args)...)
        end
        invalid()
    end
    ex = try Meta.parse(s) catch; invalid() end
    d = try
        build(ex)
    catch e
        e isa ErrorException && rethrow()
        error("Invalid prior \"$s\": $(sprint(showerror, e))")
    end
    d isa Union{Real, Distributions.ContinuousUnivariateDistribution} || invalid()
    return d
end

"""
    TransformedPrior(base, f, finv, logdxdy)
Prior on `y = f(x)` given the prior `base` on `x`, for a monotonically increasing `f` with inverse `finv` and
`logdxdy(y) = log(dx/dy)`. Used for priors on age (Gyr) and distance (pc) when sampling logAge and distance modulus.
"""
struct TransformedPrior{D <: Distributions.ContinuousUnivariateDistribution, F, G, J} <: Distributions.ContinuousUnivariateDistribution
    base::D
    f::F
    finv::G
    logdxdy::J
end
function Distributions.logpdf(d::TransformedPrior, y::Real)
    x = d.finv(y)
    return Distributions.insupport(d.base, x) ? Distributions.logpdf(d.base, x) + d.logdxdy(y) : oftype(float(y), -Inf)
end
Distributions.pdf(d::TransformedPrior, y::Real) = exp(Distributions.logpdf(d, y))
Distributions.cdf(d::TransformedPrior, y::Real) = Distributions.cdf(d.base, d.finv(y))
Distributions.quantile(d::TransformedPrior, q::Real) = d.f(Distributions.quantile(d.base, q))
Distributions.insupport(d::TransformedPrior, y::Real) = Distributions.insupport(d.base, d.finv(y))
Base.minimum(d::TransformedPrior) = d.f(minimum(d.base))
Base.maximum(d::TransformedPrior) = d.f(maximum(d.base))
# Half the 16th to 84th percentile range, as the moments of the transformed variable have no closed form
Distributions.std(d::TransformedPrior) = (Distributions.quantile(d, 0.8413447460685429) - Distributions.quantile(d, 0.15865525393145707)) / 2
# Priors on age (Gyr) as logAge = log10(age) + 9, and on distance (pc) as distance modulus = 5 log10(distance) - 5;
# ages and distances are positive, so priors extending below 0 are truncated there
positive(d) = minimum(d) < 0 ? Distributions.truncated(d, 0, maximum(d)) : d
age_prior(d) = TransformedPrior(positive(d), x -> log10(x) + 9, y -> exp10(y - 9), y -> (y - 9) * log(10) + log(log(10)))
distance_prior(d) = TransformedPrior(positive(d), x -> 5 * log10(x) - 5, y -> exp10((y + 5) / 5), y -> (y + 5) / 5 * log(10) + log(log(10) / 5))
age_prior(x::Real) = log10(x) + 9
distance_prior(x::Real) = 5 * log10(x) - 5

# Restricts the parameter `p` (a fixed value or a prior) to [lo, hi]: fixed values outside are an error, and priors are
# truncated, in the coordinates of `base` for a TransformedPrior so that the truncation is exact
function restrict(p::Real, lo, hi, name)
    lo <= p <= hi || error("Invalid configuration: parameters.$name = $p is outside the valid range [$lo, $hi].")
    return p
end
function restrict(p::Distributions.ContinuousUnivariateDistribution, lo, hi, name)
    minimum(p) >= lo && maximum(p) <= hi && return p
    @info "Truncating the prior on parameters.$name to the valid range [$lo, $hi]."
    if p isa TransformedPrior
        return TransformedPrior(Distributions.truncated(p.base, max(minimum(p.base), p.finv(lo)), min(maximum(p.base), p.finv(hi))), p.f, p.finv, p.logdxdy)
    end
    return Distributions.truncated(p, max(minimum(p), lo), min(maximum(p), hi))
end

"""
    parse_parameters(dict, binary_model; background::Bool=true)
Parses the `parameters` section of a fit_ssp configuration into a `NamedTuple` with a fixed value or prior for each of
`logAge`, `MH`, `dmod`, `Av`, `binary_fraction`, and `background_fraction`. Age may be given as `logAge` or `age` (Gyr) and
distance as `distance_modulus` or `distance` (pc). Without a `background` (`data.background: none`), the background
fraction is fixed to 0.
"""
function parse_parameters(dict, binary_model; background::Bool=true)
    p = dict["parameters"]
    function one_of(a, b, fb)
        haskey(p, a) == haskey(p, b) && error("Invalid configuration: give exactly one of parameters.$a or parameters.$b.")
        return haskey(p, a) ? parse_prior(p[a]) : fb(parse_prior(p[b]))
    end
    required(k) = haskey(p, k) ? parse_prior(p[k]) : error("Invalid configuration: parameters.$k is required.")
    logAge = one_of("logAge", "age", age_prior)
    dmod = one_of("distance_modulus", "distance", distance_prior)
    binary_fraction = if binary_model isa NoBinaries
        f = parse_prior(get(p, "binary_fraction", 0.0))
        f isa Real && iszero(f) || error("Invalid configuration: parameters.binary_fraction requires binaries.model RandomBinaryPairs or BinaryMassRatio.")
        f
    else
        haskey(p, "binary_fraction") || error("Invalid configuration: parameters.binary_fraction is required for binaries.model $(nameof(typeof(binary_model))).")
        restrict(parse_prior(p["binary_fraction"]), 0, 1, "binary_fraction")
    end
    background_fraction = if background
        restrict(parse_prior(get(p, "background_fraction", "Uniform(0, 1)")), 0, 1, "background_fraction")
    else
        f = parse_prior(get(p, "background_fraction", 0.0))
        f isa Real && iszero(f) || error("Invalid configuration: data.background: none fixes the background fraction to 0, so parameters.background_fraction must be omitted or 0.")
        f
    end
    return (; logAge, MH=required("MH"), dmod, Av=required("Av"), binary_fraction, background_fraction)
end

# Binary model for fit_ssp, whose fraction is set by parameters.binary_fraction
function parse_ssp_binaries(dict)
    b = dict["binaries"]
    haskey(b, "binary_fraction") && error("Invalid configuration: for fit_ssp, give the binary fraction as parameters.binary_fraction rather than binaries.binary_fraction.")
    model = strip_whitespace(b["model"])
    model == "NoBinaries" && return NoBinaries()
    model == "RandomBinaryPairs" && return RandomBinaryPairs(0.0)
    model == "BinaryMassRatio" && return BinaryMassRatio(0.0) # ponytail: default Uniform(0.1, 1) mass-ratio distribution only, as in parse_binaries
    error("Binary model $model unrecognized; valid options are NoBinaries, RandomBinaryPairs, and BinaryMassRatio.")
end

# Optional background Hess diagram source in `data.background`: a field photometry file with the same filter columns as
# the photometry, or a Hess diagram file written by write_histogram, and the fraction of the background spread uniformly.
# Without the section the background is uniform; `background: none` disables it.
function parse_background(dict)
    b = get(dict["data"], "background", nothing)
    isnothing(b) && return (photometry_file=nothing, hess_file=nothing, floor=0.05, enabled=true)
    b == "none" && return (photometry_file=nothing, hess_file=nothing, floor=0.05, enabled=false)
    b isa AbstractDict || error("Invalid configuration: data.background must be `none` or a section with photometry_file or hess_file.")
    haskey(b, "photometry_file") == haskey(b, "hess_file") && error("Invalid configuration: data.background must give exactly one of photometry_file or hess_file.")
    path(k) = haskey(b, k) ? joinpath(dict["data"]["path"], b[k]) : nothing
    floor = Float64(get(b, "floor", 0.05))
    0 <= floor <= 1 || error("Invalid configuration: data.background.floor must be between 0 and 1.")
    return (photometry_file=path("photometry_file"), hess_file=path("hess_file"), floor, enabled=true)
end

"""
    parse_ssp_config(file::AbstractString)
    parse_ssp_config(config::AbstractDict)
Parses a fit_ssp configuration file or dictionary into a `NamedTuple` of inputs for `fit_ssp`.
"""
function parse_ssp_config(file::AbstractString)
    @info "Parsing config"
    isfile(file) || throw(ArgumentError("Config file $file not found."))
    return parse_ssp_config(YAML.load_file(file; dicttype=OrderedDict{String, Any}))
end
function parse_ssp_config(config::AbstractDict)
    data = parse_data(config)
    background = parse_background(config)
    binary_model = parse_ssp_binaries(config)
    parameters = parse_parameters(config, binary_model; background=background.enabled)
    sampling = get(config, "sampling", OrderedDict{String, Any}())
    nfree = count(v -> !(v isa Real), parameters)
    nwalkers = Int(get(sampling, "nwalkers", max(16, 4 * nfree)))
    (iseven(nwalkers) && nwalkers >= nfree + 2) || error("Invalid configuration: sampling.nwalkers must be even and at least the number of free parameters plus 2 ($(nfree + 2)).")
    fit = get(config, "fit", OrderedDict{String, Any}())
    @info "Loading stellar tracks"
    stellar_tracks = parse_tracks(config)
    @info "Loading bolometric corrections"
    bcs = parse_bcs(config)
    return (; data..., background, imf=parse_imf(config), binary_model, parameters, stellar_tracks, bcs,
            ngrid=Tuple(Int.(get(fit, "ngrid", [16, 12]))), restarts=Int(get(fit, "restarts", 3)), run_sampling=Bool(get(sampling, "run", true)),
            nsteps=Int(get(sampling, "nsteps", 3000)), nburnin=Int(get(sampling, "nburnin", 1000)), nwalkers,
            seed=get(sampling, "seed", nothing), plot_diagnostics=Bool(get(get(config, "plotting", Dict()), "diagnostics", true)),
            output_path=get(config["output"], "path", "."), output_filename=config["output"]["filename"], config)
end

end # module
