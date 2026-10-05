"""Module containing code to parse input configuration file and construct input for Systematics module."""
module Parsing

export parse_config

import YAML
using OrderedCollections: OrderedDict
using InitialMassFunctions
using StarFormationHistories: NoBinaries, RandomBinaryPairs, BinaryMassRatio, MH_from_Z, dMH_dZ, PowerLawMZR, LinearAMR, LogarithmicAMR, GaussianDispersion, snr_magerr
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

    return (phot_file=phot_file, ast_file=ast_file, ast_filters=ast_filters, filter_models=filter_models, filters=filters, badval=badval, maxerr=maxerr, minerr=minerr, xbins=xbins, ybins=ybins, gates=gates, plot_diagnostics=config["plotting"]["diagnostics"], imf=imf, binary_model=binary_model, Av=config["properties"]["Av"], dmod=config["properties"]["distance_modulus"], Mstar=config["properties"]["Mstar"], stellar_tracks=stellar_tracks, bcs=bcs, MH_model0=MH_model0, disp_model0=disp_model0, output_path=output_path, output_filename=config["output"]["filename"], ystring=ystring, xstrings=xstrings, logAge=logAge, MH=MH, T_max=T_max)
end

end # module
