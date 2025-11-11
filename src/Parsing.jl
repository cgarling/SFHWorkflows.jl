"""Module containing code to parse input configuration file and construct input for Systematics module."""
module Parsing

export parse_config

import YAML
using InitialMassFunctions
using StarFormationHistories: NoBinaries, RandomBinaryPairs, MH_from_Z, dMH_dZ, PowerLawMZR, LinearAMR, LogarithmicAMR, GaussianDispersion
using StellarTracks: PARSECLibrary, MISTLibrary, BaSTIv1Library, BaSTIv2Library
using BolometricCorrections: YBCGrid, MISTBCGrid

strip_whitespace(s::AbstractString) = replace(s, r"\s+" => "")

"""
    parse_to_vector(input::AbstractString)
Parse a string like `"F475W, F606W, F814W" to a `Vector{String}`.
"""
function parse_to_vector(input::AbstractString)
    return String.(split(strip_whitespace(input), ","))
end

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
    else
        error("Binary model $binary_model unrecognized; valid options are NoBinaries and RandomBinaryPairs.")
    end
    return binary_model
end

function parse_tracks(dict)
    # Recursion: If this is top-level dict, call again with sub-dictionaries as argument 
    if "stellartracks" ∈ keys(dict)
        # return [parse_tracks(track) for track in dict["stellartracks"]]
        st = dict["stellartracks"]
        return [parse_tracks(st[key]) for key in keys(st)]
    end
    name = lowercase(dict["name"])
    if name == "parsec"
        return PARSECLibrary()
    elseif name == "mist"
        vvcrit = get(dict, "vvcrit", 0.0)
        return MISTLibrary(vvcrit)
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
        return [parse_tracks(bc[key]) for key in keys(bc)]
    end
    name = lowercase(dict["name"])
    filterset = dict["filterset"]
    return if name == "ybc"
        try
            YBCGrid(filterset)
        catch e
            println("Failed to initialize YBC bolometric correction grid with filterset $filterset; error shown below")
            rethrow(e)
        end
    elseif name == "mist"
        try
            MISTBCGrid(filterset)
        catch e
            println("Failed to initialize MIST bolometric correction grid with filterset $filterset; error shown below")
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

function parse_metallicity(dict)
    # Recursion: If this is top-level dict, call again with sub-dictionaries as argument 
    if "metallicity" ∈ keys(dict)
        d = dict["metallicity"]
        return parse_metallicity(d)
    end
    valid_models = ("PowerLawMZR", "LinearAMR", "LogarithmicAMR") # ("powerlawmzr", "linearamr", "logarithmicamr")
    name = lowercase(dict["name"])
    if !(name ∈ lowercase.(valid_models))
        error("Metallicity model $name invalid; valid metallicity models are $(join(valid_models, ", ")).")
    end
    MH_model0 = if name == "powerlawmzr"
        PowerLawMZR(dict["alpha"]["x0"], dict["beta"]["x0"], dict["mstar0"], (dict["alpha"]["free"], dict["beta"]["free"]))
    elseif name == "linearamr"
        LinearAMR(dict["alpha"]["x0"], dict["beta"]["x0"], dict["T_max"], (dict["alpha"]["free"], dict["beta"]["free"]))
    elseif name == "logarithmicamr"
        LogarithmicAMR(dict["alpha"]["x0"], dict["beta"]["x0"], dict["T_max"], MH_from_Z, dMH_dZ, (dict["alpha"]["free"], dict["beta"]["free"]))
    end
    disp_model0 = GaussianDispersion(dict["std"], (true,))
    return MH_model0, disp_model0
end

function parse_config(file::AbstractString)
    @info "Parsing config"
    if !isfile(file)
        throw(ArgumentError("Config file $file not found."))
    end

    config = try
        YAML.load_file(file)
    catch e
        println("Failed to parse configuration YAML file $file with error: ")
        rethrow(e)
    end

    output_path = get(config["output"], "path", ".")
    if !isdir(output_path)
        try
            mkdir(output_path)
        catch e
            "Requested output_path $output_path does not exist and attempt to create directory failed. Ensure you have write permissions to this path. Full error: "
            rethrow(e)
        end
    end

    data_path = config["data"]["path"]
    phot_file = joinpath(data_path, config["data"]["photometry"]["photometry_file"])
    ast_file = joinpath(data_path, config["data"]["ASTs"]["ast_file"])
    filters = parse_to_vector(config["data"]["photometry"]["filters"])
    ybins = eval(Meta.parse(config["data"]["binning"]["ybins"]))
    xbins = eval(Meta.parse(config["data"]["binning"]["xbins"]))
    maxerr = get(config["data"]["ASTs"], "maxerr", Inf)::Float64 # If maxerr not provided, use Inf
    minerr = get(config["data"]["ASTs"], "minerr", 0.0)::Float64
    imf = parse_imf(config)
    binary_model = parse_binaries(config)
    @info "Loading stellar tracks"
    stellar_tracks = parse_tracks(config)
    @info "Loading bolometric corrections"
    bcs = parse_bcs(config)
    MH_model0, disp_model0 = parse_metallicity(config)

    return (phot_file=phot_file, ast_file=ast_file, filters=filters, badval=config["data"]["ASTs"]["badval"], maxerr=maxerr, minerr=minerr, xbins=xbins, ybins=ybins, plot_diagnostics=config["plotting"]["diagnostics"], imf=imf, binary_model=binary_model, Av=config["properties"]["Av"], dmod=config["properties"]["distance_modulus"], Mstar=config["properties"]["Mstar"], stellar_tracks=stellar_tracks, bcs=bcs, MH_model0=MH_model0, disp_model0=disp_model0,output_path=output_path)
end

end # module
