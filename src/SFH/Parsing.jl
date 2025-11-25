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
parse_to_vector(input::AbstractString) = String.(split(strip_whitespace(input), ","))


"""
    parse_range(input::AbstractString)::Vector{Float64}
Parse a range specified as a string (e.g., input = "range(0, 10; step=0.1)").
"""
# parse_range(input::AbstractString) = eval(Meta.parse(input))::StepRangeLen
parse_range(::Type{Float64}, input::AbstractString) = eval(Meta.parse(input))::StepRangeLen
function parse_range(::Type{Float32}, input::AbstractString)
    r = parse_range(Float64, input)
    return range(Float32(first(r)), Float32(last(r)); step=Float32(step(r)))::StepRangeLen
end

"""
    parse_imf(dict)
Parses the IMF portion of the configuration dictionary and returns an instantiated IMF object.
"""
function parse_imf(dtype::Union{Type{Float32}, Type{Float64}}, dict)
    imf = dict["imf"]
    model = strip_whitespace(imf["model"])
    valid_models = ("Kroupa2001", "Chabrier2001BPL", "Chabrier2001LogNormal", "Chabrier2003", "Chabrier2003System", "Salpeter1955")
    if !(model ∈ valid_models)
        error("IMF model $model invalid; valid IMF models are $(join(valid_models, ", ")).")
    end
    imf_type = eval(Meta.parse(model))
    default = imf_type() # Create instance with default limits
    mmin = dtype(get(imf, "mmin", minimum(default)))
    mmax = dtype(get(imf, "mmax", maximum(default)))
    return imf_type(mmin, mmax)
end

function parse_binaries(dtype::Union{Type{Float32}, Type{Float64}}, dict)
    # binary_model = eval(Meta.parse(dict["binaries"]["model"]))
    binary_model = strip_whitespace(dict["binaries"]["model"])
    binary_model = if binary_model == "NoBinaries"
        NoBinaries()
    elseif binary_model == "RandomBinaryPairs"
        RandomBinaryPairs(dtype(dict["binaries"]["binary_fraction"]))
    else
        error("Binary model $binary_model unrecognized; valid options are NoBinaries and RandomBinaryPairs.")
    end
    return binary_model
end

function parse_tracks(dtype::Union{Type{Float32}, Type{Float64}}, dict)
    # Recursion: If this is top-level dict, call again with sub-dictionaries as argument 
    if "stellartracks" ∈ keys(dict)
        # return [parse_tracks(track) for track in dict["stellartracks"]]
        st = dict["stellartracks"]
        # Filter out only keys that contain "track" and sort so that order is track1, track2, track3, etc.
        goodkeys = sort([key for key in keys(st) if occursin("track", key)])
        return [parse_tracks(dtype, st[key]) for key in goodkeys]
    end
    valid_models = ("parsec", "mist", "bastiv1", "bastiv2")
    name = lowercase(dict["name"])
    if !(name ∈ valid_models)
        error("Stellar track model $name invalid; valid stellar track models are $(join(valid_models, ", ")).")
    end
    if name == "parsec"
        return PARSECLibrary()
    elseif name == "mist"
        vvcrit = dtype(get(dict, "vvcrit", 0.0))
        return MISTLibrary(vvcrit)
    elseif name == "bastiv1"
        α_fe = dtype(get(dict, "alpha_fe", 0.0))
        canonical = get(dict, "canonical", false)::Bool
        agb = get(dict, "agb", false)::Bool
        eta = dtype(get(dict, "eta", 0.4))
        return BaSTIv1Library(α_fe, canonical, agb, eta)
    elseif name == "bastiv2"
        α_fe = dtype(get(dict, "alpha_fe", 0.0))
        canonical = get(dict, "canonical", false)::Bool
        diffusion = get(dict, "diffusion", true)::Bool
        yp = dtype(get(dict, "yp", 0.247))
        eta = dtype(get(dict, "eta", 0.3))
        return BaSTIv2Library(α_fe, canonical, diffusion, yp, eta)
    end
end

function parse_bcs(dict)
    # Recursion: If this is top-level dict, call again with sub-dictionaries as argument 
    if "bolometriccorrections" ∈ keys(dict)
        bc = dict["bolometriccorrections"]
        return [parse_bcs(bc[key]) for key in keys(bc)]
    end
    valid_models = ("ybc", "mist")
    name = lowercase(dict["name"])
    if !(name ∈ valid_models)
        error("Bolometric correction grid $name invalid; valid BC grids are are $(join(valid_models, ", ")).")
    end
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

function parse_metallicity(dtype::Union{Type{Float32}, Type{Float64}}, dict)
    # Recursion: If this is top-level dict, call again with sub-dictionaries as argument 
    if "metallicity" ∈ keys(dict)
        d = dict["metallicity"]
        return parse_metallicity(dtype, d)
    end
    valid_models = ("PowerLawMZR", "LinearAMR", "LogarithmicAMR") # ("powerlawmzr", "linearamr", "logarithmicamr")
    name = lowercase(dict["name"])
    if !(name ∈ lowercase.(valid_models))
        error("Metallicity model $name invalid; valid metallicity models are $(join(valid_models, ", ")).")
    end
    MH_model0 = if name == "powerlawmzr"
        PowerLawMZR(dtype(dict["alpha"]["x0"]), dtype(dict["beta"]["x0"]), dtype(log10(dict["mstar0"])), (dict["alpha"]["free"], dict["beta"]["free"]))
    elseif name == "linearamr"
        LinearAMR(dtype(dict["alpha"]["x0"]), dtype(dict["beta"]["x0"]), dtype(dict["T_max"]), (dict["alpha"]["free"], dict["beta"]["free"]))
    elseif name == "logarithmicamr"
        LogarithmicAMR(dtype(dict["alpha"]["x0"]), dtype(dict["beta"]["x0"]), dtype(dict["T_max"]), MH_from_Z, dMH_dZ, (dict["alpha"]["free"], dict["beta"]["free"]))
    end
    disp_model0 = GaussianDispersion(dtype(dict["std"]), (false,))
    return MH_model0, disp_model0
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
        YAML.load_file(file)
    catch e
        println("Failed to parse configuration YAML file $file with error: ")
        rethrow(e)
    end
    return parse_config(config)
end

# first parse data type; with that known, then parse the rest
function parse_config(config::AbstractDict)
    dtype = config["data"]["fp"]
    if dtype == "Float32"
        dtype = Float32
    elseif dtype == "Float64"
        dtype = Float64
    else
        throw(ArgumentError("`data.fp` must be either Float32 or Float64."))
    end
    return parse_config(dtype, config)
end
function parse_config(dtype::Union{Type{Float32}, Type{Float64}}, config::AbstractDict)
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
    ybins = parse_range(dtype, config["data"]["binning"]["ybins"])
    xbins = parse_range(dtype, config["data"]["binning"]["xbins"])
    maxerr = dtype(get(config["data"]["ASTs"], "maxerr", Inf)) # If maxerr not provided, use Inf
    minerr = dtype(get(config["data"]["ASTs"], "minerr", 0.0))
    imf = parse_imf(dtype, config)
    binary_model = parse_binaries(dtype, config)
    @info "Loading stellar tracks"
    stellar_tracks = parse_tracks(dtype, config)
    @info "Loading bolometric corrections"
    bcs = parse_bcs(config)
    MH_model0, disp_model0 = parse_metallicity(dtype, config)
    logAge = dtype.(eval(Meta.parse(config["stellartracks"]["logAge"])))
    MH = dtype.(eval(Meta.parse(config["stellartracks"]["MH"])))

    return (phot_file=phot_file, ast_file=ast_file, filters=filters, badval=dtype(config["data"]["ASTs"]["badval"]), maxerr=maxerr, minerr=minerr, xbins=xbins, ybins=ybins, plot_diagnostics=config["plotting"]["diagnostics"], imf=imf, binary_model=binary_model, Av=dtype(config["properties"]["Av"]), dmod=dtype(config["properties"]["distance_modulus"]), Mstar=dtype(config["properties"]["Mstar"]), stellar_tracks=stellar_tracks, bcs=bcs, MH_model0=MH_model0, disp_model0=disp_model0, output_path=output_path, output_filename=config["output"]["filename"], ystring=config["data"]["binning"]["yfilter"], xstrings=string.(split(strip_whitespace(config["data"]["binning"]["xcolor"]), ",")), logAge=logAge, MH=MH, dtype=dtype)
end

end # module
