import YAML
using InitialMassFunctions
using StarFormationHistories: NoBinaries, RandomBinaryPairs

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

function parse_config(file::AbstractString)
    file
    if !isfile(file)
        throw(ArgumentError("Config file $file not found."))
    end

    config = try
        YAML.load_file(file)
    catch e
        println("Failed to parse configuration YAML file $file with error: ")
        rethrow(e)
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

    return (phot_file=phot_file, ast_file=ast_file, filters=filters, badval=config["data"]["ASTs"]["badval"], maxerr=maxerr, minerr=minerr, xbins=xbins, ybins=ybins, plot_diagnostics=config["plotting"]["diagnostics"], imf=imf, binary_model=binary_model, Av=config["properties"]["Av"], dmod=config["properties"]["distance_modulus"], Mstar=config["properties"]["Mstar"])
end

# parse_config("sfh_scripting.jl")
# parse_config("config.yml")