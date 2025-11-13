"""Module with code to measure resolved star formation histories from provided photometry. Main function is `fit_sfh` which takes in a properly formatted YAML configuration file."""
module SFHFitting

export fit_sfh

import StarFormationHistories as SFH
using ArgCheck: @argcheck, @check
using DelimitedFiles: readdlm

# Includes
include("Parsing.jl") # Code to parse configuration file
using .Parsing
include("ASTs.jl") # Code to analyze artificial star tests
using .ASTs
include("Systematics.jl")
using .Systematics


# Top-level functions
function fit_sfh(obsfile::AbstractString, astfile::AbstractString, filters, xstrings, ystring, edges,
                 MH_model0::SFH.AbstractMetallicityModel, disp_model0::SFH.AbstractDispersionModel, Mstar::Number, stellar_tracks, bcs,
                 dmod::Number, Av::Number, imf, MH, logAge, binary_model::SFH.AbstractBinaryModel, output_filename::AbstractString; 
                 badval::Number=99.999, minerr::Number=0.0, maxerr::Number=0.2, plot_diagnostics::Bool=true, output_path::AbstractString=".") # filters=("mag1", "mag2")
    @argcheck length(xstrings) == 2
    @argcheck Mstar > 0
    filters = string.(filters)
    completeness, bias, err = process_ast_file(astfile, filters, badval, minerr, maxerr, plot_diagnostics, output_path)
    data = readdlm(obsfile, ' ', Float64)
    yidx = findfirst(==(ystring), filters)
    xidxs = [findfirst(==(x), filters) for x in xstrings]
    h = SFH.bin_cmd(view(data, :, xidxs[1]) .- view(data, :, xidxs[2]), view(data, :, yidx); edges=edges)
    s = systematics(MH_model0, disp_model0, Mstar, vec(h.weights), stellar_tracks, bcs, xstrings, ystring, dmod, Av, err, completeness, bias, imf, MH, logAge, edges; binary_model=binary_model, output=joinpath(output_path, output_filename))
    return nothing
end

fit_sfh(@nospecialize(config::NamedTuple)) = fit_sfh(config.phot_file, config.ast_file, config.filters, config.xstrings, config.ystring, (config.xbins, config.ybins), config.MH_model0, config.disp_model0, config.Mstar, config.stellar_tracks, config.bcs, config.dmod, config.Av, config.imf, config.MH, config.logAge, config.binary_model, config.output_filename; badval=config.badval, minerr=config.minerr, maxerr=config.maxerr, plot_diagnostics=config.plot_diagnostics, output_path=config.output_path)
fit_sfh(config_file::AbstractString) = fit_sfh(parse_config(config_file))


end # module