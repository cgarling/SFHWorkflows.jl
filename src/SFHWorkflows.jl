module SFHWorkflows

export fit_sfh

import StarFormationHistories as SFH
using ArgCheck: @argcheck, @check
using DelimitedFiles: readdlm
using StatsBase: Histogram
using CairoMakie: Makie, set_theme!, theme_latexfonts # CairoMakie re-exports Makie
set_theme!(theme_latexfonts(); 
        fontsize=20,
        Axis = (xticks = Makie.LinearTicks(5),
                yticks = Makie.LinearTicks(7),
                # xticks=Makie.WilkinsonTicks(10; k_min=5, k_max=5),
                # yticks=Makie.WilkinsonTicks(5; k_min=5, k_max=5)))
                # xminorticks=Makie.IntervalsBetween(5),
                xminorticksvisible=true),
            Scatter = (strokecolor=:black, strokewidth=1))
                # xminorgridvisible=true, 
                # yminorgridvisible=true))

function Makie.convert_arguments(P::Makie.CellGrid, h::Histogram)
    return Makie.convert_arguments(P, h.edges[1], h.edges[2], h.weights)
end
Makie.plottype(::Histogram) = Makie.Heatmap

# Includes
include("Parsing.jl") # Code to parse configuration file
using .Parsing
include("ASTs.jl") # Code to analyze artificial star tests
using .ASTs
include("Systematics.jl")
using .Systematics


# systematics(MH_model0::SFH.AbstractMetallicityModel,
#                      disp_model0::SFH.AbstractDispersionModel,
#                      mstar::Number, # Estimate of stellar mass of galaxy
#                      data::Union{AbstractVector{<:Number}, AbstractMatrix{<:Number}},
#                      tracklibs, bclibs,
#                      xstrings, ystring, # strings or symbols that allow us to select the correct filters from isochrone
#                      dmod, Av, err_funcs, complete_funcs, bias_funcs, imf,
#                      unique_MH, unique_logAge, edges; 
#                      normalize_value::Number=1, binary_model::SFH.AbstractBinaryModel=SFH.NoBinaries(),
#                      imf_mean::Number=SFH.mean(imf), T_max::Number=13.7, sfr_floor::Number=1e-10,
#                      output::Union{AbstractString, Nothing}=nothing)

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

fit_sfh(config::NamedTuple) = fit_sfh(config.phot_file, config.ast_file, config.filters, config.xstrings, config.ystring, (config.xbins, config.ybins), config.MH_model0, config.disp_model0, config.Mstar, config.stellar_tracks, config.bcs, config.dmod, config.Av, config.imf, config.MH, config.logAge, config.binary_model, config.output_filename; badval=config.badval, minerr=config.minerr, maxerr=config.maxerr, plot_diagnostics=config.plot_diagnostics, output_path=config.output_path)
fit_sfh(config_file::AbstractString) = fit_sfh(parse_config(config_file))

end # module