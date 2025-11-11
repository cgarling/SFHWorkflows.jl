module SFHWorkflows

export fit_sfh

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

# Top-level functions
function fit_sfh(obsfile::AbstractString, astfile::AbstractString, filters; 
                 badval::Number=99.999, minerr::Number=0.0, maxerr::Number=0.2, plot_diagnostics::Bool=true, output_path::AbstractString=".") # filters=("mag1", "mag2")
    filters = string.(filters)
    completeness, bias, err = process_ast_file(astfile, filters, badval, minerr, maxerr, plot_diagnostics, output_path)
end

fit_sfh(config::NamedTuple) = fit_sfh(config.phot_file, config.ast_file, config.filters; badval=config.badval, minerr=config.minerr, maxerr=config.maxerr, plot_diagnostics=config.plot_diagnostics, output_path=config.output_path)
fit_sfh(config_file::AbstractString) = fit_sfh(parse_config(config_file))

end # module