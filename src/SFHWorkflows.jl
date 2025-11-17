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
        Scatter = (strokecolor=:black, strokewidth=1),
        Lines   = (linewidth=3,))
                # xminorgridvisible=true, 
                # yminorgridvisible=true))

function Makie.convert_arguments(P::Makie.CellGrid, h::Histogram)
    return Makie.convert_arguments(P, h.edges[1], h.edges[2], h.weights)
end
Makie.plottype(::Histogram) = Makie.Heatmap

include(joinpath("SFH", "SFHFitting.jl"))
using .SFHFitting

end # module