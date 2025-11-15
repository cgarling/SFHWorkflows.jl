"""Module containing functions to plot the results of SFH fits."""
module Plotting

export plot_cmd_residuals

import StarFormationHistories as SFH
using CairoMakie
using StatsBase: fit, Histogram

parse_xcolor(xcolor) = join(xcolor, " - ")
parse_xcolor(xcolor::AbstractString) = string(xcolor)

"""
    reposition_rect(rect_obs::Observable{HyperRectangle{2, Int64}}, dx::Int, dy::Int)

Returns a new Observable{HyperRectangle{2, Int64}} that is the original rectangle
shifted by (dx, dy) relative to `rect_obs`.
"""
function reposition_bbox(rect_obs::Observable{Makie.HyperRectangle{2, Int64}}, dx::Int, dy::Int)
    lift(rect_obs) do r
        Makie.HyperRectangle(r.origin .+ (dx, dy), r.widths)
    end
end

function plot_cmd_residuals(data::Histogram, result, xcolor, yfilter::AbstractString, 
                            galaxy_name::AbstractString, output_file::AbstractString; 
                            idx::Int=1, normalize_value::Number=1)
    xcolor = parse_xcolor(xcolor)
    coeffs = SFH.calculate_coeffs(result.results[idx], result.logAge[idx], result.MH[idx])
    model_hess = sum(coeffs .* result.templates[idx] ./ normalize_value)
    # Significance = Residual / σ; sometimes called Pearson residual
    signif = (data.weights .- model_hess) ./
            sqrt.(model_hess)
    signif[data.weights .== 0] .= NaN

    xlims = extrema(data.edges[1])
    ylims = reverse(extrema(data.edges[2]))  # reverse y-axis

    # Create figure and axes
    figsize = (750, 750)
    fig = Figure(size = figsize)
    axs = [Axis(fig[i, j], xlabel = xcolor, ylabel = yfilter, xgridvisible = false, ygridvisible = false, xtickalign=1, xminortickalign=1, ytickalign=1, yminortickalign=1, xticksmirrored = true, yticksmirrored = true) for (i, j) in ((1, 1), (1, 2), (2, 1), (2, 2))]
    colgap!(fig.layout, 0)
    rowgap!(fig.layout, 0)

    tl_kws = Dict(:space => :relative, :text_align => (:left, :top), :background_color => (:white, 1.0), :text_color => :black,
                  :strokewidth => 0)
    # Panel a: Data
    hm1 = heatmap!(axs[1], data.edges[1], data.edges[2], data.weights;
                colormap = Reverse(:gray1), colorrange = (1.0, maximum(data.weights)), colorscale = log10)
    textlabel!(axs[1], 0.02, 0.82; text = "a) $galaxy_name", tl_kws...)

    # Panel b: Model
    hm2 = heatmap!(axs[2], data.edges[1], data.edges[2], model_hess;
                colormap = Reverse(:gray1), colorrange = (1.0, maximum(data.weights)), colorscale = log10)
    textlabel!(axs[2], 0.02, 0.82; text = "b) $galaxy_name Model", tl_kws...) # , align = (:left, :top)

    # Panel c: Data - Model
    hm3 = heatmap!(axs[3], data.edges[1], data.edges[2], data.weights .- model_hess;
                colormap = :seismic, colorrange = (-50, 50))
    textlabel!(axs[3], 0.02, 0.82; text = "c) Data - Model", tl_kws...) # , align = (:left, :top)

    # Panel d: Significance
    hm4 = heatmap!(axs[4], data.edges[1], data.edges[2], signif;
                colormap = :seismic, colorrange = (-15, 15))
    textlabel!(axs[4], 0.02, 0.82; text = "d) Residual\nSignificance", tl_kws...) # , align = (:left, :top)

    # Format axes
    hidexdecorations!(axs[1], ticks=false, minorticks=false)
    hidedecorations!(axs[2], ticks=false, minorticks=false)
    hideydecorations!(axs[4], ticks=false, minorticks=false)
    # Colorbars
    hm1_ticks = 0:floor(Int, log10(maximum(data.weights)))
    hm1_ticks = (exp10.(hm1_ticks), [L"10^%$i" for i in hm1_ticks])
    # Colorbar(fig, hm1, vertical = false, minorticksvisible = true, bbox=reposition_bbox(axs[1].scene.viewport, -20, 0),
    #          ticks = hm1_ticks, width=225,
    #          alignmode = Outside(10), halign = :center, valign = :top, flipaxis = false, size = 15)
    # Colorbar(fig, hm2, vertical = false, minorticksvisible = true, bbox=reposition_bbox(axs[2].scene.viewport, -20, 0), 
    #          ticks = hm1_ticks, width=225,
    #          alignmode = Outside(10), halign = :center, valign = :top, flipaxis = false, size = 15)
    # Colorbar(fig, hm3, vertical = false, minorticksvisible = true, bbox=reposition_bbox(axs[3].scene.viewport, -20, 0), width=225,
    #          alignmode = Outside(10), halign = :center, valign = :top, flipaxis = false, size = 15)
    # Colorbar(fig, hm4, vertical = false, minorticksvisible = true, bbox=reposition_bbox(axs[4].scene.viewport, -20, 0), width=225,
    #          alignmode = Outside(10), halign = :center, valign = :top, flipaxis = false, size = 15)
    hm_vec = (hm1, hm2, hm3, hm4)
    for i in eachindex(hm_vec)
        kws = Dict(:vertical => false, :minorticksvisible => true, :bbox => reposition_bbox(axs[i].scene.viewport, -20, 0), 
                   :width => 225, :alignmode => Outside(10), :halign => :center, :valign => :top, :flipaxis => false, :size => 15)
        cb = if i == 1 || i == 2
            Colorbar(fig, hm_vec[i]; ticks = hm1_ticks, kws...)
        else
            Colorbar(fig, hm_vec[i]; kws...)
        end
        cb.axis.attributes[:background_color] = :white
    end

    # Shared axis limits and ticks
    for ax in axs
        limits!(ax, xlims, ylims)
        # hidespines!(ax, :t, :r)
    end
    save(output_file, fig)
end

end # module