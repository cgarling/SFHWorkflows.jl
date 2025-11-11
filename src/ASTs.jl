module ASTs

export process_ast_file

import StarFormationHistories as SFH
# import CSV
using ArgCheck: @argcheck, @check
using DelimitedFiles: readdlm
using PDFmerger: merge_pdfs
using DataInterpolations: CubicSpline, ExtrapolationType.Constant
using Logging: ConsoleLogger, with_logger, Error
using StatsBase: fit, Histogram
using CairoMakie
using PlotUtils: zscale

# Try both TypedTables.Table and DataFrames.DataFrame; works
using TypedTables: Table
# using DataFrames: DataFrame

function process_asts(input, output, badval::Number, minerr::Number, maxerr::Number)
    ast_table = Table(input = input, output = output)
    ast_bin_edges = range(round(minimum(input), RoundUp; digits=1),
                          round(maximum(input), RoundDown; digits=1); step=0.3)
    # Don't care about warnings issued by process_ASTs, but do want to see
    # any errors
    with_logger(ConsoleLogger(Error)) do
        r = SFH.process_ASTs(ast_table, :input, :output, ast_bin_edges, 
                             x -> !isapprox(x.output, badval))
    end
    # (bin_centers1, completeness1, bias1, error1) = r
    # Filter out NaN bins
    good = .!isnan.(r[3])
    r = [rr[good] for rr in r]
    # Find when/if maxerr is exceeded and truncate results
    err = r[4]
    min_idx = findmin(err)[2] # Find minimum of error and start from there
    trunc_idx = lastindex(err)
    for i in eachindex(err)[min_idx:end]
        if err[i] >= maxerr
            trunc_idx = i - 1
            break
        end
    end
    complete_itp = CubicSpline(r[2], r[1]; extrapolation=Constant)
    bias_itp = CubicSpline(r[3], r[1]; extrapolation=Constant)
    err_itp = CubicSpline(max.(minerr, r[4][begin:trunc_idx]), r[1][begin:trunc_idx]; extrapolation=Constant)
    return (complete_itp, bias_itp, err_itp)
end

function plot_ast_residual(input, output, itps, mag_label, err_range, badval::Number)
        # Plot residuals, bias
        f = Figure()
        ax = Axis(f[1, 1], xlabel="Magnitude", ylabel=L"$\langle$ Output - Input $\rangle$", title=mag_label * " residuals")

        good = findall(!≈(badval), output)
        h = fit(Histogram, (input[good], output[good] .- input[good]),
                (range(first(itps[1].t), last(itps[1].t); length=100),
                 range(err_range...; length=100)))
        # Colorscale is now in Makie, see PR https://github.com/MakieOrg/Makie.jl/pull/5166
        # Added in version 0.24.5
        transform = Makie.LuptonAsinhScale(0.1, 0.01, 1)
        p = heatmap!(ax, h.edges[1], h.edges[2], h.weights; interpolate=false, colormap=:cividis,
                     # colorrange = zscale(h.weights; contrast=0.02, k_rej=2),
                     colorrange = zscale(transform.(h.weights)), # ; contrast=0.02, k_rej=2.5),
                     # colorscale = Makie.ReversibleScale(x -> asinh(x / 2) / log(10), x -> 2sinh(log(10) * x)))
                     colorscale = transform)
        hlines!(ax, [0], color=:black)
        # Overplot bias
        # scatter!(ax, itps[2].t, itps[2].u, color=:red)
        errorbars!(ax, itps[2].t, itps[2].u, itps[3].(itps[2].t); color=:red)
        lines!(ax, itps[2].t, itps[2].(itps[2].t), color=:red, linestyle=:solid)
        xlims!(ax, extrema(h.edges[1])...)
        ylims!(ax, extrema(h.edges[2])...)
        return f
end

# Processes AST file and returns completeness, error, bias results
function process_ast_file(astfile::AbstractString, filters, badval::Number, minerr::Number, maxerr::Number, plot_diagnostics::Bool)
    astmags = readdlm(astfile, ' ', Float64)
    # @check iseven(size(astmags, 2))
    # nfilters = size(astmags, 2) ÷ 2
    # @check nfilters == length(filters) "Mismatch between `length(filters)` and number of columns in AST file $astfile."
    @check size(astmags, 2) == 4 "ast_file must have 4 columns; `(inputmag1, inputmag2, (outmag1 - inputmag1), (outmag2 - inputmag2)`."

    @info "Processing ASTs"
    # MATCH convention for the AST file is
    # (inputmag1, inputmag2, (outmag1 - inputmag1), (outmag2 - inputmag2)
    input1, input2 = view(astmags, :, 1), view(astmags, :, 2)
    output1 = view(astmags, :, 3) .+ input1
    output2 = view(astmags, :, 4) .+ input2
    # Add badval's back into output
    bad1 = findall(≈(badval), view(astmags, :, 3))
    # output1[bad1] .= badval 
    bad2 = findall(≈(badval), view(astmags, :, 4))
    # output2[bad2] .= badval
    # Require detection in both bands
    bad = bad1 .| bad2
    output1[bad] .= badval
    output2[bad] .= badval

    r1 = process_asts(input1, output1, badval, minerr, maxerr)
    r2 = process_asts(input2, output2, badval, minerr, maxerr)

    if plot_diagnostics
        # errmax1 = abs(maximum(input1 .- output1))
        # errlim1 = (max(-1.5*maxerr, -errmax1), 
        #            min(1.5*maxerr, errmax1))
        # errmax2 = abs(maximum(input2 .- output2))
        # errlim2 = (max(-1.5*maxerr, -errmax2), 
        #            min(1.5*maxerr, errmax2))
        errlim = (max(-1.5*maxerr, -0.4), min(1.5*maxerr, 0.4))
        f = plot_ast_residual(input1, output1, r1, filters[1], errlim, badval)
        # display(f)
        save("residuals1.pdf", f)
        f = plot_ast_residual(input2, output2, r2, filters[2], errlim, badval)
        # display(f)
        save("residuals2.pdf", f)

        # Plot completeness
        f = Figure()
        ax = Axis(f[1, 1], xlabel="Magnitude", ylabel="Completeness")
        # Plot mag1
        scatter!(ax, r1[1].t, r1[1].u, color=:black, label=filters[1])
        lines!(ax, r1[1].t, r1[1].(r1[1].t), color=:black, linestyle=:dash, label=filters[1])
        # Plot mag2
        scatter!(ax, r2[1].t, r2[1].u, color=:red, label=filters[2])
        lines!(ax, r2[1].t, r2[1].(r2[1].t), color=:red, linestyle=:dash, label=filters[2])
        axislegend(ax, merge=true, unique=true)
        # display(f)
        save("completeness.pdf", f)

        # # Plot bias
        # f = Figure()
        # ax = Axis(f[1, 1], xlabel="Magnitude", ylabel=L"Bias $\langle m_\text{out} - m_\text{in} \rangle$")
        # ylims!(ax, -0.1, 0.1)
        # # Plot mag1
        # scatter!(ax, r1[2].t, r1[2].u, color=:black, label=filters[1])
        # lines!(ax, r1[2].t, r1[2].(r1[2].t), color=:black, linestyle=:dash, label=filters[1])
        # # Plot mag2
        # scatter!(ax, r2[2].t, r2[2].u, color=:red, label=filters[2])
        # lines!(ax, r2[2].t, r2[2].(r2[2].t), color=:red, linestyle=:dash, label=filters[2])
        # axislegend(ax, merge=true, unique=true)
        # display(f)
        # save("bias.pdf", f)

        # Plot error
        f = Figure()
        ax = Axis(f[1, 1], xlabel="Magnitude", ylabel="Median Photometric Error")
        # Plot mag1
        scatter!(ax, r1[3].t, r1[3].u, color=:black, label=filters[1])
        lines!(ax, r1[3].t, r1[3].(r1[3].t), color=:black, linestyle=:dash, label=filters[1])
        # Plot mag2
        scatter!(ax, r2[3].t, r2[3].u, color=:red, label=filters[2])
        lines!(ax, r2[3].t, r2[3].(r2[3].t), color=:red, linestyle=:dash, label=filters[2])
        ylims!(0.0, maxerr*1.5)
        axislegend(ax, merge=true, unique=true, position=:lt)
        # display(f)
        save("error.pdf", f)

        # When finished, merge pdfs into one
        merge_pdfs(["residuals1.pdf", "residuals2.pdf", "error.pdf", "completeness.pdf"], "diagnostics.pdf"; cleanup=true)
    end
    completeness = [r1[1], r2[1]]
    bias = [r1[2], r2[2]]
    err = [r1[3], r2[3]]
    return (completeness = completeness, bias = bias, err = err)
end

end # module
