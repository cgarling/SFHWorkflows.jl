"""Module containing code to process artificial star test results into input appropriate for use in Systematics module."""
module ASTs

export process_ast_file, snr_model

import StarFormationHistories as SFH
# import CSV
using ArgCheck: @argcheck, @check
using DelimitedFiles: readdlm
using PDFmerger: merge_pdfs
using DataInterpolations: CubicSpline, LinearInterpolation, ExtrapolationType, ExtrapolationType.Constant
using Interpolations: interpolate, extrapolate, Gridded, Linear, Flat
using SpecialFunctions: erf
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

# Replaces NaN entries with the nearest non-NaN value along the first dimension, then along the second.
# `StarFormationHistories.process_ASTs` returns NaN for rows or columns of the grid without (detected) artificial stars.
function fill_nan(A::AbstractMatrix)
    A = copy(A)
    for v in Iterators.flatten((eachcol(A), eachrow(A)))
        good = findall(!isnan, v)
        isempty(good) && continue
        for i in eachindex(v)
            isnan(v[i]) && (v[i] = v[good[argmin(abs.(good .- i))]])
        end
    end
    return A
end

# Processes AST file and returns joint completeness, bias, and error models; each is a function of the input
# magnitudes `(m1, m2)` in the two AST filters, and bias and error return one value per filter
function process_ast_file(astfile::AbstractString, filters, badval::Number, minerr::Number, maxerr::Number, plot_diagnostics::Bool, output_path::AbstractString)
    astmags = readdlm(astfile, Float64)
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
    output1[isapprox.(badval, view(astmags, :, 3))] .= badval
    output2[isapprox.(badval, view(astmags, :, 4))] .= badval

    # The catalog requires detection in both filters, which is not separable into per-filter completeness, so the models
    # are measured on a 2-D grid of input magnitudes; sparse cells fall back to 1-D values inside process_ASTs
    bins = Tuple(range(round(minimum(x), RoundUp; digits=1), round(maximum(x), RoundDown; digits=1); step=0.3) for x in (input1, input2))
    table = Table(in1=input1, in2=input2, out1=output1, out2=output2)
    centers, C, B, E = SFH.process_ASTs(table, (:in1, :in2), (:out1, :out2), bins, r -> !isapprox(r.out1, badval) && !isapprox(r.out2, badval))
    itp(A) = extrapolate(interpolate(centers, fill_nan(A), Gridded(Linear())), Flat())
    completeness = itp(C)
    b = itp.(B)
    e = itp.(map(x -> clamp.(x, minerr, maxerr), E))

    if plot_diagnostics
        # Per-filter (1-D) models, used only to summarize the ASTs in the diagnostic plots
        r1 = process_asts(input1, output1, badval, minerr, maxerr)
        r2 = process_asts(input2, output2, badval, minerr, maxerr)
        # errmax1 = abs(maximum(input1 .- output1))
        # errlim1 = (max(-1.5*maxerr, -errmax1), 
        #            min(1.5*maxerr, errmax1))
        # errmax2 = abs(maximum(input2 .- output2))
        # errlim2 = (max(-1.5*maxerr, -errmax2), 
        #            min(1.5*maxerr, errmax2))
        errlim = (max(-1.5*maxerr, -0.4), min(1.5*maxerr, 0.4))
        f = plot_ast_residual(input1, output1, r1, filters[1], errlim, badval)
        # display(f)
        save(joinpath(output_path, "residuals1.pdf"), f)
        f = plot_ast_residual(input2, output2, r2, filters[2], errlim, badval)
        # display(f)
        save(joinpath(output_path, "residuals2.pdf"), f)

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
        save(joinpath(output_path, "completeness.pdf"), f)

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
        save(joinpath(output_path, "error.pdf"), f)

        # Plot joint completeness used in the fit
        f = Figure()
        ax = Axis(f[1, 1], xlabel=filters[1], ylabel=filters[2], title="Joint completeness")
        hm = heatmap!(ax, centers..., C; colorrange=(0, 1))
        Colorbar(f[1, 2], hm)
        save(joinpath(output_path, "completeness2d.pdf"), f)

        # When finished, merge pdfs into one
        merge_pdfs(map(Base.Fix1(joinpath, output_path), ["residuals1.pdf", "residuals2.pdf", "error.pdf", "completeness.pdf", "completeness2d.pdf"]), joinpath(output_path, "diagnostics.pdf"); cleanup=true)
    end
    return (completeness = completeness, bias = (m1, m2) -> (b[1](m1, m2), b[2](m1, m2)), err = (m1, m2) -> (e[1](m1, m2), e[2](m1, m2)))
end

"""
    (completeness, bias, err) = snr_model(mag, snr; bias=nothing, minerr=0, snr50=5, width=1)
Returns completeness, bias, and photometric error functions of apparent magnitude for a filter whose signal-to-noise
ratio `snr` is tabulated at apparent magnitudes `mag`. `log10(snr)` is interpolated linearly in magnitude and
extrapolated linearly beyond the table, so the SNR keeps falling past the faint end rather than flattening.
The error is `max(minerr, StarFormationHistories.magerr_snr(SNR))`. The bias is zero unless `bias` gives its value at
each of `mag`, in which case it is interpolated linearly and held constant beyond the table. The completeness is
`0.5 * (1 + erf((SNR - snr50) / (sqrt(2) * width)))`.
"""
function snr_model(mag, snr; bias=nothing, minerr::Number=0, snr50::Number=5, width::Number=1)
    @argcheck length(mag) == length(snr) >= 2
    @argcheck isnothing(bias) || length(bias) == length(mag) "`bias` must have one value per entry of `mag`."
    @argcheck all(>(0), snr) "All tabulated SNR values must be positive."
    @argcheck width > 0
    p = sortperm(mag)
    # Linear extrapolation past the faint end must keep SNR falling, or completeness would rise again
    @argcheck issorted(snr[p]; rev=true) "Tabulated SNR must not increase toward fainter magnitudes."
    @argcheck snr[p[end]] < snr[p[end-1]] "Tabulated SNR must decrease between the two faintest magnitudes."
    logsnr = LinearInterpolation(log10.(snr[p]), mag[p]; extrapolation=ExtrapolationType.Linear)
    completeness(m) = (1 + erf((exp10(logsnr(m)) - snr50) / (sqrt(2) * width))) / 2
    err(m) = max(minerr, SFH.magerr_snr(exp10(logsnr(m))))
    biasfunc = isnothing(bias) ? zero : LinearInterpolation(bias[p], mag[p]; extrapolation=Constant)
    return (completeness = completeness, bias = biasfunc, err = err)
end

end # module
