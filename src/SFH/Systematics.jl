"""Module containing code that will run SFH fitting given prepared inputs."""
module Systematics

export systematics, mag_select

import StarFormationHistories as SFH
import CSV
using TypedTables: Table, columnnames, getproperties
using BolometricCorrections: chemistry, MH, Z, filternames, AbstractBCGrid, AbstractBCTable, MISTBCGrid, gridname
using StellarTracks: AbstractTrackLibrary, isochrone
using ArgCheck: @argcheck
using Printf: @printf, @sprintf, Format, format
using LinearAlgebra: BLAS
using StatsBase: quantile
# using Interpolations: interpolate, Gridded, Linear #, Throw, extrapolate, Flat
# using Roots: find_zero # For taus

"""
    convert_MH(mh, tracks, bcs)
Given an [M/H] value `mh` defined for the chemistry of the track object `tracks` (which can be any type from StellarTracks.jl that supports the `chemistry` API), returns the equivalent [M/H] value for the chemistry of the bolometric correction object `bcs` (which can be any type from BolometricCorrections.jl that supports the `chemistry` API).
"""
convert_MH(mh, tracks, bcs) = MH(chemistry(bcs), Z(chemistry(tracks), mh))

"""
    (xidx, yidxs) = mag_select(bcgrid::AbstractBCGrid, ystring::Union{AbstractString, Symbol}, xstrings)
Returns the indices into `filternames(bcgrid)` that correspond to the provided `ystring` and `xstrings` arguments. 
Used to programmatically determine the proper columns to select from an `AbstractBCGrid` using user-specific filter names.
```jldoctest
julia> using BolometricCorrections: MISTBCGrid, filternames

julia> m = MISTBCGrid("hst_acs_wfc")

julia> mag_select(m, "F150W", (:F090W, :F150W))
(6, [2, 6])
```
"""
function mag_select(bcgrid::Union{AbstractBCGrid, AbstractBCTable}, ystring::Union{AbstractString, Symbol}, xstrings)
    @argcheck length(xstrings) == 2
    Base.require_one_based_indexing(xstrings)
    # ystring = lowercase(string(ystring))
    # xstrings = [lowercase(string(xs)) for xs in xstrings]
    # filters = [lowercase(string(f)) for f in filternames(bcgrid)]
    ystring = string(ystring)
    xstrings = [string(xs) for xs in xstrings]
    filters = [string(f) for f in filternames(bcgrid)]
    ymatches = findfirst(==(ystring), filters)
    if !isnothing(ymatches) # Look for exact match first
        yidx = ymatches
    else # Look for partial match second
        ymatches = findall(occursin.(ystring, filters))
        if length(ymatches) == 0
            throw(ArgumentError("No filter in $bcgrid contains the provided `ystring` $ystring. The full filter set is $filters."))
        elseif length(ymatches) > 1
            throw(ArgumentError("More than one filter in $bcgrid matched with provided `ystring` $ystring; a more specific `ystring` is required. The list of matching filters is $(filters[ymatches])."))
        else
            yidx = ymatches[1]
        end
    end

    xidxs = Vector{Int}(undef, 2)
    for i in eachindex(xstrings, xidxs)
        xs = xstrings[i]
        xmatches = findfirst(==(xs), filters)
        if !isnothing(xmatches) # Look for exact match first
            xidxs[i] = xmatches
        else # Look for partial match second
            xmatches = findall(occursin.(xs, filters))
            if length(xmatches) == 0
                throw(ArgumentError("No filter in $bcgrid contains the provided filter $xs. The full filter set is $filters."))
            elseif length(xmatches) > 1
                throw(ArgumentError("More than one filter in $bcgrid matched with provided filter $xs; a more specific filter name is required. The list of matching filters is $(filters[xmatches])."))
            else
                xidxs[i] = xmatches[1]
            end
        end      
    end
    # return yidx, xidxs
    return Symbol(filters[yidx]), (Symbol(filters[xidxs[1]]), Symbol(filters[xidxs[2]]))
end

function templates(tracklib::AbstractTrackLibrary, bclib::AbstractBCGrid,
                   xstrings, ystring,
                   dmod, Av, err_funcs, complete_funcs, bias_funcs, imf,
                   unique_MH, unique_logAge, edges;
                   normalize_value::Number=1, binary_model::SFH.AbstractBinaryModel = SFH.NoBinaries(),
                   imf_mean::Number = SFH.mean(imf))

    # Figure out the filters we want to use
    ysymb, xsymbs = mag_select(bclib, ystring, xstrings)
    if ysymb in xsymbs
        iso_symb = xsymbs # All the symbols to pull from the isochrone
        yidx = findfirst(==(ysymb), iso_symb) # index into iso_symb of y filter
        xidxs = eachindex(xsymbs) # index into iso_symb of the x filters
    else
        iso_symb = (ysymb, xsymbs...)
        yidx = 1 # index into iso_symb of y filter
        xidxs = [2, 3] # index into iso_symb of the x filters
    end

    # Convert unique_MH (defined for the tracklib chemistry) to the bclib chemistry
    unique_bc_MH = convert_MH.(unique_MH, Ref(tracklib), Ref(bclib))

    templates = Vector{Matrix{Float64}}(undef, length(unique_logAge) * length(unique_MH))
    template_logAge = Vector{Float64}(undef, length(templates))
    template_MH = similar(template_logAge)
    # @threads for (i, (mh, logage)) in collect(enumerate(Iterators.product(unique_MH, unique_logAge))) # issorted(mdf_template_logAge) == true
    Threads.@threads for i in eachindex(unique_MH)
    # for i in eachindex(unique_MH)
        mh = unique_MH[i]
        # Deal with MH outside range of tracklib
        if (mh < minimum(MH(tracklib)))
            tracklib_mh = minimum(MH(tracklib))
        elseif (mh > maximum(MH(tracklib)))
            tracklib_mh = maximum(MH(tracklib))
        else
            tracklib_mh = mh
        end
        # # Deal with MH outside range of bclib
        # if (mh < minimum(MH(bclib)))
        #     bclib_mh = minimum(MH(bclib))
        # elseif (mh > maximum(MH(bclib)))
        #     bclib_mh = maximum(MH(bclib))
        # else
        #     bclib_mh = mh
        # end
        # bct = bclib(bclib_mh) 
        # Interpolate bclib at correct MH, Av
        bct = bclib(unique_bc_MH[i], Av)

        Threads.@threads for j in eachindex(unique_logAge)
        # for j in eachindex(unique_logAge)
            logage = unique_logAge[j]
            ind = j + ((i-1) * length(unique_logAge)) # index into templates and other buffers for (i,j)
            iso = isochrone(tracklib, bct, logage, tracklib_mh)
            iso_mags = [getproperty(iso, k) for k in iso_symb] # (xsymb, ysymbs...)]
            m_ini = iso.m_ini
            templates[ind] = SFH.partial_cmd_smooth(m_ini, iso_mags, err_funcs, yidx, xidxs, imf, 
                                                    complete_funcs, bias_funcs; 
                                                    dmod=dmod, normalize_value=normalize_value, edges=edges, 
                                                    mean_mass=imf_mean, binary_model=binary_model).weights
            template_logAge[ind] = logage
            template_MH[ind] = mh
        end
    end
    # Sort the template_logAge and template_MH so we have a guarantee for later
    permidx = sortperm(template_logAge)
    templates = templates[permidx]
    template_logAge = template_logAge[permidx]
    template_MH = template_MH[permidx]
    # @printf "Number of templates constructed: %i" length(templates)
    return (templates = templates, logAge = template_logAge, MH = template_MH)
end

# If `fname` has a leading directory and it does not exist, create it
function _mkpath(fname::AbstractString)
    dir = dirname(fname)
    if !isempty(dir) && !isdir(dir)
        mkpath(dir)
    end
end
# Set up for writing output
function write_systable(fname::AbstractString, table)
    _mkpath(fname) # Ensure the directory exists
    # Define transformations for CSV formatting
    transform(col, val) = val
    function transform(col, val::Number)
        if (col == 1) | (col == 2)
            return @sprintf("%.3f", val)
        # SFR parameters are in columns with # multiples of 7, 8, 9
        # and can be low, so format as exponential for better precision
        elseif any(==(0), map(Base.Fix1(rem, col), (7, 8, 9)))
            return @sprintf("%.7e", val)
        else 
            return @sprintf("%.7f", val)
        end
    end
    CSV.write(fname, table; delim=' ', transform=transform)
end
write_systable(::Nothing, ::Any) = nothing

function write_masstable(fname::AbstractString, table)
    _mkpath(fname) # Ensure the directory exists
    # Define transformations for CSV formatting
    transform(col, val) = val
    transform(col, val::Number) = @sprintf("%.7e", val)
    CSV.write(fname, table; delim=' ', transform=transform)
end
write_masstable(::Nothing, ::Any) = nothing

function fit_sfh(MH_model0::SFH.AbstractMetallicityModel,
                 disp_model0::SFH.AbstractDispersionModel,
                 mstar::Number, # Estimate of stellar mass of galaxy
                 data::Union{AbstractVector{<:Number}, AbstractMatrix{<:Number}},
                 tracklib::AbstractTrackLibrary, bclib::AbstractBCGrid,
                 xstrings, ystring, # strings or symbols that allow us to select the correct filters from isochrone
                 dmod, Av, err_funcs, complete_funcs, bias_funcs, imf,
                 unique_MH, unique_logAge, edges; 
                 normalize_value::Number=1, binary_model::SFH.AbstractBinaryModel=SFH.NoBinaries(),
                 imf_mean::Number=SFH.mean(imf), T_max::Number=13.7)

    @argcheck mstar > 0
    # Construct templates
    all_templates = templates(tracklib, bclib, xstrings, ystring, dmod, Av, err_funcs, complete_funcs, bias_funcs,
                              imf, unique_MH, unique_logAge, edges;
                              normalize_value=normalize_value, binary_model=binary_model, imf_mean=imf_mean)
    result = SFH.fit_sfh(MH_model0, disp_model0, SFH.stack_models(all_templates.templates), vec(data), all_templates.logAge, all_templates.MH;
                         x0=SFH.construct_x0_mdf(all_templates.logAge, T_max; normalize_value=mstar / normalize_value))
    return merge((result=result,), all_templates)
end

function systematics(MH_model0::SFH.AbstractMetallicityModel,
                     disp_model0::SFH.AbstractDispersionModel,
                     mstar::Number, # Estimate of stellar mass of galaxy
                     data::Union{AbstractVector{<:Number}, AbstractMatrix{<:Number}},
                     tracklibs, bclibs,
                     xstrings, ystring, # strings or symbols that allow us to select the correct filters from isochrone
                     dmod, Av, err_funcs, complete_funcs, bias_funcs, imf,
                     unique_MH, unique_logAge, edges; 
                     normalize_value::Number=1, binary_model::SFH.AbstractBinaryModel=SFH.NoBinaries(),
                     imf_mean::Number=SFH.mean(imf), T_max::Number=13.7, sfr_floor::Number=1e-10,
                     output::Union{AbstractString, Nothing}=nothing)

    @argcheck length(tracklibs) >= 1
    @argcheck length(bclibs) >= 1
    data = vec(data)
    logAge_l = sort(unique_logAge)
    logAge_u = vcat(logAge_l[begin+1:end], log10(T_max) + 9)

    blas_threads = BLAS.get_num_threads() # Log number of BLAS threads so we can revert back to this at the end
    BLAS.set_num_threads(1) # Set BLAS threads to 1 for more efficient solves

    try # try-finally to make sure BLAS threads revert after solves
        println(BLAS.get_num_threads())
        nsolutions = length(tracklibs) * length(bclibs)
        results = Vector{SFH.CompositeBFGSResult}(undef, nsolutions)
        templates = Vector{Vector{Matrix{Float64}}}(undef, nsolutions)
        logAge = Vector{Vector{Float64}}(undef, nsolutions) # These should be the same for every solution
        MH = Vector{Vector{Float64}}(undef, nsolutions)     # These should be the same for every solution
        birth_masses = Matrix{Float64}(undef, nsolutions, 3)
        # present_masses = Matrix{Float64}(undef, nsolutions, 3)
        tables = Vector{Table}(undef, nsolutions)

        @info "Entering threaded SFH loop"
        Threads.@threads for i in eachindex(tracklibs)
            tracklib = tracklibs[i]
            Threads.@threads for j in eachindex(bclibs)
                bclib = bclibs[j]
                ind = j + ((i-1) * length(bclibs)) # index into result for (i,j)
                fit_result = fit_sfh(MH_model0, disp_model0, mstar, data, tracklib, bclib, xstrings,
                                    ystring, dmod, Av, err_funcs,
                                    complete_funcs, bias_funcs, imf, unique_MH, unique_logAge, edges; normalize_value=normalize_value,
                                    binary_model=binary_model, imf_mean=imf_mean, T_max=T_max)
                results[ind] = fit_result[1]
                templates[ind] = fit_result.templates
                logAge[ind] = fit_result.logAge
                MH[ind] = fit_result.MH

                coeffs = SFH.calculate_coeffs(fit_result[1], fit_result.logAge, fit_result.MH) .* normalize_value
                ul, cum_sfh, sfr, mean_MH = SFH.calculate_cum_sfr(coeffs, fit_result.logAge, fit_result.MH, T_max; sorted=true)
                
                # Calculate quantiles for uncertainties
                quantile_results_map = SFH.cum_sfr_quantiles(fit_result[1].map, fit_result.logAge, fit_result.MH, T_max, 10_000, (0.16, 0.5, 0.84))
                # quantile_results_map = try
                #     SFH.cum_sfr_quantiles(fit_result[1].map, fit_result.logAge, fit_result.MH, T_max, 10_000, (0.16, 0.5, 0.84))
                # catch e
                #     println("Failed to compute MAP quantiles for tracklib=$(gridname(tracklib)), bclib=$(gridname(bclib)). Thrown error was:")
                #     showerror(stderr, e)
                # end
                quantile_results_map[2] .*= normalize_value # Correct SFR normalization
                quantile_results = SFH.cum_sfr_quantiles(fit_result[1], fit_result.logAge, fit_result.MH, T_max, 10_000, (0.16, 0.5, 0.84))
                # quantile_results = try
                #     SFH.cum_sfr_quantiles(fit_result[1], fit_result.logAge, fit_result.MH, T_max, 10_000, (0.16, 0.5, 0.84))
                # catch e
                #     println("Failed to compute MLE quantiles for tracklib=$(gridname(tracklib)), bclib=$(gridname(bclib)) -- using MAP quantiles instead. Thrown error was:")
                #     showerror(stderr, e)
                #     quantile_results_map
                # end
                quantile_results[2] .*= normalize_value # Correct SFR normalization
                # For SFRs that are ~0, make them =0, and use the MAP uncertainty estimate
                # for the upper uncertainty limit.
                for k in eachindex(sfr)
                    if sfr[k] <= sfr_floor
                        sfr[k] = 0
                        quantile_results[2][k,begin] = 0
                        quantile_results[2][k,end] = quantile_results_map[2][k,end]
                    end
                end

                # Write birth masses into output array
                birth_masses[ind, :] .= (SFH.integrate_sfr(logAge_l, logAge_u, quantile_results[2][:, i]) for i in 1:3)

                # The complicated @eval is necessary because TypedTables.Table cannot be constructed
                # from separate data and column names
                # tables[ind] = @eval Table($(Symbol(:sfr_lower, ind))=$(quantile_results[2][:,begin]),
                #                           $(Symbol(:sfr, ind))=$sfr,
                #                           $(Symbol(:sfr_upper, ind))=$(quantile_results[2][:,end]),
                #                           $(Symbol(:cum_sfh_lower, ind))=$(quantile_results[1][:,begin]),
                #                           $(Symbol(:cum_sfh, ind))=$cum_sfh,
                #                           $(Symbol(:cum_sfh_upper, ind))=$(quantile_results[1][:,end]),
                #                           $(Symbol(:MH_lower, ind))=$(quantile_results[3][:,begin]),
                #                           $(Symbol(:MH, ind))=$mean_MH,
                #                           $(Symbol(:MH_upper, ind))=$(quantile_results[3][:,end]))
                
                # colnames = Tuple(Symbol(x, ind) for x in (:sfr_lower, :sfr, :sfr_upper, :cum_sfh_lower, :cum_sfh, :cum_sfh_upper, :MH_lower, :MH, :MH_upper))
                colnames = Tuple(Symbol(x * "_" * gridname(tracklib) * "_" * gridname(bclib)) for x in ("sfr_lower", "sfr", "sfr_upper", "cum_sfh_lower", "cum_sfh", "cum_sfh_upper", "MH_lower", "MH", "MH_upper"))
                datacolumns = (quantile_results[2][:,begin], sfr, quantile_results[2][:,end], quantile_results[1][:,begin], cum_sfh, quantile_results[1][:,end], quantile_results[3][:,begin], mean_MH, quantile_results[3][:,end])
                tables[ind] = Table(NamedTuple{colnames}(datacolumns))

                # Try Tables.MatrixTable for more efficient storage
                # doesn't seem any more efficient
                # tables[ind] = Tables.table([quantile_results[2][:,begin] sfr quantile_results[2][:,end] quantile_results[1][:,begin] cum_sfh quantile_results[1][:,end] quantile_results[3][:,begin] mean_MH quantile_results[3][:,end]]; header=[:sfr_lower, :sfr, :sfr_upper, :cum_sfh_lower, :cum_sfh, :cum_sfh_upper, :MH_lower, :MH, :MH_upper])

            end
        end
        @info "SFH fits complete; measuring statistics"
        masstable = Table(name = vec([gridname(i) * "_" * gridname(j) for i=tracklibs, j=bclibs]),
                        mstar_lower = birth_masses[:,1], mstar = birth_masses[:,2], mstar_upper = birth_masses[:,3])
        write_masstable(splitext(output)[1]*"_mass"*splitext(output)[2], masstable)

        # Derive systematic uncertainty on cum_sfh, mean_MH by simply taking the extrema
        # of the results for each combination of inputs
        rtable = Table(tables...) # concatenate tables from results
        cnames = columnnames(rtable)
        cum_sfh_lower_sys = map(minimum,
                                getproperties(rtable, cnames[findall(Base.Fix1(occursin, "cum_sfh_lower"), string.(cnames))]))
        cum_sfh_upper_sys = map(maximum,
                                getproperties(rtable, cnames[findall(Base.Fix1(occursin, "cum_sfh_upper"), string.(cnames))]))
        mean_MH_lower_sys = map(minimum,
                                getproperties(rtable, cnames[findall(Base.Fix1(occursin, "MH_lower"), string.(cnames))]))
        mean_MH_upper_sys = map(maximum,
                                getproperties(rtable, cnames[findall(Base.Fix1(occursin, "MH_upper"), string.(cnames))]))

        final_table = Table(Table(logAge_lower=logAge_l, 
                                logAge_upper=logAge_u,
                                cum_sfh_lower_sys=cum_sfh_lower_sys,
                                cum_sfh_upper_sys=cum_sfh_upper_sys,
                                MH_lower_sys=mean_MH_lower_sys,
                                MH_upper_sys=mean_MH_upper_sys),
                            rtable)
        # If output::AbstractString, will try to write final_table to output
        write_systable(output, final_table)
        return (results=results, templates=templates, logAge=logAge, MH=MH, 
                table=final_table)
    finally
        BLAS.set_num_threads(blas_threads)
        println(BLAS.get_num_threads())
    end
end

end # module