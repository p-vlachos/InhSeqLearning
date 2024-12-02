# using Pkg
# Pkg.activate(".")
# using InhSequences
# using Statistics
# using HDF5

using ColorSchemes
using LaTeXStrings
using CairoMakie
# using CurveFit
using Colors


function plotNetworkActivity(times::Matrix{Float64}, popmembers::Matrix{Int64}, ipopmembers::Matrix{Int64}; interval::AbstractVector, name::AbstractString, output_dir::AbstractString="./output_analysis/")
    Ne::Int64 = 4000
    Ni2::Int64 = 250
    Ncells::Int64 = size(times)[1]
    Npop::Int64 = size(popmembers, 2)
    Nmembers_max::Int64 = size(popmembers, 1)
    Ni_members::Int64 = 27
    plot_interval::Vector{Int64} = collect((minimum(interval)):1000:maximum(interval))

    labelsize::Int64 = 24
    ticklabelsize::Int64 = 22
    linewidth::Float64 = 2.
    markersize::Float64 = 1.

    restcells::Vector{Int64} = deleteat!(map(*, ones(Int, 4000), range(1,stop=4000)), sort(unique(popmembers))[2:end])
    ylim_max::Int64 = count(i->(i>0), popmembers) + length(restcells) + (Ncells-Ni2-Ne) + length(ipopmembers)

    ns::Vector{Int64} = zeros(Ncells)
    for cc = 1:Ncells
        ns[cc] = count(i->i>0, times[cc, :])
    end

    # Plot params
    ytick_seq::Vector{Int64} = zeros(10)
    rowcount::Int64 = 1
    ytickcount::Int64 = 1

    fig = Figure(resolution=(1080, 720))
    ax = CairoMakie.Axis(fig[1, 1], xlabel=L"\text{simulation time (s)}", ylabel=L"\text{sequences (neurons)}", xlabelsize=labelsize, ylabelsize=labelsize,
                xticks=(plot_interval, [L"%$(x)" for x = 1:length(plot_interval)]), xticklabelsize=ticklabelsize, xgridvisible=false,
                yticklabelsize=ticklabelsize, ygridvisible=false,
                limits=(minimum(interval), maximum(interval), 1, ylim_max))
    # Excitatoy assembly members
    for pp = 1:Npop
        for cc = 1:Nmembers_max
            (popmembers[cc, pp] < 1) && (break)
            indx = round(Int, popmembers[cc, pp])
            x_vals = times[indx, 1:round(Int, ns[indx])]
            y_vals = rowcount * ones(length(x_vals))
            CairoMakie.scatter!(ax, x_vals, y_vals, color=ColorSchemes.Blues[7], markersize=markersize)
            rowcount += 1
        end
        if mod(pp, 4) == 0
            hlines!(ax, rowcount, linewidth=linewidth*.8, color=ColorSchemes.Blues[7])
        else
            hlines!(ax, rowcount, linewidth=linewidth*.3, color=ColorSchemes.Blues[7])
        end
        if mod(pp+2, 4) == 0
            ytick_seq[ytickcount] = rowcount
            ytickcount += 1
        end
    end
    hlines!(ax, rowcount, linewidth=linewidth*.5, color=:gray20)
    # Excitatory non-members 
    for cc in restcells
        x_vals = times[cc, 1:round(Int, ns[cc])]
        y_vals = rowcount * ones(length(x_vals))
        CairoMakie.scatter!(ax, x_vals, y_vals, color=ColorSchemes.Greys[7], markersize=markersize)
        rowcount += 1
    end
    hlines!(ax, rowcount, linewidth=linewidth*.5, color=:gray20)
    # Inhibitory I₁
    for cc = (Ne+1):(Ncells-Ni2+1)
        x_vals = times[cc, 1:round(Int, ns[cc])]
        y_vals = rowcount * ones(length(x_vals))
        CairoMakie.scatter!(ax, x_vals, y_vals, color=ColorSchemes.Reds[9], markersize=markersize)
        rowcount += 1
    end
    hlines!(ax, rowcount, linewidth=linewidth*.5, color=:gray20)
    # I₂ assembly members
    for pp = 1:Npop
        for cc = 1:Ni_members
            indx = round(Int, popmembers[cc, pp])
            x_vals = times[indx, 1:round(Int, ns[indx])]
            y_vals = rowcount * ones(length(x_vals))
            CairoMakie.scatter!(ax, x_vals, y_vals, color=ColorSchemes.Reds[6], markersize=markersize)
            rowcount += 1
        end
        if mod(pp, 4) == 0
            hlines!(ax, rowcount, linewidth=linewidth*.3, color=ColorSchemes.Reds[9])
        end
        if mod(pp+2, 4) == 0
            ytick_seq[ytickcount] = rowcount
            ytickcount += 1
        end
    end
    ax.yticks = (ytick_seq, [L"\text{A}", L"\text{B}", L"\text{C}", L"\text{D}", L"\text{E}", L"\text{A}", L"\text{B  }", L"\text{C}", L"\text{D  }", L"\text{E}"])
    
    (!isdir(output_dir)) && (mkpath(output_dir))
    save(joinpath(output_dir, string("network_activity", name, ".png")), fig)
end

function plotWeightsEE(weightsEE::Matrix{Float64}; name::AbstractString, output_dir::AbstractString="./output_analysis/")
    cl1 = ColorScheme(range(colorant"gray5", colorant"gray80", length=100))
    cl_exc = ColorScheme(range(colorant"gray80", colorant"dodgerblue4", length=100))
    excitation_cs = vcat(get(cl1, LinRange(0, 1, 100)), get(cl_exc, LinRange(0, 1, 100)))

    labelsize::Int64 = 24
    ticklabelsize::Int64 = 22
    linewidth::Float64 = 2.
    markersize::Float64 = 1.

    fig = Figure(resolution=(920, 720))
    ax = CairoMakie.Axis(fig[1, 1], xlabel=L"\text{presynaptic }E\text{-assembly}", ylabel=L"\text{postsynaptic }E\text{-assembly}", xlabelsize=labelsize, ylabelsize=labelsize,
                xticks=(1:4:20, [L"%$(x)" for x=1:4:20]), xticklabelsize=ticklabelsize, xgridvisible=false,
                yticks=(1:4:20, [L"%$(x)" for x=1:4:20]), yticklabelsize=ticklabelsize, ygridvisible=false)
    ht = heatmap!(ax, weightsEE, colormap=excitation_cs)
    hlines!([4.5, 8.5, 12.5, 16.5], linewidth=linewidth, color=:firebrick, linestyle=:dot)
    vlines!([4.5, 8.5, 12.5, 16.5], linewidth=linewidth, color=:firebrick, linestyle=:dot)

    vlines!(4.5, ymin=0., ymax=.4, linewidth=linewidth, color=:firebrick)
    vlines!(8.5, ymin=.2, ymax=.6, linewidth=linewidth, color=:firebrick)
    vlines!(12.5, ymin=.4, ymax=.8, linewidth=linewidth, color=:firebrick)
    vlines!(16.5, ymin=.6, ymax=1., linewidth=linewidth, color=:firebrick)

    hlines!(4.5, xmin=0., xmax=.4, linewidth=linewidth, color=:firebrick)
    hlines!(8.5, xmin=.2, xmax=.6, linewidth=linewidth, color=:firebrick)
    hlines!(12.5, xmin=.4, xmax=.8, linewidth=linewidth, color=:firebrick)
    hlines!(16.5, xmin=.6, xmax=1., linewidth=linewidth, color=:firebrick)

    Colorbar(fig[1, 2], ht, label=L"\text{coupling strength (pF)}", height=Relative(1.), labelsize=labelsize, ticklabelsize=ticklabelsize)

    (!isdir(output_dir)) && (mkpath(output_dir))
    save(joinpath(output_dir, string("weights_EE", name, ".png")), fig)
end


function plotWeightsIE(weightsIE::Matrix{Float64}; name::AbstractString, output_dir::AbstractString="./output_analysis/")
    cl1 = ColorScheme(range(colorant"gray5", colorant"gray80", length=100))
    cl_inh = ColorScheme(range(colorant"gray80", colorant"firebrick", length=100))
    inhibition_cs = vcat(get(cl1, LinRange(0, 1, 100)), get(cl_inh, LinRange(0, 1, 100)))

    labelsize::Int64 = 24
    ticklabelsize::Int64 = 22
    linewidth::Float64 = 2.
    markersize::Float64 = 1.

    fig = Figure(resolution=(920, 720))
    ax = CairoMakie.Axis(fig[1, 1], xlabel=L"\text{presynaptic }E \text{-assembly}", ylabel=L"\text{postsynaptic }I_2 \text{-assembly}", xlabelsize=labelsize, ylabelsize=labelsize,
                xticks=(1:4:20, [L"%$(x)" for x=1:4:20]), xticklabelsize=ticklabelsize, xgridvisible=false,
                yticks=(1:4:20, [L"%$(x)" for x=1:4:20]), yticklabelsize=ticklabelsize, ygridvisible=false)
    ht = heatmap!(ax, weightsIE, colormap=inhibition_cs)
    hlines!([4.5, 8.5, 12.5, 16.5], linewidth=linewidth, color=:dodgerblue4, linestyle=:dot)
    vlines!([4.5, 8.5, 12.5, 16.5], linewidth=linewidth, color=:dodgerblue4, linestyle=:dot)

    vlines!(4.5, ymin=0., ymax=.4, linewidth=linewidth, color=:dodgerblue4)
    vlines!(8.5, ymin=.2, ymax=.6, linewidth=linewidth, color=:dodgerblue4)
    vlines!(12.5, ymin=.4, ymax=.8, linewidth=linewidth, color=:dodgerblue4)
    vlines!(16.5, ymin=.6, ymax=1., linewidth=linewidth, color=:dodgerblue4)

    hlines!(4.5, xmin=0., xmax=.4, linewidth=linewidth, color=:dodgerblue4)
    hlines!(8.5, xmin=.2, xmax=.6, linewidth=linewidth, color=:dodgerblue4)
    hlines!(12.5, xmin=.4, xmax=.8, linewidth=linewidth, color=:dodgerblue4)
    hlines!(16.5, xmin=.6, xmax=1., linewidth=linewidth, color=:dodgerblue4)

    Colorbar(fig[1, 2], ht, label=L"\text{coupling strength (pF)}", height=Relative(1.), labelsize=labelsize, ticklabelsize=ticklabelsize)

    (!isdir(output_dir)) && (mkpath(output_dir))
    save(joinpath(output_dir, string("weights_IE", name, ".png")), fig)
end


function plotWeightsEI(weightsEI::Matrix{Float64}; name::AbstractString, output_dir::AbstractString="./output_analysis/")
    cl1 = ColorScheme(range(colorant"gray5", colorant"gray80", length=100))
    cl_inh = ColorScheme(range(colorant"gray80", colorant"dodgerblue2", length=100))
    inhibition_cs = vcat(get(cl1, LinRange(0, 1, 100)), get(cl_inh, LinRange(0, 1, 100)))

    labelsize::Int64 = 24
    ticklabelsize::Int64 = 22
    linewidth::Float64 = 2.
    markersize::Float64 = 1.

    fig = Figure(resolution=(920, 720))
    ax = CairoMakie.Axis(fig[1, 1], xlabel=L"\text{presynaptic }I_2 \text{-assembly}", ylabel=L"\text{postsynaptic }E \text{-assembly}", xlabelsize=labelsize, ylabelsize=labelsize,
                xticks=(1:4:20, [L"%$(x)" for x=1:4:20]), xticklabelsize=ticklabelsize, xgridvisible=false,
                yticks=(1:4:20, [L"%$(x)" for x=1:4:20]), yticklabelsize=ticklabelsize, ygridvisible=false)
    ht = heatmap!(ax, weightsEI, colormap=inhibition_cs)
    hlines!([4.5, 8.5, 12.5, 16.5], linewidth=linewidth, color=:firebrick, linestyle=:dot)
    vlines!([4.5, 8.5, 12.5, 16.5], linewidth=linewidth, color=:firebrick, linestyle=:dot)

    vlines!(4.5, ymin=0., ymax=.4, linewidth=linewidth, color=:firebrick)
    vlines!(8.5, ymin=.2, ymax=.6, linewidth=linewidth, color=:firebrick)
    vlines!(12.5, ymin=.4, ymax=.8, linewidth=linewidth, color=:firebrick)
    vlines!(16.5, ymin=.6, ymax=1., linewidth=linewidth, color=:firebrick)

    hlines!(4.5, xmin=0., xmax=.4, linewidth=linewidth, color=:firebrick)
    hlines!(8.5, xmin=.2, xmax=.6, linewidth=linewidth, color=:firebrick)
    hlines!(12.5, xmin=.4, xmax=.8, linewidth=linewidth, color=:firebrick)
    hlines!(16.5, xmin=.6, xmax=1., linewidth=linewidth, color=:firebrick)

    Colorbar(fig[1, 2], ht, label=L"\text{coupling strength (pF)}", height=Relative(1.), labelsize=labelsize, ticklabelsize=ticklabelsize)

    (!isdir(output_dir)) && (mkpath(output_dir))
    save(joinpath(output_dir, string("weights_EI", name, ".png")), fig)
end