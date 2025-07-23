using ColorSchemes
using LaTeXStrings
using CairoMakie
using Colors


function plotNetworkActivity(times::Matrix{Float64}, popmembers::Matrix{Int64}, ipopmembers::Matrix{Int64}; seq_length::Int64=4, interval::AbstractVector, name::AbstractString, output_dir::AbstractString="./output_analysis/")
    Ncells::Int64 = size(times)[1]
    Ne::Int64 = round(Int, Ncells*.8)
    Ni2::Int64 = round(Int, (Ncells*0.2)/2)
    Npop::Int64 = size(popmembers, 2)
    Nmembers_max::Int64 = size(popmembers, 1)
    Ni_members::Int64 = 20
    plot_interval::Vector{Int64} = collect((minimum(interval)):1000:maximum(interval))

    labelsize::Int64 = 24
    ticklabelsize::Int64 = 22
    linewidth::Float64 = 2.
    markersize::Float64 = 1.

    restcells::Vector{Int64} = deleteat!(map(*, ones(Int, Ne), range(1, stop=Ne)), sort(unique(popmembers))[2:end])
    ylim_max::Int64 = count(i->(i>0), popmembers) + length(restcells) + (Ncells-Ni2-Ne) + length(ipopmembers)

    ns::Vector{Int64} = zeros(Ncells)
    for cc = 1:Ncells
        ns[cc] = count(i->i>0, times[cc, :])
    end

    # Plot params
    ytick_seq::Vector{Int64} = zeros(round(Int, Npop/seq_length)*2)
    rowcount::Int64 = 1
    ytickcount::Int64 = 1

    fig = Figure(size=(1080, 720))
    ax = Axis(fig[1, 1], xlabel=L"\text{simulation time (s)}", ylabel=L"\text{sequences (neurons)}", xlabelsize=labelsize, ylabelsize=labelsize,
                xticks=(plot_interval, [L"%$(x)" for x = 0:length(plot_interval)-1]), xticklabelsize=ticklabelsize, xgridvisible=false,
                yticklabelsize=ticklabelsize, ygridvisible=false,
                limits=(minimum(interval), maximum(interval), 1, ylim_max))
    # Excitatoy assembly members
    for pp = 1:Npop
        for cc = 1:Nmembers_max
            (popmembers[cc, pp] < 1) && (break)
            indx = round(Int, popmembers[cc, pp])
            x_vals = times[indx, 1:round(Int, ns[indx])]
            y_vals = rowcount * ones(length(x_vals))
            scatter!(ax, x_vals, y_vals, color=ColorSchemes.Blues[7], markersize=markersize)
            rowcount += 1
        end
        if mod(pp, seq_length) == 0
            hlines!(ax, rowcount, linewidth=linewidth*.8, color=ColorSchemes.Blues[7])
        else
            hlines!(ax, rowcount, linewidth=linewidth*.3, color=ColorSchemes.Blues[7])
        end
        if mod(pp+2, seq_length) == 0
            ytick_seq[ytickcount] = rowcount
            ytickcount += 1
        end
    end
    hlines!(ax, rowcount, linewidth=linewidth*.5, color=:gray20)
    # Excitatory non-members 
    for cc in restcells
        x_vals = times[cc, 1:round(Int, ns[cc])]
        y_vals = rowcount * ones(length(x_vals))
        scatter!(ax, x_vals, y_vals, color=ColorSchemes.Greys[7], markersize=markersize)
        rowcount += 1
    end
    hlines!(ax, rowcount, linewidth=linewidth*.5, color=:gray20)
    # Inhibitory I₁
    for cc = (Ne+1):(Ncells-Ni2+1)
        x_vals = times[cc, 1:round(Int, ns[cc])]
        y_vals = rowcount * ones(length(x_vals))
        scatter!(ax, x_vals, y_vals, color=ColorSchemes.Reds[9], markersize=markersize)
        rowcount += 1
    end
    hlines!(ax, rowcount, linewidth=linewidth*.5, color=:gray20)
    # I₂ assembly members
    for pp = 1:Npop
        for cc = 1:Ni_members
            indx = round(Int, popmembers[cc, pp])
            x_vals = times[indx, 1:round(Int, ns[indx])]
            y_vals = rowcount * ones(length(x_vals))
            scatter!(ax, x_vals, y_vals, color=ColorSchemes.Reds[6], markersize=markersize)
            rowcount += 1
        end
        if mod(pp, seq_length) == 0
            hlines!(ax, rowcount, linewidth=linewidth*.3, color=ColorSchemes.Reds[9])
        end
        if mod(pp+2, seq_length) == 0
            ytick_seq[ytickcount] = rowcount
            ytickcount += 1
        end
    end

    ytick_labels::Vector{Char} = [Char(i+64) for i in 1:(round(Int, Npop/seq_length))]
    ytick_labels = vcat(ytick_labels, ytick_labels)
    ax.yticks = (ytick_seq, [L"\text{%$(x)}" for x in ytick_labels])

    (!isdir(output_dir)) && (mkpath(output_dir))
    save(joinpath(output_dir, string("network_activity", name, ".png")), fig)
end

function plotWeightsEE(weightsEE::Matrix{Float64}; seq_length::Int64=3, name::AbstractString, output_dir::AbstractString="./output_analysis/")
    cl1 = ColorScheme(range(colorant"gray5", colorant"gray80", length=100))
    cl_exc = ColorScheme(range(colorant"gray80", colorant"dodgerblue4", length=100))
    excitation_cs = vcat(get(cl1, LinRange(0, 1, 100)), get(cl_exc, LinRange(0, 1, 100)))

    Npop::Int64 = size(weightsEE)[1]
    line_positions::Vector{Float64} = collect(seq_length+.5:seq_length:Npop+.5)

    labelsize::Int64 = 24
    ticklabelsize::Int64 = 22
    linewidth::Float64 = 2.

    fig = Figure(size=(920, 720))
    ax = Axis(fig[1, 1], xlabel=L"\text{presynaptic }E\text{-assembly}", ylabel=L"\text{postsynaptic }E\text{-assembly}", xlabelsize=labelsize, ylabelsize=labelsize,
                xticks=(1:seq_length:Npop, [L"%$(x)" for x=1:seq_length:Npop]), xticklabelsize=ticklabelsize, xgridvisible=false,
                yticks=(1:seq_length:Npop, [L"%$(x)" for x=1:seq_length:Npop]), yticklabelsize=ticklabelsize, ygridvisible=false)
    ht = heatmap!(ax, weightsEE, colormap=excitation_cs)
    hlines!(line_positions, linewidth=linewidth, color=:firebrick, linestyle=:dot)
    vlines!(line_positions, linewidth=linewidth, color=:firebrick, linestyle=:dot)

    for (ind, ln) in enumerate(line_positions)
        vlines!(ln, ymin=(ind-1)*seq_length/Npop, ymax=ind*seq_length/Npop, linewidth=linewidth, color=:firebrick)
        hlines!(ln, xmin=(ind-1)*seq_length/Npop, xmax=ind*seq_length/Npop, linewidth=linewidth, color=:firebrick)
    end

    Colorbar(fig[1, 2], ht, label=L"\text{coupling strength (pF)}", height=Relative(1.), labelsize=labelsize, ticklabelsize=ticklabelsize)

    (!isdir(output_dir)) && (mkpath(output_dir))
    save(joinpath(output_dir, string("weights_EE", name, ".png")), fig)
end


function plotWeightsIE(weightsIE::Matrix{Float64}; seq_length::Int64=3, name::AbstractString, output_dir::AbstractString="./output_analysis/")
    cl1 = ColorScheme(range(colorant"gray5", colorant"gray80", length=100))
    cl_inh = ColorScheme(range(colorant"gray80", colorant"firebrick", length=100))
    inhibition_cs = vcat(get(cl1, LinRange(0, 1, 100)), get(cl_inh, LinRange(0, 1, 100)))

    Npop::Int64 = size(weightsIE)[1]
    line_positions::Vector{Float64} = collect(seq_length+.5:seq_length:Npop+.5)

    labelsize::Int64 = 24
    ticklabelsize::Int64 = 22
    linewidth::Float64 = 2.
    
    fig = Figure(size=(920, 720))
    ax = Axis(fig[1, 1], xlabel=L"\text{presynaptic }E \text{-assembly}", ylabel=L"\text{postsynaptic }I_2 \text{-assembly}", xlabelsize=labelsize, ylabelsize=labelsize,
                xticks=(1:seq_length:Npop, [L"%$(x)" for x=1:seq_length:Npop]), xticklabelsize=ticklabelsize, xgridvisible=false,
                yticks=(1:seq_length:Npop, [L"%$(x)" for x=1:seq_length:Npop]), yticklabelsize=ticklabelsize, ygridvisible=false)
    ht = heatmap!(ax, weightsIE, colormap=inhibition_cs)
    hlines!(line_positions, linewidth=linewidth, color=:dodgerblue4, linestyle=:dot)
    vlines!(line_positions, linewidth=linewidth, color=:dodgerblue4, linestyle=:dot)

    for (ind, ln) in enumerate(line_positions)
        vlines!(ln, ymin=(ind-1)*seq_length/Npop, ymax=ind*seq_length/Npop, linewidth=linewidth, color=:dodgerblue4)
        hlines!(ln, xmin=(ind-1)*seq_length/Npop, xmax=ind*seq_length/Npop, linewidth=linewidth, color=:dodgerblue4)
    end

    Colorbar(fig[1, 2], ht, label=L"\text{coupling strength (pF)}", height=Relative(1.), labelsize=labelsize, ticklabelsize=ticklabelsize)

    (!isdir(output_dir)) && (mkpath(output_dir))
    save(joinpath(output_dir, string("weights_IE", name, ".png")), fig)
end


function plotWeightsEI(weightsEI::Matrix{Float64}; seq_length::Int64=3, name::AbstractString, output_dir::AbstractString="./output_analysis/")
    cl1 = ColorScheme(range(colorant"gray5", colorant"gray80", length=100))
    cl_inh = ColorScheme(range(colorant"gray80", colorant"dodgerblue2", length=100))
    inhibition_cs = vcat(get(cl1, LinRange(0, 1, 100)), get(cl_inh, LinRange(0, 1, 100)))

    Npop::Int64 = size(weightsEI)[1]
    line_positions::Vector{Float64} = collect(seq_length+.5:seq_length:Npop+.5)

    labelsize::Int64 = 24
    ticklabelsize::Int64 = 22
    linewidth::Float64 = 2.

    fig = Figure(size=(920, 720))
    ax = Axis(fig[1, 1], xlabel=L"\text{presynaptic }E \text{-assembly}", ylabel=L"\text{postsynaptic }I_2 \text{-assembly}", xlabelsize=labelsize, ylabelsize=labelsize,
                xticks=(1:seq_length:Npop, [L"%$(x)" for x=1:seq_length:Npop]), xticklabelsize=ticklabelsize, xgridvisible=false,
                yticks=(1:seq_length:Npop, [L"%$(x)" for x=1:seq_length:Npop]), yticklabelsize=ticklabelsize, ygridvisible=false)
    ht = heatmap!(ax, weightsEI, colormap=inhibition_cs)
    hlines!(line_positions, linewidth=linewidth, color=:firebrick, linestyle=:dot)
    vlines!(line_positions, linewidth=linewidth, color=:firebrick, linestyle=:dot)

    for (ind, ln) in enumerate(line_positions)
        vlines!(ln, ymin=(ind-1)*seq_length/Npop, ymax=ind*seq_length/Npop, linewidth=linewidth, color=:firebrick)
        hlines!(ln, xmin=(ind-1)*seq_length/Npop, xmax=ind*seq_length/Npop, linewidth=linewidth, color=:firebrick)
    end

    Colorbar(fig[1, 2], ht, label=L"\text{coupling strength (pF)}", height=Relative(1.), labelsize=labelsize, ticklabelsize=ticklabelsize)

    (!isdir(output_dir)) && (mkpath(output_dir))
    save(joinpath(output_dir, string("weights_EI", name, ".png")), fig)
end