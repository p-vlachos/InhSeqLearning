using Pkg
Pkg.activate(".")
using InhSequences
using Statistics
using HDF5

# include("quantification.jl")

using LaTeXStrings
using CairoMakie
using Colors

############## THIS FOR TESTING ########################
mode = "spontaneous" # "spontaneous" "stimulation"
sim_name = string("network_7_", mode, ".h5")
sim_savedpath = string("./networks_trained_", mode, "/")
output_dir = "./output_testing/"

fid = h5open(joinpath(sim_savedpath, sim_name), "r")
popmembers = read(fid["data"]["popmembers"])
weights = read(fid["data"]["weights"])
# weightsEE = read(fid["data"]["weightsEE"])
# weightsEI = read(fid["data"]["weightsEI"])
# weightsIE = read(fid["data"]["weightsIE"])
times = read(fid["data"]["times"])
close(fid)
########################################################


function plot_activity_network(times::Matrix{Float64})
    # This function creates a raster plot of the network activity

    fig = Figure(resolution=(720, 480))
    ax = Axis(fig[1, 1], xlabel="", ylabel="", xlabelsize=20, ylabelsize=20,
                xticks=([], []), xticklabelsize=20,
                yticks=([], []), yticklabelsize=20,
                limits=(0, 5000, 1, 5000))
    CairoMakie.scatter!()

end

# Parameters
Ne = 4000
Ni2 = 250
Ncells = size(times)[1]
Npop = size(popmembers, 2)
Nmembers_max = size(popmembers, 1)
Ni_members = 27
interval = 2001:5000

ipopmembers = findI2populations(weights, Npop, popmembers, iipop_len=Ni_members)
restcells = deleteat!(map(*, ones(Int, 4000), range(1,stop=4000)), sort(unique(popmembers))[2:end])
ylim_max = count(i->(i>0), popmembers) + length(restcells) + (Ncells-Ni2-Ne) + length(ipopmembers)

rates = convolveSpikes(times, interval=interval)
emembers = unique(filter(i->i>0, popmembers))

ns = zeros(Ncells)
for cc = 1:Ncells
    ns[cc] = count(i->i>0, times[cc, :])
end

# Plot params
mrksize = 2.

rowcount = 1
fig = CairoMakie.Figure(resolution=(720, 480))
g = fig[1, 1] = GridLayout()
ax = CairoMakie.Axis(g[2, 1], xlabel=L"\text{Simulation time (s)}", ylabel=L"\text{Neuron #}", xlabelsize=24, ylabelsize=24,
            xticks=(collect((minimum(interval)-1):500:maximum(interval)), [L"0", "", L"1", "", L"2", "", L"3"]), xticklabelsize=18, xgridvisible=false,
            yticks=([], []), yticklabelsize=20, ygridvisible=false,
            limits=(minimum(interval)-1, maximum(interval), 1, ylim_max))
# Excitatoy assembly members
for pp = 1:Npop
    for cc = 1:Nmembers_max
        (popmembers[cc, pp] < 1) && (break)
        indx = round(Int, popmembers[cc, pp])
        x_vals = times[indx, 1:round(Int, ns[indx])]
        y_vals = rowcount * ones(length(x_vals))
        CairoMakie.scatter!(x_vals, y_vals, color=RGBA(14/255, 134/255, 212/255), markersize=mrksize)
        rowcount += 1
    end
end
hlines!(rowcount, linewidth=0.5, color=:gray20)
# Excitatory non-members 
for cc in restcells
    x_vals = times[cc, 1:round(Int, ns[cc])]
    y_vals = rowcount * ones(length(x_vals))
    CairoMakie.scatter!(x_vals, y_vals, color=RGBA(103/255, 136/255, 150/255), markersize=mrksize)
    rowcount += 1
end
hlines!(rowcount, linewidth=0.5, color=:gray20)
# Inhibitory I₁
for cc = (Ne+1):(Ncells-Ni2+1)
    x_vals = times[cc, 1:round(Int, ns[cc])]
    y_vals = rowcount * ones(length(x_vals))
    CairoMakie.scatter!(x_vals, y_vals, color=RGBA(131/255, 30/255, 19/255), markersize=mrksize)
    rowcount += 1
end
hlines!(rowcount, linewidth=0.5, color=:gray20)
# I₂ assembly members
for pp = 1:Npop
    for cc = 1:Ni_members
        indx = round(Int, popmembers[cc, pp])
        x_vals = times[indx, 1:round(Int, ns[indx])]
        y_vals = rowcount * ones(length(x_vals))
        CairoMakie.scatter!(x_vals, y_vals, color=RGBA(220/255, 50/255, 30/255), markersize=mrksize)
        rowcount += 1
    end
end



ax_top = CairoMakie.Axis(g[1, 1], xlabel="", ylabel=L"\text{Firing rate (Hz)}", xlabelsize=20, ylabelsize=20,
                xticks=([], []), xticklabelsize=20,
                yticks=([1, 3, 5], [L"1", L"3", L"5"]), yticklabelsize=18,
                limits=(1, 3_000, 0., 6.6))
# Plot mean firing rates
CairoMakie.lines!(ax_top, vec(mean(rates[emembers, :], dims=1)), color=RGBA(14/255, 134/255, 212/255), linewidth=1.)
CairoMakie.lines!(ax_top, vec(mean(rates[restcells, :], dims=1)), color=RGBA(103/255, 136/255, 150/255), linewidth=1.)
CairoMakie.lines!(ax_top, vec(mean(rates[(Ne+1):(Ncells-Ni2), :], dims=1)), color=RGBA(131/255, 30/255, 19/255), linewidth=1.)
CairoMakie.lines!(ax_top, vec(mean(rates[(Ncells-Ni2+1):Ncells, :], dims=1)), color=RGBA(220/255, 50/255, 30/255), linewidth=1.)

# linkxaxes!(ax, ax_top); #xlims!(ax, low=1, high=5000)
# hidedecorations!(ax, grid=true)
rowsize!(g, 1, Relative(1/5))

fig






(!ispath(output_dir)) && (mkpath(output_dir))
savefig(string(output_dir, "test_plot_1_sim_", sim_num,"_spontaneous.png"), dpi=150)
PyPlot.clf()









if doplot
    (!ispath(output_dir)) && (mkpath(output_dir))
    plot_eeWeights(new_weights, popmembers, Npop, sim_outDir)
    plot_eiWeights(new_weights, popmembers, ipopmembers, Npop, sim_outDir)
    plot_ieWeights(new_weights, popmembers, ipopmembers, Npop, sim_outDir)
end


mean(rates[emembers, :], dims=1)

rates[emembers, :]

rates[[2, 3, 4], 1]

minimum(emembers)

