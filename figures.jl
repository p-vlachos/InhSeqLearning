using Pkg
Pkg.activate(".")
using InhSequences
using Statistics
using HDF5

using ColorSchemes
using LaTeXStrings
using CairoMakie
using CurveFit
using Colors
using FileIO

output_dir = "./output_analysis/Figures"

makefig_1 = true
makefig_2 = true
makefig_3 = true

##############################################################################################
######                               --- FIGURE 1 ---                                   ###### 
##############################################################################################
if makefig_1

fig = Figure(size=(4800, 3280))

g_model = fig[1, 1:4] = GridLayout()
g_raster = fig[1, 5:12] = GridLayout()
g_react = fig[2, 1:5] = GridLayout()
g_corr = fig[2, 6:9] = GridLayout()
g_seq = fig[2, 10:12] = GridLayout()


labelsize = 66
ticklabelsize = 54
legendlabelsize = 54
textsize = 66
subplotlabelsize = 66

dodge = 1
dodge_gap = .1
v_width = 1.
b_width = 1.
b_strokewidth = 1.

linewidth = 4
markersize = 5

##############################################################################################
###################################  --- Network Model --- ###################################
##############################################################################################

img_model = load(assetpath(joinpath(pwd(), "./analysis_data/network.png")))

ax_model = Axis(g_model[1, 1], aspect=DataAspect())
img = image!(ax_model, rotr90(img_model))

hidespines!(ax_model)
hidedecorations!(ax_model)

##############################################################################################
####################################  --- Raster plot --- ####################################
##############################################################################################

########################################################
################ --- Pre-Training --- ##################
########################################################
# --- Load Data ---
sim_name = string("network_pre-train.h5")
sim_savedpath = string("./networks_trained/")
fid = h5open(joinpath(sim_savedpath, sim_name), "r")
popmembers = read(fid["data"]["popmembers"])
weights = read(fid["data"]["weights"])
times = read(fid["data"]["times"])
close(fid)

# Parameters
Ne = 4000
Ni2 = 400
seq_length = 4
Ni_members = 25

Ncells = size(times)[1]
Npop = size(popmembers, 2)
Nmaxmembers = size(popmembers, 1)
interval = collect(12_000:17_000)

# Plot params
ytick_seq = zeros(round(Int, Npop/seq_length)*2)
x_ticks = minimum(interval):1000:maximum(interval)
rowcount = 1
ytickcount = 1

ipopmembers = findI2populations(weights, popmembers, iipop_len=Ni_members)
restcells = deleteat!(map(*, ones(Int, Ne), range(1, stop=Ne)), sort(unique(popmembers))[2:end])
ylim_max = count(i->(i>0), popmembers) + length(restcells) + (Ncells-Ni2-Ne) + length(ipopmembers)

# rates = convolveSpikes(times, interval=interval, gaussian=false, sigma=10.)
rates = binRates(times, interval=interval, window=100)
emembers = unique(filter(i->i>0, popmembers))

ns = zeros(Ncells)
for cc = 1:Ncells
    ns[cc] = count(i->i>0, times[cc, :])
end

# Plot params
ytick_seq = zeros(4)

rowcount = 1
ytickcount = 1
ytick_prev = 0

ax = Axis(g_raster[2, 1], xlabel=L"\text{simulation time (s)}", ylabel=L"\text{neurons}", xlabelsize=labelsize, ylabelsize=labelsize,
            xticks=(x_ticks, [L"%$(x)" for x in 1:length(x_ticks)]), xticklabelsize=ticklabelsize, xgridvisible=false,
            yticklabelsize=ticklabelsize, ygridvisible=false,
            limits=(minimum(interval), maximum(interval), 1, ylim_max))
# Excitatoy assembly members
for pp = 1:Npop
    for cc = 1:Nmaxmembers
        (popmembers[cc, pp] < 1) && (break)
        indx = round(Int, popmembers[cc, pp])
        x_vals = times[indx, 1:round(Int, ns[indx])]
        y_vals = rowcount * ones(length(x_vals))
        scatter!(ax, x_vals, y_vals, color=ColorSchemes.Blues[7], markersize=markersize)
        rowcount += 1
    end
end
hlines!(ax, rowcount, linewidth=linewidth*.5, color=:gray20)
ytick_seq[ytickcount] = rowcount/2
ytick_prev = rowcount
ytickcount += 1
# Excitatory non-members 
for cc in restcells
    x_vals = times[cc, 1:round(Int, ns[cc])]
    y_vals = rowcount * ones(length(x_vals))
    scatter!(ax, x_vals, y_vals, color=ColorSchemes.Greys[7], markersize=markersize)
    rowcount += 1
end
hlines!(ax, rowcount, linewidth=linewidth*.5, color=:gray20)
ytick_seq[ytickcount] = ytick_prev+(rowcount-ytick_prev)/2
ytick_prev = rowcount
ytickcount += 1
# Inhibitory I₁
for cc = (Ne+1):(Ncells-Ni2+1)
    x_vals = times[cc, 1:round(Int, ns[cc])]
    y_vals = rowcount * ones(length(x_vals))
    scatter!(ax, x_vals, y_vals, color=ColorSchemes.Reds[9], markersize=markersize)
    rowcount += 1
end
hlines!(ax, rowcount, linewidth=linewidth*.5, color=:gray20)
ytick_seq[ytickcount] = ytick_prev+(rowcount-ytick_prev)/2
ytick_prev = rowcount
ytickcount += 1
# I₂ assembly members
for pp = 1:Npop
    for cc = 1:Ni_members
        indx = round(Int, popmembers[cc, pp])
        x_vals = times[indx, 1:round(Int, ns[indx])]
        y_vals = rowcount * ones(length(x_vals))
        scatter!(ax, x_vals, y_vals, color=ColorSchemes.Reds[6], markersize=markersize)
        rowcount += 1
    end
end
ytick_seq[ytickcount] = ytick_prev+(rowcount-ytick_prev)/2
ax.yticks = (ytick_seq, [L"E\text{-members}", L"E\text{-rest}", L"I_1", L"I_2"])
ax.yticklabelrotation = pi/2

ax_top = Axis(g_raster[1, 1], xlabel="", ylabel=L"\text{firing rate (Hz)}", xlabelsize=labelsize, ylabelsize=labelsize,
                title=L"\textbf{pre stimulation}", titlesize=labelsize,
                ylabelpadding=10.,
                xticks=([], []), xticklabelsize=ticklabelsize,
                yticks=([0, 3, 6, 9], [L"0", L"3", L"6", L"9"]), yticklabelsize=ticklabelsize,
                limits=(minimum(interval), maximum(interval), -0.2, 10.1))
# Plot mean firing rates
lines!(ax_top, interval, vec(mean(rates[emembers, :], dims=1)), color=ColorSchemes.Blues[7], linewidth=linewidth)
lines!(ax_top, interval, vec(mean(rates[restcells, :], dims=1)), color=ColorSchemes.Greys[7], linewidth=linewidth)
lines!(ax_top, interval, vec(mean(rates[(Ne+1):(Ncells-Ni2), :], dims=1)), color=ColorSchemes.Reds[9], linewidth=linewidth)
lines!(ax_top, interval, vec(mean(rates[(Ncells-Ni2+1):Ncells, :], dims=1)), color=ColorSchemes.Reds[6], linewidth=linewidth)

linkxaxes!(ax, ax_top); #xlims!(ax, low=minimum(interval)+500, high=maximum(interval))
hidedecorations!(ax, grid=true, ticks=false, ticklabels=false, label=false)
hideydecorations!(ax, grid=true, ticks=true, ticklabels=false, label=false)
rowsize!(g_raster, 1, Relative(1/5))

########################################################
################### --- Training --- ###################
########################################################
# --- Load Data ---
sim_name = string("network_train.h5")
sim_savedpath = string("./networks_trained/")
fid = h5open(joinpath(sim_savedpath, sim_name), "r")
popmembers = read(fid["data"]["popmembers"])
weights = read(fid["data"]["weights"])
times = read(fid["data"]["times"])
close(fid)

interval = collect(10_000:25_000)

restcells = deleteat!(map(*, ones(Int, Ne), range(1, stop=Ne)), sort(unique(popmembers))[2:end])
ylim_max = count(i->(i>0), popmembers) + length(restcells) + (Ncells-Ni2-Ne) + length(ipopmembers)

ns = zeros(Ncells)
for cc = 1:Ncells
    ns[cc] = count(i->i>0, times[cc, :])
end

# Plot params
ytick_seq = zeros(round(Int, Npop/seq_length)*2)
rowcount = 1
ytickcount = 1
x_ticks = minimum(interval):1000:maximum(interval)

ipopmembers = findI2populations(weights, popmembers, iipop_len=Ni_members)
restcells = deleteat!(map(*, ones(Int, Ne), range(1, stop=Ne)), sort(unique(popmembers))[2:end])
ylim_max = count(i->(i>0), popmembers) + length(restcells) + (Ncells-Ni2-Ne) + length(ipopmembers)

# rates = convolveSpikes(times, interval=interval, gaussian=false, sigma=10.)
rates = binRates(times, interval=interval, window=100)
emembers = unique(filter(i->i>0, popmembers))

ns = zeros(Ncells)
for cc = 1:Ncells
    ns[cc] = count(i->i>0, times[cc, :])
end

ax = Axis(g_raster[2, 2], xlabel=L"\text{simulation time (s)}", ylabel=L"\text{sequences (neurons)}", xlabelsize=labelsize, ylabelsize=labelsize,
            xticks=(x_ticks, [L"%$(x)" for x in 1:length(x_ticks)]), xticklabelsize=ticklabelsize, xgridvisible=false,
            yticklabelsize=ticklabelsize, ygridvisible=false, yaxisposition=:right,
            limits=(minimum(interval), maximum(interval), 1, ylim_max))
# Excitatoy assembly members
for pp = 1:Npop
    for cc = 1:Nmaxmembers
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
    if mod(pp+1, seq_length) == 0
        ytick_seq[ytickcount] = rowcount + (Nmaxmembers/2)
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
    if mod(pp+1, seq_length) == 0
        ytick_seq[ytickcount] = rowcount + (Ni_members/2)
        ytickcount += 1
    end
end

ax_top = Axis(g_raster[1, 2], title=L"\textbf{stimulation}",  titlesize=labelsize,
                limits=(minimum(interval), maximum(interval), -0.2, 10.1))
# Plot mean firing rates
lines!(ax_top, interval, vec(mean(rates[emembers, :], dims=1)), color=ColorSchemes.Blues[7], linewidth=linewidth)
lines!(ax_top, interval, vec(mean(rates[restcells, :], dims=1)), color=ColorSchemes.Greys[7], linewidth=linewidth)
lines!(ax_top, interval, vec(mean(rates[(Ne+1):(Ncells-Ni2), :], dims=1)), color=ColorSchemes.Reds[9], linewidth=linewidth)
lines!(ax_top, interval, vec(mean(rates[(Ncells-Ni2+1):Ncells, :], dims=1)), color=ColorSchemes.Reds[6], linewidth=linewidth)

linkxaxes!(ax, ax_top); #xlims!(ax, low=minimum(interval)+500, high=maximum(interval))
hidedecorations!(ax, grid=true, ticks=false, ticklabels=false, label=false)
hideydecorations!(ax, grid=true, ticks=true, ticklabels=true, label=true)
hideydecorations!(ax_top, grid=false, ticks=true, ticklabels=true, label=true)
hidexdecorations!(ax_top, grid=false, ticks=true, ticklabels=true, label=true)
rowsize!(g_raster, 1, Relative(1/5))

########################################################
################ --- Post-Training --- #################
########################################################
# --- Load Data ---
net_num = 1
# seed_num = 1
mode = "spontaneous"
# sim_name = string("network_1_", mode, "_seed_", seed_num, ".h5")
sim_name = string("network_1_", mode, ".h5")
# sim_savedpath = string("./networks_trained_", mode, "/network_", net_num, "/")
sim_savedpath = string("./networks_trained_", mode, "/")
fid = h5open(joinpath(sim_savedpath, sim_name), "r")
popmembers = read(fid["data"]["popmembers"])
weights = read(fid["data"]["weights"])
times = read(fid["data"]["times"])
close(fid)

interval = collect(12_500:17_500)

restcells = deleteat!(map(*, ones(Int, Ne), range(1, stop=Ne)), sort(unique(popmembers))[2:end])
ylim_max = count(i->(i>0), popmembers) + length(restcells) + (Ncells-Ni2-Ne) + length(ipopmembers)

ns = zeros(Ncells)
for cc = 1:Ncells
    ns[cc] = count(i->i>0, times[cc, :])
end

# Plot params
ytick_seq = zeros(round(Int, Npop/seq_length)*2)
x_ticks = minimum(interval):1000:maximum(interval)
rowcount = 1
ytickcount = 1

ipopmembers = findI2populations(weights, popmembers, iipop_len=Ni_members)
restcells = deleteat!(map(*, ones(Int, Ne), range(1, stop=Ne)), sort(unique(popmembers))[2:end])
ylim_max = count(i->(i>0), popmembers) + length(restcells) + (Ncells-Ni2-Ne) + length(ipopmembers)

# rates = convolveSpikes(times, interval=interval, gaussian=false, sigma=10.)
rates = binRates(times, interval=interval)
emembers = unique(filter(i->i>0, popmembers))

ns = zeros(Ncells)
for cc = 1:Ncells
    ns[cc] = count(i->i>0, times[cc, :])
end

ax = Axis(g_raster[2, 3], xlabel=L"\text{simulation time (s)}", ylabel=L"\text{sequences}", xlabelsize=labelsize, ylabelsize=labelsize,
            xticks=(x_ticks, [L"%$(x)" for x in 1:length(x_ticks)]), xticklabelsize=ticklabelsize, xgridvisible=false,
            yticklabelsize=ticklabelsize, ygridvisible=false, yaxisposition=:right,
            limits=(minimum(interval), maximum(interval), 1, ylim_max))
# Excitatoy assembly members
for pp = 1:Npop
    for cc = 1:Nmaxmembers
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
    if mod(pp, seq_length) == 0
        ytick_seq[ytickcount] = rowcount - (Nmaxmembers*seq_length/2)
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
    if mod(pp, seq_length) == 0
        ytick_seq[ytickcount] = rowcount - (Ni_members*seq_length/2)
        ytickcount += 1
    end
end

ytick_labels = [Char(i+64) for i in 1:(round(Int, Npop/seq_length))]
ytick_labels = vcat(ytick_labels, ytick_labels)
ax.yticks = (ytick_seq, [L"\text{%$(x)}" for x in ytick_labels])


ax_top = Axis(g_raster[1, 3], xlabel="", ylabel="", xlabelsize=labelsize, ylabelsize=labelsize,
                title=L"\textbf{after stimulation}",  titlesize=labelsize,
                xticks=([], []), xticklabelsize=ticklabelsize,
                yticks=([0, 1, 2], [L"0", L"1", L"2"]), yticklabelsize=ticklabelsize, yaxisposition=:right,
                limits=(minimum(interval), maximum(interval), -0.2, 10.1))
# Plot mean firing rates
lines!(ax_top, interval, vec(mean(rates[emembers, :], dims=1)), color=ColorSchemes.Blues[7], linewidth=linewidth)
lines!(ax_top, interval, vec(mean(rates[restcells, :], dims=1)), color=ColorSchemes.Greys[7], linewidth=linewidth)
lines!(ax_top, interval, vec(mean(rates[(Ne+1):(Ncells-Ni2), :], dims=1)), color=ColorSchemes.Reds[9], linewidth=linewidth)
lines!(ax_top, interval, vec(mean(rates[(Ncells-Ni2+1):Ncells, :], dims=1)), color=ColorSchemes.Reds[6], linewidth=linewidth)

linkxaxes!(ax, ax_top); #xlims!(ax, low=minimum(interval)+500, high=maximum(interval))
hidedecorations!(ax, grid=true, ticks=false, ticklabels=false, label=false)
hideydecorations!(ax, grid=true, ticks=true, ticklabels=false, label=false)
hideydecorations!(ax_top, grid=true, ticks=true, ticklabels=true, label=true)
rowsize!(g_raster, 1, Relative(1/5))

##############################################################################################
#################################  --- Reactivation plot --- #################################
##############################################################################################

# Parameters
seq_number = 1
seq_length = 4
interval = 15_400:15_900

rates = binRates(times, interval=interval, window=100)
ylim_max = (Nmaxmembers*3) + (Ni_members*3)

ns = zeros(Ncells)
for cc = 1:Ncells
    ns[cc] = count(i->i>0, times[cc, :])
end

# Plot params
markersize = 8.
linewidth = 3.

rowcount = 1
ytick_seq = zeros(round(Int, seq_length*2))
x_ticks = minimum(interval):100:maximum(interval)
ytick_prev = 0.
ytickcount = 1
ax = Axis(g_react[1:2, 1], xlabel=L"\text{simulation time (ms)}", ylabel=L"\text{sequence %$(Char(64+seq_number))}",
            xlabelsize=labelsize, ylabelsize=labelsize,
            xticks=(x_ticks, [L"%$(x)" for x in 1:length(x_ticks)]),
            xticklabelsize=ticklabelsize, xgridvisible=false, yticklabelsize=ticklabelsize, ygridvisible=false,
            limits=(minimum(interval), maximum(interval), 1, ylim_max))
# Excitatoy assembly members
for pp = (((seq_number-1)*seq_length)+1):(seq_number*seq_length)
    for indx in filter(i->i>0, popmembers[:, pp])
        x_vals = times[indx, 1:round(Int, ns[indx])]
        y_vals = rowcount * ones(length(x_vals))
        scatter!(ax, x_vals, y_vals, color=ColorSchemes.Blues[5+ytickcount], markersize=markersize)
        rowcount += 1
    end
    hlines!(ax, rowcount, linewidth=linewidth*.5, color=ColorSchemes.Blues[9])
    if ytickcount == 1
        ytick_seq[ytickcount] = rowcount/2; ytickcount += 1; ytick_prev = rowcount
    else
        ytick_seq[ytickcount] = ytick_prev + (rowcount-ytick_prev)/2; ytickcount += 1; ytick_prev = rowcount
    end
end
hlines!(ax, rowcount, linewidth=linewidth*.5, color=:gray20)
# I₂ assembly members
for pp = (((seq_number-1)*seq_length)+1):(seq_number*seq_length)
    for indx in ipopmembers[:, pp]
        x_vals = times[indx, 1:round(Int, ns[indx])]
        y_vals = rowcount * ones(length(x_vals))
        scatter!(ax, x_vals, y_vals, color=ColorSchemes.Reds[1+ytickcount], markersize=markersize)
        rowcount += 1
    end
    hlines!(ax, rowcount, linewidth=linewidth*.5, color=ColorSchemes.Reds[9])
    ytick_seq[ytickcount] = ytick_prev + (rowcount-ytick_prev)/2; ytickcount += 1; ytick_prev = rowcount
end

ax.yticks = (ytick_seq, [L"%$(x)" for x in 1:length(ytick_seq)])

ax_exc = Axis(g_react[1, 2], xlabel="", ylabel=L"E-\text{assembly firing rate (Hz)}", xlabelsize=labelsize, ylabelsize=labelsize,
                xticks=(x_ticks, [L"%$(x)" for x in 1:length(x_ticks)]), xticklabelsize=ticklabelsize,
                yticks=([0, 4, 8, 12, 16, 20], [L"0", L"4", L"8", L"12", L"16", L"20"]), yticklabelsize=ticklabelsize)#,
                # limits=(minimum(interval), maximum(interval), -0.2, 18.3))
# Plot mean firing rates
for (indx, pp) in enumerate((((seq_number-1)*seq_length)+1):(seq_number*seq_length))
    lines!(ax_exc, interval, vec(mean(rates[filter(i->i>0, popmembers[:, pp]), :], dims=1)), color=ColorSchemes.Blues[5+indx], linewidth=linewidth, label=L"\text{E}%$(indx)")
end
axislegend(ax_exc, position=(0.99, 0.99), labelsize=legendlabelsize)

ax_inh = Axis(g_react[2, 2], xlabel=L"\text{simulation time (ms)}", ylabel=L"I_2-\text{assembly firing rate (Hz)}", xlabelsize=labelsize, ylabelsize=labelsize,
                xticks=(x_ticks, [L"%$(x)" for x in 1:length(x_ticks)]), xticklabelsize=ticklabelsize,
                yticks=([0, 4, 8, 12], [L"0", L"4", L"8", L"12"]), yticklabelsize=ticklabelsize)#,
                # limits=(minimum(interval), maximum(interval), -0.2, 12.8))
# Plot mean firing rates
for (indx, pp) in enumerate((((seq_number-1)*seq_length)+1):(seq_number*seq_length))
    lines!(ax_inh, interval, vec(mean(rates[ipopmembers[:, pp], :], dims=1)), color=ColorSchemes.Reds[5+indx], linewidth=linewidth, label=L"\text{I}_2%$(indx)")
end
axislegend(ax_inh, position=(0.99, 0.99), labelsize=legendlabelsize)

linkxaxes!(ax_exc, ax_inh) #; xlims!(ax, low=minimum(interval)+500, high=maximum(interval))
hidexdecorations!(ax_exc, grid=false, ticks=true, ticklabels=true, label=true)

colgap!(g_react, 1, Relative(.05))

# ##############################################################################################
# ############################# --- Cross-correlation & Delays --- #############################
# ##############################################################################################

# activity_type = "spontaneous"   # Choose between "stimulation" or "spontaneous"
# sim_name = string("cross_corr_", activity_type, ".h5")
# sim_savedpath = "./analysis_data/"

# fid = h5open(joinpath(sim_savedpath, sim_name), "r")
# crossEE = read(fid["data"]["crossEE"])
# crossEI = read(fid["data"]["crossEI"])
# crossIE = read(fid["data"]["crossIE"])
# close(fid)


# Nseq = Npop/seq_length
# # Average over all simulations
# avg_cross = mean(crossEE, dims=4)[:, :, :, 1]
# # Average over all sequences (for 4 seq. of length 3)
# cross0 = mean(avg_cross[:, 1, 1], dims=2)
# cross0 += mean(avg_cross[:, 4, 4], dims=2)
# cross0 += mean(avg_cross[:, 7, 7], dims=2)
# cross0 += mean(avg_cross[:, 10, 10], dims=2)

# cross1 = mean(avg_cross[:, 1, 2], dims=2)
# cross1 += mean(avg_cross[:, 4, 5], dims=2)
# cross1 += mean(avg_cross[:, 7, 8], dims=2)
# cross1 += mean(avg_cross[:, 10, 11], dims=2)

# cross2 = mean(avg_cross[:, 1, 3], dims=2)
# cross2 += mean(avg_cross[:, 4, 6], dims=2)
# cross2 += mean(avg_cross[:, 7, 9], dims=2)
# cross2 += mean(avg_cross[:, 10, 12], dims=2)

# ax_mean = Axis(g_corr[1, 1], 
#             xlabel=L"\text{lag (ms)}", xlabelsize=labelsize,
#             ylabel=L"\text{cross-correlation}", ylabelsize=labelsize,
#             xticks=([-200., -100., 0., 100., 200.], [L"-200", L"-100", L"0", L"100", L"200"]), xticklabelsize=ticklabelsize,
#             yticks=([0., .5, 1.], [L"0", L"0.5", L"1"]), yticklabelsize=ticklabelsize,
#             limits=(-280, 280, -0.15, 1.1))
# lines!(ax_mean, -400:400, cross0./Nseq, color=ColorSchemes.Blues[6], linewidth=linewidth, label=L"E1-E1")
# lines!(ax_mean, -400:400, cross1./Nseq, color=ColorSchemes.Blues[7], linewidth=linewidth, label=L"E1-E2")
# lines!(ax_mean, -400:400, cross2./Nseq, color=ColorSchemes.Blues[8], linewidth=linewidth, label=L"E1-E3")
# axislegend(ax_mean, position=(0.95, 0.95), labelsize=legendlabelsize)

# # ____________ --- Delays --- ____________
# # Average over all simulations
# avg_crossEE = mean(crossEE, dims=4)[:, :, :, 1]
# avg_crossEI = mean(crossEI, dims=4)[:, :, :, 1]
# avg_crossIE = mean(crossIE, dims=4)[:, :, :, 1]

# # Calculate the peak cross-correlation for each assembly
# max_delayEE = zeros(Npop, Npop)
# max_delayEI = zeros(Npop, Npop)
# max_delayIE = zeros(Npop, Npop)
# for ipop = 1:Npop
#     for iipop = 1:Npop
#         max_delayEE[ipop, iipop] = findmax(avg_crossEE[:, ipop, iipop])[2] - 400
#         max_delayEI[ipop, iipop] = findmax(avg_crossEI[:, ipop, iipop])[2] - 400
#         max_delayIE[ipop, iipop] = findmax(avg_crossIE[:, ipop, iipop])[2] - 400
#     end
# end

# # Extract the relative cross-correlation delays
# delaysEE = zeros(4, 3)
# delaysEI = zeros(4, 3)
# delaysIE = zeros(4, 2)
# for ipop = 1:seq_length:Npop
#     delaysEE[div(ipop, seq_length) + 1, 1] = max_delayEE[ipop, ipop]
#     delaysEE[div(ipop, seq_length) + 1, 2] = max_delayEE[ipop, ipop + 1]
#     delaysEE[div(ipop, seq_length) + 1, 3] = max_delayEE[ipop, ipop + 2]

#     delaysEI[div(ipop, seq_length) + 1, 1] = max_delayEI[ipop, ipop]
#     delaysEI[div(ipop, seq_length) + 1, 2] = max_delayEI[ipop, ipop + 1]
#     delaysEI[div(ipop, seq_length) + 1, 3] = max_delayEI[ipop, ipop + 2]
    
#     delaysIE[div(ipop, seq_length) + 1, 1] = max_delayIE[ipop, ipop + 1]
#     delaysIE[div(ipop, seq_length) + 1, 2] = max_delayIE[ipop, ipop + 2]
# end

# # Fit the mean values to get the recall speed
# fitted_mean = linear_fit(mean(delaysEE, dims=1)[:], [1, 2, 3])

# # Compute the error bars
# low_errorsEE = mean(delaysEE, dims=1)[:] - minimum(delaysEE, dims=1)[:]
# high_errorsEE = maximum(delaysEE, dims=1)[:] - mean(delaysEE, dims=1)[:]

# low_errorsEI = mean(delaysEI, dims=1)[:] - minimum(delaysEI, dims=1)[:]
# high_errorsEI = maximum(delaysEI, dims=1)[:] - mean(delaysEI, dims=1)[:]

# low_errorsIE = mean(delaysIE, dims=1)[:] - minimum(delaysIE, dims=1)[:]
# high_errorsIE = maximum(delaysIE, dims=1)[:] - mean(delaysIE, dims=1)[:]

# # Plot the data
# ax = Axis(g_corr[2, 1], xlabel=L"\text{delay (ms)}", ylabel=L"\text{activity trajectory}",
#                     yticks=([1, 1.25, 1.5, 2, 2.25, 2.5, 3, 3.25], [L"E1→E1", L"E1→I_21", L"I_21→E2", L"E2→E2", L"E2→I_22", L"I_22→E3", L"E3→E3", L"E3→I_23"]),
#                     # xticks=([0, 25, 50, 75, 100, 125], [L"0", L"25", L"50", L"75", L"100", L"125"]),
#                     xlabelsize=labelsize, ylabelsize=labelsize, xticklabelsize=ticklabelsize, yticklabelsize=ticklabelsize)#, limits=(-0.5, 146, 0.8, 4.5))

# scatter!(ax, mean(delaysEE, dims=1)[:], [1, 2, 3], markersize=markersize, color=ColorSchemes.Blues[6])
# errorbars!(ax, mean(delaysEE, dims=1)[:], [1, 2, 3], low_errorsEE, high_errorsEE, whiskerwidth=20,  direction=:x, linewidth=linewidth, color=ColorSchemes.Blues[6])

# scatter!(ax, mean(delaysEI, dims=1)[:], [1.25, 2.25, 3.25], markersize=markersize, color=ColorSchemes.Reds[6])
# errorbars!(ax, mean(delaysEI, dims=1)[:], [1.25, 2.25, 3.25], low_errorsEI, high_errorsEI, whiskerwidth=20,  direction=:x, linewidth=linewidth, color=ColorSchemes.Reds[6])

# scatter!(ax, mean(delaysIE, dims=1)[:], [1.5, 2.5], markersize=markersize, color=ColorSchemes.BuPu[6])
# errorbars!(ax, mean(delaysIE, dims=1)[:], [1.5, 2.5], low_errorsIE, high_errorsIE, whiskerwidth=20,  direction=:x, linewidth=linewidth, color=ColorSchemes.BuPu[6])


# ##############################################################################################
# ################################### --- Sequentiality --- ####################################
# ##############################################################################################

# sim_savedpath = "./analysis_data/sequentiality/"
# Nsimulations = 10
# Nseeds = 10
# seq_scores = zeros(Nsimulations, Nseeds)

# for sim = 1:Nsimulations
#     sim_name = string("sequentiality_network_", sim,"_", activity_type, ".h5")
#     fid = h5open(joinpath(sim_savedpath, sim_name), "r")
#     seq_scores[sim, :] = read(fid["data"]["scores"])
#     close(fid)
# end

# mean_scores = mean(seq_scores, dims=2)[:, 1]
# seq_scores = vec(seq_scores)

# categories = ones(length(seq_scores)+length(mean_scores))
# categories[end-length(mean_scores):end] .= 2

# ax = Axis(g_seq[1, 1], ylabel=L"\text{sequentiality (%)}", ylabelsize=labelsize,
#             xticks=([1, 2], [L"\text{overall}", L"\text{mean / model seed}"]), xticklabelsize=ticklabelsize,
#             yticks=([.9, .95, 1.], [L"90", L"95", L"100"]), yticklabelsize=ticklabelsize,
#             limits=(0, 3, .88, 1.02))

# violin!(ax, categories, vcat(seq_scores, mean_scores), dodge=dodge, dodge_gap=dodge_gap, width=v_width)
# boxplot!(ax, categories, vcat(seq_scores, mean_scores), color=(:white, .0), strokewidth=b_strokewidth, width=b_width, strokecolor=(:black, 1.))



Label(g_model[1, 0, TopLeft()], L"\textbf{a}",
    fontsize=subplotlabelsize,
    font=:bold,
    padding=(0, 5, 5, 0),
    halign=:right)

Label(g_raster[1, 0, TopLeft()], L"\textbf{b}",
    fontsize=subplotlabelsize,
    font=:bold,
    padding=(0, 5, 5, 0),
    halign=:right)

Label(g_react[1, 0, TopLeft()], L"\textbf{c}",
    fontsize=subplotlabelsize,
    font=:bold,
    padding=(0, 5, 5, 0),
    halign=:right)

# Label(g_corr[1, 0, TopLeft()], L"\textbf{d}",
#     fontsize=subplotlabelsize,
#     font=:bold,
#     padding=(0, 5, 5, 0),
#     halign=:right)

# Label(g_seq[1, 0, TopLeft()], L"\textbf{e}",
#     fontsize=subplotlabelsize,
#     font=:bold,
#     padding=(0, 5, 5, 0),
#     halign=:right)
    
colsize!(g_model, 0, Relative(0))
colsize!(g_raster, 0, Relative(0))
colsize!(g_react, 0, Relative(0))
# colsize!(g_corr, 0, Relative(0))
# colsize!(g_seq, 0, Relative(0))

rowgap!(fig.layout, 1, Relative(.05))
colgap!(fig.layout, 4, Relative(.05))
colgap!(fig.layout, 5, Relative(.05))
colgap!(fig.layout, 9, Relative(.05))

(!isdir(output_dir)) && (mkpath(output_dir))
save(joinpath(output_dir, string("figure_1.png")), fig)
end



##############################################################################################
######                               --- FIGURE 2 ---                                   ###### 
##############################################################################################

if makefig_2

fig = Figure(size=(4800, 3800))

g_model = fig[1, 1:2] = GridLayout()
g_hist = fig[1, 3:4] = GridLayout()
g_struct = fig[2, 1:4] = GridLayout()
g_learn = fig[3, 1:2] = GridLayout()
g_rest = fig[3, 3:4] = GridLayout()


labelsize = 66
ticklabelsize = 54
legendlabelsize = 54
textsize = 66
subplotlabelsize = 66

dodge = 1
dodge_gap = .1
v_width = 1.
b_width = 1.
b_strokewidth = 1.

linewidth = 6
markersize = 5

colorbarsize = 40

##############################################################################################
###################################  --- Network Model --- ###################################
##############################################################################################

img_model = load(assetpath(joinpath(pwd(), "./analysis_data/network.png")))

ax_model = Axis(g_model[1, 1], aspect=DataAspect())
img = image!(ax_model, rotr90(img_model))

hidespines!(ax_model)
hidedecorations!(ax_model)

##############################################################################################
############################# --- Histogram synaptic weights --- #############################
##############################################################################################
sim_savedpath = string("./networks_trained/")

# Choose the parameters
Ncells = 5000
Ni2 = 400
Npop = 20
Nsims = 1

x_vals = 1:Ni2
ie_weights = zeros(Ni2, Npop, Nsims)
sorted_ie_weights = zeros(Ni2, Npop, Nsims)
ie_perms = zeros(Ni2, Npop, Nsims)

for sim_num = 1:Nsims
    # Load data
    sim_name = string("network_", sim_num,".h5")
    fid = h5open(joinpath(sim_savedpath, sim_name), "r")
    popmembers = read(fid["data"]["popmembers"])
    weights = read(fid["data"]["weights"])
    close(fid)

    # Compute the power laws that fit the distributions
    for ipop = 1:Npop
        members = filter(i->(i>0), popmembers[:, ipop])    			# Get members of excitatory assembly
        ie_weights[:, ipop, sim_num] .= vec(sum(weights[members, Ncells-Ni2+1:Ncells], dims=1))  # Get sum of all weights projected from each E-assemble to each 2nd ipopulation neuron
        ie_perms[:, ipop, sim_num] .= sortperm(ie_weights[:, ipop, sim_num], rev=true)
        sorted_ie_weights[:, ipop, sim_num] .= ie_weights[round.(Int, ie_perms[:, ipop, sim_num]), ipop, sim_num]
    end
end

plot_data_band_min = minimum(minimum(sorted_ie_weights, dims=3)[:, :, 1], dims=2)[:]
plot_data_avg = mean(mean(sorted_ie_weights, dims=3)[:, :, 1], dims=2)[:]
plot_data_band_max = maximum(maximum(sorted_ie_weights, dims=3)[:, :, 1], dims=2)[:]

ax = Axis(g_hist[1, 1], xlabel=L"I_2 \, \text{neuron index (sorted)}", ylabel=L"\text{mean synaptic strength (pF)}",
                xticks=([1, 25, 50, 100, 150, 200, 250], [L"1", L"\mathbf{25}", L"50", L"100", L"150", L"200", L"250"]),
                yticks=([50, 100, 150, 200, 250], [L"50", L"100", L"150", L"200", L"250"]),
                xlabelsize=labelsize, ylabelsize=labelsize, xticklabelsize=ticklabelsize, yticklabelsize=ticklabelsize,
                limits=(0., 255.5, 0.5, 255.5))
band!(x_vals, plot_data_band_min, plot_data_band_max, label=L"\text{extremum values}")
lines!(x_vals, plot_data_avg, linewidth=linewidth*2, color=:gold3, label=L"\text{average values}")
vlines!(25, ymin=0, ymax=251, linewidth=linewidth, linestyle=(:dash, :loose), color=:gray20)
# text!(25, 240, text=L"\text{max slope change}", rotation=3*π/2, fontsize=textsize)
axislegend(ax, position=(.95, .95), labelsize=labelsize, titlesize=textsize)

##############################################################################################
############################# --- Histogram synaptic weights --- #############################
##############################################################################################
# --- Load Data ---
net_num = 1
seed_num = 1
mode = "spontaneous"
# sim_name = string("network_", net_num,"_", mode, "_seed_", seed_num, ".h5")
# sim_savedpath = string("./networks_trained_", mode, "/network_", net_num, "/")
sim_name = string("network_", net_num,"_", mode, ".h5")
sim_savedpath = string("./networks_trained_", mode, "/")
fid = h5open(joinpath(sim_savedpath, sim_name), "r")
popmembers = read(fid["data"]["popmembers"])
weights = read(fid["data"]["weights"])
close(fid)

# Parameters
Ne = 4000
Ni2 = 400
seq_length = 4
Ni_members = 25

Npop = size(popmembers, 2)
Nmaxmembers = size(popmembers, 1)

ipopmembers = findI2populations(weights, popmembers, iipop_len=Ni_members)

weightsEE = sum(weights[filter(i->i>0, popmembers[:, :]), filter(i->i>0, popmembers[:, :])], dims=(1, 3))[1, :, 1, :] ./ count(i->i>0, weights[filter(i->i>0, popmembers[:, :]), filter(i->i>0, popmembers[:, :])], dims=(1, 3))[1, :, 1, :]
weightsEI = sum(weights[filter(i->i>0, popmembers[:, :]), ipopmembers[:, :]], dims=(1, 3))[1, :, 1, :] ./ count(i->i>0, weights[filter(i->i>0, popmembers[:, :]), ipopmembers[:, :]], dims=(1, 3))[1, :, 1, :]
weightsIE = sum(weights[ipopmembers[:, :], filter(i->i>0, popmembers[:, :])], dims=(1, 3))[1, :, 1, :] ./ count(i->i>0, weights[ipopmembers[:, :], filter(i->i>0, popmembers[:, :])], dims=(1, 3))[1, :, 1, :]


########################################################
################ --- E-to-E weights --- ################
########################################################

cl1 = ColorScheme(range(colorant"gray5", colorant"gray80", length=100))
cl_exc = ColorScheme(range(colorant"gray80", colorant"dodgerblue4", length=100))
excitation_cs = vcat(get(cl1, LinRange(0, 1, 100)), get(cl_exc, LinRange(0, 1, 100)))

Npop = size(weightsEE)[1]
line_positions = collect(seq_length+.5:seq_length:Npop+.5)

ax = Axis(g_struct[1, 1], xlabel=L"\text{presynaptic }E\text{-assembly}", ylabel=L"\text{postsynaptic }E\text{-assembly}", xlabelsize=labelsize, ylabelsize=labelsize,
            xticks=(1:seq_length:Npop, [L"%$(x)" for x=1:seq_length:Npop]), xticklabelsize=ticklabelsize, xgridvisible=false,
            yticks=(1:seq_length:Npop, [L"%$(x)" for x=1:seq_length:Npop]), yticklabelsize=ticklabelsize, ygridvisible=false,
            aspect=1)
ht = heatmap!(ax, weightsEE, colormap=excitation_cs)
hlines!(ax, line_positions, linewidth=linewidth, color=:firebrick, linestyle=:dot)
vlines!(ax, line_positions, linewidth=linewidth, color=:firebrick, linestyle=:dot)

for (ind, ln) in enumerate(line_positions)
    vlines!(ax, ln, ymin=(ind-1)*seq_length/Npop, ymax=(ind+1)*seq_length/Npop, linewidth=linewidth, color=:firebrick)
    hlines!(ax, ln, xmin=(ind-1)*seq_length/Npop, xmax=(ind+1)*seq_length/Npop, linewidth=linewidth, color=:firebrick)
end

Colorbar(g_struct[1, 2], ht, size=colorbarsize, label=L"\text{mean coupling strength (pF)}", height=Relative(.6), labelsize=labelsize, ticklabelsize=ticklabelsize)

########################################################
################ --- I₂-to-E weights --- ###############
########################################################

cl1 = ColorScheme(range(colorant"gray5", colorant"gray80", length=100))
cl_inh = ColorScheme(range(colorant"gray80", colorant"firebrick", length=100))
inhibition_cs = vcat(get(cl1, LinRange(0, 1, 100)), get(cl_inh, LinRange(0, 1, 100)))
    
ax = Axis(g_struct[1, 3], xlabel=L"\text{presynaptic }I_2 \text{-assembly}", ylabel=L"\text{postsynaptic }E \text{-assembly}", xlabelsize=labelsize, ylabelsize=labelsize,
            xticks=(1:seq_length:Npop, [L"%$(x)" for x=1:seq_length:Npop]), xticklabelsize=ticklabelsize, xgridvisible=false,
            yticks=(1:seq_length:Npop, [L"%$(x)" for x=1:seq_length:Npop]), yticklabelsize=ticklabelsize, ygridvisible=false,
            aspect=1)
ht = heatmap!(ax, weightsIE, colormap=inhibition_cs)
hlines!(line_positions, linewidth=linewidth, color=:dodgerblue4, linestyle=:dot)
vlines!(line_positions, linewidth=linewidth, color=:dodgerblue4, linestyle=:dot)

for (ind, ln) in enumerate(line_positions)
    vlines!(ln, ymin=(ind-1)*seq_length/Npop, ymax=(ind+1)*seq_length/Npop, linewidth=linewidth, color=:dodgerblue4)
    hlines!(ln, xmin=(ind-1)*seq_length/Npop, xmax=(ind+1)*seq_length/Npop, linewidth=linewidth, color=:dodgerblue4)
end

Colorbar(g_struct[1, 4], ht, size=colorbarsize, label=L"\text{mean coupling strength (pF)}", height=Relative(.6), labelsize=labelsize, ticklabelsize=ticklabelsize)

########################################################
################ --- I₂-to-E weights --- ###############
########################################################

cl1 = ColorScheme(range(colorant"gray5", colorant"gray80", length=100))
cl_inh = ColorScheme(range(colorant"gray80", colorant"dodgerblue2", length=100))
inhibition_cs = vcat(get(cl1, LinRange(0, 1, 100)), get(cl_inh, LinRange(0, 1, 100)))

ax = Axis(g_struct[1, 5], xlabel=L"\text{presynaptic }E \text{-assembly}", ylabel=L"\text{postsynaptic }I_2 \text{-assembly}", xlabelsize=labelsize, ylabelsize=labelsize,
            xticks=(1:seq_length:Npop, [L"%$(x)" for x=1:seq_length:Npop]), xticklabelsize=ticklabelsize, xgridvisible=false,
            yticks=(1:seq_length:Npop, [L"%$(x)" for x=1:seq_length:Npop]), yticklabelsize=ticklabelsize, ygridvisible=false,
            aspect=1)
ht = heatmap!(ax, weightsEI, colormap=inhibition_cs)
hlines!(line_positions, linewidth=linewidth, color=:firebrick, linestyle=:dot)
vlines!(line_positions, linewidth=linewidth, color=:firebrick, linestyle=:dot)

for (ind, ln) in enumerate(line_positions)
    vlines!(ln, ymin=(ind-1)*seq_length/Npop, ymax=(ind+1)*seq_length/Npop, linewidth=linewidth, color=:firebrick)
    hlines!(ln, xmin=(ind-1)*seq_length/Npop, xmax=(ind+1)*seq_length/Npop, linewidth=linewidth, color=:firebrick)
end

Colorbar(g_struct[1, 6], ht, size=colorbarsize, label=L"\text{mean coupling strength (pF)}", height=Relative(.6), labelsize=labelsize, ticklabelsize=ticklabelsize)

##############################################################################################
############################# --- Structure during learning --- ##############################
##############################################################################################
# --- Load Data ---
sim_name = string("network_1.h5")
sim_savedpath = string("./networks_trained/")
fid = h5open(joinpath(sim_savedpath, sim_name), "r")
popmembers = read(fid["data"]["popmembers"])
weights = read(fid["data"]["weights"])
weightsEE = read(fid["data"]["weightsEE"])
weightsEI = read(fid["data"]["weightsEI"])
weightsIE = read(fid["data"]["weightsIE"])
close(fid)

ipopmembers = findI2populations(weights, popmembers, iipop_len=Ni_members)

########################################################
################ --- E-to-E weights --- ################
########################################################
ax_ee = Axis(g_learn[1, 2], #xlabel=L"\text{training time (%)}", ylabel=L"\text{mean synaptic strength (pF)}",
            xticks=([0, 125, 250, 375, 500], [L"0", L"25", L"50", L"75", L"100"]),
            yticks=([5, 10, 15], [L"5", L"10", L"15"]),
            xlabelsize=labelsize, ylabelsize=labelsize, xticklabelsize=ticklabelsize, yticklabelsize=ticklabelsize,
            limits=(-2.5, 505.5, 1., 19.5))

for ipop = 1:Npop
    lines!(ax_ee, weightsEE[ipop, ipop, 1:500], linewidth=linewidth, color=ColorSchemes.Blues[7], label=L"\text{E-assembly}")
end
axislegend(ax_ee, position=(.9, .7), labelsize=labelsize, titlesize=textsize, merge=true)

########################################################
################ --- E-to-I₂ weights --- ###############
########################################################
ax_ei = Axis(g_learn[1, 3], #xlabel=L"\text{training time (%)}", ylabel=L"\text{mean synaptic strength (pF)}",
            xticks=([0, 125, 250, 375, 500], [L"0", L"25", L"50", L"75", L"100"]),
            yticks=([1.5, 2.5, 3.5], [L"1.5", L"2.5", L"3.5"]),
            xlabelsize=labelsize, ylabelsize=labelsize, xticklabelsize=ticklabelsize, yticklabelsize=ticklabelsize,
            limits=(-2.5, 505.5, 1.1, 3.9))
for ipop = 1:Npop
    lines!(ax_ei, mean(weightsEI[ipopmembers[:, ipop] .- (Ncells-Ni2), ipop, 1:500], dims=1)[:], linewidth=linewidth, color=ColorSchemes.Blues[5], label=L"E\text{-to-}I_2")
end
axislegend(ax_ei, position=(.1, .9), labelsize=labelsize, titlesize=textsize, merge=true)

########################################################
################ --- I₁-to-E weights --- ###############
########################################################
ax_i1 = Axis(g_learn[2, 2], xlabel=L"\text{training time (%)}", #ylabel=L"\text{mean synaptic strength (pF)}",
            xticks=([0, 125, 250, 375, 500], [L"0", L"25", L"50", L"75", L"100"]),
            yticks=([60, 80, 100], [L"60", L"80", L"100"]),
            xlabelsize=labelsize, ylabelsize=labelsize, xticklabelsize=ticklabelsize, yticklabelsize=ticklabelsize,
            limits=(-2.5, 505.5, 48.1, 115))
for ipop = 2:Npop
    (mod(ipop, seq_length) == 1) && (continue)
    lines!(ax_i1, mean(weightsIE[1:(Ncells-Ne-Ni2), ipop, 1:500], dims=1)[:], linewidth=linewidth, color=ColorSchemes.Reds[9], label=L"I_1\text{-to-}E")
end
axislegend(ax_i1, position=(.9, .3), labelsize=labelsize, titlesize=textsize, merge=true)

########################################################
################ --- I₂-to-E weights --- ###############
########################################################
ax_i2 = Axis(g_learn[2, 3], xlabel=L"\text{training time (%)}", #ylabel=L"\text{mean synaptic strength (pF)}",
            xticks=([0, 125, 250, 375, 500], [L"0", L"25", L"50", L"75", L"100"]),
            yticks=([30, 60, 90], [L"30", L"60", L"90"]),
            xlabelsize=labelsize, ylabelsize=labelsize, xticklabelsize=ticklabelsize, yticklabelsize=ticklabelsize,
            limits=(-2.5, 505.5, 15.1, 106.))
for ipop = 1:Npop
    lines!(ax_i2, mean(weightsIE[ipopmembers[:, ipop] .- (Ne), ipop, 1:500], dims=1)[:], linewidth=linewidth, color=:gray70, label=L"I_2-\text{to}-E \text{ recurrent}")
end

for ipop = 2:Npop
    (mod(ipop, seq_length) == 1) && (continue)
    lines!(ax_i2, mean(weightsIE[ipopmembers[:, ipop] .- (Ne), ipop-1, 1:500], dims=1)[:], linewidth=linewidth, color=:black, label=L"I_2\text{-to-}E \text{ preceding}")
end

for ipop = 2:Npop
    (mod(ipop, seq_length) == 1) && (continue)
    lines!(ax_i2, mean(weightsIE[ipopmembers[:, ipop-1] .- (Ne), ipop, 1:500], dims=1)[:], linewidth=linewidth, color=ColorSchemes.Reds[6], label=L"I_2\text{-to-}E \text{ succeeding}")
end
axislegend(ax_i2, position=(.1, .9), labelsize=labelsize, titlesize=textsize, merge=true)

linkxaxes!(ax_ee, ax_i1); linkxaxes!(ax_ei, ax_i2)


Label(g_learn[1:2, 1], L"\text{mean synaptic strength (pF)}",
    fontsize=labelsize, rotation=π/2)

##############################################################################################
############################### --- Structure during rest --- ################################
##############################################################################################
# --- Load Data ---
sim_name = string("network_permit_rest.h5")
sim_savedpath = string("./networks_trained/")
fid = h5open(joinpath(sim_savedpath, sim_name), "r")
popmembers = read(fid["data"]["popmembers"])
weights = read(fid["data"]["weights"])
weightsEE = read(fid["data"]["weightsEE"])
weightsEI = read(fid["data"]["weightsEI"])
weightsIE = read(fid["data"]["weightsIE"])
close(fid)

ipopmembers = findI2populations(weights, popmembers, iipop_len=Ni_members)

########################################################
################ --- E-to-E weights --- ################
########################################################
ax_ee = Axis(g_rest[1, 2], #xlabel=L"\text{training time (%)}", ylabel=L"\text{mean synaptic strength (pF)}",
            xticks=([0, 125, 250, 375, 500], [L"0", L"25", L"50", L"75", L"100"]),
            # yticks=([5, 10, 15], [L"5", L"10", L"15"]),
            xlabelsize=labelsize, ylabelsize=labelsize, xticklabelsize=ticklabelsize, yticklabelsize=ticklabelsize)#,
            # limits=(-2.5, 505.5, 1., 19.5))

for ipop = 1:Npop
    lines!(ax_ee, weightsEE[ipop, ipop, 500:end], linewidth=linewidth, color=ColorSchemes.Blues[7], label=L"\text{E-assembly}")
end
axislegend(ax_ee, position=(.9, .7), labelsize=labelsize, titlesize=textsize, merge=true)

########################################################
################ --- E-to-I₂ weights --- ###############
########################################################
ax_ei = Axis(g_rest[1, 3], #xlabel=L"\text{training time (%)}", ylabel=L"\text{mean synaptic strength (pF)}",
            xticks=([0, 125, 250, 375, 500], [L"0", L"25", L"50", L"75", L"100"]),
            # yticks=([1.5, 2.5, 3.5], [L"1.5", L"2.5", L"3.5"]),
            xlabelsize=labelsize, ylabelsize=labelsize, xticklabelsize=ticklabelsize, yticklabelsize=ticklabelsize)#,
            # limits=(-2.5, 505.5, 1.1, 3.9))
for ipop = 1:Npop
    lines!(ax_ei, mean(weightsEI[ipopmembers[:, ipop] .- (Ncells-Ni2), ipop, 500:end], dims=1)[:], linewidth=linewidth, color=ColorSchemes.Blues[5], label=L"E\text{-to-}I_2")
end
axislegend(ax_ei, position=(.1, .9), labelsize=labelsize, titlesize=textsize, merge=true)

########################################################
################ --- I₁-to-E weights --- ###############
########################################################
ax_i1 = Axis(g_rest[2, 2], xlabel=L"\text{training time (%)}", #ylabel=L"\text{mean synaptic strength (pF)}",
            xticks=([0, 125, 250, 375, 500], [L"0", L"25", L"50", L"75", L"100"]),
            # yticks=([60, 80, 100], [L"60", L"80", L"100"]),
            xlabelsize=labelsize, ylabelsize=labelsize, xticklabelsize=ticklabelsize, yticklabelsize=ticklabelsize)#,
            # limits=(-2.5, 505.5, 48.1, 115))
for ipop = 2:Npop
    (mod(ipop, seq_length) == 1) && (continue)
    lines!(ax_i1, mean(weightsIE[1:(Ncells-Ne-Ni2), ipop, 500:end], dims=1)[:], linewidth=linewidth, color=ColorSchemes.Reds[9], label=L"I_1\text{-to-}E")
end
axislegend(ax_i1, position=(.9, .3), labelsize=labelsize, titlesize=textsize, merge=true)

########################################################
################ --- I₂-to-E weights --- ###############
########################################################
ax_i2 = Axis(g_rest[2, 3], xlabel=L"\text{training time (%)}", #ylabel=L"\text{mean synaptic strength (pF)}",
            xticks=([0, 125, 250, 375, 500], [L"0", L"25", L"50", L"75", L"100"]),
            # yticks=([30, 60, 90], [L"30", L"60", L"90"]),
            xlabelsize=labelsize, ylabelsize=labelsize, xticklabelsize=ticklabelsize, yticklabelsize=ticklabelsize)#),
            # limits=(-2.5, 505.5, 15.1, 106.))
for ipop = 1:Npop
    lines!(ax_i2, mean(weightsIE[ipopmembers[:, ipop] .- (Ne), ipop, 500:end], dims=1)[:], linewidth=linewidth, color=:gray70, label=L"I_2-\text{to}-E \text{ recurrent}")
end

for ipop = 2:Npop
    (mod(ipop, seq_length) == 1) && (continue)
    lines!(ax_i2, mean(weightsIE[ipopmembers[:, ipop] .- (Ne), ipop-1, 500:end], dims=1)[:], linewidth=linewidth, color=:black, label=L"I_2\text{-to-}E \text{ preceding}")
end

for ipop = 2:Npop
    (mod(ipop, seq_length) == 1) && (continue)
    lines!(ax_i2, mean(weightsIE[ipopmembers[:, ipop-1] .- (Ne), ipop, 500:end], dims=1)[:], linewidth=linewidth, color=ColorSchemes.Reds[6], label=L"I_2\text{-to-}E \text{ succeeding}")
end
axislegend(ax_i2, position=(.1, .9), labelsize=labelsize, titlesize=textsize, merge=true)

linkxaxes!(ax_ee, ax_i1); linkxaxes!(ax_ei, ax_i2)


Label(g_rest[1:2, 1], L"\text{mean synaptic strength (pF)}",
    fontsize=labelsize, rotation=π/2)
#

Label(g_model[1, 0, TopLeft()], L"\textbf{a}",
    fontsize=subplotlabelsize,
    font=:bold,
    padding=(0, 5, 5, 0),
    halign=:right)

Label(g_hist[1, 0, TopLeft()], L"\textbf{b}",
    fontsize=subplotlabelsize,
    font=:bold,
    padding=(0, 5, 5, 0),
    halign=:right)

Label(g_struct[1, 0, TopLeft()], L"\textbf{c}",
    fontsize=subplotlabelsize,
    font=:bold,
    padding=(0, 5, 5, 0),
    halign=:right)

Label(g_learn[1, 0, TopLeft()], L"\textbf{d}",
    fontsize=subplotlabelsize,
    font=:bold,
    padding=(0, 5, 5, 0),
    halign=:right)

Label(g_rest[1, 0, TopLeft()], L"\textbf{e}",
    fontsize=subplotlabelsize,
    font=:bold,
    padding=(0, 5, 5, 0),
    halign=:right)
    
colsize!(g_model, 0, Relative(0))
colsize!(g_hist, 0, Relative(0))
colsize!(g_struct, 0, Relative(0))
colsize!(g_learn, 0, Relative(0))
colsize!(g_rest, 0, Relative(0))

# colgap!(g_struct, 2, Relative(.1))
colgap!(g_struct, 1, Relative(0))
colgap!(g_struct, 3, Relative(0))
colgap!(g_struct, 5, Relative(0))

rowgap!(fig.layout, 1, Relative(.05))
rowgap!(fig.layout, 2, Relative(.05))
colgap!(fig.layout, 2, Relative(.05))

rowsize!(fig.layout, 1, Relative(.3))
rowsize!(fig.layout, 2, Relative(.25))


(!isdir(output_dir)) && (mkpath(output_dir))
save(joinpath(output_dir, string("figure_2.png")), fig)
end


##############################################################################################
######                               --- FIGURE 3 ---                                   ###### 
##############################################################################################

if makefig_3

fig = Figure(size=(4800, 3800))

g_i1_act = fig[1, 1] = GridLayout()
g_i1_anal = fig[1, 2:3] = GridLayout()
g_i2_act = fig[2, 1] = GridLayout()
g_i2_anal = fig[2, 2:3] = GridLayout()
g_ei_act = fig[3, 1] = GridLayout()
g_ei_anal = fig[3, 2:3] = GridLayout()


labelsize = 66
ticklabelsize = 54
legendlabelsize = 54
# legendticklabelsize = 54
# colorbarlabelsize = 66
# colorbarticklabelsize = 54
textsize = 66
subplotlabelsize = 66

dodge = 1
dodge_gap = .1
v_width = 1.
b_width = 1.
b_strokewidth = 1.

linewidth = 6
markersize = 5

colorbarsize = 40

##############################################################################################
#############################  --- iSTDP₁ knockout activity --- ##############################
##############################################################################################
# --- Load Data ---
mode = "spontaneous"
sim_name = string("network_i1STDP-knockout_spontaneous_seed_", seed_num, ".h5")
sim_savedpath = string("./networks_trained_", mode, "/network_i1STDP-knockout/")
fid = h5open(joinpath(sim_savedpath, sim_name), "r")
popmembers = read(fid["data"]["popmembers"])
weights = read(fid["data"]["weights"])
times = read(fid["data"]["times"])
close(fid)

# Parameters
Ne = 4000
Ni2 = 400
seq_length = 4
Ni_members = 25

Ncells = size(times)[1]
Npop = size(popmembers, 2)
Nmaxmembers = size(popmembers, 1)
interval = collect(12_500:17_500)

restcells = deleteat!(map(*, ones(Int, Ne), range(1, stop=Ne)), sort(unique(popmembers))[2:end])
ylim_max = count(i->(i>0), popmembers) + length(restcells) + (Ncells-Ni2-Ne) + length(ipopmembers)

ns = zeros(Ncells)
for cc = 1:Ncells
    ns[cc] = count(i->i>0, times[cc, :])
end

# Plot params
ytick_seq = zeros(round(Int, Npop/seq_length)*2)
x_ticks = minimum(interval):1000:maximum(interval)
rowcount = 1
ytickcount = 1

ipopmembers = findI2populations(weights, popmembers, iipop_len=Ni_members, Ni2=Ni2)
restcells = deleteat!(map(*, ones(Int, Ne), range(1, stop=Ne)), sort(unique(popmembers))[2:end])
ylim_max = count(i->(i>0), popmembers) + length(restcells) + (Ncells-Ni2-Ne) + length(ipopmembers)

# rates = convolveSpikes(times, interval=interval, gaussian=false, sigma=10.)
rates = binRates(times, interval=interval)
emembers = unique(filter(i->i>0, popmembers))

ns = zeros(Ncells)
for cc = 1:Ncells
    ns[cc] = count(i->i>0, times[cc, :])
end

ax = Axis(g_i1_act[2, 1], xlabel=L"\text{simulation time (s)}", ylabel=L"\text{sequences}", xlabelsize=labelsize, ylabelsize=labelsize,
            xticks=(x_ticks, [L"%$(x)" for x in 1:length(x_ticks)]), xticklabelsize=ticklabelsize, xgridvisible=false,
            yticklabelsize=ticklabelsize, ygridvisible=false,
            limits=(minimum(interval), maximum(interval), 1, ylim_max))
# Excitatoy assembly members
for pp = 1:Npop
    for cc = 1:Nmaxmembers
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
    if mod(pp, seq_length) == 0
        ytick_seq[ytickcount] = rowcount - (Nmaxmembers*seq_length/2)
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
    if mod(pp, seq_length) == 0
        ytick_seq[ytickcount] = rowcount - (Ni_members*seq_length/2)
        ytickcount += 1
    end
end

ytick_labels = [Char(i+64) for i in 1:(round(Int, Npop/seq_length))]
ytick_labels = vcat(ytick_labels, ytick_labels)
ax.yticks = (ytick_seq, [L"\text{%$(x)}" for x in ytick_labels])


ax_top = Axis(g_i1_act[1, 1], xlabel="", ylabel="", xlabelsize=labelsize, ylabelsize=labelsize,
                # title=L"\textbf{after stimulation}",  titlesize=labelsize,
                xticks=([], []), xticklabelsize=ticklabelsize,
                yticks=([0, 1, 2], [L"0", L"1", L"2"]), yticklabelsize=ticklabelsize, yaxisposition=:right,
                limits=(minimum(interval), maximum(interval), -0.2, 10.1))
# Plot mean firing rates
lines!(ax_top, interval, vec(mean(rates[emembers, :], dims=1)), color=ColorSchemes.Blues[7], linewidth=linewidth)
lines!(ax_top, interval, vec(mean(rates[restcells, :], dims=1)), color=ColorSchemes.Greys[7], linewidth=linewidth)
lines!(ax_top, interval, vec(mean(rates[(Ne+1):(Ncells-Ni2), :], dims=1)), color=ColorSchemes.Reds[9], linewidth=linewidth)
lines!(ax_top, interval, vec(mean(rates[(Ncells-Ni2+1):Ncells, :], dims=1)), color=ColorSchemes.Reds[6], linewidth=linewidth)

linkxaxes!(ax, ax_top); #xlims!(ax, low=minimum(interval)+500, high=maximum(interval))
hidedecorations!(ax, grid=true, ticks=false, ticklabels=false, label=false)
hideydecorations!(ax, grid=true, ticks=true, ticklabels=false, label=false)
hideydecorations!(ax_top, grid=true, ticks=true, ticklabels=true, label=true)
rowsize!(g_i1_act, 1, Relative(1/5))    


##############################################################################################
#############################  --- iSTDP₁ knockout analysis --- ##############################
##############################################################################################

# --- Load Data ---
sim_name = string("network_i1STDP-knockout.h5")
sim_savedpath = string("./networks_trained/")
fid = h5open(joinpath(sim_savedpath, sim_name), "r")
popmembers = read(fid["data"]["popmembers"])
weights = read(fid["data"]["weights"])
times = read(fid["data"]["times"])
weightsEE_tr = read(fid["data"]["weightsEE"])
weightsEI_tr = read(fid["data"]["weightsEI"])
weightsIE_tr = read(fid["data"]["weightsIE"])
close(fid)


weightsEE = sum(weights[popmembers[:, :], popmembers[:, :]], dims=(1, 3))[1, :, 1, :] ./ count(i->i>0, weights[popmembers[:, :], popmembers[:, :]], dims=(1, 3))[1, :, 1, :]
weightsEI = sum(weights[popmembers[:, :], ipopmembers[:, :]], dims=(1, 3))[1, :, 1, :] ./ count(i->i>0, weights[popmembers[:, :], ipopmembers[:, :]], dims=(1, 3))[1, :, 1, :]
weightsIE = sum(weights[ipopmembers[:, :], popmembers[:, :]], dims=(1, 3))[1, :, 1, :] ./ count(i->i>0, weights[ipopmembers[:, :], popmembers[:, :]], dims=(1, 3))[1, :, 1, :]


########################################################
################ --- E-to-E weights --- ################
########################################################

cl1 = ColorScheme(range(colorant"gray5", colorant"gray80", length=100))
cl_exc = ColorScheme(range(colorant"gray80", colorant"dodgerblue4", length=100))
excitation_cs = vcat(get(cl1, LinRange(0, 1, 100)), get(cl_exc, LinRange(0, 1, 100)))

Npop = size(weightsEE)[1]
line_positions = collect(seq_length+.5:seq_length:Npop+.5)

ax = Axis(g_i1_anal[1:3, 1:2], xlabel=L"\text{pre }E\text{-assembly}", ylabel=L"\text{post }E\text{-assembly}", xlabelsize=labelsize, ylabelsize=labelsize,
            xticks=(1:seq_length:Npop, [L"%$(x)" for x=1:seq_length:Npop]), xticklabelsize=ticklabelsize, xgridvisible=false,
            yticks=(1:seq_length:Npop, [L"%$(x)" for x=1:seq_length:Npop]), yticklabelsize=ticklabelsize, ygridvisible=false,
            aspect=1)
ht = heatmap!(ax, weightsEE, colormap=excitation_cs)
hlines!(ax, line_positions, linewidth=linewidth, color=:firebrick, linestyle=:dot)
vlines!(ax, line_positions, linewidth=linewidth, color=:firebrick, linestyle=:dot)

for (ind, ln) in enumerate(line_positions)
    vlines!(ax, ln, ymin=(ind-1)*seq_length/Npop, ymax=(ind+1)*seq_length/Npop, linewidth=linewidth, color=:firebrick)
    hlines!(ax, ln, xmin=(ind-1)*seq_length/Npop, xmax=(ind+1)*seq_length/Npop, linewidth=linewidth, color=:firebrick)
end

Colorbar(g_i1_anal[1:3, 3], ht, size=colorbarsize, label=L"\text{mean syn. str. (pF)}", height=Relative(.6), labelsize=labelsize, ticklabelsize=ticklabelsize)

########################################################
################ --- I₂-to-E weights --- ###############
########################################################

cl1 = ColorScheme(range(colorant"gray5", colorant"gray80", length=100))
cl_inh = ColorScheme(range(colorant"gray80", colorant"firebrick", length=100))
inhibition_cs = vcat(get(cl1, LinRange(0, 1, 100)), get(cl_inh, LinRange(0, 1, 100)))
    
ax = Axis(g_i1_anal[4:6, 1:2], xlabel=L"\text{pre }I_2 \text{-assembly}", ylabel=L"\text{post }E \text{-assembly}", xlabelsize=labelsize, ylabelsize=labelsize,
            xticks=(1:seq_length:Npop, [L"%$(x)" for x=1:seq_length:Npop]), xticklabelsize=ticklabelsize, xgridvisible=false,
            yticks=(1:seq_length:Npop, [L"%$(x)" for x=1:seq_length:Npop]), yticklabelsize=ticklabelsize, ygridvisible=false,
            aspect=1)
ht = heatmap!(ax, weightsIE, colormap=inhibition_cs)
hlines!(line_positions, linewidth=linewidth, color=:dodgerblue4, linestyle=:dot)
vlines!(line_positions, linewidth=linewidth, color=:dodgerblue4, linestyle=:dot)

for (ind, ln) in enumerate(line_positions)
    vlines!(ln, ymin=(ind-1)*seq_length/Npop, ymax=(ind+1)*seq_length/Npop, linewidth=linewidth, color=:dodgerblue4)
    hlines!(ln, xmin=(ind-1)*seq_length/Npop, xmax=(ind+1)*seq_length/Npop, linewidth=linewidth, color=:dodgerblue4)
end

Colorbar(g_i1_anal[4:6, 3], ht, size=colorbarsize, label=L"\text{mean syn. str. (pF)}", height=Relative(.6), labelsize=labelsize, ticklabelsize=ticklabelsize)

########################################################
################ --- I₂-to-E weights --- ###############
########################################################

cl1 = ColorScheme(range(colorant"gray5", colorant"gray80", length=100))
cl_inh = ColorScheme(range(colorant"gray80", colorant"dodgerblue2", length=100))
inhibition_cs = vcat(get(cl1, LinRange(0, 1, 100)), get(cl_inh, LinRange(0, 1, 100)))

ax = Axis(g_i1_anal[7:9, 1:2], xlabel=L"\text{pre }E \text{-assembly}", ylabel=L"\text{post }I_2 \text{-assembly}", xlabelsize=labelsize, ylabelsize=labelsize,
            xticks=(1:seq_length:Npop, [L"%$(x)" for x=1:seq_length:Npop]), xticklabelsize=ticklabelsize, xgridvisible=false,
            yticks=(1:seq_length:Npop, [L"%$(x)" for x=1:seq_length:Npop]), yticklabelsize=ticklabelsize, ygridvisible=false,
            aspect=1)
ht = heatmap!(ax, weightsEI, colormap=inhibition_cs)
hlines!(line_positions, linewidth=linewidth, color=:firebrick, linestyle=:dot)
vlines!(line_positions, linewidth=linewidth, color=:firebrick, linestyle=:dot)

for (ind, ln) in enumerate(line_positions)
    vlines!(ln, ymin=(ind-1)*seq_length/Npop, ymax=(ind+1)*seq_length/Npop, linewidth=linewidth, color=:firebrick)
    hlines!(ln, xmin=(ind-1)*seq_length/Npop, xmax=(ind+1)*seq_length/Npop, linewidth=linewidth, color=:firebrick)
end

Colorbar(g_i1_anal[7:9, 3], ht, size=colorbarsize, label=L"\text{mean syn. str. (pF)}", height=Relative(.6), labelsize=labelsize, ticklabelsize=ticklabelsize)

########################################################
################ --- Extra plot here --- ###############
########################################################

########################################################
################ --- E-to-E weights --- ################
########################################################
ax_ee = Axis(g_i1_anal[1:3, 4:5], #xlabel=L"\text{training time (%)}", ylabel=L"\text{mean synaptic strength (pF)}",
            # xticks=([0, 125, 250, 375, 500], [L"0", L"25", L"50", L"75", L"100"]),
            yticks=([5, 10, 15], [L"5", L"10", L"15"]),
            xlabelsize=labelsize, ylabelsize=labelsize, xticklabelsize=ticklabelsize, yticklabelsize=ticklabelsize)#,
            # limits=(-2.5, 505.5, 1., 19.5))
for ipop = 1:Npop
    lines!(ax_ee, weightsEE_tr[ipop, ipop, 1:500], linewidth=linewidth, color=ColorSchemes.Blues[7], label=L"\text{E-assembly}")
end
# axislegend(ax_ee, position=(.9, .7), labelsize=labelsize, titlesize=textsize, merge=true)

########################################################
################ --- E-to-I₂ weights --- ###############
########################################################
ax_ei = Axis(g_i1_anal[4:6, 4:5], #xlabel=L"\text{training time (%)}", 
            ylabel=L"\text{mean synaptic strength (pF)}",
            # xticks=([0, 125, 250, 375, 500], [L"0", L"25", L"50", L"75", L"100"]),
            yticks=([1.5, 2.5, 3.5], [L"1.5", L"2.5", L"3.5"]),
            xlabelsize=labelsize, ylabelsize=labelsize, xticklabelsize=ticklabelsize, yticklabelsize=ticklabelsize)#,
            # limits=(-2.5, 505.5, 1.1, 3.9))
for ipop = 1:Npop
    lines!(ax_ei, mean(weightsEI_tr[ipopmembers[:, ipop] .- (Ncells-Ni2), ipop, 1:500], dims=1)[:], linewidth=linewidth, color=ColorSchemes.Blues[5], label=L"E\text{-to-}I_2")
end
# axislegend(ax_ei, position=(.1, .9), labelsize=labelsize, titlesize=textsize, merge=true)

########################################################
################ --- I₂-to-E weights --- ###############
########################################################
ax_i2 = Axis(g_i1_anal[7:9, 4:5], xlabel=L"\text{training time (%)}", #ylabel=L"\text{mean synaptic strength (pF)}",
            xticks=([0, 125, 250, 375, 500], [L"0", L"25", L"50", L"75", L"100"]),
            yticks=([30, 60, 90], [L"30", L"60", L"90"]),
            xlabelsize=labelsize, ylabelsize=labelsize, xticklabelsize=ticklabelsize, yticklabelsize=ticklabelsize)#,
            # limits=(-2.5, 505.5, 15.1, 106.))
for ipop = 1:Npop
    lines!(ax_i2, mean(weightsIE_tr[ipopmembers[:, ipop] .- (Ne), ipop, 1:500], dims=1)[:], linewidth=linewidth, color=:gray70, label=L"I_2-\text{to}-E \text{ recurrent}")
end

for ipop = 2:Npop
    (mod(ipop, seq_length) == 1) && (continue)
    lines!(ax_i2, mean(weightsIE_tr[ipopmembers[:, ipop] .- (Ne), ipop-1, 1:500], dims=1)[:], linewidth=linewidth, color=:black, label=L"I_2\text{-to-}E \text{ preceding}")
end

for ipop = 2:Npop
    (mod(ipop, seq_length) == 1) && (continue)
    lines!(ax_i2, mean(weightsIE_tr[ipopmembers[:, ipop-1] .- (Ne), ipop, 1:500], dims=1)[:], linewidth=linewidth, color=ColorSchemes.Reds[6], label=L"I_2\text{-to-}E \text{ succeeding}")
end
# axislegend(ax_i2, position=(.1, .9), labelsize=labelsize, titlesize=textsize, merge=true)

linkxaxes!(ax_ee, ax_ei, ax_i2)#; linkxaxes!(ax_ei, ax_i2); linkxaxes!(ax_ee, ax_i2)
hidexdecorations!(ax_ee, grid=false, ticks=true, ticklabels=true, label=true)
hidexdecorations!(ax_ei, grid=false, ticks=true, ticklabels=true, label=true)

# Label(g_i1_anal[4:5, 4], L"\text{mean synaptic strength (pF)}",
#     fontsize=labelsize, rotation=π/2)
#


##############################################################################################
#############################  --- iSTDP₂ knockout activity --- ##############################
##############################################################################################
# --- Load Data ---
mode = "spontaneous"
sim_name = string("network_i2STDP-knockout_spontaneous_seed_", seed_num, ".h5")
sim_savedpath = string("./networks_trained_", mode, "/network_i2STDP-knockout/")
fid = h5open(joinpath(sim_savedpath, sim_name), "r")
popmembers = read(fid["data"]["popmembers"])
weights = read(fid["data"]["weights"])
times = read(fid["data"]["times"])
close(fid)

# Parameters
Ne = 4000
Ni2 = 400
seq_length = 4
Ni_members = 25

Ncells = size(times)[1]
Npop = size(popmembers, 2)
Nmaxmembers = size(popmembers, 1)
interval = collect(12_500:17_500)

restcells = deleteat!(map(*, ones(Int, Ne), range(1, stop=Ne)), sort(unique(popmembers))[2:end])
ylim_max = count(i->(i>0), popmembers) + length(restcells) + (Ncells-Ni2-Ne) + length(ipopmembers)

ns = zeros(Ncells)
for cc = 1:Ncells
    ns[cc] = count(i->i>0, times[cc, :])
end

# Plot params
ytick_seq = zeros(round(Int, Npop/seq_length)*2)
x_ticks = minimum(interval):1000:maximum(interval)
rowcount = 1
ytickcount = 1

ipopmembers = findI2populations(weights, popmembers, iipop_len=Ni_members, Ni2=Ni2)
restcells = deleteat!(map(*, ones(Int, Ne), range(1, stop=Ne)), sort(unique(popmembers))[2:end])
ylim_max = count(i->(i>0), popmembers) + length(restcells) + (Ncells-Ni2-Ne) + length(ipopmembers)

# rates = convolveSpikes(times, interval=interval, gaussian=false, sigma=10.)
rates = binRates(times, interval=interval)
emembers = unique(filter(i->i>0, popmembers))

ns = zeros(Ncells)
for cc = 1:Ncells
    ns[cc] = count(i->i>0, times[cc, :])
end

ax = Axis(g_i2_act[2, 1], xlabel=L"\text{simulation time (s)}", ylabel=L"\text{sequences}", xlabelsize=labelsize, ylabelsize=labelsize,
            xticks=(x_ticks, [L"%$(x)" for x in 1:length(x_ticks)]), xticklabelsize=ticklabelsize, xgridvisible=false,
            yticklabelsize=ticklabelsize, ygridvisible=false,
            limits=(minimum(interval), maximum(interval), 1, ylim_max))
# Excitatoy assembly members
for pp = 1:Npop
    for cc = 1:Nmaxmembers
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
    if mod(pp, seq_length) == 0
        ytick_seq[ytickcount] = rowcount - (Nmaxmembers*seq_length/2)
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
    if mod(pp, seq_length) == 0
        ytick_seq[ytickcount] = rowcount - (Ni_members*seq_length/2)
        ytickcount += 1
    end
end

ytick_labels = [Char(i+64) for i in 1:(round(Int, Npop/seq_length))]
ytick_labels = vcat(ytick_labels, ytick_labels)
ax.yticks = (ytick_seq, [L"\text{%$(x)}" for x in ytick_labels])


ax_top = Axis(g_i2_act[1, 1], xlabel="", ylabel="", xlabelsize=labelsize, ylabelsize=labelsize,
                # title=L"\textbf{after stimulation}",  titlesize=labelsize,
                xticks=([], []), xticklabelsize=ticklabelsize,
                yticks=([0, 1, 2], [L"0", L"1", L"2"]), yticklabelsize=ticklabelsize, yaxisposition=:right,
                limits=(minimum(interval), maximum(interval), -0.2, 10.1))
# Plot mean firing rates
lines!(ax_top, interval, vec(mean(rates[emembers, :], dims=1)), color=ColorSchemes.Blues[7], linewidth=linewidth)
lines!(ax_top, interval, vec(mean(rates[restcells, :], dims=1)), color=ColorSchemes.Greys[7], linewidth=linewidth)
lines!(ax_top, interval, vec(mean(rates[(Ne+1):(Ncells-Ni2), :], dims=1)), color=ColorSchemes.Reds[9], linewidth=linewidth)
lines!(ax_top, interval, vec(mean(rates[(Ncells-Ni2+1):Ncells, :], dims=1)), color=ColorSchemes.Reds[6], linewidth=linewidth)

linkxaxes!(ax, ax_top); #xlims!(ax, low=minimum(interval)+500, high=maximum(interval))
hidedecorations!(ax, grid=true, ticks=false, ticklabels=false, label=false)
hideydecorations!(ax, grid=true, ticks=true, ticklabels=false, label=false)
hideydecorations!(ax_top, grid=true, ticks=true, ticklabels=true, label=true)
rowsize!(g_i2_act, 1, Relative(1/5))    


##############################################################################################
#############################  --- iSTDP₂ knockout analysis --- ##############################
##############################################################################################
# --- Load Data ---
sim_name = string("network_i2STDP-knockout.h5")
sim_savedpath = string("./networks_trained/")
fid = h5open(joinpath(sim_savedpath, sim_name), "r")
popmembers = read(fid["data"]["popmembers"])
weights = read(fid["data"]["weights"])
times = read(fid["data"]["times"])
weightsEE_tr = read(fid["data"]["weightsEE"])
weightsEI_tr = read(fid["data"]["weightsEI"])
weightsIE_tr = read(fid["data"]["weightsIE"])
close(fid)

weightsEE = sum(weights[popmembers[:, :], popmembers[:, :]], dims=(1, 3))[1, :, 1, :] ./ count(i->i>0, weights[popmembers[:, :], popmembers[:, :]], dims=(1, 3))[1, :, 1, :]
weightsEI = sum(weights[popmembers[:, :], ipopmembers[:, :]], dims=(1, 3))[1, :, 1, :] ./ count(i->i>0, weights[popmembers[:, :], ipopmembers[:, :]], dims=(1, 3))[1, :, 1, :]
weightsIE = sum(weights[ipopmembers[:, :], popmembers[:, :]], dims=(1, 3))[1, :, 1, :] ./ count(i->i>0, weights[ipopmembers[:, :], popmembers[:, :]], dims=(1, 3))[1, :, 1, :]


########################################################
################ --- E-to-E weights --- ################
########################################################

cl1 = ColorScheme(range(colorant"gray5", colorant"gray80", length=100))
cl_exc = ColorScheme(range(colorant"gray80", colorant"dodgerblue4", length=100))
excitation_cs = vcat(get(cl1, LinRange(0, 1, 100)), get(cl_exc, LinRange(0, 1, 100)))

Npop = size(weightsEE)[1]
line_positions = collect(seq_length+.5:seq_length:Npop+.5)

ax = Axis(g_i2_anal[1:3, 1:2], xlabel=L"\text{pre }E\text{-assembly}", ylabel=L"\text{post }E\text{-assembly}", xlabelsize=labelsize, ylabelsize=labelsize,
            xticks=(1:seq_length:Npop, [L"%$(x)" for x=1:seq_length:Npop]), xticklabelsize=ticklabelsize, xgridvisible=false,
            yticks=(1:seq_length:Npop, [L"%$(x)" for x=1:seq_length:Npop]), yticklabelsize=ticklabelsize, ygridvisible=false,
            aspect=1)
ht = heatmap!(ax, weightsEE, colormap=excitation_cs)
hlines!(ax, line_positions, linewidth=linewidth, color=:firebrick, linestyle=:dot)
vlines!(ax, line_positions, linewidth=linewidth, color=:firebrick, linestyle=:dot)

for (ind, ln) in enumerate(line_positions)
    vlines!(ax, ln, ymin=(ind-1)*seq_length/Npop, ymax=(ind+1)*seq_length/Npop, linewidth=linewidth, color=:firebrick)
    hlines!(ax, ln, xmin=(ind-1)*seq_length/Npop, xmax=(ind+1)*seq_length/Npop, linewidth=linewidth, color=:firebrick)
end

Colorbar(g_i2_anal[1:3, 3], ht, size=colorbarsize, label=L"\text{mean syn. str. (pF)}", height=Relative(.6), labelsize=labelsize, ticklabelsize=ticklabelsize)

########################################################
################ --- I₂-to-E weights --- ###############
########################################################

cl1 = ColorScheme(range(colorant"gray5", colorant"gray80", length=100))
cl_inh = ColorScheme(range(colorant"gray80", colorant"dodgerblue2", length=100))
inhibition_cs = vcat(get(cl1, LinRange(0, 1, 100)), get(cl_inh, LinRange(0, 1, 100)))

ax = Axis(g_i2_anal[4:6, 1:2], xlabel=L"\text{pre }E \text{-assembly}", ylabel=L"\text{post }I_2 \text{-assembly}", xlabelsize=labelsize, ylabelsize=labelsize,
            xticks=(1:seq_length:Npop, [L"%$(x)" for x=1:seq_length:Npop]), xticklabelsize=ticklabelsize, xgridvisible=false,
            yticks=(1:seq_length:Npop, [L"%$(x)" for x=1:seq_length:Npop]), yticklabelsize=ticklabelsize, ygridvisible=false,
            aspect=1)
ht = heatmap!(ax, weightsEI, colormap=inhibition_cs)
hlines!(line_positions, linewidth=linewidth, color=:firebrick, linestyle=:dot)
vlines!(line_positions, linewidth=linewidth, color=:firebrick, linestyle=:dot)

for (ind, ln) in enumerate(line_positions)
    vlines!(ln, ymin=(ind-1)*seq_length/Npop, ymax=(ind+1)*seq_length/Npop, linewidth=linewidth, color=:firebrick)
    hlines!(ln, xmin=(ind-1)*seq_length/Npop, xmax=(ind+1)*seq_length/Npop, linewidth=linewidth, color=:firebrick)
end

Colorbar(g_i2_anal[4:6, 3], ht, size=colorbarsize, label=L"\text{mean syn. str. (pF)}", height=Relative(.6), labelsize=labelsize, ticklabelsize=ticklabelsize)

########################################################
################ --- Extra plot here --- ###############
########################################################


########################################################
################ --- E-to-E weights --- ################
########################################################
ax_ee = Axis(g_i2_anal[1:2, 4:5], #xlabel=L"\text{training time (%)}", ylabel=L"\text{mean synaptic strength (pF)}",
            # xticks=([0, 125, 250, 375, 500], [L"0", L"25", L"50", L"75", L"100"]),
            yticks=([5, 10, 15], [L"5", L"10", L"15"]),
            xlabelsize=labelsize, ylabelsize=labelsize, xticklabelsize=ticklabelsize, yticklabelsize=ticklabelsize)#,
            # limits=(-2.5, 505.5, 1., 19.5))
for ipop = 1:Npop
    lines!(ax_ee, weightsEE_tr[ipop, ipop, 1:500], linewidth=linewidth, color=ColorSchemes.Blues[7], label=L"\text{E-assembly}")
end
# axislegend(ax_ee, position=(.9, .7), labelsize=labelsize, titlesize=textsize, merge=true)

########################################################
################ --- E-to-I₂ weights --- ###############
########################################################
ax_ei = Axis(g_i2_anal[3:4, 4:5], #xlabel=L"\text{training time (%)}", ylabel=L"\text{mean synaptic strength (pF)}",
            # xticks=([0, 125, 250, 375, 500], [L"0", L"25", L"50", L"75", L"100"]),
            yticks=([1.5, 2.5, 3.5], [L"1.5", L"2.5", L"3.5"]),
            xlabelsize=labelsize, ylabelsize=labelsize, xticklabelsize=ticklabelsize, yticklabelsize=ticklabelsize)#,
            # limits=(-2.5, 505.5, 1.1, 3.9))
for ipop = 1:Npop
    lines!(ax_ei, mean(weightsEI_tr[ipopmembers[:, ipop] .- (Ncells-Ni2), ipop, 1:500], dims=1)[:], linewidth=linewidth, color=ColorSchemes.Blues[5], label=L"E\text{-to-}I_2")
end
# axislegend(ax_ei, position=(.1, .9), labelsize=labelsize, titlesize=textsize, merge=true)

########################################################
################ --- I₁-to-E weights --- ###############
########################################################
ax_i1 = Axis(g_i2_anal[5:6, 4:5], xlabel=L"\text{training time (%)}", #ylabel=L"\text{mean synaptic strength (pF)}",
            xticks=([0, 125, 250, 375, 500], [L"0", L"25", L"50", L"75", L"100"]),
            yticks=([60, 80, 100], [L"60", L"80", L"100"]),
            xlabelsize=labelsize, ylabelsize=labelsize, xticklabelsize=ticklabelsize, yticklabelsize=ticklabelsize)#,
            # limits=(-2.5, 505.5, 48.1, 115))
for ipop = 2:Npop
    (mod(ipop, seq_length) == 1) && (continue)
    lines!(ax_i1, mean(weightsIE_tr[1:(Ncells-Ne-Ni2), ipop, 1:500], dims=1)[:], linewidth=linewidth, color=ColorSchemes.Reds[9], label=L"I_1\text{-to-}E")
end
# axislegend(ax_i1, position=(.9, .3), labelsize=labelsize, titlesize=textsize, merge=true)

linkxaxes!(ax_ee, ax_i1, ax_ei)#; linkxaxes!(ax_ee, ax_ei)
hidexdecorations!(ax_ee, grid=false, ticks=true, ticklabels=true, label=true)
hidexdecorations!(ax_ei, grid=false, ticks=true, ticklabels=true, label=true)

# Label(g_i1_anal[4:5, 4], L"\text{mean synaptic strength (pF)}",
#     fontsize=labelsize, rotation=π/2)
#


##############################################################################################
#############################  --- eiSTDP knockout activity --- ##############################
##############################################################################################
# --- Load Data ---
mode = "spontaneous"
sim_name = string("network_eiSTDP-knockout_spontaneous_seed_", seed_num, ".h5")
sim_savedpath = string("./networks_trained_", mode, "/network_eiSTDP-knockout/")
fid = h5open(joinpath(sim_savedpath, sim_name), "r")
popmembers = read(fid["data"]["popmembers"])
weights = read(fid["data"]["weights"])
times = read(fid["data"]["times"])
close(fid)

# Parameters
Ne = 4000
Ni2 = 400
seq_length = 4
Ni_members = 25

Ncells = size(times)[1]
Npop = size(popmembers, 2)
Nmaxmembers = size(popmembers, 1)
interval = collect(12_500:17_500)

restcells = deleteat!(map(*, ones(Int, Ne), range(1, stop=Ne)), sort(unique(popmembers))[2:end])
ylim_max = count(i->(i>0), popmembers) + length(restcells) + (Ncells-Ni2-Ne) + length(ipopmembers)

ns = zeros(Ncells)
for cc = 1:Ncells
    ns[cc] = count(i->i>0, times[cc, :])
end

# Plot params
ytick_seq = zeros(round(Int, Npop/seq_length)*2)
x_ticks = minimum(interval):1000:maximum(interval)
rowcount = 1
ytickcount = 1

ipopmembers = findI2populations(weights, popmembers, iipop_len=Ni_members, Ni2=Ni2)
restcells = deleteat!(map(*, ones(Int, Ne), range(1, stop=Ne)), sort(unique(popmembers))[2:end])
ylim_max = count(i->(i>0), popmembers) + length(restcells) + (Ncells-Ni2-Ne) + length(ipopmembers)

# rates = convolveSpikes(times, interval=interval, gaussian=false, sigma=10.)
rates = binRates(times, interval=interval)
emembers = unique(filter(i->i>0, popmembers))

ns = zeros(Ncells)
for cc = 1:Ncells
    ns[cc] = count(i->i>0, times[cc, :])
end

ax = Axis(g_ei_act[2, 1], xlabel=L"\text{simulation time (s)}", ylabel=L"\text{sequences}", xlabelsize=labelsize, ylabelsize=labelsize,
            xticks=(x_ticks, [L"%$(x)" for x in 1:length(x_ticks)]), xticklabelsize=ticklabelsize, xgridvisible=false,
            yticklabelsize=ticklabelsize, ygridvisible=false,
            limits=(minimum(interval), maximum(interval), 1, ylim_max))
# Excitatoy assembly members
for pp = 1:Npop
    for cc = 1:Nmaxmembers
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
    if mod(pp, seq_length) == 0
        ytick_seq[ytickcount] = rowcount - (Nmaxmembers*seq_length/2)
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
    if mod(pp, seq_length) == 0
        ytick_seq[ytickcount] = rowcount - (Ni_members*seq_length/2)
        ytickcount += 1
    end
end

ytick_labels = [Char(i+64) for i in 1:(round(Int, Npop/seq_length))]
ytick_labels = vcat(ytick_labels, ytick_labels)
ax.yticks = (ytick_seq, [L"\text{%$(x)}" for x in ytick_labels])


ax_top = Axis(g_ei_act[1, 1], xlabel="", ylabel="", xlabelsize=labelsize, ylabelsize=labelsize,
                # title=L"\textbf{after stimulation}",  titlesize=labelsize,
                xticks=([], []), xticklabelsize=ticklabelsize,
                yticks=([0, 1, 2], [L"0", L"1", L"2"]), yticklabelsize=ticklabelsize, yaxisposition=:right,
                limits=(minimum(interval), maximum(interval), -0.2, 10.1))
# Plot mean firing rates
lines!(ax_top, interval, vec(mean(rates[emembers, :], dims=1)), color=ColorSchemes.Blues[7], linewidth=linewidth)
lines!(ax_top, interval, vec(mean(rates[restcells, :], dims=1)), color=ColorSchemes.Greys[7], linewidth=linewidth)
lines!(ax_top, interval, vec(mean(rates[(Ne+1):(Ncells-Ni2), :], dims=1)), color=ColorSchemes.Reds[9], linewidth=linewidth)
lines!(ax_top, interval, vec(mean(rates[(Ncells-Ni2+1):Ncells, :], dims=1)), color=ColorSchemes.Reds[6], linewidth=linewidth)

linkxaxes!(ax, ax_top); #xlims!(ax, low=minimum(interval)+500, high=maximum(interval))
hidedecorations!(ax, grid=true, ticks=false, ticklabels=false, label=false)
hideydecorations!(ax, grid=true, ticks=true, ticklabels=false, label=false)
hideydecorations!(ax_top, grid=true, ticks=true, ticklabels=true, label=true)
rowsize!(g_ei_act, 1, Relative(1/5))    


##############################################################################################
#############################  --- eiSTDP knockout analysis --- ##############################
##############################################################################################
# --- Load Data ---
sim_name = string("network_eiSTDP-knockout.h5")
sim_savedpath = string("./networks_trained/")
fid = h5open(joinpath(sim_savedpath, sim_name), "r")
popmembers = read(fid["data"]["popmembers"])
weights = read(fid["data"]["weights"])
times = read(fid["data"]["times"])
weightsEE_tr = read(fid["data"]["weightsEE"])
weightsEI_tr = read(fid["data"]["weightsEI"])
weightsIE_tr = read(fid["data"]["weightsIE"])
close(fid)

weightsEE = sum(weights[popmembers[:, :], popmembers[:, :]], dims=(1, 3))[1, :, 1, :] ./ count(i->i>0, weights[popmembers[:, :], popmembers[:, :]], dims=(1, 3))[1, :, 1, :]
weightsEI = sum(weights[popmembers[:, :], ipopmembers[:, :]], dims=(1, 3))[1, :, 1, :] ./ count(i->i>0, weights[popmembers[:, :], ipopmembers[:, :]], dims=(1, 3))[1, :, 1, :]
weightsIE = sum(weights[ipopmembers[:, :], popmembers[:, :]], dims=(1, 3))[1, :, 1, :] ./ count(i->i>0, weights[ipopmembers[:, :], popmembers[:, :]], dims=(1, 3))[1, :, 1, :]


########################################################
################ --- E-to-E weights --- ################
########################################################

cl1 = ColorScheme(range(colorant"gray5", colorant"gray80", length=100))
cl_exc = ColorScheme(range(colorant"gray80", colorant"dodgerblue4", length=100))
excitation_cs = vcat(get(cl1, LinRange(0, 1, 100)), get(cl_exc, LinRange(0, 1, 100)))

Npop = size(weightsEE)[1]
line_positions = collect(seq_length+.5:seq_length:Npop+.5)

ax = Axis(g_ei_anal[1:3, 1:2], xlabel=L"\text{pre }E\text{-assembly}", ylabel=L"\text{post }E\text{-assembly}", xlabelsize=labelsize, ylabelsize=labelsize,
            xticks=(1:seq_length:Npop, [L"%$(x)" for x=1:seq_length:Npop]), xticklabelsize=ticklabelsize, xgridvisible=false,
            yticks=(1:seq_length:Npop, [L"%$(x)" for x=1:seq_length:Npop]), yticklabelsize=ticklabelsize, ygridvisible=false,
            aspect=1)
ht = heatmap!(ax, weightsEE, colormap=excitation_cs)
hlines!(ax, line_positions, linewidth=linewidth, color=:firebrick, linestyle=:dot)
vlines!(ax, line_positions, linewidth=linewidth, color=:firebrick, linestyle=:dot)

for (ind, ln) in enumerate(line_positions)
    vlines!(ax, ln, ymin=(ind-1)*seq_length/Npop, ymax=(ind+1)*seq_length/Npop, linewidth=linewidth, color=:firebrick)
    hlines!(ax, ln, xmin=(ind-1)*seq_length/Npop, xmax=(ind+1)*seq_length/Npop, linewidth=linewidth, color=:firebrick)
end

Colorbar(g_ei_anal[1:3, 3], ht, size=colorbarsize, label=L"\text{mean syn. str. (pF)}", height=Relative(.6), labelsize=labelsize, ticklabelsize=ticklabelsize)

########################################################
################ --- I₂-to-E weights --- ###############
########################################################

cl1 = ColorScheme(range(colorant"gray5", colorant"gray80", length=100))
cl_inh = ColorScheme(range(colorant"gray80", colorant"firebrick", length=100))
inhibition_cs = vcat(get(cl1, LinRange(0, 1, 100)), get(cl_inh, LinRange(0, 1, 100)))
    
ax = Axis(g_ei_anal[4:6, 1:2], xlabel=L"\text{pre }I_2 \text{-assembly}", ylabel=L"\text{post }E \text{-assembly}", xlabelsize=labelsize, ylabelsize=labelsize,
            xticks=(1:seq_length:Npop, [L"%$(x)" for x=1:seq_length:Npop]), xticklabelsize=ticklabelsize, xgridvisible=false,
            yticks=(1:seq_length:Npop, [L"%$(x)" for x=1:seq_length:Npop]), yticklabelsize=ticklabelsize, ygridvisible=false,
            aspect=1)
ht = heatmap!(ax, weightsIE, colormap=inhibition_cs)
hlines!(line_positions, linewidth=linewidth, color=:dodgerblue4, linestyle=:dot)
vlines!(line_positions, linewidth=linewidth, color=:dodgerblue4, linestyle=:dot)

for (ind, ln) in enumerate(line_positions)
    vlines!(ln, ymin=(ind-1)*seq_length/Npop, ymax=(ind+1)*seq_length/Npop, linewidth=linewidth, color=:dodgerblue4)
    hlines!(ln, xmin=(ind-1)*seq_length/Npop, xmax=(ind+1)*seq_length/Npop, linewidth=linewidth, color=:dodgerblue4)
end

Colorbar(g_ei_anal[4:6, 3], ht, size=colorbarsize, label=L"\text{mean syn. str. (pF)}", height=Relative(.6), labelsize=labelsize, ticklabelsize=ticklabelsize)

########################################################
################ --- Extra plot here --- ###############
########################################################


########################################################
################ --- E-to-E weights --- ################
########################################################
ax_ee = Axis(g_ei_anal[1:2, 4:5], #xlabel=L"\text{training time (%)}", ylabel=L"\text{mean synaptic strength (pF)}",
            # xticks=([0, 125, 250, 375, 500], [L"0", L"25", L"50", L"75", L"100"]),
            yticks=([5, 10, 15], [L"5", L"10", L"15"]),
            xlabelsize=labelsize, ylabelsize=labelsize, xticklabelsize=ticklabelsize, yticklabelsize=ticklabelsize)#,
            # limits=(-2.5, 505.5, 1., 19.5))
for ipop = 1:Npop
    lines!(ax_ee, weightsEE_tr[ipop, ipop, 1:500], linewidth=linewidth, color=ColorSchemes.Blues[7], label=L"\text{E-assembly}")
end
# axislegend(ax_ee, position=(.9, .7), labelsize=labelsize, titlesize=textsize, merge=true)

########################################################
################ --- I₁-to-E weights --- ###############
########################################################
ax_i1 = Axis(g_ei_anal[3:4, 4:5], xlabel=L"\text{training time (%)}", #ylabel=L"\text{mean synaptic strength (pF)}",
            # xticks=([0, 125, 250, 375, 500], [L"0", L"25", L"50", L"75", L"100"]),
            yticks=([60, 80, 100], [L"60", L"80", L"100"]),
            xlabelsize=labelsize, ylabelsize=labelsize, xticklabelsize=ticklabelsize, yticklabelsize=ticklabelsize)#,
            # limits=(-2.5, 505.5, 48.1, 115))
for ipop = 2:Npop
    (mod(ipop, seq_length) == 1) && (continue)
    lines!(ax_i1, mean(weightsIE_tr[1:(Ncells-Ne-Ni2), ipop, 1:500], dims=1)[:], linewidth=linewidth, color=ColorSchemes.Reds[9], label=L"I_1\text{-to-}E")
end
# axislegend(ax_i1, position=(.9, .3), labelsize=labelsize, titlesize=textsize, merge=true)

########################################################
################ --- I₂-to-E weights --- ###############
########################################################
ax_i2 = Axis(g_ei_anal[5:6, 4:5], xlabel=L"\text{training time (%)}", #ylabel=L"\text{mean synaptic strength (pF)}",
            xticks=([0, 125, 250, 375, 500], [L"0", L"25", L"50", L"75", L"100"]),
            yticks=([30, 60, 90], [L"30", L"60", L"90"]),
            xlabelsize=labelsize, ylabelsize=labelsize, xticklabelsize=ticklabelsize, yticklabelsize=ticklabelsize)#,
            # limits=(-2.5, 505.5, 15.1, 106.))
for ipop = 1:Npop
    lines!(ax_i2, mean(weightsIE_tr[ipopmembers[:, ipop] .- (Ne), ipop, 1:500], dims=1)[:], linewidth=linewidth, color=:gray70, label=L"I_2-\text{to}-E \text{ recurrent}")
end

for ipop = 2:Npop
    (mod(ipop, seq_length) == 1) && (continue)
    lines!(ax_i2, mean(weightsIE_tr[ipopmembers[:, ipop] .- (Ne), ipop-1, 1:500], dims=1)[:], linewidth=linewidth, color=:black, label=L"I_2\text{-to-}E \text{ preceding}")
end

for ipop = 2:Npop
    (mod(ipop, seq_length) == 1) && (continue)
    lines!(ax_i2, mean(weightsIE_tr[ipopmembers[:, ipop-1] .- (Ne), ipop, 1:500], dims=1)[:], linewidth=linewidth, color=ColorSchemes.Reds[6], label=L"I_2\text{-to-}E \text{ succeeding}")
end
# axislegend(ax_i2, position=(.1, .9), labelsize=labelsize, titlesize=textsize, merge=true)

linkxaxes!(ax_ee, ax_i1, ax_i2)
hidexdecorations!(ax_ee, grid=false, ticks=true, ticklabels=true, label=true)
hidexdecorations!(ax_i1, grid=false, ticks=true, ticklabels=true, label=true)

# Label(g_i1_anal[4:5, 4], L"\text{mean synaptic strength (pF)}",
#     fontsize=labelsize, rotation=π/2)
#


Label(g_i1_act[1, 0, TopLeft()], L"\textbf{a}",
    fontsize=subplotlabelsize,
    font=:bold,
    padding=(0, 5, 5, 0),
    halign=:right)

Label(g_i2_act[1, 0, TopLeft()], L"\textbf{b}",
    fontsize=subplotlabelsize,
    font=:bold,
    padding=(0, 5, 5, 0),
    halign=:right)

Label(g_ei_act[1, 0, TopLeft()], L"\textbf{c}",
    fontsize=subplotlabelsize,
    font=:bold,
    padding=(0, 5, 5, 0),
    halign=:right)
#

colsize!(g_i1_act, 0, Relative(0))
colsize!(g_i2_act, 0, Relative(0))
colsize!(g_ei_act, 0, Relative(0))

rowgap!(fig.layout, 1, Relative(.05))
rowgap!(fig.layout, 2, Relative(.05))

(!isdir(output_dir)) && (mkpath(output_dir))
save(joinpath(output_dir, string("figure_3.png")), fig)
end