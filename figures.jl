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
makefig_2 = false

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

linewidth = 4
markersize = 5

# colorbarsize = 40


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
Ne = 3000
Ni2 = 250
seq_length=3
Ncells = size(times)[1]
Npop = size(popmembers, 2)
Nmaxmembers = size(popmembers, 1)
Ni_members = 27
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
mode = "spontaneous"
sim_name = string("network_1_", mode, ".h5")
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
seq_length = 3
interval = 14_500:15_000

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
    for indx in popmembers[:, pp]
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
    lines!(ax_exc, interval, vec(mean(rates[popmembers[:, pp], :], dims=1)), color=ColorSchemes.Blues[5+indx], linewidth=linewidth, label=L"\text{E}%$(indx)")
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

##############################################################################################
############################# --- Cross-correlation & Delays --- #############################
##############################################################################################

activity_type = "spontaneous"   # Choose between "stimulation" or "spontaneous"
sim_name = string("cross_corr_", activity_type, ".h5")
sim_savedpath = "./analysis_data/"

fid = h5open(joinpath(sim_savedpath, sim_name), "r")
crossEE = read(fid["data"]["crossEE"])
crossEI = read(fid["data"]["crossEI"])
crossIE = read(fid["data"]["crossIE"])
close(fid)


Nseq = Npop/seq_length
# Average over all simulations
avg_cross = mean(crossEE, dims=4)[:, :, :, 1]
# Average over all sequences (for 4 seq. of length 3)
cross0 = mean(avg_cross[:, 1, 1], dims=2)
cross0 += mean(avg_cross[:, 4, 4], dims=2)
cross0 += mean(avg_cross[:, 7, 7], dims=2)
cross0 += mean(avg_cross[:, 10, 10], dims=2)

cross1 = mean(avg_cross[:, 1, 2], dims=2)
cross1 += mean(avg_cross[:, 4, 5], dims=2)
cross1 += mean(avg_cross[:, 7, 8], dims=2)
cross1 += mean(avg_cross[:, 10, 11], dims=2)

cross2 = mean(avg_cross[:, 1, 3], dims=2)
cross2 += mean(avg_cross[:, 4, 6], dims=2)
cross2 += mean(avg_cross[:, 7, 9], dims=2)
cross2 += mean(avg_cross[:, 10, 12], dims=2)

ax_mean = Axis(g_corr[1, 1], 
            xlabel=L"\text{lag (ms)}", xlabelsize=labelsize,
            ylabel=L"\text{cross-correlation}", ylabelsize=labelsize,
            xticks=([-200., -100., 0., 100., 200.], [L"-200", L"-100", L"0", L"100", L"200"]), xticklabelsize=ticklabelsize,
            yticks=([0., .5, 1.], [L"0", L"0.5", L"1"]), yticklabelsize=ticklabelsize,
            limits=(-280, 280, -0.15, 1.1))
lines!(ax_mean, -400:400, cross0./Nseq, color=ColorSchemes.Blues[6], linewidth=linewidth, label=L"E1-E1")
lines!(ax_mean, -400:400, cross1./Nseq, color=ColorSchemes.Blues[7], linewidth=linewidth, label=L"E1-E2")
lines!(ax_mean, -400:400, cross2./Nseq, color=ColorSchemes.Blues[8], linewidth=linewidth, label=L"E1-E3")
axislegend(ax_mean, position=(0.95, 0.95), labelsize=legendlabelsize)

# ____________ --- Delays --- ____________
# Average over all simulations
avg_crossEE = mean(crossEE, dims=4)[:, :, :, 1]
avg_crossEI = mean(crossEI, dims=4)[:, :, :, 1]
avg_crossIE = mean(crossIE, dims=4)[:, :, :, 1]

# Calculate the peak cross-correlation for each assembly
max_delayEE = zeros(Npop, Npop)
max_delayEI = zeros(Npop, Npop)
max_delayIE = zeros(Npop, Npop)
for ipop = 1:Npop
    for iipop = 1:Npop
        max_delayEE[ipop, iipop] = findmax(avg_crossEE[:, ipop, iipop])[2] - 400
        max_delayEI[ipop, iipop] = findmax(avg_crossEI[:, ipop, iipop])[2] - 400
        max_delayIE[ipop, iipop] = findmax(avg_crossIE[:, ipop, iipop])[2] - 400
    end
end

# Extract the relative cross-correlation delays
delaysEE = zeros(4, 3)
delaysEI = zeros(4, 3)
delaysIE = zeros(4, 2)
for ipop = 1:seq_length:Npop
    delaysEE[div(ipop, seq_length) + 1, 1] = max_delayEE[ipop, ipop]
    delaysEE[div(ipop, seq_length) + 1, 2] = max_delayEE[ipop, ipop + 1]
    delaysEE[div(ipop, seq_length) + 1, 3] = max_delayEE[ipop, ipop + 2]

    delaysEI[div(ipop, seq_length) + 1, 1] = max_delayEI[ipop, ipop]
    delaysEI[div(ipop, seq_length) + 1, 2] = max_delayEI[ipop, ipop + 1]
    delaysEI[div(ipop, seq_length) + 1, 3] = max_delayEI[ipop, ipop + 2]
    
    delaysIE[div(ipop, seq_length) + 1, 1] = max_delayIE[ipop, ipop + 1]
    delaysIE[div(ipop, seq_length) + 1, 2] = max_delayIE[ipop, ipop + 2]
end

# Fit the mean values to get the recall speed
fitted_mean = linear_fit(mean(delaysEE, dims=1)[:], [1, 2, 3])

# Compute the error bars
low_errorsEE = mean(delaysEE, dims=1)[:] - minimum(delaysEE, dims=1)[:]
high_errorsEE = maximum(delaysEE, dims=1)[:] - mean(delaysEE, dims=1)[:]

low_errorsEI = mean(delaysEI, dims=1)[:] - minimum(delaysEI, dims=1)[:]
high_errorsEI = maximum(delaysEI, dims=1)[:] - mean(delaysEI, dims=1)[:]

low_errorsIE = mean(delaysIE, dims=1)[:] - minimum(delaysIE, dims=1)[:]
high_errorsIE = maximum(delaysIE, dims=1)[:] - mean(delaysIE, dims=1)[:]

# Plot the data
ax = Axis(g_corr[2, 1], xlabel=L"\text{delay (ms)}", ylabel=L"\text{activity trajectory}",
                    yticks=([1, 1.25, 1.5, 2, 2.25, 2.5, 3, 3.25], [L"E1→E1", L"E1→I_21", L"I_21→E2", L"E2→E2", L"E2→I_22", L"I_22→E3", L"E3→E3", L"E3→I_23"]),
                    # xticks=([0, 25, 50, 75, 100, 125], [L"0", L"25", L"50", L"75", L"100", L"125"]),
                    xlabelsize=labelsize, ylabelsize=labelsize, xticklabelsize=ticklabelsize, yticklabelsize=ticklabelsize)#, limits=(-0.5, 146, 0.8, 4.5))

scatter!(ax, mean(delaysEE, dims=1)[:], [1, 2, 3], markersize=markersize, color=ColorSchemes.Blues[6])
errorbars!(ax, mean(delaysEE, dims=1)[:], [1, 2, 3], low_errorsEE, high_errorsEE, whiskerwidth=20,  direction=:x, linewidth=linewidth, color=ColorSchemes.Blues[6])

scatter!(ax, mean(delaysEI, dims=1)[:], [1.25, 2.25, 3.25], markersize=markersize, color=ColorSchemes.Reds[6])
errorbars!(ax, mean(delaysEI, dims=1)[:], [1.25, 2.25, 3.25], low_errorsEI, high_errorsEI, whiskerwidth=20,  direction=:x, linewidth=linewidth, color=ColorSchemes.Reds[6])

scatter!(ax, mean(delaysIE, dims=1)[:], [1.5, 2.5], markersize=markersize, color=ColorSchemes.BuPu[6])
errorbars!(ax, mean(delaysIE, dims=1)[:], [1.5, 2.5], low_errorsIE, high_errorsIE, whiskerwidth=20,  direction=:x, linewidth=linewidth, color=ColorSchemes.BuPu[6])


##############################################################################################
################################### --- Sequentiality --- ####################################
##############################################################################################

sim_savedpath = "./analysis_data/sequentiality/"
Nsimulations = 4
Nseeds = 10
seq_scores = zeros(Nsimulations, Nseeds)

for sim = 1:Nsimulations
    sim_name = string("sequentiality_network_", sim,"_", activity_type, ".h5")
    fid = h5open(joinpath(sim_savedpath, sim_name), "r")
    seq_scores[sim, :] = read(fid["data"]["scores"])
    close(fid)
end

mean_scores = mean(seq_scores, dims=2)[:, 1]
seq_scores = vec(seq_scores)

categories = ones(length(seq_scores)+length(mean_scores))
categories[end-length(mean_scores):end] .= 2

ax = Axis(g_seq[1, 1], ylabel=L"\text{sequentiality (%)}", ylabelsize=labelsize,
            xticks=([1, 2], [L"\text{overall}", L"\text{mean / model seed}"]), xticklabelsize=ticklabelsize,
            yticks=([.9, .95, 1.], [L"90", L"95", L"100"]), yticklabelsize=ticklabelsize,
            limits=(0, 3, .88, 1.02))

violin!(ax, categories, vcat(seq_scores, mean_scores), dodge=dodge, dodge_gap=dodge_gap, width=v_width)
boxplot!(ax, categories, vcat(seq_scores, mean_scores), color=(:white, .0), strokewidth=b_strokewidth, width=b_width, strokecolor=(:black, 1.))



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

Label(g_corr[1, 0, TopLeft()], L"\textbf{d}",
    fontsize=subplotlabelsize,
    font=:bold,
    padding=(0, 5, 5, 0),
    halign=:right)

Label(g_seq[1, 0, TopLeft()], L"\textbf{e}",
    fontsize=subplotlabelsize,
    font=:bold,
    padding=(0, 5, 5, 0),
    halign=:right)
    
colsize!(g_model, 0, Relative(0))
colsize!(g_raster, 0, Relative(0))
colsize!(g_react, 0, Relative(0))
colsize!(g_corr, 0, Relative(0))
colsize!(g_seq, 0, Relative(0))

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

fig = Figure(size=(4800, 3280))

g_structure = fig[1, 1] = GridLayout()


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

linewidth = 4
markersize = 5

# colorbarsize = 40





(!isdir(output_dir)) && (mkpath(output_dir))
save(joinpath(output_dir, string("figure_2.png")), fig)
end