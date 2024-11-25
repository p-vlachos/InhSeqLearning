using Pkg
Pkg.activate(".")
using InhSequences
using Statistics
using HDF5

# include("quantification.jl")

using ColorSchemes
using LaTeXStrings
using CairoMakie
using CurveFit
using Colors


makefig_1 = true


##############################################################################################
######                               --- FIGURE 1 ---                                   ###### 
##############################################################################################
if makefig_1

fig = Figure(resolution=(3240, 2280))

g_raster = fig[1, 1] = GridLayout()
g_react = fig[1, 2] = GridLayout()
g_corr = fig[2, 1] = GridLayout()
g_delay = fig[2, 2] = GridLayout()


labelsize = 66
ticklabelsize = 54
legendlabelsize = 54
# legendticklabelsize = 54
# colorbarlabelsize = 66
# colorbarticklabelsize = 54
textsize = 66
subplotlabelsize = 66

# colorbarsize = 40


mode = "spontaneous" # "spontaneous" "stimulation" "plot"
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

##############################################################################################
######  --- Raster plot --- 
##############################################################################################

# Parameters
Ne = 4000
Ni2 = 250
Ncells = size(times)[1]
Npop = size(popmembers, 2)
Nmembers_max = size(popmembers, 1)
Ni_members = 27
interval = 15_000:32_000

ipopmembers = findI2populations(weights, Npop, popmembers, iipop_len=Ni_members)
restcells = deleteat!(map(*, ones(Int, 4000), range(1,stop=4000)), sort(unique(popmembers))[2:end])
ylim_max = count(i->(i>0), popmembers) + length(restcells) + (Ncells-Ni2-Ne) + length(ipopmembers)

rates = convolveSpikes(times, interval=interval, gaussian=false, sigma=10.)
emembers = unique(filter(i->i>0, popmembers))

ns = zeros(Ncells)
for cc = 1:Ncells
    ns[cc] = count(i->i>0, times[cc, :])
end

# Plot params
mrksize = 1.
ytick_seq = zeros(10)

rowcount = 1
ytickcount = 1
fig = CairoMakie.Figure(resolution=(1080, 720))
g = fig[1, 1] = GridLayout()
ax = CairoMakie.Axis(g[2, 1], xlabel=L"\text{simulation time (s)}", ylabel=L"\text{sequences (neurons)}", xlabelsize=24, ylabelsize=24,
            xticks=(collect((minimum(interval)+500):5000:maximum(interval)), [L"0", L"5", L"10", L"15"]), xticklabelsize=18, xgridvisible=false,
            yticklabelsize=16, ygridvisible=false,
            limits=(minimum(interval)+500, maximum(interval), 1, ylim_max))
# Excitatoy assembly members
for pp = 1:Npop
    for cc = 1:Nmembers_max
        (popmembers[cc, pp] < 1) && (break)
        indx = round(Int, popmembers[cc, pp])
        x_vals = times[indx, 1:round(Int, ns[indx])]
        y_vals = rowcount * ones(length(x_vals))
        CairoMakie.scatter!(ax, x_vals, y_vals, color=RGBA(14/255, 134/255, 212/255), markersize=mrksize)
        rowcount += 1
    end
    if mod(pp, 4) == 0
        hlines!(ax, rowcount, linewidth=0.8, color=RGBA(14/255, 134/255, 212/255))
    else
        hlines!(ax, rowcount, linewidth=0.3, color=RGBA(14/255, 134/255, 212/255))
    end
    if mod(pp+2, 4) == 0
        ytick_seq[ytickcount] = rowcount
        ytickcount += 1
    end
end
hlines!(ax, rowcount, linewidth=0.5, color=:gray20)
# Excitatory non-members 
for cc in restcells
    x_vals = times[cc, 1:round(Int, ns[cc])]
    y_vals = rowcount * ones(length(x_vals))
    CairoMakie.scatter!(ax, x_vals, y_vals, color=RGBA(103/255, 136/255, 150/255), markersize=mrksize)
    rowcount += 1
end
hlines!(ax, rowcount, linewidth=0.5, color=:gray20)
# Inhibitory I₁
for cc = (Ne+1):(Ncells-Ni2+1)
    x_vals = times[cc, 1:round(Int, ns[cc])]
    y_vals = rowcount * ones(length(x_vals))
    CairoMakie.scatter!(ax, x_vals, y_vals, color=RGBA(131/255, 30/255, 19/255), markersize=mrksize)
    rowcount += 1
end
hlines!(ax, rowcount, linewidth=0.5, color=:gray20)
# I₂ assembly members
for pp = 1:Npop
    for cc = 1:Ni_members
        indx = round(Int, popmembers[cc, pp])
        x_vals = times[indx, 1:round(Int, ns[indx])]
        y_vals = rowcount * ones(length(x_vals))
        CairoMakie.scatter!(ax, x_vals, y_vals, color=RGBA(220/255, 50/255, 30/255), markersize=mrksize)
        rowcount += 1
    end
    if mod(pp, 4) == 0
        hlines!(ax, rowcount, linewidth=0.3, color=RGBA(220/255, 50/255, 30/255))
    end
    if mod(pp+2, 4) == 0
        ytick_seq[ytickcount] = rowcount
        ytickcount += 1
    end
end
ax.yticks = (ytick_seq, [L"\text{A}", L"\text{B}", L"\text{C}", L"\text{D}", L"\text{E}", L"\text{A}", L"\text{B  }", L"\text{C}", L"\text{D  }", L"\text{E}"])

ax_top = CairoMakie.Axis(g[1, 1], xlabel="", ylabel=L"\text{firing rate (Hz)}", xlabelsize=20, ylabelsize=20,
                xticks=([], []), xticklabelsize=20,
                yticks=([0, 1, 2], [L"0", L"1", L"2"]), yticklabelsize=18,
                limits=(minimum(interval)+500, maximum(interval), -0.2, 2.5))
# Plot mean firing rates
CairoMakie.lines!(ax_top, interval, vec(mean(rates[emembers, :], dims=1)), color=RGBA(14/255, 134/255, 212/255), linewidth=1.)
CairoMakie.lines!(ax_top, interval, vec(mean(rates[restcells, :], dims=1)), color=RGBA(103/255, 136/255, 150/255), linewidth=1.)
CairoMakie.lines!(ax_top, interval, vec(mean(rates[(Ne+1):(Ncells-Ni2), :], dims=1)), color=RGBA(131/255, 30/255, 19/255), linewidth=1.)
CairoMakie.lines!(ax_top, interval, vec(mean(rates[(Ncells-Ni2+1):Ncells, :], dims=1)), color=RGBA(220/255, 50/255, 30/255), linewidth=1.)

linkxaxes!(ax, ax_top); #xlims!(ax, low=minimum(interval)+500, high=maximum(interval))
hidedecorations!(ax, grid=true, ticks=false, ticklabels=false, label=false)
hideydecorations!(ax, grid=true, ticks=true, ticklabels=false, label=false)
rowsize!(g, 1, Relative(1/5))

fig

(!ispath(output_dir)) && (mkpath(output_dir))
save(string(output_dir, "test_plot_sim_spontaneous_alpha.png"), fig, dpi=150)

##############################################################################################
######    --- Raster plot (zoom) ---
##############################################################################################

# Parameters
Ne = 4000
Ni2 = 250
Ncells = size(times)[1]
Npop = size(popmembers, 2)
Nmembers_max = size(popmembers, 1)
Ni_members = 27
seq_number = 3
seq_length = 4
interval = 25_950:26_350

ipopmembers = findI2populations(weights, Npop, popmembers, iipop_len=Ni_members)
ylim_max = count(i->(i>0), popmembers[:, (((seq_number-1)*seq_length)+1):(seq_number*seq_length)]) + length(ipopmembers[:, (((seq_number-1)*seq_length)+1):(seq_number*seq_length)])

rates = convolveSpikes(times, interval=interval, gaussian=false, sigma=10.)

emembers1 = unique(filter(i->i>0, popmembers[:, (((seq_number-1)*seq_length)+1)]))
emembers2 = unique(filter(i->i>0, popmembers[:, (((seq_number-1)*seq_length)+1)+1]))
emembers3 = unique(filter(i->i>0, popmembers[:, (((seq_number-1)*seq_length)+1)+2]))
emembers4 = unique(filter(i->i>0, popmembers[:, (((seq_number-1)*seq_length)+1)+3]))

ns = zeros(Ncells)
for cc = 1:Ncells
    ns[cc] = count(i->i>0, times[cc, :])
end

# Plot params
mrksize = 8.
linewidth = 3.
legendlabelsize = 16.


rowcount = 1
ytick_seq = zeros(8)
ytick_prev = 0.
ytickcount = 1
fig = CairoMakie.Figure(resolution=(920, 720))
g = fig[1, 1] = GridLayout()
ax = CairoMakie.Axis(g[1:2, 1], xlabel=L"\text{simulation time (ms)}", ylabel=L"\text{sequence C}", xlabelsize=24, ylabelsize=24,
            xticks=(collect((minimum(interval)):100:maximum(interval)), [L"0", L"100", L"200", L"300", L"400"]), xticklabelsize=18, xgridvisible=false,
            yticklabelsize=16, ygridvisible=false,
            limits=(minimum(interval), maximum(interval), 1, ylim_max))
# Excitatoy assembly members
for pp = (((seq_number-1)*seq_length)+1):(seq_number*seq_length)
    for cc = 1:Nmembers_max
        (popmembers[cc, pp] < 1) && (break)
        indx = round(Int, popmembers[cc, pp])
        x_vals = times[indx, 1:round(Int, ns[indx])]
        y_vals = rowcount * ones(length(x_vals))
        CairoMakie.scatter!(ax, x_vals, y_vals, color=ColorSchemes.Blues[5+ytickcount], markersize=mrksize)
        rowcount += 1
    end
    hlines!(ax, rowcount, linewidth=0.5, color=ColorSchemes.Blues[9])
    if ytickcount == 1
        ytick_seq[ytickcount] = rowcount/2; ytickcount += 1; ytick_prev = rowcount
    else
        ytick_seq[ytickcount] = ytick_prev + (rowcount-ytick_prev)/2; ytickcount += 1; ytick_prev = rowcount
    end
end
hlines!(ax, rowcount, linewidth=0.5, color=:gray20)
# I₂ assembly members
for pp = (((seq_number-1)*seq_length)+1):(seq_number*seq_length)
    for cc = 1:Ni_members
        indx = round(Int, popmembers[cc, pp])
        x_vals = times[indx, 1:round(Int, ns[indx])]
        y_vals = rowcount * ones(length(x_vals))
        CairoMakie.scatter!(ax, x_vals, y_vals, color=ColorSchemes.Reds[1+ytickcount], markersize=mrksize)
        rowcount += 1
    end
    hlines!(ax, rowcount, linewidth=0.5, color=ColorSchemes.Reds[9])
    ytick_seq[ytickcount] = ytick_prev + (rowcount-ytick_prev)/2; ytickcount += 1; ytick_prev = rowcount
end
fig

ax.yticks = (ytick_seq, [L"1", L"2", L"3", L"4", L"1", L"2", L"3", L"4"])

ax_exc = CairoMakie.Axis(g[1, 2], xlabel="", ylabel=L"E-\text{assembly firing rate (Hz)}", xlabelsize=20, ylabelsize=20,
xticks=(collect((minimum(interval)):100:maximum(interval)), [L"0", L"100", L"200", L"300", L"400"]), xticklabelsize=20,
                yticks=([0, 4, 8, 12, 16, 20], [L"0", L"4", L"8", L"12", L"16", L"20"]), yticklabelsize=18,
                limits=(minimum(interval), maximum(interval), -0.2, 18.3))
# Plot mean firing rates
CairoMakie.lines!(ax_exc, interval, vec(mean(rates[emembers1, :], dims=1)), color=ColorSchemes.Blues[6], linewidth=linewidth, label=L"E1")#L"E\text{-assembly }1")
CairoMakie.lines!(ax_exc, interval, vec(mean(rates[emembers2, :], dims=1)), color=ColorSchemes.Blues[7], linewidth=linewidth, label=L"E2")#L"E\text{-assembly }2")
CairoMakie.lines!(ax_exc, interval, vec(mean(rates[emembers3, :], dims=1)), color=ColorSchemes.Blues[8], linewidth=linewidth, label=L"E3")#L"E\text{-assembly }3")
CairoMakie.lines!(ax_exc, interval, vec(mean(rates[emembers4, :], dims=1)), color=ColorSchemes.Blues[9], linewidth=linewidth, label=L"E4")#L"E\text{-assembly }4")
axislegend(ax_exc, position=(0.99, 0.99), labelsize=legendlabelsize)
fig

ax_inh = CairoMakie.Axis(g[2, 2], xlabel=L"\text{simulation time (ms)}", ylabel=L"I_2-\text{assembly firing rate (Hz)}", xlabelsize=20, ylabelsize=20,
                xticks=(collect((minimum(interval)):100:maximum(interval)), [L"0", L"100", L"200", L"300", L"400"]), xticklabelsize=20,
                yticks=([0, 4, 8, 12], [L"0", L"4", L"8", L"12"]), yticklabelsize=18,
                limits=(minimum(interval), maximum(interval), -0.2, 12.8))
# Plot mean firing rates
CairoMakie.lines!(ax_inh, interval, vec(mean(rates[ipopmembers[:, (((seq_number-1)*seq_length)+1)], :], dims=1)), color=ColorSchemes.Reds[6], linewidth=linewidth, label=L"I_21")#L"I_2\text{-assembly }1")
CairoMakie.lines!(ax_inh, interval, vec(mean(rates[ipopmembers[:, (((seq_number-1)*seq_length)+1)+1], :], dims=1)), color=ColorSchemes.Reds[7], linewidth=linewidth, label=L"I_22")#L"I_2\text{-assembly }2")
CairoMakie.lines!(ax_inh, interval, vec(mean(rates[ipopmembers[:, (((seq_number-1)*seq_length)+1)+2], :], dims=1)), color=ColorSchemes.Reds[8], linewidth=linewidth, label=L"I_23")#L"I_2\text{-assembly }3")
CairoMakie.lines!(ax_inh, interval, vec(mean(rates[ipopmembers[:, (((seq_number-1)*seq_length)+1)+3], :], dims=1)), color=ColorSchemes.Reds[9], linewidth=linewidth, label=L"I_24")#L"I_2\text{-assembly }4")
axislegend(ax_inh, position=(0.99, 0.99), labelsize=legendlabelsize)
fig

linkxaxes!(ax_exc, ax_inh) #; xlims!(ax, low=minimum(interval)+500, high=maximum(interval))
# hidedecorations!(ax_exc, grid=true, ticks=false, ticklabels=false, label=false)
# hidedecorations!(ax_inh, grid=true, ticks=false, ticklabels=false, label=false)

hidexdecorations!(ax_exc, grid=false, ticks=true, ticklabels=true, label=true)
# hideydecorations!(ax_exc, grid=false, ticks=true, ticklabels=false, label=true)
# rowsize!(g, 1, Relative(1/5))

fig

(!ispath(output_dir)) && (mkpath(output_dir))
save(string(output_dir, "test_plot_sim_spontaneous.png"), fig, dpi=150)


##############################################################################################
######    --- Cross-correlation ---
##############################################################################################
## _________________________________

sim_name = string("average_cross_corr.h5")
# sim_name = string("cross_corr.h5")
sim_name = string("cross_corr_stimulation.h5")
sim_name = string("cross_corr_spontaneous.h5")
sim_savedpath = "./analysis_data/"
output_dir = "./output_analysis/"

fid = h5open(joinpath(sim_savedpath, sim_name), "r")
# popmembers = read(fid["data"]["popmembers"])
# weights = read(fid["data"]["weights"])
crossEE = read(fid["data"]["crossEE"])
crossEI = read(fid["data"]["crossEI"])
crossIE = read(fid["data"]["crossIE"])
# times = read(fid["data"]["times"])
close(fid)


# Average over all simulations
avg_cross = mean(crossEE, dims=4)
# Average over all 1st sequences
cross0 = mean(avg_cross[:, 1, 1, 1], dims=2)
cross0 += mean(avg_cross[:, 5, 5, 1], dims=2)
cross0 += mean(avg_cross[:, 9, 9, 1], dims=2)
cross0 += mean(avg_cross[:, 13, 13, 1], dims=2)
cross0 += mean(avg_cross[:, 17, 17, 1], dims=2)

cross1 = mean(avg_cross[:, 1, 2, 1], dims=2)
cross1 += mean(avg_cross[:, 5, 6, 1], dims=2)
cross1 += mean(avg_cross[:, 9, 10, 1], dims=2)
cross1 += mean(avg_cross[:, 13, 14, 1], dims=2)
cross1 += mean(avg_cross[:, 17, 18, 1], dims=2)

cross2 = mean(avg_cross[:, 1, 3, 1], dims=2)
cross2 += mean(avg_cross[:, 5, 7, 1], dims=2)
cross2 += mean(avg_cross[:, 9, 11, 1], dims=2)
cross2 += mean(avg_cross[:, 13, 15, 1], dims=2)
cross2 += mean(avg_cross[:, 17, 19, 1], dims=2)

cross3 = mean(avg_cross[:, 1, 4, 1], dims=2)
cross3 += mean(avg_cross[:, 5, 8, 1], dims=2)
cross3 += mean(avg_cross[:, 9, 12, 1], dims=2)
cross3 += mean(avg_cross[:, 13, 16, 1], dims=2)
cross3 += mean(avg_cross[:, 17, 20, 1], dims=2)


fig = Figure(resolution=(720, 480))
ax = Axis(fig[1, 1], 
            xlabel=L"\text{lag (ms)}", xlabelsize=24,
            ylabel=L"\text{cross-correlation}", ylabelsize=24,
            xticks=([-200., -100., 0., 100., 200.], [L"-200", L"-100", L"0", L"100", L"200"]), xticklabelsize=20,
            yticks=([0., .5, 1.], [L"0", L"0.5", L"1"]), yticklabelsize=20,
            limits=(-280, 280, -0.2, 1.2))
lines!(ax, -400:400, cross0./5, color=ColorSchemes.Blues[6], linewidth=linewidth, label=L"E1-E1")
lines!(ax, -400:400, cross1./5, color=ColorSchemes.Blues[7], linewidth=linewidth, label=L"E1-E2")
lines!(ax, -400:400, cross2./5, color=ColorSchemes.Blues[8], linewidth=linewidth, label=L"E1-E3")
lines!(ax, -400:400, cross3./5, color=ColorSchemes.Blues[9], linewidth=linewidth, label=L"E1-E4")
axislegend(ax, position=(0.95, 0.95), labelsize=legendlabelsize)
fig






fig = Figure(resolution=(720, 480))
ax = Axis(fig[1, 1], 
            xlabel=L"\text{lag (ms)}", xlabelsize=24,
            ylabel=L"\text{cross-covariance}", ylabelsize=24,
            xticks=([-200., -100., 0., 100., 200.], [L"-200", L"-100", L"0", L"100", L"200"]), xticklabelsize=20,
            yticks=([0., .5, 1.], [L"0", L"0.5", L"1"]), yticklabelsize=20,
            limits=(-280, 280, -0.2, 1.1))

for pop = 20:-1:5
    lines!(ax, -400:400, avg_cross[:, 1, pop, 1], color=ColorSchemes.Greys[5], linewidth=linewidth/2, label=L"\text{Other}")
end
lines!(ax, -400:400, avg_cross[:, 1, 1, 1], color=ColorSchemes.Blues[6], linewidth=linewidth, label=L"E1-E1")
lines!(ax, -400:400, avg_cross[:, 1, 2, 1], color=ColorSchemes.Blues[7], linewidth=linewidth, label=L"E1-E2")
lines!(ax, -400:400, avg_cross[:, 1, 3, 1], color=ColorSchemes.Blues[8], linewidth=linewidth, label=L"E1-E3")
lines!(ax, -400:400, avg_cross[:, 1, 4, 1], color=ColorSchemes.Blues[9], linewidth=linewidth, label=L"E1-E4")
axislegend(ax, position=(0.95, 0.95), labelsize=legendlabelsize, merge=true)
fig


##############################################################################################
######    --- delays ---
##############################################################################################
# Load data
activity_type = "spontaneous"   # Choose between "stimulation" or "spontaneous"
sim_name = string("cross_corr_", activity_type, ".h5")
sim_savedpath = "./analysis_data/"
output_dir = "./output_analysis/"

fid = h5open(joinpath(sim_savedpath, sim_name), "r")
crossEE = read(fid["data"]["crossEE"])
crossEI = read(fid["data"]["crossEI"])
crossIE = read(fid["data"]["crossIE"])
close(fid)

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
delaysEE = zeros(5, 4)
delaysEI = zeros(5, 4)
delaysIE = zeros(5, 3)
for ipop = 1:4:Npop
    delaysEE[div(ipop, 4) + 1, 1] = max_delayEE[ipop, ipop]
    delaysEE[div(ipop, 4) + 1, 2] = max_delayEE[ipop, ipop + 1]
    delaysEE[div(ipop, 4) + 1, 3] = max_delayEE[ipop, ipop + 2]
    delaysEE[div(ipop, 4) + 1, 4] = max_delayEE[ipop, ipop + 3]

    delaysEI[div(ipop, 4) + 1, 1] = max_delayEI[ipop, ipop]
    delaysEI[div(ipop, 4) + 1, 2] = max_delayEI[ipop, ipop + 1]
    delaysEI[div(ipop, 4) + 1, 3] = max_delayEI[ipop, ipop + 2]
    delaysEI[div(ipop, 4) + 1, 4] = max_delayEI[ipop, ipop + 3]
    
    delaysIE[div(ipop, 4) + 1, 1] = max_delayIE[ipop, ipop + 1]
    delaysIE[div(ipop, 4) + 1, 2] = max_delayIE[ipop, ipop + 2]
    delaysIE[div(ipop, 4) + 1, 3] = max_delayIE[ipop, ipop + 3]
end

# Fit the mean values to get the recall speed
fitted_mean = linear_fit(mean(delaysEE, dims=1)[:], [1, 2, 3, 4])

# Compute the error bars
low_errorsEE = mean(delaysEE, dims=1)[:] - minimum(delaysEE, dims=1)[:]
high_errorsEE = maximum(delaysEE, dims=1)[:] - mean(delaysEE, dims=1)[:]

low_errorsEI = mean(delaysEI, dims=1)[:] - minimum(delaysEI, dims=1)[:]
high_errorsEI = maximum(delaysEI, dims=1)[:] - mean(delaysEI, dims=1)[:]

low_errorsIE = mean(delaysIE, dims=1)[:] - minimum(delaysIE, dims=1)[:]
high_errorsIE = maximum(delaysIE, dims=1)[:] - mean(delaysIE, dims=1)[:]

# Plot the data
f = CairoMakie.Figure(resolution=(720, 480))
ax = CairoMakie.Axis(f[1, 1], xlabel=L"\text{delay (ms)}", ylabel=L"\text{activity trajectory}",
                    yticks=([1, 1.25, 1.5, 2, 2.25, 2.5, 3, 3.25, 3.5, 4, 4.25], [L"E1→E1", L"E1→I_21", L"I_21→E2", L"E2→E2", L"E2→I_22", L"I_22→E3", L"E3→E3", L"E3→I_23", L"I_23→E4", L"E4→E4", L"E4→I_24"]),
                    xticks=([0, 25, 50, 75, 100, 125], [L"0", L"25", L"50", L"75", L"100", L"125"]),
                    xlabelsize=24, ylabelsize=24, xticklabelsize=18, yticklabelsize=18, limits=(-0.5, 146, 0.8, 4.5))

CairoMakie.scatter!(mean(delaysEE, dims=1)[:], [1, 2, 3, 4], markersize=10, color=ColorSchemes.Blues[6])
CairoMakie.errorbars!(mean(delaysEE, dims=1)[:], [1, 2, 3, 4], low_errorsEE, high_errorsEE, whiskerwidth=20,  direction=:x, linewidth=2, color=ColorSchemes.Blues[6])

CairoMakie.scatter!(mean(delaysEI, dims=1)[:], [1.25, 2.25, 3.25, 4.25], markersize=10, color=ColorSchemes.Reds[6])
CairoMakie.errorbars!(mean(delaysEI, dims=1)[:], [1.25, 2.25, 3.25, 4.25], low_errorsEI, high_errorsEI, whiskerwidth=20,  direction=:x, linewidth=2, color=ColorSchemes.Reds[6])

CairoMakie.scatter!(mean(delaysIE, dims=1)[:], [1.5, 2.5, 3.5], markersize=10, color=ColorSchemes.BuPu[6])
CairoMakie.errorbars!(mean(delaysIE, dims=1)[:], [1.5, 2.5, 3.5], low_errorsIE, high_errorsIE, whiskerwidth=20,  direction=:x, linewidth=2, color=ColorSchemes.BuPu[6])

# CairoMakie.lines!([0, 25, 50, 75, 100, 125, 150], fitted_mean[1] .+ [0, 25, 50, 75, 100, 125, 150] .* fitted_mean[2], linestyle=:dashdot, label=L"\text{Slope: } %$(round(fitted_mean[2], digits=2))")
# axislegend(ax, position=(.9, .2), labelsize=20, titlesize=22)

f



##############################################################################################
######    --- Weights ---
##############################################################################################





