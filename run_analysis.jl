using Pkg
Pkg.activate(".")
using InhSequences
using ImageFiltering
using Statistics
using StatsBase
using HDF5

include("quantification.jl")

save_data = false
Npop = 20
ccov_range = -400:400
tt_range = 5_000:45_000

crossEE = zeros(length(ccov_range), Npop, Npop, 10)
crossEI = zeros(length(ccov_range), Npop, Npop, 10)
crossIE = zeros(length(ccov_range), Npop, Npop, 10)

for sim_num = 1:10
    # sim_name = string("network_", sim_num,"_spontaneous.h5")
	# sim_savedpath = "./networks_trained_spontaneous/"
    sim_name = string("network_", sim_num,"_stimulation.h5")
	sim_savedpath = "./networks_trained_stimulation/"
	output_dir = "./output_analysis/"

    fid = h5open(joinpath(sim_savedpath, sim_name), "r")
    popmembers = read(fid["data"]["popmembers"])
    weights = read(fid["data"]["weights"])
    weightsEE = read(fid["data"]["weightsEE"])
    weightsEI = read(fid["data"]["weightsEI"])
    weightsIE = read(fid["data"]["weightsIE"])
    times = read(fid["data"]["times"])
    close(fid)
    
    ipopmembers = findI2populations(weights, 20, popmembers, iipop_len=25)

    # Calculate population firing rates
    pop_rates = getPopulationRates(times, popmembers, interval=tt_range, sigma=5.)
    ipop_rates = getPopulationRates(times, ipopmembers, interval=tt_range, sigma=5.)

    for ipop = 1:Npop
        for iipop = 1:Npop
            crossEE[:, ipop, iipop, sim_num] = crosscor(pop_rates[:, ipop], pop_rates[:, iipop], ccov_range)
            crossEI[:, ipop, iipop, sim_num] = crosscor(pop_rates[:, ipop], ipop_rates[:, iipop], ccov_range)
            crossIE[:, ipop, iipop, sim_num] = crosscor(ipop_rates[:, ipop], pop_rates[:, iipop], ccov_range)
        end
    end
end

if save_data
    # Save the data
    sim_savepath = "./analysis_data/"
    (!ispath(sim_savepath)) && (mkpath(sim_savepath))
    cd(sim_savepath[3:end-1])
    fid = h5open("cross_corr_stimulation.h5","w")
    g = create_group(fid,"data")
    g["crossEE"] = crossEE
    g["crossEI"] = crossEI
    g["crossIE"] = crossIE
    close(fid)
    cd("..")
end



































ipopmembers = findI2populations(weights, 20, popmembers, iipop_len=50)
# ipopmembers = ipopmembers'

cnv_spikes = convolveSpikes(times[filter(i->(i>0), popmembers[:, 1]), :], interval=1:10_000, sigma=5.)
Plots.plot(vec(mean(cnv_spikes, dims=1)), legend=false)

pop_rates = getPopulationRates(times, popmembers, interval=1:10_000, sigma=5.)
ipop_rates = getPopulationRates(times, ipopmembers, interval=1:10_000, sigma=5.)
Plots.plot(pop_rates[1:5000, 5:9], legend=false)
Plots.plot(ipop_rates[1:5000, 5:9], legend=false)



x = crosscor(pop_rates[1000:end-1000, 1], ipop_rates[1000:end-1000, 4], -400:400)
# Plots.plot(-400:400, x)
Plots.plot!(-400:400, x)


Plots.plot(crossEE[:, 1, 1, 1])

# PLOTS HERE FOR CROSSCOR
Plots.plot(-400:400, crossEE[:, 1, 1, 1])
Plots.plot!(xlabel="delay", ylabel="summed crosscor")
Plots.plot!(-400:400, crossEE[:, 1, 2, 1])
Plots.plot!(-400:400, crossEE[:, 1, 3, 1])
Plots.plot!(-400:400, crossEE[:, 1, 4, 1])

# PLOTS HERE FOR CROSSCOR
Plots.plot(-400:400, avg_crossEE[:, 1, 1])
Plots.plot!(xlabel="Delay (ms?)", ylabel="Mean cross-correlation")
Plots.plot!(-400:400, avg_crossEE[:, 1, 2])
Plots.plot!(-400:400, avg_crossEE[:, 1, 3])
Plots.plot!(-400:400, avg_crossEE[:, 1, 4])

# Average all the rest (not belonging to the sequence) and plot as ribbon
Plots.plot!(-200:200, mean(crossEE[:, 1, 5:20], dims=2)[:, 1], lw=4)


Plots.plot!(-200:200, crossEI[:, 1, 1])
Plots.plot!(-200:200, crossEI[:, 2, 2])
Plots.plot!(-200:200, crossEI[:, 3, 3])
Plots.plot!(-200:200, crossEI[:, 4, 4])

Plots.plot(-200:200, crossIE[:, 1, 2])
Plots.plot!(-200:200, crossIE[:, 2, 2])
Plots.plot!(-200:200, crossIE[:, 3, 3])
Plots.plot!(-200:200, crossIE[:, 4, 4])

################################################################
# ___________ --- Make average replay delay plot --- ___________
################################################################
using LaTeXStrings
using CairoMakie
using CurveFit

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
ax = CairoMakie.Axis(f[1, 1], xlabel="Delay (ms)", ylabel="Assembly",
                    yticks=([1, 1.25, 1.5, 2, 2.25, 2.5, 3, 3.25, 3.5, 4, 4.25], ["E₁➡E₁", "E₁➡I₁", "I₁➡E₂", "E₂➡E₂", "E₂➡I₂", "I₂➡E₃", "E₃➡E₃", "E₃➡I₃", "I₃➡E₄", "E₄➡E₄", "E₄➡I₄"]),
                    xticks=[0, 25, 50, 75, 100, 125],
                    xlabelsize=24, ylabelsize=24, xticklabelsize=18, yticklabelsize=18, limits=(-0.5, 146, 0.8, 4.5))

CairoMakie.scatter!(mean(delaysEE, dims=1)[:], [1, 2, 3, 4], markersize=10)
CairoMakie.errorbars!(mean(delaysEE, dims=1)[:], [1, 2, 3, 4], low_errorsEE, high_errorsEE, whiskerwidth=20,  direction=:x, linewidth=2)

CairoMakie.scatter!(mean(delaysEI, dims=1)[:], [1.25, 2.25, 3.25, 4.25], markersize=10)
CairoMakie.errorbars!(mean(delaysEI, dims=1)[:], [1.25, 2.25, 3.25, 4.25], low_errorsEI, high_errorsEI, whiskerwidth=20,  direction=:x, linewidth=2)

CairoMakie.scatter!(mean(delaysIE, dims=1)[:], [1.5, 2.5, 3.5], markersize=10)
CairoMakie.errorbars!(mean(delaysIE, dims=1)[:], [1.5, 2.5, 3.5], low_errorsIE, high_errorsIE, whiskerwidth=20,  direction=:x, linewidth=2)

CairoMakie.lines!([0, 25, 50, 75, 100, 125, 150], fitted_mean[1] .+ [0, 25, 50, 75, 100, 125, 150] .* fitted_mean[2], linestyle=:dashdot, label=L"\text{Slope: } %$(round(fitted_mean[2], digits=2))")
axislegend(ax, position=(.9, .2), labelsize=20, titlesize=22)

f
(!isdir(output_dir)) && (mkpath(output_dir))
save(joinpath(output_dir, string("mean_replay_delay_", activity_type,".png")), f)

#####################################################################
# ___________ --- Make average cross-correlation plot --- ___________
#####################################################################











#####################################################################
# ___________ --- Make weight distribution plot --- ___________
#####################################################################
sim_savedpath = "./networks_trained/"
output_dir = "./output_analysis/"

# Choose the parameters
Ni2 = 250
Npop = 20
Nsims = 10

x_vals = 1:Ni2
ie_weights = zeros(Ni2, Npop, Nsims)
sorted_ie_weights = zeros(Ni2, Npop, Nsims)
ie_perms = zeros(Ni2, Npop, Nsims)
fitted_params_p = zeros(2, Npop, Nsims)
fitted_params = zeros(8, Npop, Nsims)

for sim_num = 1:Nsims
    # Load data
    sim_name = string("network_", sim_num,".h5")
    fid = h5open(joinpath(sim_savedpath, sim_name), "r")
    popmembers = read(fid["data"]["popmembers"])
    weights = read(fid["data"]["weights"])
    close(fid)

    # Compute the power laws that fit the distributions
    for ipop = 1:Npop
        members = filter(i->(i>0), popmembers[ipop, :])    			# Get members of excitatory assembly
        ie_weights[:, ipop, sim_num] .= vec(sum(weights[members, 5000-Ni2+1:5000], dims=1))  # Get sum of all weights projected from each E-assemble to each 2nd ipopulation neuron
        ie_perms[:, ipop, sim_num] .= sortperm(ie_weights[:, ipop, sim_num], rev=true)
        sorted_ie_weights[:, ipop, sim_num] .= ie_weights[round.(Int, ie_perms[:, ipop, sim_num]), ipop, sim_num]
        # fitted_params_p[:, ipop, sim_num] .= power_fit(x_vals, sorted_ie_weights[:, ipop, sim_num])
        fitted_params[:, ipop, sim_num] .= poly_fit(x_vals, sorted_ie_weights[:, ipop, sim_num], 7)
    end
end

plot_data_band_min = minimum(minimum(sorted_ie_weights, dims=3)[:, :, 1], dims=2)[:]
plot_data_avg = mean(mean(sorted_ie_weights, dims=3)[:, :, 1], dims=2)[:]
plot_data_band_max = maximum(maximum(sorted_ie_weights, dims=3)[:, :, 1], dims=2)[:]

# f = CairoMakie.Figure(resolution=(720, 480))
# ax = CairoMakie.Axis(f[1, 1], xlabel=L"I_2 \, \text{neuron index}", ylabel=L"\text{Synaptic strength (sum; pF)}",
#                     yticks=([50, 100, 150, 200], [L"50", L"100", L"150", L"200"]),
#                     xticks=([1, 27, 50, 100, 150, 200, 250], [L"1", L"\mathbf{27}", L"50", L"100", L"150", L"200", L"250"]),
#                     xlabelsize=24, ylabelsize=24, xticklabelsize=18, yticklabelsize=18, limits=(0.5, 250.5, 0.5, 249))
# CairoMakie.band!(x_vals, plot_data_band_min, plot_data_band_max)
# # CairoMakie.lines!(x_vals, mean(fitted_params_p[1, :, :]) .* (x_vals .^ mean(fitted_params_p[2, :, :])), linewidth=4, color=:navajowhite, label=L"y = %$(round(mean(fitted_params[1, :, :]), digits=2))\,(\pm %$(round(std(fitted_params[1, :, :]), digits=2)))\,x^{%$(round(mean(fitted_params[2, :, :]), digits=2))\,(\pm %$(round(std(fitted_params[2, :, :]), digits=2)))}")
# CairoMakie.lines!(x_vals, mean(fitted_params[1, :, :]) .+ (x_vals .* mean(fitted_params[2, :, :])) .+ ((x_vals .^ 2) .* mean(fitted_params[3, :, :])) .+ ((x_vals .^ 3) .* mean(fitted_params[4, :, :])) .+ ((x_vals .^ 4) .* mean(fitted_params[5, :, :])) .+ ((x_vals .^ 5) .* mean(fitted_params[6, :, :])) .+ ((x_vals .^ 6) .* mean(fitted_params[7, :, :])) .+ ((x_vals .^ 7) .* mean(fitted_params[8, :, :])), linewidth=4, color=:gold3, label=L"7 \text{th degree polynomial}")
# CairoMakie.vlines!(27, ymin=0, ymax=251, linestyle=(:dash, :loose), color=:gray20)
# CairoMakie.text!(27, 240, text=L"\text{max slope change}", rotation=3*π/2, fontsize=22)
# axislegend(ax, position=(.95, .95), labelsize=20, titlesize=18, L"\mathbf{\text{Fitted power law:}}", titlehalign=:left)
# f
# (!isdir(output_dir)) && (mkpath(output_dir))
# save(joinpath(output_dir, string("power_fit_E_to_I2.png")), f)

f = CairoMakie.Figure(resolution=(720, 480))
ax = CairoMakie.Axis(f[1, 1], xlabel=L"I_2 \, \text{neuron index}", ylabel=L"\text{Synaptic strength (sum; pF)}",
                    yticks=([50, 100, 150, 200], [L"50", L"100", L"150", L"200"]),
                    xticks=([1, 27, 50, 100, 150, 200, 250], [L"1", L"\mathbf{27}", L"50", L"100", L"150", L"200", L"250"]),
                    xlabelsize=24, ylabelsize=24, xticklabelsize=18, yticklabelsize=18, limits=(0.5, 250.5, 0.5, 249))
CairoMakie.band!(x_vals, plot_data_band_min, plot_data_band_max)
# CairoMakie.errorbars!(x_vals, plot_data_avg, plot_data_avg .- plot_data_band_min, plot_data_band_max .- plot_data_avg)
CairoMakie.lines!(x_vals, plot_data_avg, linewidth=4, color=:gold3, label=L"7 \text{th degree polynomial}")
CairoMakie.vlines!(27, ymin=0, ymax=251, linestyle=(:dash, :loose), color=:gray20)
CairoMakie.text!(27, 240, text=L"\text{max slope change}", rotation=3*π/2, fontsize=22)
axislegend(ax, position=(.95, .95), labelsize=20, titlesize=18, L"\mathbf{\text{Fitted power law:}}", titlehalign=:left)
f
(!isdir(output_dir)) && (mkpath(output_dir))
save(joinpath(output_dir, string("power_fit_E_to_I2.png")), f)

#####################################################################
#####################################################################

mean(fitted_params[1, :, :])
mean(fitted_params[2, :, :])
mean(fitted_params[3, :, :])
mean(fitted_params[4, :, :])
mean(fitted_params[5, :, :])
mean(fitted_params[6, :, :])
mean(fitted_params[7, :, :])
mean(fitted_params[8, :, :])

CairoMakie.lines(x_vals, mean(fitted_params[1, :, :]) .+ (x_vals .* mean(fitted_params[2, :, :])) .+ ((x_vals .^ 2) .* mean(fitted_params[3, :, :])) .+ ((x_vals .^ 3) .* mean(fitted_params[4, :, :])) .+ ((x_vals .^ 4) .* mean(fitted_params[5, :, :])) .+ ((x_vals .^ 5) .* mean(fitted_params[6, :, :])) .+ ((x_vals .^ 6) .* mean(fitted_params[7, :, :])) .+ ((x_vals .^ 7) .* mean(fitted_params[8, :, :])), linewidth=4, color=:gold3, label=L"y = %$(round(mean(fitted_params[1, :, :]), digits=2))\,(\pm %$(round(std(fitted_params[1, :, :]), digits=2)))\,x^{%$(round(mean(fitted_params[2, :, :]), digits=2))\,(\pm %$(round(std(fitted_params[2, :, :]), digits=2)))}")


#####################################################################
#####################################################################
sim_name = string("network_1.h5")
sim_savedpath = "./networks_trained/"
output_dir = "./output_analysis/"

fid = h5open(joinpath(sim_savedpath, sim_name), "r")
popmembers = read(fid["data"]["popmembers"])
weights = read(fid["data"]["weights"])
weightsEE = read(fid["data"]["weightsEE"])
weightsEI = read(fid["data"]["weightsEI"])
weightsIE = read(fid["data"]["weightsIE"])
# times = read(fid["data"]["times"])
close(fid)

ipopmembers = findI2populations(weights, 20, popmembers', iipop_len=50)
popmembers = transpose(popmembers)

# Create EI plot_data
plot_data = zeros(10_000, 20, 20)
for ipop = 1:20
    for iipop = 1:20
        plot_data[:, ipop, iipop] = mean(weightsEI[ipopmembers[:, iipop] .- 4750, ipop, :], dims=1)
    end
end

pre = 1
Plots.plot()
for post = 1:20
    if post == pre
        Plots.plot!(plot_data[:, pre, post], lw=4)
    else
        Plots.plot!(plot_data[:, pre, post])
    end
end
Plots.plot!(legend=false)

# Create IE plot_data
plot_data = zeros(10_000, 20, 20)
for ipop = 1:20
    for iipop = 1:20
        plot_data[:, ipop, iipop] = mean(weightsIE[ipopmembers[:, iipop] .- 4000, ipop, :], dims=1)
    end
end


plot_data ./= 1250 
pre = 1
Plots.plot()
for post = 1:20
    if post == pre-1
        Plots.plot!(plot_data[:, pre, post], lw=4)
    elseif post == pre
        Plots.plot!(plot_data[:, pre, post], lw=6)
    else
        Plots.plot!(plot_data[:, pre, post], c=:grey)
    end
end
Plots.plot!(legend=false)












weightsEE_data = (10_000, Npop, Npop)
weightsEI_data = (10_000, Npop, Npop)
weightsIE_data = (10_000, Npop, Npop)

for sim = 1:10
    sim_name = string("network_", sim, ".h5")
    sim_savedpath = "./networks_trained/"
    output_dir = "./output_analysis/"

    fid = h5open(joinpath(sim_savedpath, sim_name), "r")
    popmembers = read(fid["data"]["popmembers"])
    weights = read(fid["data"]["weights"])
    weightsEE = read(fid["data"]["weightsEE"])
    weightsEI = read(fid["data"]["weightsEI"])
    weightsIE = read(fid["data"]["weightsIE"])
    # times = read(fid["data"]["times"])
    close(fid)

    ipopmembers = findI2populations(weights, 20, popmembers', iipop_len=50)
    popmembers = transpose(popmembers)

    for ipop = 1:20
        for iipop = 1:20
            plot_data[:, ipop, iipop] += mean(weightsEE[ipopmembers[:, iipop] .- 4750, ipop, :], dims=1)
        end
    end
    for ipop = 1:20
        for iipop = 1:20
            plot_data[:, ipop, iipop] += mean(weightsEI[ipopmembers[:, iipop] .- 4750, ipop, :], dims=1)
        end
    end
    for ipop = 1:20
        for iipop = 1:20
            plot_data[:, ipop, iipop] += mean(weightsEI[ipopmembers[:, iipop] .- 4750, ipop, :], dims=1)
        end
    end
end




ipop = 19
members = filter(i->(i>0), popmembers[:, ipop])    			# Get members of excitatory assembly
ie_weights = vec(sum(weights[members, 4751:5000], dims=1))  # Get sum of all weights projected from each E-assemble to each 2nd ipopulation neuron
x = sortperm(ie_weights, rev=true) #[(end-iipop_len+1):end]  


# 
Plots.bar!(ie_weights[x], a)
Plots.vline!([25])

Plots.plot!(xscale=:log, xlims=(1, 300))
