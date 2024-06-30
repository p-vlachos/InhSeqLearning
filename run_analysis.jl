using Pkg
Pkg.activate(".")
using InhSequences
using ImageFiltering
using Statistics
using StatsBase
# using Plots
using HDF5

include("quantification.jl")

save_data = true
Npop = 20
ccov_range = -400:400
tt_range = 5_000:45_000

crossEE = zeros(length(ccov_range), Npop, Npop, 10)
crossEI = zeros(length(ccov_range), Npop, Npop, 10)
crossIE = zeros(length(ccov_range), Npop, Npop, 10)

for sim_num = 1:10
    sim_name = string("network_", sim_num,"_spontaneous.h5")
	sim_savedpath = "./networks_trained_spontaneous/"
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
    fid = h5open("cross_corr.h5","w")
    g = create_group(fid,"data")
    g["crossEE"] = crossEE
    g["crossEI"] = crossEI
    g["crossIE"] = crossIE
    close(fid)
    cd("..")
end

## _________________________________

sim_name = string("average_cross_corr.h5")
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

# PLOTS HERE FOR CROSSCOR
Plots.plot(-200:200, crossEE[:, 1, 1, 1])
Plots.plot!(xlabel="delay", ylabel="summed crosscor")
Plots.plot!(-200:200, crossEE[:, 1, 2, 1])
Plots.plot!(-200:200, crossEE[:, 1, 3, 1])
Plots.plot!(-200:200, crossEE[:, 1, 4, 1])

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


# ------------- HERE I MAKE THE NEW PLOT -------------
using CairoMakie

# Calculate the peak cross-correlation for each assembly
max_delayEE = zeros(Npop, Npop)
max_delayEI = zeros(Npop, Npop)
max_delayIE = zeros(Npop, Npop)
for ipop = 1:Npop
    for iipop = 1:Npop
        max_delayEE[ipop, iipop] = findmax(crossEE[:, ipop, iipop])[2] - 200
        max_delayEI[ipop, iipop] = findmax(crossEI[:, ipop, iipop])[2] - 200
        max_delayIE[ipop, iipop] = findmax(crossIE[:, ipop, iipop])[2] - 200
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

# Compute the error bars
low_errorsEE = mean(delaysEE, dims=1)[:] - minimum(delaysEE, dims=1)[:]
high_errorsEE = maximum(delaysEE, dims=1)[:] - mean(delaysEE, dims=1)[:]

low_errorsEI = mean(delaysEI, dims=1)[:] - minimum(delaysEI, dims=1)[:]
high_errorsEI = maximum(delaysEI, dims=1)[:] - mean(delaysEI, dims=1)[:]

low_errorsIE = mean(delaysIE, dims=1)[:] - minimum(delaysIE, dims=1)[:]
high_errorsIE = maximum(delaysIE, dims=1)[:] - mean(delaysIE, dims=1)[:]

# Plot the data
f = CairoMakie.Figure()
ax = CairoMakie.Axis(f[1, 1], xlabel="Delay (ms)", ylabel="Assembly",
                    yticks=([1, 1.33, 1.66, 2, 2.33, 2.66, 3, 3.33, 3.66, 4, 4.33], ["E₁➡E₁", "E₁➡I₁", "I₁➡E₂", "E₂➡E₂", "E₂➡I₂", "I₂➡E₃", "E₃➡E₃", "E₃➡I₃", "I₃➡E₄", "E₄➡E₄", "E₄➡I₄"]))

CairoMakie.plot!(mean(delaysEE, dims=1)[:], [1, 2, 3, 4])
CairoMakie.errorbars!(mean(delaysEE, dims=1)[:], [1, 2, 3, 4], low_errorsEE, high_errorsEE, whiskerwidth=20,  direction=:x)

CairoMakie.plot!(mean(delaysEI, dims=1)[:], [1.33, 2.33, 3.33, 4.33])
CairoMakie.errorbars!(mean(delaysEI, dims=1)[:], [1.33, 2.33, 3.33, 4.33], low_errorsEI, high_errorsEI, whiskerwidth=20,  direction=:x)

CairoMakie.plot!(mean(delaysIE, dims=1)[:], [1.66, 2.66, 3.66])
CairoMakie.errorbars!(mean(delaysIE, dims=1)[:], [1.66, 2.66, 3.66], low_errorsIE, high_errorsIE, whiskerwidth=20,  direction=:x)
f

# ------------- UNIL HERE -------------

sim_name = string("network_2.h5")
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

pre = 6
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
pre = 4
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
