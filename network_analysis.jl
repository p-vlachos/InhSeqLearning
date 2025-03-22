using Pkg
Pkg.activate(".")
using InhSequences
using Statistics
using HDF5

# include("plots.jl")

simulation = 1
# sim_name = string("network", simulation, "_.h5")
# sim_savedpaths = ["./networks_trained/", "./networks_trained_spontaneous/", "./networks_trained_stimulation/"]
sim_savedpaths = ["./networks_trained_spontaneous/", "./networks_trained_stimulation/"]

for sim_savedpath in sim_savedpaths
    for sim = simulation:simulation

        # if sim_savedpath == "./networks_trained/"
        #     sim_name = string("network_", sim,".h5")
        # else  _i1STDP-knockout
        if sim_savedpath == "./networks_trained_spontaneous/"
            # sim_name = string("network_", sim,"_spontaneous.h5")
            sim_name = string("network_", sim,"_i1STDP-knockout_spontaneous.h5")
        else
            # sim_name = string("network_", sim,"_stimulation.h5")
            sim_name = string("network_", sim,"_i1STDP-knockout_stimulation.h5")
        end

        fid = h5open(joinpath(sim_savedpath, sim_name), "r")
        popmembers = read(fid["data"]["popmembers"])
        weights = read(fid["data"]["weights"])
        weightsEE = read(fid["data"]["weightsEE"])
        weightsEI = read(fid["data"]["weightsEI"])
        weightsIE = read(fid["data"]["weightsIE"])
        times = read(fid["data"]["times"])
        close(fid)

        Ncells = size(weights)[1]
        Ne = 3000
        Ni2 = 750
        ipopmembers = findI2populations(weights, popmembers, iipop_len=27, Ni2=Ni2)
        # output_dir = string("./output_analysis/simulation_", sim, "/tests/")
        output_dir = string("./output_analysis/simulation_", sim)

        if sim_savedpath == "./networks_trained_spontaneous/"
            # seq_score = sequentialityScore(getPopulationRates(times, popmembers; interval=5_000:10_000))
            seq_score, _, sequence = findOptimalDecoder(times, popmembers, interval=5_000:20_000, dt=.125, seq_length=3)
            @info "For the ", sim_savedpath[3:end-1], " case sequentiality is: ", seq_score
            plotNetworkActivity(times, popmembers, ipopmembers; interval=5_000:20_000, seq_length=3, name=string("_testinActivitySpontaneous_", sim), output_dir=output_dir)
            plotWeightsEE(mean(weights[popmembers[:, :], popmembers[:, :]], dims=(1, 3))[1, :, 1, :], name="_testinEE", output_dir=output_dir, seq_length=3)
            plotWeightsEI(mean(weights[popmembers[:, :], ipopmembers[:, :]], dims=(1, 3))[1, :, 1, :], name="_testinEI", output_dir=output_dir, seq_length=3)
            plotWeightsIE(mean(weights[ipopmembers[:, :], popmembers[:, :]], dims=(1, 3))[1, :, 1, :], name="_testinIE", output_dir=output_dir, seq_length=3)
        else
            plotNetworkActivity(times, popmembers, ipopmembers; interval=5_000:20_000, seq_length=3, name=string("_testinActivityStimulation_", sim), output_dir=output_dir)
            # seq_score = sequentialityScore(getPopulationRates(times, popmembers; interval=5_000:10_000))
            seq_score, _, _ = findOptimalDecoder(times, popmembers, interval=5_000:20_000, dt=.125, seq_length=3)
            @info "For the ", sim_savedpath[3:end-1], " case sequentiality is: ", seq_score
        end
    end
end




##########################################################################################
############################# --- Compute sequentiality --- ##############################
##########################################################################################

Nsimulations = 4
Nseeds = 10
sim_savedpaths = ["./networks_trained_spontaneous/"]#, "./networks_trained_stimulation/"]

for sim_savedpath in sim_savedpaths
    for sim = 2:Nsimulations
        if sim_savedpath == "./networks_trained_spontaneous/"
            sim_name = string("network_", sim,"_spontaneous")
        else
            sim_name = string("network_", sim,"_stimulation.h5")
        end

        seq_scores = zeros(Nseeds)
        opt_params = zeros(2, Nseeds)
        
        for iseed = 1:Nseeds
            saved_net_path = string(sim_savedpath, "network_", sim)
            fid = h5open(joinpath(saved_net_path, string(sim_name, "_seed_", iseed, ".h5")), "r")
            popmembers = read(fid["data"]["popmembers"])
            weights = read(fid["data"]["weights"])
            times = read(fid["data"]["times"])
            close(fid)

            Ncells = size(weights)[1]
            Ne = 3000
            Ni2 = 250
            ipopmembers = findI2populations(weights, popmembers, iipop_len=27, Ni2=Ni2)
            output_dir = string("./output_analysis/simulation_", sim)

            # Test sequentiality
            seq_score, opt_param, sequence = findOptimalDecoder(times, popmembers, interval=10_000:50_000, dt=.125, seq_length=3)
            @info "For the ", sim_savedpath[3:end-1], " case ", iseed, " sequentiality is: ", seq_score
            seq_scores[iseed] = seq_score
            opt_params[:, iseed] .= opt_param 
        end

        if save_data
            # Save the data
            sim_savepath = "./analysis_data/sequentiality/"
            (!ispath(sim_savepath)) && (mkpath(sim_savepath))
            cd(sim_savepath[3:end-1])
            if sim_savedpath == "./networks_trained_spontaneous/"
                save_file_name = string("sequentiality_network_", sim, "_spontaneous.h5")
            else
                save_file_name = string("sequentiality_network_", sim, "_stimulation.h5")
            end
            fid = h5open(save_file_name,"w")
            g = create_group(fid,"data")
            g["scores"] = seq_scores
            g["params"] = opt_params
            close(fid)
            cd("../..")
        end
    end
end


##########################################################################################
##########################################################################################
##########################################################################################

using StatsBase

save_data = true
Npop = 12
ccov_range = -400:400
tt_range = 10_000:50_000

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
    
    ipopmembers = findI2populations(weights, popmembers, iipop_len=27)

    # Calculate population firing rates
    # pop_rates = getPopulationRates(times, popmembers, interval=tt_range, sigma=5.)
    # ipop_rates = getPopulationRates(times, ipopmembers, interval=tt_range, sigma=5.)

    pop_rates = getPopulationBinRates(times, popmembers; interval=tt_range)
    ipop_rates = getPopulationBinRates(times, ipopmembers; interval=tt_range)

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

############################################################################################
############################################################################################
############################################################################################

#####################################################################
# ___________ --- Make weight distribution plot --- ___________
#####################################################################
sim_savedpath = "./networks_trained/"
output_dir = "./output_analysis/"

# Choose the parameters
Ncells = 3750
Ni2 = 250
Npop = 12
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
        members = filter(i->(i>0), popmembers[:, ipop])    			# Get members of excitatory assembly
        ie_weights[:, ipop, sim_num] .= vec(sum(weights[members, Ncells-Ni2+1:Ncells], dims=1))  # Get sum of all weights projected from each E-assemble to each 2nd ipopulation neuron
        ie_perms[:, ipop, sim_num] .= sortperm(ie_weights[:, ipop, sim_num], rev=true)
        sorted_ie_weights[:, ipop, sim_num] .= ie_weights[round.(Int, ie_perms[:, ipop, sim_num]), ipop, sim_num]
        # fitted_params_p[:, ipop, sim_num] .= power_fit(x_vals, sorted_ie_weights[:, ipop, sim_num])
        fitted_params[:, ipop, sim_num] .= poly_fit(x_vals, sorted_ie_weights[:, ipop, sim_num], 7)
    end
end

plot_data_band_min = minimum(minimum(sorted_ie_weights, dims=3)[:, :, 1], dims=2)[:]
plot_data_avg = mean(mean(sorted_ie_weights, dims=3)[:, :, 1], dims=2)[:]
plot_data_band_max = maximum(maximum(sorted_ie_weights, dims=3)[:, :, 1], dims=2)[:]

# f = Figure(resolution=(720, 480))
# ax = Axis(f[1, 1], xlabel=L"I_2 \, \text{neuron index}", ylabel=L"\text{Synaptic strength (sum; pF)}",
#                     yticks=([50, 100, 150, 200], [L"50", L"100", L"150", L"200"]),
#                     xticks=([1, 27, 50, 100, 150, 200, 250], [L"1", L"\mathbf{27}", L"50", L"100", L"150", L"200", L"250"]),
#                     xlabelsize=24, ylabelsize=24, xticklabelsize=18, yticklabelsize=18, limits=(0.5, 250.5, 0.5, 249))
# band!(x_vals, plot_data_band_min, plot_data_band_max)
# lines!(x_vals, mean(fitted_params_p[1, :, :]) .* (x_vals .^ mean(fitted_params_p[2, :, :])), linewidth=4, color=:navajowhite, label=L"y = %$(round(mean(fitted_params[1, :, :]), digits=2))\,(\pm %$(round(std(fitted_params[1, :, :]), digits=2)))\,x^{%$(round(mean(fitted_params[2, :, :]), digits=2))\,(\pm %$(round(std(fitted_params[2, :, :]), digits=2)))}")
# # lines!(x_vals, mean(fitted_params[1, :, :]) .+ (x_vals .* mean(fitted_params[2, :, :])) .+ ((x_vals .^ 2) .* mean(fitted_params[3, :, :])) .+ ((x_vals .^ 3) .* mean(fitted_params[4, :, :])) .+ ((x_vals .^ 4) .* mean(fitted_params[5, :, :])) .+ ((x_vals .^ 5) .* mean(fitted_params[6, :, :])) .+ ((x_vals .^ 6) .* mean(fitted_params[7, :, :])) .+ ((x_vals .^ 7) .* mean(fitted_params[8, :, :])), linewidth=4, color=:gold3, label=L"7 \text{th degree polynomial}")
# vlines!(27, ymin=0, ymax=251, linestyle=(:dash, :loose), color=:gray20)
# text!(27, 240, text=L"\text{max slope change}", rotation=3*π/2, fontsize=22)
# axislegend(ax, position=(.95, .95), labelsize=20, titlesize=18, L"\mathbf{\text{Fitted power law:}}", titlehalign=:left)
# f
# (!isdir(output_dir)) && (mkpath(output_dir))
# save(joinpath(output_dir, string("power_fit_E_to_I2.png")), f)

f = Figure(resolution=(720, 480))
ax = Axis(f[1, 1], xlabel=L"I_2 \, \text{neuron index}", ylabel=L"\text{Synaptic strength (sum; pF)}",
                    yticks=([50, 100, 150, 200], [L"50", L"100", L"150", L"200"]),
                    xticks=([1, 27, 50, 100, 150, 200, 250], [L"1", L"\mathbf{27}", L"50", L"100", L"150", L"200", L"250"]),
                    xlabelsize=24, ylabelsize=24, xticklabelsize=18, yticklabelsize=18)#, limits=(0.5, 250.5, 0.5, 249))
band!(x_vals, plot_data_band_min, plot_data_band_max)
# errorbars!(x_vals, plot_data_avg, plot_data_avg .- plot_data_band_min, plot_data_band_max .- plot_data_avg)
lines!(x_vals, plot_data_avg, linewidth=4, color=:gold3, label=L"7 \text{th degree polynomial}")
vlines!(25, ymin=0, ymax=251, linestyle=(:dash, :loose), color=:gray20)
text!(25, 240, text=L"\text{max slope change}", rotation=3*π/2, fontsize=22)
axislegend(ax, position=(.95, .95), labelsize=20, titlesize=18, L"\mathbf{\text{Fitted power law:}}", titlehalign=:left)
f
(!isdir(output_dir)) && (mkpath(output_dir))
save(joinpath(output_dir, string("power_fit_E_to_I2.png")), f)


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

ipopmembers = findI2populations(weights, popmembers, iipop_len=25)
popmembers = transpose(popmembers)

# Create EI plot_data
plot_data = zeros(1_000, Npop, Npop)
for ipop = 1:Npop
    for iipop = 1:Npop
        plot_data[:, ipop, iipop] = mean(weightsEI[ipopmembers[:, iipop] .- (Ncells-Ni2), ipop, :], dims=1)
    end
end


fig = Figure()
ax = Axis(fig[1, 1])

pre = 1
for post = 1:Npop
    if post == pre
        lines!(ax, plot_data[:, pre, post])
    else
        lines!(ax, plot_data[:, pre, post])
    end
end
fig

# Create IE plot_data
plot_data = zeros(1_000, Npop, Npop)
for ipop = 1:Npop
    for iipop = 1:Npop
        plot_data[:, ipop, iipop] = mean(weightsIE[ipopmembers[:, iipop] .- Ne, ipop, :], dims=1)
    end
end


plot_data ./= 1250 

fig = Figure()
ax = Axis(fig[1, 1])
pre = 1
for post = 1:Npop
    if post == pre-1
        lines!(ax, plot_data[:, pre, post], color=:blue)
    elseif post == pre
        lines!(ax, plot_data[:, pre, post], color=:red)
    else
        lines!(ax, plot_data[:, pre, post], color=:gray)
    end
end
fig
