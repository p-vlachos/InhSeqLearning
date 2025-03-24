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

Nsimulations = 10
Nseeds = 10
sim_savedpaths = ["./networks_trained_spontaneous/", "./networks_trained_stimulation/"]

for sim_savedpath in sim_savedpaths
    for sim = 1:Nsimulations
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



# --- Load Data ---
sim_name = string("network_4.h5")
sim_savedpath = string("./networks_trained/")
fid = h5open(joinpath(sim_savedpath, sim_name), "r")
popmembers = read(fid["data"]["popmembers"])
weights = read(fid["data"]["weights"])
# times = read(fid["data"]["times"])
weightsEE = read(fid["data"]["weightsEE"])
weightsEI = read(fid["data"]["weightsEI"])
weightsIE = read(fid["data"]["weightsIE"])
close(fid)


ipopmembers = findI2populations(weights, popmembers, iipop_len=Ni_members)

# E
fig = Figure()
ax = Axis(fig[1, 1])
for ipop = 1:Npop
    lines!(ax, weightsEE[ipop, ipop, 1:500], color=ColorSchemes.Blues[7])
end
fig

# EI
fig = Figure()
ax = Axis(fig[1, 1])
for ipop = 1:Npop
    lines!(ax, mean(weightsEI[ipopmembers[:, ipop] .- (Ncells-Ni2), ipop, 1:500], dims=1)[:], color=ColorSchemes.Blues[5])
end
fig

# I2
fig = Figure()
ax = Axis(fig[1, 1])
for ipop = 2:Npop
    (mod(ipop, seq_length) == 1) && (continue)
    lines!(ax, mean(weightsIE[ipopmembers[:, ipop-1] .- (Ne), ipop, 1:500], dims=1)[:], color=ColorSchemes.Reds[6])
end
fig

for ipop = 2:Npop
    (mod(ipop, seq_length) == 1) && (continue)
    lines!(ax, mean(weightsIE[ipopmembers[:, ipop] .- (Ne), ipop-1, 1:500], dims=1)[:], color=:black)
end
fig

for ipop = 1:Npop
    # (mod(ipop, seq_length) == 1) && (continue)
    lines!(ax, mean(weightsIE[ipopmembers[:, ipop] .- (Ne), ipop, 1:500], dims=1)[:], color=:gray70)
end
fig

# I1
fig = Figure()
ax = Axis(fig[1, 1])
for ipop = 2:Npop
    (mod(ipop, seq_length) == 1) && (continue)
    lines!(ax, mean(weightsIE[1:(Ncells-Ne-Ni2), ipop, 1:500], dims=1)[:], color=ColorSchemes.Reds[9])
end
fig