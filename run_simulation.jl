using Pkg
Pkg.activate(".")
using InhSequences
using Statistics
using HDF5

# include("quantification.jl")
sim_savedpath = "./networks_trained/"

# Set the parameters for simulation
loadnet = false				# Choose between simulating an existing or a novel network
savenet = true				# Save the network after stimulation

Nsimulations = 10			# Number of simulation to run
stim_mode = "spontaneous"	# Choose between spontaneous or brief stimulation (only for loaded networks)
random_seeds = [2061, 5987, 3642, 9465, 1837, 6487, 3791, 6482, 6485, 9316]

for sim_num = 1:Nsimulations
	if loadnet
		T = 50_000
		if stim_mode == "spontaneous"
			stim = zeros(4, 20)						# Spontaneous activity (no stimulation)
		else
			stim = makeStimSeq_brief(T)				# Brief stimulation of 1st assembly in sequences (random order)
		end
	else
		T = 1_500_000
		stim = makeStimSeq(5, 20, seq_len=4)		# Full train 5 sequences
	end

	sim_name = string("network_", sim_num,".h5")

	if loadnet
		sim_savepath = string("./networks_trained_", sim_mode, "/")
		fid = h5open(joinpath(sim_savedpath, sim_name), "r")
		popmembers_old = read(fid["data"]["popmembers"])
		weights_old = read(fid["data"]["weights"])
		close(fid)

		times, weights, popmembers, weightsEE, weightsIE, weightsEI = sim(stim, weights_old, popmembers, T, random_seed=random_seeds[sim_num])
		# times, weights, popmembers = sim(stim, weights_old, popmembers, T)
	else
		sim_savepath = string("./networks_trained/")
		times, weights, popmembers, weightsEE, weightsIE, weightsEI = simnew(stim, T, random_seed=random_seeds[sim_num])
		# times, weights, popmembers = simnew(stim, T)
	end

	if savenet
		(!ispath(sim_savepath)) && (mkpath(sim_savepath))
		if loadnet
			cd(string("networks_trained_", stim_mode))
			fid = h5open(string(sim_name, "_", stim_mode, ".h5"), "w")
		else
			cd("networks_trained")
			fid = h5open(sim_name,"w")
		end
		g = create_group(fid,"data")
		g["popmembers"] = popmembers
		g["weights"] = weights
		g["weightsEE"] = weightsEE
		g["weightsEI"] = weightsEI
		g["weightsIE"] = weightsIE
		g["times"] = times
		close(fid)
		cd("..")
	end
end
