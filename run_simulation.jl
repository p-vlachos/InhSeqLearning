using Pkg
Pkg.activate(".")
using InhSequences
using Statistics
using HDF5

sim_savedpath = "./networks_trained/"

# Set the parameters for simulation
loadnet = false				# Choose between simulating an existing or a novel network
savenet = true				# Save the network after stimulation

random_seeds = [2061, 5987, 3642, 9465, 1837, 8174, 9729, 8537, 1835, 6935]	# Train/test seeds
# random_seeds = [6487, 3791, 6482, 6485, 9316]	# Knockout seeds

Nsimulations = length(random_seeds)			# Number of simulation to run

stim_modes = ["spontaneous"]#, "stimulation"]	# Choose between spontaneous or brief stimulation (only for loaded networks)

for sim_num = 1:1
	# sim_name = string("network_", sim_num,".h5")
	# sim_name = string("network_", sim_num,"_i1STDP-knockout.h5")	
	# sim_name = string("network_", sim_num,"_pre-train.h5")
	sim_name = string("network_", sim_num,"_train.h5")
	if loadnet
		T = 10_000
		if stim_mode == "stimulation"
			# Brief stimulation of 1st assembly in sequences
			stim = makeStimSeq_brief(T, Npop=12, seq_len=3, seq_num=4, stim_rate=4., randomize=false)
		else
			# Spontaneous activity (no stimulation)
			stim = zeros(4, 12)
		end
	else
		# Full train 4 sequences x 3 assemblies
		# T = 1_000_000
		T = 50_000
		stim = makeStimSeq(4, 20, seq_len=3)
	end

	if loadnet
		sim_savepath = string("./networks_trained_", stim_mode, "/")
		fid = h5open(joinpath(sim_savedpath, sim_name), "r")
		weights_old = read(fid["data"]["weights"])
		popmembers = read(fid["data"]["popmembers"])
		close(fid)
		times, weights, popmembers, weightsEE, weightsIE, weightsEI = sim(stim, weights_old, popmembers, T, random_seed=random_seeds[sim_num])
	else
		sim_savepath = string("./networks_trained/")
		times, weights, popmembers, weightsEE, weightsIE, weightsEI = simnew(stim, T, random_seed=random_seeds[sim_num])
	end

	if savenet
		(!ispath(sim_savepath)) && (mkpath(sim_savepath))
		if loadnet
			cd(string("networks_trained_", stim_mode))
			fid = h5open(string(sim_name[1:end-3], "_", stim_mode, ".h5"), "w")
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

# ______________________________ --- END --- ______________________________
for stim_mode in stim_modes
	for sim_num = 2:Nsimulations
		for iseed in eachindex(random_seeds)
			T = 50_000
			if stim_mode == "stimulation"
				stim = makeStimSeq_brief(T, Npop=12, seq_len=3, seq_num=4, stim_rate=4., randomize=false)	# Brief stimulation of 1st assembly in sequences
			else
				stim = zeros(4, 12)		# Spontaneous activity (no stimulation)
			end

			sim_name = string("network_", sim_num,".h5")

			# Load data
			sim_savepath = string("./networks_trained_", stim_mode, "/network_", sim_num)
			fid = h5open(joinpath(sim_savedpath, sim_name), "r")
			weights_old = read(fid["data"]["weights"])
			popmembers = read(fid["data"]["popmembers"])
			close(fid)

			times, weights, popmembers, _, _, _ = sim(stim, weights_old, popmembers, T, random_seed=random_seeds[iseed])

			if savenet
				(!ispath(sim_savepath)) && (mkpath(sim_savepath))
				cd(sim_savepath)
				fid = h5open(string(sim_name[1:end-3], "_", stim_mode, "_seed_", iseed, ".h5"), "w")
				g = create_group(fid,"data")
				g["popmembers"] = popmembers
				g["weights"] = weights
				g["times"] = times
				close(fid)
				cd("../..")
			end
		end
	end
end