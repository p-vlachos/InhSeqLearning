using Pkg
Pkg.activate(".")
using InhSequences
using Statistics
using HDF5

# include("quantification.jl")
sim_savedpath = "./networks_trained/"

# Set the parameters for simulation
loadnet = true				# Choose between simulating an existing or a novel network
savenet = true				# Save the network after stimulation

stim_mode = "stimulation"	# Choose between spontaneous or brief stimulation (only for loaded networks)
random_seeds = 3791 # [2061, 5987, 3642, 9465, 1837, 6487, 3791, 6482, 6485, 9316]
Nsimulations = length(random_seeds)			# Number of simulation to run


for sim_num = 7:7#Nsimulations
	if loadnet
		# T = 50_000
		T = 20_000
		if stim_mode == "spontaneous"
			stim = zeros(4, 20)						# Spontaneous activity (no stimulation)
		elseif stim_mode == "stimulation"
			stim = makeStimSeq_brief(T, stim_rate=8., randomize=false)				# Brief stimulation of 1st assembly in sequences (random order)
		else
			stim = [1. 5. 9. 13. 17. 1. 5. 9. 13. 17. 1. 5. 9. 13. 17. 1. 5. 9. 13. 17.;
					2_001. 2_201. 2_401. 2_601. 2_801. 7_001. 7_201. 7_401. 7_601. 7_801. 12_001. 12_201. 12_401. 12_601. 12_801. 17_001. 17_201. 17_401. 17_601. 17_801.;
					2_030. 2_230. 2_430. 2_630. 2_830. 7_030. 7_230. 7_430. 7_630. 7_830. 12_030. 12_230. 12_430. 12_630. 12_830. 17_030. 17_230. 17_430. 17_630. 17_830.;
					8. 8. 8. 8. 8. 8. 8. 8. 8. 8. 8. 8. 8. 8. 8. 8. 8. 8. 8. 8.]
		end
	else
		# T = 1_500_000
		stim = makeStimSeq(5, 20, seq_len=4)		# Full train 5 sequences
		T = 20_000
		# stim = zeros(4, 20)		# Full train 5 sequences
	end

	# sim_name = string("network_pre_train.h5")
	sim_name = string("network_", sim_num,".h5")
	# sim_name = string("network_i2STDP_knockout.h5")

	if loadnet
		sim_savepath = string("./networks_trained_", stim_mode, "/")
		fid = h5open(joinpath(sim_savedpath, sim_name), "r")
		weights_old = read(fid["data"]["weights"])
		popmembers = read(fid["data"]["popmembers"])
		close(fid)

		# times, weights, popmembers, weightsEE, weightsIE, weightsEI = sim(stim, weights_old, popmembers, T, random_seed=random_seeds[sim_num])
		times, weights, popmembers = sim(stim, weights_old, popmembers, T)
	else
		sim_savepath = string("./networks_trained/")
		# times, weights, popmembers, weightsEE, weightsIE, weightsEI = simnew(stim, T, random_seed=random_seeds[sim_num])
		times, weights, popmembers = simnew(stim, T)
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
		# g["weightsEE"] = weightsEE
		# g["weightsEI"] = weightsEI
		# g["weightsIE"] = weightsIE
		g["times"] = times
		close(fid)
		cd("..")
	end
end
