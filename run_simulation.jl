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

stim_modes = ["spontaneous", "stimulation"]	# Choose between spontaneous or brief stimulation (only for loaded networks)
random_seeds = [2061, 5987, 3642, 9465, 1837, 6487, 3791, 6482, 6485, 9316, 8174, 9729, 8537, 1835, 6935, 7164, 9386] #[3791, 2061, 5987, 3642, 9465, 1837, 6487, 6482, 6485, 9316] #3791 # [2061, 5987, 3642, 9465, 1837, 6487, 3791, 6482, 6485, 9316]
Nsimulations = length(random_seeds)			# Number of simulation to run


for stim_mode in stim_modes
	for sim_num = 6:6
		if loadnet
			# T = 50_000
			T = 10_000
			if stim_mode == "spontaneous"
				stim = zeros(4, 16)						# Spontaneous activity (no stimulation)
			elseif stim_mode == "stimulation"
				# stim = makeStimSeq_brief(T, stim_rate=8., randomize=false)				# Brief stimulation of 1st assembly in sequences (random order)
				stim = makeStimSeq_brief(T, Npop=12, seq_len=4, seq_num=3, stim_rate=8., randomize=false)
			else
				stim = [1. 5. 9. 13. 17. 1. 5. 9. 13. 17. 1. 5. 9. 13. 17. 1. 5. 9. 13. 17.;
						2_001. 2_201. 2_401. 2_601. 2_801. 7_001. 7_201. 7_401. 7_601. 7_801. 12_001. 12_201. 12_401. 12_601. 12_801. 17_001. 17_201. 17_401. 17_601. 17_801.;
						2_030. 2_230. 2_430. 2_630. 2_830. 7_030. 7_230. 7_430. 7_630. 7_830. 12_030. 12_230. 12_430. 12_630. 12_830. 17_030. 17_230. 17_430. 17_630. 17_830.;
						8. 8. 8. 8. 8. 8. 8. 8. 8. 8. 8. 8. 8. 8. 8. 8. 8. 8. 8. 8.]
			end
		else
			T = 1_000_000
			# stim = makeStimSeq(5, 20, seq_len=4)		# Full train 5 sequences
			stim = makeStimSeq(3, 20, seq_len=4)		# Full train 5 sequences	#200 disjoint
			# T = 1_000_000
			# stim = makeStimSeq(3, 20, seq_len=4)		# Full train 5 sequences
			# T = 20_000
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

			# Define here manipulations
			seq_len=4
			ipopmembers = findI2populations(weights_old, popmembers, iipop_len=27)
			Npop = size(popmembers)[2]
			Ncells = size(weights_old)[1]
			Ne = round(Int, Ncells*0.8) 
			Ni2 = 250
			last_elements = seq_len:seq_len:Npop
			# @. weights_old[(Ncells-Ni2+1):Ncells, 1:Ne] *= 2.		# I2-to-E strength
			# @. weights_old[(Ne+1):(Ncells-Ni2), 1:Ne] *= 0.5		# I1-to-E strength
			for pop = 1:Npop
				# @. weights_old[ipopmembers[:, pop], popmembers[:, pop]] *= .5		# I2pop-to-Epop strength (within)
				# @. weights_old[(Ne+1):(Ncells-Ni2), pop] *= .5		# I1-to-Epop strength
				for ipop = 1:Npop
					# if ipop==pop
						# continue
					# end
					# @. weights_old[ipopmembers[:, ipop], popmembers[:, pop]] *= .5		# I2pop-to-Epop strength (between all)
					# if ipop == pop-1 && !(ipop in last_elements)			# This targets the most weakly inhibited (I2-to-E)
					# 	@. weights_old[ipopmembers[:, ipop], popmembers[:, pop]] *= .5
					# end
					# if ipop-1 == pop && !(pop in last_elements)			# This targets the most strongly inhibited (I2-to-E)
					# 	@. weights_old[ipopmembers[:, ipop], popmembers[:, pop]] *= .5
					# end
					# if ipop == pop			
					# 	# @. weights_old[popmembers[:, pop], ipopmembers[:, ipop]] *= 2.		# This targets the E-to-I2
					# 	@. weights_old[ipopmembers[:, ipop], popmembers[:, pop]] *= .5		# This targets the I2-to-E
					# end
				end
			end

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
end
