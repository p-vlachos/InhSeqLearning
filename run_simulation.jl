using Pkg
Pkg.activate(".")
using InhSequences
using Statistics
using HDF5

sim_savedpath = "./networks_trained/"

# Set the parameters for simulation
loadnet = true				# Choose between simulating an existing or a novel network
savenet = true				# Save the network after stimulation

# random_seeds = [2061, 5987, 3642, 9465, 1837, 8174, 9729, 8537, 1835, 6935]	# Train seeds
random_seeds = [6487, 3791, 6482, 6485, 9316]	# Knockout seeds

Nsimulations = length(random_seeds)			# Number of simulation to run

stim_modes = ["spontaneous", "stimulation"]	# Choose between spontaneous or brief stimulation (only for loaded networks)

for sim_num = 1:1
	# sim_name = string("network_", sim_num,".h5")
	# sim_name = string("network_", sim_num,"_i1STDP-knockout.h5")	
	sim_name = string("network_", sim_num,"_pre-train.h5")
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
	for sim_num = 1:1
		if loadnet
			T = 20_000
			if stim_mode == "spontaneous"
				stim = zeros(4, 12)						# Spontaneous activity (no stimulation)
			elseif stim_mode == "stimulation"
				stim = makeStimSeq_brief(T, Npop=12, seq_len=3, seq_num=4, stim_rate=4., randomize=false)	# Brief stimulation of 1st assembly in sequences
			else
				stim = [1. 5. 9. 13. 17. 1. 5. 9. 13. 17. 1. 5. 9. 13. 17. 1. 5. 9. 13. 17.;
						2_001. 2_201. 2_401. 2_601. 2_801. 7_001. 7_201. 7_401. 7_601. 7_801. 12_001. 12_201. 12_401. 12_601. 12_801. 17_001. 17_201. 17_401. 17_601. 17_801.;
						2_030. 2_230. 2_430. 2_630. 2_830. 7_030. 7_230. 7_430. 7_630. 7_830. 12_030. 12_230. 12_430. 12_630. 12_830. 17_030. 17_230. 17_430. 17_630. 17_830.;
						8. 8. 8. 8. 8. 8. 8. 8. 8. 8. 8. 8. 8. 8. 8. 8. 8. 8. 8. 8.]
			end
		else
			T = 1_000_000
			stim = makeStimSeq(4, 20, seq_len=3)		# Full train 5 sequences	#200 disjoint
		end

		sim_name = string("network_", sim_num,".h5")
		# sim_name = string("network_", sim_num,"_i1STDP-knockout.h5")

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
end




# # Define here manipulations
# seq_len=4
# ipopmembers = findI2populations(weights_old, popmembers, iipop_len=27)
# Npop = size(popmembers)[2]
# Ncells = size(weights_old)[1]
# Ne = round(Int, Ncells*0.8) 
# Ni2 = 250
# last_elements = seq_len:seq_len:Npop
# # @. weights_old[(Ncells-Ni2+1):Ncells, 1:Ne] *= .8		# I2-to-E strength
# # @. weights_old[(Ne+1):(Ncells-Ni2), 1:Ne] *= 1.2		# I1-to-E strength
# for pop = 1:Npop
# 	# @. weights_old[ipopmembers[:, pop], popmembers[:, pop]] *= .5		# I2pop-to-Epop strength (within)
# 	# @. weights_old[(Ne+1):(Ncells-Ni2), pop] *= .5		# I1-to-Epop strength
# 	for ipop = 1:Npop
# 		# if ipop==pop
# 			# continue
# 		# end
# 		# @. weights_old[ipopmembers[:, ipop], popmembers[:, pop]] *= .5		# I2pop-to-Epop strength (between all)
# 		# if ipop == pop-1 && !(ipop in last_elements)			# This targets the most weakly inhibited (I2-to-E)
# 		# 	@. weights_old[ipopmembers[:, ipop], popmembers[:, pop]] *= 2.
# 		# end
# 		# if ipop-1 == pop && !(pop in last_elements)			# This targets the most strongly inhibited (I2-to-E)
# 		# 	@. weights_old[ipopmembers[:, ipop], popmembers[:, pop]] *= 2.
# 		# end
# 		# if ipop == pop			
# 		# 	@. weights_old[popmembers[:, pop], ipopmembers[:, ipop]] *= 1.2		# This targets the E-to-I2
# 		# 	# @. weights_old[ipopmembers[:, ipop], popmembers[:, pop]] *= 2.		# This targets the I2-to-E
# 		# end
# 	end
# end