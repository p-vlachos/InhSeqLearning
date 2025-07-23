using Pkg
Pkg.activate(".")
using InhSequences
using Statistics
using HDF5

sim_savedpath = "./networks_trained_original/"

# Set the parameters for simulation
loadnet = true				# Choose between simulating an existing or a novel network
savenet = true				# Save the network after stimulation

random_seeds = [2061, 5987, 3642, 9465, 1837, 8174, 9729, 8537, 1835, 6935, 9316, 8425, 7167, 9921, 2785, 6309, 1482]	# Train/test seeds	(no. 11 is for showing decay)
# random_seeds = [6487, 3791, 6482, 6485]	# Knockout seeds

Nsimulations = length(random_seeds)			# Number of simulation to run

stim_modes = ["spontaneous", "stimulation"]	# Choose between spontaneous or brief stimulation (only for loaded networks)

# sim_names = ["network_eiSTDP-knockout.h5", "network_i2STDP-knockout.h5", "network_i1STDP-knockout.h5"]
for stim_mode in stim_modes
for sim_num = 1:1
	# sim_name = string("network_", sim_num,".h5")
	# sim_name = string("network_", sim_num,"_eiSTDP-knockout.h5")
	# sim_name = string("network_", sim_num,"_i2STDP-knockout.h5")
	# sim_name = string("network_", sim_num,"_i1STDP-knockout.h5")
	sim_name = string("network_", sim_num,"_pre-train.h5")
	# sim_name = string("network_", sim_num,"_train.h5")
	if loadnet
		T = 20_000
		if stim_mode == "stimulation"
			# Brief stimulation of 1st assembly in sequences
			# stim = makeStimSeq_brief(T, Npop=12, seq_len=3, seq_num=4, stim_rate=4., randomize=false)
			stim = makeStimSeq_brief(T, Npop=20, seq_len=4, seq_num=5, stim_rate=8., randomize=false)
		else
			# Spontaneous activity (no stimulation)
			stim = zeros(4, 12)
		end
	else
		# # Full train 4 sequences x 3 assemblies
		# T = 1_000_000
		# stim = makeStimSeq(4, 20, seq_len=3)
		# Full train 5 sequences x 4 assemblies
		T = 1_500_000
		# T = 40_000
		stim = makeStimSeq(5, 20, seq_len=4)
	end

	if loadnet
		sim_savepath = string("./networks_trained_original_", stim_mode, "/")
		fid = h5open(joinpath(sim_savedpath, sim_name), "r")
		weights_old = read(fid["data"]["weights"])
		popmembers = read(fid["data"]["popmembers"])
		close(fid)
		times, weights, popmembers, weightsEE, weightsIE, weightsEI = sim(stim, weights_old, popmembers, T, random_seed=random_seeds[sim_num])
	else
		sim_savepath = string("./networks_trained_original/")
		times, weights, popmembers, weightsEE, weightsIE, weightsEI = simnew(stim, T, random_seed=random_seeds[sim_num])
	end

	if savenet
		(!ispath(sim_savepath)) && (mkpath(sim_savepath))
		if loadnet
			cd(string("networks_trained_original_", stim_mode))
			fid = h5open(string(sim_name[1:end-3], "_", stim_mode, ".h5"), "w")
		else
			cd("networks_trained_original")
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

####################################################################################
#### _________________________ --- Manipulations --- __________________________ ####
####################################################################################
include("plots.jl")

stim_modes = ["stimulation"]#, "stimulation"]	# Choose between spontaneous or brief stimulation (only for loaded networks)

manipulation_type = ["Eass-Eass_within", "Eass-Eass_between", "Eass-Erest", "Erest-Eass", # E-to-E manipulations
					"Eass-I1", "Erest-I1", "Eass-I2ass(within)", "Eass-I2ass(next)", "Eass-I2ass(pre)", "Eass-I2", "Erest-I2", "E-I2", # E-to-I manipulations
					"I1-Eass", "I1-Erest", "I2ass-Eass(within)", "I2ass-Eass(next)", "I2ass-Eass(pre)", "I2-Eass", "I2-Erest", "I2-E", # I-to-E manipulations
					"I1-I1", "I1-I2", "I2-I2", "I2-I1"] # I-to-I manipulations
manipulation_value = [0.7, 0.85, 1.15, 1.3]

# Define parameters
T = 20_000
Npop = 20
Ncells = 5000
Ne = 4000
Ni2 = 500
seq_length = 4

sim_name = string("network_1.h5")
fid = h5open(joinpath(sim_savedpath, sim_name), "r")
popmembers = read(fid["data"]["popmembers"])
weights_old = read(fid["data"]["weights"])
close(fid)
ipopmembers = findI2populations(weights_old, popmembers, iipop_len=50)
restcells = deleteat!(map(*, ones(Int, Ne), range(1, stop=Ne)), sort(unique(popmembers))[2:end])

for stim_mode in stim_modes
	if stim_mode == "stimulation"
		stim = makeStimSeq_brief(T, Npop=20, seq_len=4, seq_num=5, stim_rate=8., randomize=false)	# Brief stimulation of 1st assembly in sequences
	else
		stim = zeros(4, 20)		# Spontaneous activity (no stimulation)
	end
	sim_savepath = string("./networks_trained_original_testing_", stim_mode, "/")
	for (ind, mant) in enumerate(manipulation_type)
		if ind < 24
			continue
		end
		for manv in manipulation_value
			fid = h5open(joinpath(sim_savedpath, sim_name), "r")
			weights_old = read(fid["data"]["weights"])
			close(fid)
			# Apply manipulation
			if ind == 1			# Eass-to-Eass (within)
				for ipop = 1:Npop
					emembers = filter(i->i>0, popmembers[:, ipop])
					weights_old[emembers, emembers] *= manv
				end
			elseif ind == 2		# Eass-to-Eass (between)
				for ipop = 1:Npop
					for iipop = 1:Npop
						if ipop != iipop
							emembers = filter(i->i>0, popmembers[:, ipop])
							iemembers = filter(i->i>0, popmembers[:, iipop])
							weights_old[emembers, iemembers] *= manv
						end
					end
				end
			elseif ind == 3		# Eass-to-Erest
				for ipop = 1:Npop
					emembers = filter(i->i>0, popmembers[:, ipop])
					weights_old[emembers, restcells] *= manv
				end
			elseif ind == 4		# Erest-to-Eass
				for ipop = 1:Npop
					emembers = filter(i->i>0, popmembers[:, ipop])
					weights_old[restcells, emembers] *= manv
				end
			elseif ind == 5		# Eass-to-I1
				for ipop = 1:Npop
					emembers = filter(i->i>0, popmembers[:, ipop])
					weights_old[emembers, (Ne+1):(Ncells-Ni2)] *= manv
				end
			elseif ind == 6		# Erest-to-I1
				for ipop = 1:Npop
					weights_old[restcells, (Ne+1):(Ncells-Ni2)] *= manv
				end
			elseif ind == 7		# Eass-to-I2ass (within)
				for ipop = 1:Npop
					emembers = filter(i->i>0, popmembers[:, ipop])
					weights_old[emembers, ipopmembers[:, ipop]] *= manv
				end
			elseif ind == 8		# Eass-to-I2ass (next)
				for ipop = 1:Npop
					if ipop in 4:seq_length:Npop
						continue
					end
					emembers = filter(i->i>0, popmembers[:, ipop])
					weights_old[emembers, ipopmembers[:, (ipop+1)]] *= manv
				end
			elseif ind == 9		# Eass-to-I2ass (pre)
				for ipop = 1:Npop
					if ipop in 1:seq_length:Npop
						continue
					end
					emembers = filter(i->i>0, popmembers[:, ipop])
					weights_old[emembers, ipopmembers[:, (ipop-1)]] *= manv
				end
			elseif ind == 10	# Eass-to-I2
				for ipop = 1:Npop
					emembers = filter(i->i>0, popmembers[:, ipop])
					weights_old[emembers, (Ncells-Ni2+1):Ncells] *= manv
				end
			elseif ind == 11	# Erest-to-I2
				for ipop = 1:Npop
					weights_old[restcells, (Ncells-Ni2+1):Ncells] *= manv
				end
			elseif ind == 12	# E-to-I2
				weights_old[1:Ne, (Ncells-Ni2+1):Ncells] *= manv
			elseif ind == 13	# I1-to-Eass
				for ipop = 1:Npop
					emembers = filter(i->i>0, popmembers[:, ipop])
					weights_old[(Ne+1):(Ncells-Ni2), emembers] *= manv
				end
			elseif ind == 14	# I1-to-Erest
				weights_old[(Ne+1):(Ncells-Ni2), restcells] *= manv
			elseif ind == 15	# I2ass-to-Eass (within)
				for ipop = 1:Npop
					emembers = filter(i->i>0, popmembers[:, ipop])
					weights_old[ipopmembers[:, ipop], emembers] *= manv
				end
			elseif ind == 16	# I2(post)-to-Eass
				for ipop = 1:Npop
					if ipop in 4:seq_length:Npop
						continue
					end
					emembers = filter(i->i>0, popmembers[:, ipop])
					weights_old[ipopmembers[:, (ipop+1)], emembers] *= manv
				end
			elseif ind == 17	# I2(pre)-to-Eass
				for ipop = 1:Npop
					if ipop in 1:seq_length:Npop
						continue
					end
					emembers = filter(i->i>0, popmembers[:, ipop])
					weights_old[ipopmembers[:, (ipop-1)], emembers] *= manv
				end
			elseif ind == 18	# I2-to-Eass
				for ipop = 1:Npop
					emembers = filter(i->i>0, popmembers[:, ipop])
					weights_old[(Ncells-Ni2+1):Ncells, emembers] *= manv
				end
			elseif ind == 19	# I2-to-Erest
				weights_old[(Ncells-Ni2+1):Ncells, restcells] *= manv
			elseif ind == 20	# I2-to-E
				weights_old[(Ncells-Ni2+1):Ncells, 1:Ne] *= manv
			elseif ind == 21	# I1-to-I1
				weights_old[(Ne+1):(Ncells-Ni2), (Ne+1):(Ncells-Ni2)] *= manv
			elseif ind == 22	# I1-to-I2
				weights_old[(Ne+1):(Ncells-Ni2), (Ncells-Ni2+1):Ncells] *= manv
			elseif ind == 23	# I2-to-I2
				weights_old[(Ncells-Ni2+1):Ncells, (Ncells-Ni2+1):Ncells] *= manv
			elseif ind == 24	# I2-to-I1
				weights_old[(Ncells-Ni2+1):Ncells, (Ne+1):(Ncells-Ni2)] *= manv
			end

			# Run simulation
			times, weights, popmembers, weightsEE, weightsIE, weightsEI = sim(stim, weights_old, popmembers, T, random_seed=random_seeds[1])

			if savenet
				(!ispath(sim_savepath)) && (mkpath(sim_savepath))
				cd(string("networks_trained_original_testing_", stim_mode))
				fid = h5open(string(sim_name[1:end-5], "_", mant, "_", manv, ".h5"), "w")
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

			output_dir = string("./output_analysis_original/simulation_original_", stim_mode, "_", ind, "/", mant)
			(!ispath(output_dir)) && (mkpath(output_dir))

			if stim_mode == "spontaneous"
				plotNetworkActivity(times, popmembers, ipopmembers; interval=5_000:20_000, seq_length=4, name=string("_testinActivitySpontaneous_", manv), output_dir=output_dir)
				avgweightsEE = zeros(Npop, Npop)
				avgweightsEI = zeros(Npop, Npop)
				avgweightsIE = zeros(Npop, Npop)
				for ipop = 1:Npop
					for iipop = 1:Npop
						avgweightsEE[ipop, iipop] = sum(weights[filter(i->i>0, popmembers[:, ipop]), filter(i->i>0, popmembers[:, iipop])]) / count(i->i>0, weights[filter(i->i>0, popmembers[:, ipop]), filter(i->i>0, popmembers[:, iipop])])
						avgweightsEI[ipop, iipop] = sum(weights[filter(i->i>0, popmembers[:, ipop]), ipopmembers[:, iipop]]) / count(i->i>0, weights[filter(i->i>0, popmembers[:, ipop]), ipopmembers[:, iipop]])
						avgweightsIE[ipop, iipop] = sum(weights[ipopmembers[:, ipop], filter(i->i>0, popmembers[:, iipop])]) / count(i->i>0, weights[ipopmembers[:, ipop], filter(i->i>0, popmembers[:, iipop])])
					end
				end
				plotWeightsEE(avgweightsEE, name=string("_testinEE_", manv), output_dir=output_dir, seq_length=4)
				plotWeightsEI(avgweightsEI, name=string("_testinEI_", manv), output_dir=output_dir, seq_length=4)
				plotWeightsIE(avgweightsIE, name=string("_testinIE_", manv), output_dir=output_dir, seq_length=4)
			else
				plotNetworkActivity(times, popmembers, ipopmembers; interval=5_000:20_000, seq_length=4, name=string("_testinActivityStimulation_", manv), output_dir=output_dir)
			end
		end
	end
end






# # ______________________________ --- END --- ______________________________
# for stim_mode in stim_modes
# 	for sim_num = 11:11#Nsimulations
# 		for iseed in eachindex(random_seeds)
# 			T = 50_000
# 			if stim_mode == "stimulation"
# 				stim = makeStimSeq_brief(T, Npop=12, seq_len=3, seq_num=4, stim_rate=4., randomize=false)	# Brief stimulation of 1st assembly in sequences
# 			else
# 				stim = zeros(4, 12)		# Spontaneous activity (no stimulation)
# 			end

# 			sim_name = sim_names[sim_num]
# 			# sim_name = string("network_", sim_num,".h5")

# 			# Load data
# 			sim_savepath = string("./networks_trained_", stim_mode, "/network_", sim_name[9:end-3])
# 			fid = h5open(joinpath(sim_savedpath, sim_name), "r")
# 			weights_old = read(fid["data"]["weights"])
# 			popmembers = read(fid["data"]["popmembers"])
# 			close(fid)

# 			times, weights, popmembers, _, _, _ = sim(stim, weights_old, popmembers, T, random_seed=random_seeds[iseed])

# 			if savenet
# 				(!ispath(sim_savepath)) && (mkpath(sim_savepath))
# 				cd(sim_savepath)
# 				fid = h5open(string(sim_name[1:end-3], "_", stim_mode, "_seed_", iseed, ".h5"), "w")
# 				g = create_group(fid,"data")
# 				g["popmembers"] = popmembers
# 				g["weights"] = weights
# 				g["times"] = times
# 				close(fid)
# 				cd("../..")
# 			end
# 		end
# 	end
# end