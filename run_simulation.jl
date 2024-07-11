using Pkg
Pkg.activate(".")
using InhSequences
using Statistics
using HDF5

include("quantification.jl")

load = true
doplot = false
loadtrained = load
# savenet = !load
savenet = true

# load = true
# doplot = false
# loadtrained = false
# savenet = true

# Stimulus matrix is in the following format:
# column 1: index of stimulated population
# columns 2/3: start/stop time of stimulus
# column 4: rate of external drive (in kHz)

# for sim_num = 1:20

for sim_num = 7:10 #20
	if !load
		T = 2_000_000
		stim = make_seq(5, 20, seq_len=4)		# Full train 5 sequences
	else
		T = 50_000
		stim = zeros(4, 20)						# Spontaneous activity (no stimulation)
		# stim .= [1 1001 1030 8; 1 1201 1230 8; 1 1401 1430 8; 1 1601 1630 8; 1 1801 1830 8;
		# 		5 2001 2030 8; 5 2201 2230 8; 5 2401 2430 8; 5 2601 2630 8; 5 2801 2830 8;
		# 		9 3001 3030 8; 9 3201 3230 8; 9 3401 3430 8; 9 3601 3630 8; 9 3801 3830 8;
		# 		13 4001 4030 8; 13 4201 4230 8; 13 4401 4430 8; 13 4601 4630 8; 13 4801 4830 8;
		# 		17 5001 5030 8; 17 5201 5230 8; 17 5401 5430 8; 17 5601 5630 8; 17 5801 5830 8]'
		# stim .= transpose(stim)
		# stim = makeSeqStim(T)
	end
	# stim = make_seq(1, 20, seq_len=20)	# Single sequence

	sim_name = string("network_", sim_num,".h5")
	sim_savedpath = "./networks_trained/"
	# sim_savepath = "./networks_trained_spontaneous/"
	sim_savepath = "./networks_trained_stimulation/"
	output_dir = "./output_trained/"

	if loadtrained
		fid = h5open(joinpath(sim_savedpath, sim_name), "r")
		popmembers_old = read(fid["data"]["popmembers"])
		weights_old = read(fid["data"]["weights"])
		close(fid)

		popmembers = zeros(Int, size(popmembers_old'))
		popmembers .= round.(Int, popmembers_old')
		times, ns, Ne, Ncells, T, new_weights, weightsEE, weightsIE, weightsEI, popmembers = sim(stim, weights_old, popmembers, T)

		Nmaxmembers = 300
		Ni2 = 250
		dt = .1
		Npop = size(popmembers)[1]
		Npop == 9 ? Nseq = 3 : Nseq = 5
		Nsteps = round(Int, T/dt)
		fr_threshold = 4   # Spiking threshold for assemblies (Hz; Chenkov:30 spikes/sec)
		fRate_window = 5
		time_limit = 200
		overlap_limit = 20
		Npop == 9 ? seq_length = 3 : seq_length = 4

		sim_outDir = string(output_dir, "/sim_", sim_num,"/")
	else
		((times, ns, Ne, Ncells, T, new_weights, weightsEE, weightsIE, weightsEI), popmembers) = simnew(stim, T)
	end

	if savenet
		(!ispath(sim_savepath)) && (mkpath(sim_savepath))
		cd("networks_trained_stimulation")
		# cd("networks_trained_spontaneous")
		# cd("networks")
		if loadtrained
			# fid = h5open(string(sim_name, "_spontaneous.h5"), "w")
			fid = h5open(string(sim_name, "_stimulation.h5"), "w")
		else
			fid = h5open(sim_name,"w")
		end
		g = create_group(fid,"data")
		g["popmembers"] = popmembers
		g["weights"] = new_weights
		g["weightsEE"] = weightsEE
		g["weightsEI"] = weightsEI
		g["weightsIE"] = weightsIE
		g["times"] = times
		close(fid)
		cd("..")
	end

	@info string("mean excitatory firing rate: ", mean(1000 * ns[1:Ne] / T), " Hz")
	@info string("mean inhibitory firing rate: ", mean(1000 * ns[(Ne + 1):Ncells] / T), " Hz")
end