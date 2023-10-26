using Pkg
Pkg.activate(".")
using InhSequences
using Statistics
using PyPlot
using HDF5

include("quantification.jl")

load = true
doplot = load
loadtrained = load
savenet = !load

# load = true
# doplot = false
# loadtrained = false
# savenet = true

# Stimulus matrix is in the following format:
# column 1: index of stimulated population
# columns 2/3: start/stop time of stimulus
# column 4: rate of external drive (in kHz)

for sim_num = 1:20
	# stim = zeros(1, 4)						# Spontaneous activity (no stimulation)
	if !load
		T = 1_500_000
		stim = make_seq(5, 20, seq_len=4)		# Full train 5 sequences
	else
		T = 8_000
		stim = [1. 1001 1030 8; 1 1201 1230 8; 1 1401 1430 8; 1 1601 1630 8; 1 1801 1830 8;
				5 2001 2030 8; 5 2201 2230 8; 5 2401 2430 8; 5 2601 2630 8; 5 2801 2830 8;
				9 3001 3030 8; 9 3201 3230 8; 9 3401 3430 8; 9 3601 3630 8; 9 3801 3830 8;
				13 4001 4030 8; 13 4201 4230 8; 13 4401 4430 8; 13 4601 4630 8; 13 4801 4830 8;
				17 5001 5030 8; 17 5201 5230 8; 17 5401 5430 8; 17 5601 5630 8; 17 5801 5830 8]
	end
	# stim = make_seq(1, 20, seq_len=20)	# Single sequence

	sim_name = string("network_", sim_num,".h5")
	sim_savepath = "./networks/"
	output_dir = "./output/"

	if loadtrained
		fid = h5open(joinpath(sim_savepath, sim_name), "r")
		popmembers = read(fid["data"]["popmembers"])
		weights_old = read(fid["data"]["weights"])
		close(fid)
		times, ns, Ne, Ncells, T, new_weights = sim(stim, weights_old, popmembers, T)

		Nmaxmembers = 300
		Ni2 = 500
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

		if doplot
			Npop = size(popmembers, 1)
			Nmaxmembers = size(popmembers, 2)
			ipopsize = 25
			ipopmembers = findI2populations(new_weights, Npop, popmembers, iipop_len=ipopsize)
			restcells = deleteat!(map(*, ones(Int, 4000), range(1,stop=4000)), sort(unique(popmembers))[2:end])
			println("creating plot")
			figure(figsize=(4, 4))
			xlim(0, T)
			ylim(0, sum(popmembers .> 0)+length(restcells)+500+sum(ipopmembers .> 0))
			ylabel("Neuron")
			xlabel("Simulation Time (ms)")
			tight_layout()
			# Plot raster with the order of rows determined by population membership
			rowcount = 0
			for pp = 1:Npop
				print("\rpopulation ", pp)
				for cc = 1:Nmaxmembers
					(popmembers[pp, cc] < 1) && break
					rowcount += 1
					ind = popmembers[pp, cc]
					vals = times[ind, 1:ns[ind]]
					y = rowcount * ones(length(vals))
					scatter(vals, y, s = .3, c="blue", marker="o", linewidths=0)
				end
			end
			# Plot raster of cells not belonging to any population
			for cc in restcells
				rowcount += 1
				vals = times[cc, 1:ns[cc]]
				y = rowcount * ones(length(vals))
				scatter(vals, y, s = .3, c="green", marker="o", linewidths=0)
			end
			# Plot raster for interneurons (1st i-population)
			for cc = 4001:4500
				rowcount += 1
				vals = times[cc, 1:ns[cc]]
				y = rowcount * ones(length(vals))
				scatter(vals, y, s = .3, c="red", marker="o", linewidths=0)
			end
			# Plot raster for interneurons (2nd i-population)
			for pp = 1:Npop
				for cc = 1:ipopsize
					rowcount += 1
					ind = ipopmembers[pp, cc]
					vals = times[ind, 1:ns[ind]]
					y = rowcount * ones(length(vals))
					scatter(vals, y, s = .3, c="red", marker="o", linewidths=0)
				end
			end

			(!ispath(sim_outDir)) && (mkpath(sim_outDir))
			savefig(string(sim_outDir, "network_original_", sim_num,".png"), dpi=150)
			print("\rDone creating plot\n")
			PyPlot.clf()
		end

		if doplot
			(!ispath(sim_outDir)) && (mkpath(sim_outDir))
			plot_eeWeights(new_weights, popmembers, Npop, sim_outDir)
			plot_eiWeights(new_weights, popmembers, ipopmembers, Npop, sim_outDir)
			plot_ieWeights(new_weights, popmembers, ipopmembers, Npop, sim_outDir)
		end
	else
		((times, ns, Ne, Ncells, T, new_weights), popmembers) = simnew(stim, T)
	end

	if savenet
		(!ispath(sim_savepath)) && (mkpath(sim_savepath))
		cd("networks")
		if loadtrained
			fid = h5open(string(sim_name, "_rest.h5"), "w")
		else
			fid = h5open(sim_name,"w")
		end
		g = create_group(fid,"data")
		g["popmembers"] = popmembers
		g["weights"] = new_weights
		close(fid)
		cd("..")
	end

	println("mean excitatory firing rate: ", mean(1000 * ns[1:Ne] / T), " Hz")
	println("mean inhibitory firing rate: ", mean(1000 * ns[(Ne + 1):Ncells] / T), " Hz")
end