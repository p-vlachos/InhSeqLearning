function makeStim(Npop::Int64, Npres::Int64; stim_rate::Float64=8.)
	# Create stimuli following protocol from LKD paper
	Ntotal::Int64 = Npop * Npres
	stim::Matrix{Float64} = zeros(4, Ntotal)
	stim_delay::Float64 = 10_000.	# Should be greater than or equal to 'stdpdelay' (ms)
	active::Float64 = 1_000. 		# Stimulation period for each assembly (ms)
	inactive::Float64 = 3_000.		# Gap between each assembly stimulation (ms)
    stim[:, 1] = [1., stim_delay, stim_delay+active, stim_rate]
	for ipop = 2:Ntotal
		stim[1, ipop] = mod(ipop - 1, Npop) + 1
		stim[2, ipop] = stim_delay + ((ipop - 1) * (active + inactive))
		stim[3, ipop] = stim[2, ipop] + active
		stim[4, ipop] = stim_rate
	end
	stim[1, :] = stim[1, shuffle(1:Ntotal)]
	return stim
end

function makeStimSeq(Nseq::Int64, Npres::Int64; seq_len::Int64=4, stim_rate::Float64=8., overlap::Float64=0., randomize::Bool=true)
	# Create sequences of stimuli
	# NOTE: "overlap" parameter might not work properly (haven't tested it)
	Nstim::Int64 = Nseq * seq_len
	Ntotal::Int64 = Nstim * Npres
	stim::Matrix{Float64} = zeros(4, Ntotal)
	stim_delay::Float64 = 10_000.  	# Should be greater than or equal to 'stdpdelay' (ms)
	active::Float64 = 1_000. 		# Stimulation period for each assembly (ms)
	inactive::Float64 = 3_000.		# Gap between each assembly stimulation (ms)
	stim[:, 1] = [1.0, stim_delay, stim_delay + active, stim_rate]

	for ipop = 2:Ntotal
		stim[1, ipop] = mod(ipop - 1, Nstim) + 1
		if Nseq == 1
			stim[2, ipop] = stim[3, ipop-1] - overlap
		elseif mod(ipop, seq_len) == 1
			stim[2, ipop] = (stim[3, ipop-1] + inactive)
		else
			stim[2, ipop] = stim[3, ipop-1] - overlap
		end
		stim[3, ipop] = stim[2, ipop] + active
		stim[4, ipop] = stim_rate
	end

	if randomize
		ind = shuffle(filter(i->(mod(i, seq_len) == 1), stim[1, :]))
		for ipop = 1:Ntotal
			(mod(stim[1, ipop], seq_len) == 1) ? (stim[1, ipop] = pop!(ind)) : (stim[1, ipop] = stim[1, ipop-1] + 1)
		end
	end
	return stim
end

function makeStimSeq_brief(T::Int64; Npop::Int64=20, stim_rate::Float64=8., seq_len::Int64=4, seq_num::Int64=5, randomize=true)
	# Create sequences of stimuli
	stim_delay::Float64 = 1_000. 	# Should be greater than or equal to 'stdpdelay' (ms)
	stim_duration::Float64 = 30. 	# Stimulation period for each assembly (ms)
	stim_interval::Float64 = 200.	# Interval between each assembly stimulation (ms)

	Nstim::Int64 = seq_len * seq_num
	Ntotal::Int64 = round(Int, (T - stim_delay) / stim_interval)
	stim::Matrix{Float64} = zeros(4, Ntotal)

	# Create the stimulus
	stim[:, 1] = [Npop, 0., 0., 0.]		# This is to ensure that the maximum number of assemblies is known by the sim
	stim[:, 2] = [1., stim_delay, stim_delay+stim_duration, stim_rate]
	for ipop = 3:Ntotal
		stim[1, ipop] = mod(stim[ipop-1, 1] + seq_len, Nstim)
		stim[2, ipop] = stim[2, ipop-1] + stim_interval
		stim[3, ipop] = stim[2, ipop] + stim_duration
		stim[4, ipop] = stim_rate
	end

	# Randomize the order of stimulation
	if randomize
		stim = stim[:, randperm(Ntotal)]
	end
	return stim
end

function findI2populations(weights::Matrix{Float64}, Npop::Int64, popmembers; iipop_len::Int64=27)
	#	Finds the most highly connected assemblies from E to 2nd i-population (fixed length)
	ipopmembers = zeros(iipop_len, Npop)
	for ipop = 1:Npop
	    members = filter(i->(i>0), popmembers[:, ipop])    			# Get members of excitatory assembly
	    ie_weights = vec(sum(weights[members, 4751:5000], dims=1))  # Get sum of all weights projected from each E-assemble to each 2nd ipopulation neuron
		x = sortperm(ie_weights)[(end-iipop_len+1):end]             # Get a permutation for the shorted summed weights
	    ind = 1
	    for ii in x
            ipopmembers[ind, ipop] = ii # Find and store the iipop_len neurons with the lowest summed values
            ind += 1
	    end
	end
	ipopmembers = convert(Array{Int,2}, ipopmembers .+ 4750)
	return ipopmembers
end

Θ(x::Float64) = x > 0. ? x : 0.
function alpha_fun(t; t0::Float64, tau::Float64=100.)
    (abs(t - t0)/ tau > 10) && (return 0.)
    return (t-t0) / tau * exp(1 - (t - t0) / tau) * Θ(1. * (t - t0))
end

function gaussian_fun(interval::AbstractVector; t0::Float64, sigma::Float64=5.)
    rate::Vector{Float64} = zeros(length(interval))
    sum::Float64 = 0.
    for t = 1:length(interval)
        rate[t] = exp(-0.5 * ((t - t0) / sigma)^2)
        sum += rate[t]
    end
    return rate ./= sum  # Normalize
end

function convolveSpikes(spikeTimes::Matrix{Float64}; interval::AbstractVector, sigma::Float64=5., gaussian::Bool=true)
	Ncells::Int64 = size(spikeTimes)[1]
	rate::Matrix = zeros(Ncells, length(interval))
	for cc = 1:Ncells
		for t0 in filter(i->(i>0), spikeTimes[cc, :])
			(gaussian) ? (x = gaussian_fun(interval, t0=t0, sigma=sigma)) : (x = alpha_fun.(interval, t0=t0, tau=sigma))
			x[isnan.(x)] .= 0.
			rate[cc, :] .+= x
		end
	end
	return rate .* 1000 # Convert to Hz
end

function getPopulationRates(spikeTimes::Matrix{Float64}, popmembers::Matrix{Int64}; interval::AbstractVector, sigma::Float64=5., gaussian=true)
	Npop::Int64 = size(popmembers)[2]
	rates::Matrix{Float64} = zeros(length(interval), Npop)
	for ipop = 1:Npop
		rates[:, ipop] = mean(convolveSpikes(spikeTimes[filter(i->(i>0), popmembers[:, ipop]), :], interval=interval, sigma=sigma, gaussian=gaussian), dims=1)
	end
	return rates
end


# ---------- OLD FUNCTIONS -----------
# function times2spikes(times::Matrix{Float64})
# 	Ncells::Int64 = size(times)[1]
# 	spikes::Matrix{Int64} = zeros(Ncells, round(Int, maximum(times)))
# 	for cc = 1:Ncells
# 		for tt in times[1, :]
# 			(tt == 0) && (continue)
# 			t = round(Int, tt)
# 			spikes[cc, t] += 1
# 		end
# 	end
# 	return spikes	
# end