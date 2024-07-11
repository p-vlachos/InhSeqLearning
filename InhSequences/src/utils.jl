function make_stim(Nstim, Npres; stim_rate=8)
	# Create stimuli following protocol from LKD paper
	Ntotal = round(Int, Nstim * Npres)
	stim = zeros(Ntotal, 4)
	stim_delay = 10000  # Should be greater than or equal to 'stdpdelay' (ms)
	active = 1000 		# Stimulation period for each assembly (ms)
	inactive = 3000		# Gap between each assembly stimulation (ms)
    stim[1, :] = [1 , stim_delay, stim_delay+active, stim_rate]

	for i = 2:Ntotal
		stim[i, 1] = mod(i - 1, Nstim) + 1
		stim[i, 2] = stim_delay + ((i - 1) * (active + inactive))
		stim[i, 3] = stim[i, 2] + active
		stim[i, 4] = stim_rate
	end
	stim[:, 1] = stim[shuffle(1:Ntotal), 1]
	return stim
end

function make_seq(Nseq, Npres; seq_len=3, stim_rate=8, overlap=0, randomize=true)
	# Create sequences of stimuli
	Nstim = convert(Int, Nseq*seq_len)
	Ntotal = round(Int, Nstim * Npres)
	stim = zeros(Ntotal, 4)
	stim_delay = 10000  # Should be greater than or equal to 'stdpdelay' (ms)
	active = 1000 		# Stimulation period for each assembly (ms)
	inactive = 3000		# Gap between each assembly stimulation (ms)
	stim[1, :] = [1.0, stim_delay, stim_delay + active, stim_rate]

	for i = 2:Ntotal
		stim[i, 1] = mod(i - 1, Nstim) + 1
		if Nseq == 1
			stim[i, 2] = stim[i-1, 3] - overlap
		elseif mod(i, seq_len) == 1
			stim[i, 2] = (stim[i-1, 3] + inactive)
		else
			stim[i, 2] = stim[i-1, 3] - overlap
		end
		stim[i, 3] = stim[i, 2] + active
		stim[i, 4] = stim_rate
	end

	if randomize
		ind = shuffle(filter(i->mod(i, seq_len)==1, stim[:, 1]))
		for i = 1:Ntotal
			(mod(stim[i, 1], seq_len) == 1) ? (stim[i, 1] = pop!(ind)) : (stim[i, 1] = stim[i-1, 1] + 1)
		end
	end
	return stim
end

function make_stim_seq(Nstim, Npres; stim_rate=8, seq_len=4, overlap=0, randomize=true)
	# Create stimuli following protocol from LKD paper
	Ntotal = round(Int, Nstim * Npres)
	Ntotal_wSeq = round(Int, Ntotal * 2)
	stim = zeros(Ntotal_wSeq, 4)
	stim_delay = 10000  # Should be equal to 'stdpdelay'
	active = 1000 		# Presentation period for each repetition (ms)
	inactive = 3000		# Gap between stimulation (ms)
    stim[1, :] = [1.0, stim_delay, stim_delay+active, stim_rate]

	for i = 2:Ntotal
		stim[i, 1] = mod(i - 1, Nstim) + 1
		stim[i, 2] = stim_delay + ((i - 1) * (active + inactive))
		stim[i, 3] = stim[i, 2] + active
		stim[i, 4] = stim_rate
	end
	stim[1:Ntotal, 1] = stim[shuffle(1:Ntotal), 1]

	# Add sequences of stimuli
	seq_delay = (stim[Ntotal, 3] + 50000)  # Resting period post non-sequential training
	active = 1000 		# Presentation period for each repetition (ms)
	inactive = 3000		# Gap between stimulation (ms)
	stim[Ntotal+1, :] = [1.0, seq_delay, seq_delay + active, stim_rate]

	for i = (Ntotal + 2):Ntotal_wSeq
		stim[i, 1] = mod(i - 1, Nstim) + 1
		if mod(i, seq_len) == 1
			stim[i, 2] = (stim[i-1, 3] + inactive)
		else
			stim[i, 2] = stim[i-1, 3] - overlap
		end
		stim[i, 3] = stim[i, 2] + active
		stim[i, 4] = stim_rate
	end

	if randomize
		ind = shuffle(filter(i->mod(i, seq_len)==1, stim[:, 1]))
		for i = Ntotal:Ntotal_wSeq
			if mod(stim[i, 1], seq_len) == 1
				stim[i, 1] = pop!(ind)
			else
				stim[i, 1] = stim[i-1, 1] + 1
			end
		end
	end
	return stim
end

function makeSeqStim(T::Int64; stim_rate::Float64=8., seq_len::Int64=4, seq_num::Int64=5, randomize=true)
	# Create sequences of stimuli
	stim_delay::Float64 = 1_000 	# Should be greater than or equal to 'stdpdelay' (ms)
	stim_duration::Float64 = 30 	# Stimulation period for each assembly (ms)
	stim_interval::Float64 = 200	# Interval between each assembly stimulation (ms)

	Nstim::Int64 = convert(Int, seq_len*seq_num)
	Ntotal::Int64 = round(Int, (T-stim_delay)/stim_interval)
	stim::Matrix{Float64} = zeros(Ntotal, 4)

	# Create the stimulus
	stim[1, :] = [1.0, stim_delay, stim_delay+stim_duration, stim_rate]
	for i = 2:Ntotal
		stim[i, 1] = mod(stim[i-1, 1] + seq_len, Nstim)
		stim[i, 2] = stim[i-1, 2] + stim_interval
		stim[i, 3] = stim[i, 2] + stim_duration
		stim[i, 4] = stim_rate
	end

	# Randomize the order of stimulation
	if randomize
		stim = stim[randperm(Ntotal), :]
	end
	return stim
end

function findI2populations(weights::Matrix{Float64}, Npop::Int64, popmembers; iipop_len=25)
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