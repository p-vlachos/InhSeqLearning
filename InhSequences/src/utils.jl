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

function makeStimSeq(Nseq::Int64, Npres::Int64; seq_len::Int64=4, stim_rate::Float64=8., randomize::Bool=true)
	# Create sequences of stimuli
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
			stim[2, ipop] = stim[3, ipop-1]
		elseif mod(ipop, seq_len) == 1
			stim[2, ipop] = (stim[3, ipop-1] + inactive)
		else
			stim[2, ipop] = stim[3, ipop-1]
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
	stim_duration::Float64 = 100. 	# Stimulation period for each assembly (ms)
	stim_interval::Float64 = 400.	# Interval between each assembly stimulation (ms)

	Nstim::Int64 = seq_len * seq_num
	Ntotal::Int64 = round(Int, (T - stim_delay) / stim_interval)
	stim::Matrix{Float64} = zeros(4, Ntotal)

	# Create the stimulus
	stim[:, 1] = [Npop, 0., 0., 0.]		# This is to ensure that the maximum number of assemblies is known by the sim
	stim[:, 2] = [1., stim_delay, stim_delay+stim_duration, stim_rate]
	for ipop = 3:Ntotal
		stim[1, ipop] = mod(stim[1, ipop-1] + seq_len, Nstim)
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

function findI2populations(weights::Matrix{Float64}, popmembers::Matrix{Int64}; iipop_len::Int64=27, Ni2::Int64=250)
	#	Finds the most highly connected assemblies from E to 2nd i-population (fixed length)
	Npop::Int64 = size(popmembers)[2]
	Ncells::Int64 = size(weights)[1]
	ipopmembers::Matrix{Int64} = zeros(iipop_len, Npop)
	for ipop = 1:Npop
	    members = filter(i->(i>0), popmembers[:, ipop])    			# Get members of excitatory assembly
	    ie_weights = vec(sum(weights[members, (Ncells-Ni2+1):Ncells], dims=1))  # Get sum of all weights projected from each E-assemble to each 2nd ipopulation neuron
		x = sortperm(ie_weights)[(end-iipop_len+1):end]             # Get a permutation for the shorted summed weights
	    ind = 1
	    for ii in x
            ipopmembers[ind, ipop] = ii # Find and store the iipop_len neurons with the lowest summed values
            ind += 1
	    end
	end
	return ipopmembers .+ (Ncells - Ni2)
end

Θ(x::Float64) = x > 0. ? x : 0.
function alpha_fun(t; t0::Float64, tau::Float64=100.)
    (abs(t - t0) / tau > 10) && (return 0.)
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
	if gaussian
		return rate .* 1000 # Convert to Hz
	else
		return rate
	end
end

function binRates(spikeTimes::Matrix{Float64}; interval::AbstractVector, dt::Float64=.125, window::Int64=100)
	Ncells::Int64 = size(spikeTimes)[1]
	Nsteps::Int64 = round(Int, length(interval) / dt)
	rates::Matrix{Int64} = zeros(Ncells, Nsteps)
	# Compute instantaneous population rates
	for cc = 1:Ncells
		for tt in filter(i->(interval[1]<i<interval[end]), spikeTimes[cc, :])
			((tt - interval[1]) < dt) && (continue)
			rates[cc, round(Int, (tt - interval[1]) / dt)] += 1
		end
	end
	# Compute the average over each millisecond
	binSize::Int64 = round(Int, 1/dt)			# By default equal to 1 ms
	avg_rates::Matrix{Float64} = zeros(Ncells, round(Int, Nsteps/binSize))
	for (ind, tt) in enumerate(collect(binSize:binSize:Nsteps))
		avg_rates[:, ind] .= vec(mean(rates[:, tt-binSize+1:tt], dims=2))
	end
	avg_rates *= 1000	# Convert to Hz
	# Compute the smooth rates for each population
	smooth_rates::Matrix{Float64} = zeros(Ncells, round(Int, Nsteps/binSize))
	for cc = 1:Ncells
		smooth_rates[cc, :] = movavg(avg_rates[cc, :], window).x
	end
	return smooth_rates
end

function getPopulationRates(spikeTimes::Matrix{Float64}, popmembers::Matrix{Int64}; interval::AbstractVector, sigma::Float64=5., gaussian=true)
	Npop::Int64 = size(popmembers)[2]
	rates::Matrix{Float64} = zeros(length(interval), Npop)
	for ipop = 1:Npop
		rates[:, ipop] = mean(convolveSpikes(spikeTimes[filter(i->(i>0), popmembers[:, ipop]), :], interval=interval, sigma=sigma, gaussian=gaussian), dims=1)
	end
	return rates
end

function getPopulationBinRates(spikeTimes::Matrix{Float64}, popmembers::Matrix{Int64}; interval::AbstractVector, dt::Float64=.125, window::Int64=20)
	Npop::Int64 = size(popmembers)[2]
	Nsteps::Int64 = round(Int, length(interval) / dt)
	rates::Matrix{Int64} = zeros(Nsteps, Npop)
	# Compute instantaneous population rates
	for ipop = 1:Npop
		for mem in popmembers[:, ipop]
			for tt in filter(i->(interval[1]<i<interval[end]), spikeTimes[mem, :])
				((tt - interval[1]) < dt) && (continue)
				rates[round(Int, (tt - interval[1]) / dt), ipop] += 1
			end
		end
	end
	# Compute the average over each millisecond
	binSize::Int64 = round(Int, 1/dt)			# By default equal to 1 ms
	avg_rates::Matrix{Float64} = zeros(round(Int, Nsteps/binSize), Npop)
	for (ind, tt) in enumerate(collect(binSize:binSize:Nsteps))
		avg_rates[ind, :] .= vec(mean(rates[tt-binSize+1:tt, :], dims=1))
	end
	avg_rates *= 1000/size(popmembers)[1]		# Normalize over all population size and convert to Hz
	# Compute the smooth rates for each population
	smooth_rates::Matrix{Float64} = zeros(round(Int, Nsteps/binSize), Npop)
	for ipop = 1:Npop
		smooth_rates[:, ipop] = movavg(avg_rates[:, ipop], window).x
	end
	return smooth_rates
end

function sequentialityScore(rates::Matrix{Float64}; seq_length::Int64=4)
	Npop::Int64 = size(rates)[2]
	activations::Vector{Int64} = decodeActivity(rates)
	hit::Int64 = 0
	miss::Int64 = 0
	previous::Int64 = activations[1]
	expectation::Int64 = 0
	last_elements::Vector{Int64} = seq_length:seq_length:Npop
	for current in activations
		if current == previous
			continue
		else
			if current != 0
				if expectation != 0		# Check if expectations are met
					(current == expectation) ? (hit += 1) : (miss += 1)
				end
				# Set expectations
				(current in last_elements) ? (expectation = 0) : (expectation = current + 1)
			end
			previous = current
		end
	end
	return hit / (hit + miss), activations
end

function decodeActivity(rates::Matrix{Float64}; window::Int64=20)
	# T::Int64 = size(rates)[1]
	low_limit::Float64 = 1.		# low limit for accepting activation (Hz)
	# Find the most active population at each timestep
	(act, ind) = findmax(rates, dims=2)
	activations::Vector{Int64} = [act[i] > low_limit ? ind[i][2] : 0 for i in 1:size(rates)[1]]

	# Define a time window to avoid flickering (for weak assembly activations, should not matter for good models - code unchanged from H.)
    # offset = Int((activityWindow-1)/2)
    # cleaned = []
    # for i = offset+1:recordedTime-offset
	#     item = mostActive[i]
    #     item < 1 && continue
    #     c = counter(mostActive[i-offset:i+offset])              # function "counter" from package DataStructures
    #     c[item] > offset && push!(cleaned, mostActive[i])  # if populations is most active for the majority of "window_size", record
    # end
	return activations
end

function findOptimalDecoder(spikeTimes::Matrix{Float64}, popmembers::Matrix{Int64}; interval::AbstractVector, dt::Float64=.125, seq_length::Int64=4)
    smoothingWindows::Vector{Int64} = collect(10:10:100)	# Smoothing windows for computing firing rates
    activityWindows::Vector{Int64} = collect(3:2:19)		# Windows for determining most active population
    opt_params::Vector{Float64} = zeros(2)
	sequence::Vector = zeros(length(interval))
	maxSeq::Float64 = 0.
    for sW in smoothingWindows
        for aW in activityWindows
            rates = getPopulationBinRates(spikeTimes, popmembers, interval=interval, dt=dt, window=sW)
            # Compute sequentiality score
            sequentiality, sequence = sequentialityScore(rates, seq_length=seq_length)
            if sequentiality > maxSeq
                maxSeq = sequentiality
                opt_params .= [sW, aW]
            end
        end
    end
	# Compare with random baseline - NOTE: Needs to be fixed if added
    # randScore = sequentialityScore(shuffle(sequence), seq_length=seq_length)
    return maxSeq, opt_params, sequence#, ((maxSeq - randScore)/(1 - randScore))
end

