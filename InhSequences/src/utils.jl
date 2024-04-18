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


function findI2populations(weights, Npop, popmembers; iipop_len=25)
	#	Finds the most highly connected assemblies from E to 2nd i-population (fixed length)
	ipopmembers = zeros(Npop, iipop_len)
	for ipop = 1:Npop
	    members = filter(i->(i>0), popmembers[ipop, :])    			# Get members of excitatory assembly
	    ie_weights = vec(sum(weights[members, 4751:5000], dims=1))  # Get sum of all weights projected from each E-assemble to each 2nd ipopulation neuron
		x = sortperm(ie_weights)[(end-iipop_len+1):end]             # Get a permutation for the shorted summed weights
	    ind = 1
	    for ii in x
            ipopmembers[ipop, ind] = ii # Find and store the iipop_len neurons with the lowest summed values
            ind += 1
	    end
	end
	ipopmembers = convert(Array{Int,2}, ipopmembers .+ 4750)
	return ipopmembers
end