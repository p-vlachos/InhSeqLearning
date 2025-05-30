# this file is part of litwin-kumar_doiron_formation_2014
# Copyright (C) 2014 Ashok Litwin-Kumar
# see README for more information

function simnew(stim::Matrix{Float64}, T::Int64; random_seed::Int64=2817)
	#  Generates new weights and populations with unpotentiated synapses, runs simulation

	@info "Setting up weights"
	@unpack Ne, Ni, Ni2, jee0, jei0, jie, jii, jii12, jii2, p, pmembership, Nmaxmembers = InitializationParameters()
	Ncells::Int64 = Ne + Ni
	Random.seed!(random_seed)

	# --- Set up weights ---
	# Weights are set up so that w[i,j] is weight from presynaptic i to postsynaptic j
	weights::Matrix{Float64} = zeros(Ncells, Ncells)
	weights[1:Ne, 1:Ne] .= jee0
	@. weights[1:Ne, (Ne+1):Ncells] = jie
	@. weights[(Ne+1):Ncells, 1:Ne] = jei0
	@. weights[(Ne+1):Ncells, (Ne+1):Ncells] = jii
	@. weights[(Ne+1):(Ncells-Ni2), (Ncells-Ni2+1):Ncells] = jii12
	@. weights[(Ncells-Ni2+1):Ncells, (Ncells-Ni2+1):Ncells] = jii2

	# New manipulation of weights
	@. weights[(Ne+1):(Ncells-Ni2), (Ncells-Ni2+1):Ncells] = 0.		# I₁ → I₂
	@. weights[(Ncells-Ni2+1):Ncells, (Ne+1):(Ncells-Ni2)] = 0.		# I₂ → I₁

	@. weights[1:Ne, (Ncells-Ni2+1):Ncells] = 1.		# E → I₂
	
	# @. weights[(Ncells-Ni2+1):Ncells, (Ncells-Ni2+1):Ncells] = 0.		# I₂ → I₂


	weights[rand(Ncells, Ncells) .> p] .= 0.

	# for i = 1:Ne
	# 	for j = 1:Ne
	# 		(rand() > p) && (weights[i, j] = 0.)	# E → E
	# 	end
	# 	for j = (Ne+1):(Ncells-Ni2)
	# 		(rand() > p) && (weights[i, j] = 0.)	# E → I₁
	# 	end
	# 	for j = (Ncells-Ni2+1):Ncells
	# 		(rand() > p) && (weights[i, j] = 0.)	# E → I₂
	# 	end
	# end

	# for i = (Ne+1):(Ncells-Ni2)
	# 	for j = 1:Ne
	# 		(rand() > p) && (weights[i, j] = 0.)	# I₁ → E
	# 	end
	# 	for j = (Ne+1):(Ncells-Ni2)
	# 		(rand() > p) && (weights[i, j] = 0.)	# I₁ → I₁
	# 	end
	# 	for j = (Ncells-Ni2+1):Ncells
	# 		(rand() > .5) && (weights[i, j] = 0.)	# I₁ → I₂
	# 	end
	# end
	
	# for i = (Ncells-Ni2+1):Ncells
	# 	for j = 1:Ne
	# 		(rand() > .6) && (weights[i, j] = 0.)	# I₂ → E
	# 	end
	# 	for j = (Ne+1):(Ncells-Ni2)
	# 		(rand() > p) && (weights[i, j] = 0.)	# I₂ → I₁
	# 	end
	# 	for j = (Ncells-Ni2+1):Ncells
	# 		(rand() > p) && (weights[i, j] = 0.)	# I₂ → I₂
	# 	end
	# end
	
	weights[diagind(weights)] .= 0.

	# # --- Populations ---
	# # Non-overlapping assemblies of Nmaxmembers neurons each
	# Npop::Int64 = maximum(Int, stim[1, :])	# Number of assemblies
	# popmembers::Matrix{Int64} = reshape(shuffle(1:(Npop*Nmaxmembers)), (Nmaxmembers, Npop))	# Contains indexes of neurons for each population

	# --- Populations ---
	Npop::Int64 = maximum(Int, stim[1, :])	# Number of assemblies
	popmembers::Matrix{Int64} = zeros(Int, Nmaxmembers, Npop)	# Contains indexes of neurons for each population
	@simd for pp = 1:Npop
		members::Vector{Int64} = findall(rand(Ne) .< pmembership)
		popmembers[1:length(members), pp] .= members
	end

	return sim(stim, weights, popmembers, T, random_seed=random_seed)
end


function sim(stim::Matrix{Float64}, weights::Matrix{Float64}, popmembers::Matrix{Int64}, T::Int64; random_seed::Int64=2817)
	#  Runs the simulation given weight matrix and populations

	@info "Setting up parameters..."
	@unpack Ne, Ni, Ni2, Nmaxmembers = InitializationParameters()
	@unpack taue, taui, vleake, vleaki, deltathe, C, erev, irev = NeuronalParameters()
	@unpack vth0, ath, tauth, vpeak, vre, taurefrac, aw_adapt, bw_adapt, tauw_adapt = NeuronalParameters()
	@unpack tauerise, tauedecay, tauirise, tauidecay, taui2rise, taui2decay, rex, rix, ri2x, jex, jix = SynapticParameters()
	@unpack jeemin, jeemax, jeimin, jeimax, jiemin , jiemax, jei2min, jei2max = SynapticParameters()
	@unpack altd, altp, thetaltd, thetaltp, tauu, tauv, taux = PlasticityParameters()
	@unpack tauy, eta, r0, tau_i_r, tau_i_d, ilamda, tau_ie, eta_ie, Adep_ie = PlasticityParameters()
	@unpack dt, dtnormalize, stdpdelay, Nspikes = SimulationParameters()
	Random.seed!(random_seed)

	# _________________________--- Simulation ---______________________________
	Ncells::Int64 = Ne + Ni									# Total number of neurons
	times::Matrix{Float64} = zeros(Ncells, Nspikes)			# Times of neuron spikes
	ns::Vector{Int64} = zeros(Int, Ncells)					# Number of spikes per neuron
	forwardInputsE::Vector{Float64} = zeros(Ncells) 		# Summed weight of incoming E spikes
	forwardInputsI::Vector{Float64} = zeros(Ncells)			# Summed weight of incoming I₁ spikes
	forwardInputsI2::Vector{Float64} = zeros(Ncells)		# Summed weight of incoming I₂ spikes
	forwardInputsEPrev::Vector{Float64} = zeros(Ncells) 	# As above, for previous timestep
	forwardInputsIPrev::Vector{Float64} = zeros(Ncells)		# As above, for I₁ population
	forwardInputsI2Prev::Vector{Float64} = zeros(Ncells)	# As above, for I₂ population
	v::Vector{Float64} = zeros(Ncells) 						# Membrane voltage
	nextx::Vector{Float64} = zeros(Ncells) 					# Time of next external excitatory input
	sumwee0::Vector{Float64} = zeros(Ne) 					# Initial summed E weight, for normalization
	Nee::Vector{Int64} = zeros(Int, Ne) 					# Number of E->E inputs, for normalization
	
	rx::Vector{Float64} = zeros(Ncells) 							# Rate of external input
	vth::Vector{Float64} = vth0 * ones(Ncells) 						# Adaptive threshold
	wadapt::Vector{Float64} = aw_adapt * (vre - vleake) * ones(Ne) 	# Adaptation current
	lastSpike::Vector{Float64} = -100. * ones(Ncells)				# Last time the neuron spiked
	trace_istdp::Vector{Float64} = zeros(Ncells) 					# Low-pass filtered spike train for iSTDP
	u_vstdp::Vector{Float64} = vre * zeros(Ne)						# Low-pass filtered membrane voltage for LTD
	v_vstdp::Vector{Float64} = vre * zeros(Ne)						# Low-pass filtered membrane voltage for LTP
	x_vstdp::Vector{Float64} = zeros(Ne)							# Low-pass filtered spike train for LTP
	trace_eistdp::Vector{Float64} = zeros(Ncells)					# eiSTDP trace
	trace_expDecay::Vector{Float64} = zeros(Ncells)					# iSTDP₂ trace
	Nsteps::Int64 = round(Int, T/dt)
	inormalize::Int64 = round(Int, dtnormalize/dt)

	# Auxiliary variables 
	xerise::Vector{Float64} = zeros(Ncells)				# E currents - rise exponent
	xedecay::Vector{Float64} = zeros(Ncells)			# E currents - decay exponent
	xirise::Vector{Float64} = zeros(Ncells)				# I₁ currents - rise exponent
	xidecay::Vector{Float64} = zeros(Ncells)			# I₁ currents - decay exponent
	xi2rise::Vector{Float64} = zeros(Ncells)			# I₂ currents - rise exponent
	xi2decay::Vector{Float64} = zeros(Ncells)			# I₂ currents - decay exponent
	trace_expDecay_r::Vector{Float64} = zeros(Ncells)	# iSTDP₂ - rise exponent
	trace_expDecay_d::Vector{Float64} = zeros(Ncells)	# iSTDP₂ - decay exponent
	t::Float64 = 0.
	tprev::Float64 = 0.
	ipop::Int64 = 0
	sumwee::Float64 = 0.
	ge::Float64 = 0.
	gi::Float64 = 0.
	gi1::Float64 = 0.
	gi2::Float64 = 0.
	spiked::BitVector = falses(Ncells)
	stim_size::Int64 = size(stim)[2]

	# --- Dictionaries for "sparse" looping ---
	# vSTDP and normalization
	nzRowsByColEE = Dict(m => findall(weights[1:Ne, m].!=0) for m = 1:Ne)  # E-to-E, neuron m gets input from nzRowsByColEE[m]
	nzColsByRowEE = Dict(m => findall(weights[m, 1:Ne].!=0) for m = 1:Ne)  # E-to-E, neuron m sends input to nzColsByRowEE[m]

	# iSTDP₁
	nzRowsByColI1E = Dict(m => Ne .+ findall(weights[(Ne+1):(Ncells-Ni2), m].!=0) for m = 1:Ne) # I1-to-E, E fired, input from I1, neuron m gets input from nzRowsByColI1E[m]
	nzColsByRowI1E = Dict(m => findall(weights[m, 1:Ne].!=0) for m = Ne+1:Ncells-Ni2)     		# For I1 neurons lists all excitatory postsynaptic neurons

	# iSTDP₂
	nzRowsByColI2E = Dict(m => (Ne + Ni - Ni2) .+ findall(weights[(Ncells-Ni2+1):(Ncells), m].!=0) for m = 1:Ne) # I2-to-E, E fired, input from I2, neuron m gets input from nzRowsByColI2E[m]
	nzColsByRowI2E = Dict(m => findall(weights[m, 1:Ne].!=0) for m = (Ncells-Ni2)+1:Ncells) # I2-to-E, I2 fired, lists all E presynaptic
	
	# eiSTDP
	nzColsByRowEI = Dict(m =>  (Ncells - Ni2) .+ findall(weights[m, (Ncells-Ni2+1):(Ncells)].!=0) for m = 1:Ne) # E fired, output to I2
	nzRowsByColEI = Dict(m => findall(weights[1:Ne, m].!=0) for m = (Ncells-Ni2+1):Ncells) # I2 fired, get E presynaptic

	expdist = Exponential()

	# Single cell initialization
	for cc = 1:Ncells
		v[cc] = vre + (vth0 - vre) * rand()
		if cc <= Ne
			rx[cc] = rex
			# Calculate normalization parameters
			@simd for dd = 1:Ne
				sumwee0[cc] += weights[dd, cc]
				(weights[dd, cc] > 0.) && (Nee[cc] += 1)
			end
		elseif cc <= (Ncells-Ni2)
			rx[cc] = rix
		else
			rx[cc] = ri2x
		end
		nextx[cc] = rand(expdist) / rx[cc]
	end

	# --- TRACKERS ---
	tracker_indx::Int64 = 1
	Npop::Int64 = maximum(Int, stim[1, :])	# Number of assemblies
	tracker = Tracker(T=T, Ni=Ni, Ni2=Ni2, Npop=Npop, tracker_dt=round(Int, T/1_000))
	@unpack tracker_dt , weightsEE, weightsEI, weightsIE = tracker
	tracker_dt /= dt

	@info "Starting simulation"
	# --- Begin main simulation loop ---
	iterSteps = ProgressBar(1:Nsteps)
	@inbounds @fastmath for tt in iterSteps
		
		t = dt * tt
		tprev = dt * (tt - 1)

		@. forwardInputsE = 0.
		@. forwardInputsI = 0.
		@. forwardInputsI2 = 0.

		# Check if we have entered or exited a stimulation period
		if t < (stim[3, end] + 1)
			for ss = 1:stim_size
				if (tprev < stim[2, ss]) && (t >= stim[2, ss])	# Just entered stimulation period
					ipop = round(Int, stim[1, ss])
					@simd for ii = 1:Nmaxmembers
						(popmembers[ii, ipop] != 0) && (rx[popmembers[ii, ipop]] += stim[4, ss])	# Add external drive
					end
				end
				if (tprev < stim[3, ss]) && (t >= stim[3, ss]) 	# Just exited stimulation period
					ipop = round(Int, stim[1, ss])
					@simd for ii = 1:Nmaxmembers
						(popmembers[ii, ipop] != 0) && (rx[popmembers[ii, ipop]] -= stim[4, ss])	# Subtract external drive
					end
				end
			end
		end	# End loop over stimuli

		# Synaptic normalization
		if mod(tt, inormalize) == 0
			for cc = 1:Ne
				sumwee = 0.
				@simd for dd in nzRowsByColEE[cc]
					sumwee += weights[dd, cc]
				end
				@simd for dd in nzRowsByColEE[cc]
					weights[dd, cc] = clamp((weights[dd, cc] / sumwee) * sumwee0[cc], jeemin, jeemax)
				end
			end
		end  # End normalization

		# --- Update single cells ---
		for cc = 1:Ncells
			trace_istdp[cc] -= dt * trace_istdp[cc] / tauy
			trace_eistdp[cc] -= dt * trace_eistdp[cc] / tau_ie
			trace_expDecay_r[cc] -= dt * trace_expDecay_r[cc] / tau_i_r
			trace_expDecay_d[cc] -= dt * trace_expDecay_d[cc] / tau_i_d
			trace_expDecay[cc] = (trace_expDecay_d[cc] - trace_expDecay_r[cc]) / (tau_i_d - tau_i_r)

			while(t > nextx[cc]) 	# External input
				nextx[cc] += rand(expdist) / rx[cc]
				cc <= Ne ? forwardInputsEPrev[cc] += jex : forwardInputsEPrev[cc] += jix
			end

			xerise[cc] += -dt * xerise[cc] / tauerise + forwardInputsEPrev[cc]
			xedecay[cc] += -dt * xedecay[cc] / tauedecay + forwardInputsEPrev[cc]
			xirise[cc] += -dt * xirise[cc] / tauirise + forwardInputsIPrev[cc]
			xidecay[cc] += -dt * xidecay[cc] / tauidecay + forwardInputsIPrev[cc]
			xi2rise[cc] += -dt * xi2rise[cc] / taui2rise + forwardInputsI2Prev[cc]
			xi2decay[cc] += -dt * xi2decay[cc] / taui2decay + forwardInputsI2Prev[cc]


			if cc <= Ne
				vth[cc] += dt * (vth0 - vth[cc]) / tauth
				wadapt[cc] += dt * (aw_adapt * (v[cc] - vleake) - wadapt[cc]) / tauw_adapt
				u_vstdp[cc] += dt * (v[cc] - u_vstdp[cc]) / tauu
				v_vstdp[cc] += dt * (v[cc] - v_vstdp[cc]) / tauv
				x_vstdp[cc] -= dt * x_vstdp[cc] / taux
			end

			if t > (lastSpike[cc] + taurefrac)  # Not in refractory period

				# Update membrane voltage
				ge = (xedecay[cc] - xerise[cc]) / (tauedecay - tauerise)
				gi1 = (xidecay[cc] - xirise[cc]) / (tauidecay - tauirise)
				gi2 = (xi2decay[cc] - xi2rise[cc]) / (taui2decay - taui2rise)
				gi = gi1 + gi2

				if cc <= Ne  # Excitatory neuron, has adaptation
					v[cc] += dt * ((vleake - v[cc] + deltathe * exp((v[cc] - vth[cc]) / deltathe)) / taue + ge * (erev - v[cc]) / C + gi * (irev - v[cc]) / C - wadapt[cc] / C)
					if v[cc] > vpeak
						spiked[cc] = true
						wadapt[cc] += bw_adapt
					end
				else
					v[cc] += dt * ((vleaki - v[cc]) / taui + ge * (erev - v[cc]) / C + gi * (irev - v[cc]) / C)
					(v[cc] > vth0) && (spiked[cc] = true)
				end

				if spiked[cc]	# Spike occurred
					v[cc] = vre
					lastSpike[cc] = t
					trace_istdp[cc] += 1.
					trace_eistdp[cc] += 1.
					trace_expDecay_r[cc] += 1.
					trace_expDecay_d[cc] += 1.

					if ns[cc] < Nspikes
						ns[cc] += 1
						times[cc, ns[cc]] = t
					end

					if cc <= Ne
						x_vstdp[cc] += 1. / taux
						vth[cc] = vth0 + ath
					end

					# Loop over synaptic projections
					@simd for dd = 1:Ncells
						if cc <= Ne
							forwardInputsE[dd] += weights[cc, dd]
						elseif cc <= (Ncells-Ni2)
							forwardInputsI[dd] += weights[cc, dd]
						else
							forwardInputsI2[dd] += weights[cc, dd]
						end
					end
				end  # End if(spiked)
			end  # End, if(not refractory)

			# ________________________ Plasticity Rules ________________________
			# _____ iSTDP₁ _____
			if spiked[cc] && (t > stdpdelay)
				if cc <= Ne
					# Excitatory neuron fired, modify inputs from 1st i-population
					@simd for dd in nzRowsByColI1E[cc]
						weights[dd, cc] = min(weights[dd, cc] + eta * trace_istdp[dd], jeimax)
					end
				elseif cc <= (Ncells - Ni2)
					# 1st i-population neuron fired, modify outputs to excitatory neurons
					@simd for dd in nzColsByRowI1E[cc]
						weights[cc, dd] = clamp(weights[cc, dd] + eta * (trace_istdp[dd] - 2. * r0 * tauy), jeimin, jeimax)
					end
				end
			end	# End, iSTDP₁

			# _____ iSTDP₂ _____
			if spiked[cc] && (t > stdpdelay) && (t < stim[3, end])
				if cc <= Ne
					# Excitatory neuron fired, modify inputs from 2nd i-population
					@simd for dd in nzRowsByColI2E[cc]
						weights[dd, cc] = max(weights[dd, cc] - ilamda * trace_expDecay[dd], jei2min)
					end
				elseif cc > (Ncells - Ni2)
					# 2nd i-population neuron fired, modify outputs to excitatory neurons
					@simd for dd in nzColsByRowI2E[cc]
						weights[cc, dd] = min(weights[cc, dd] + ilamda * trace_expDecay[dd], jei2max)
					end
				end
			end # End, iSTDP₂

			#  _____ eiSTDP _____
			if spiked[cc] && (t > stdpdelay) && (t < stim[3, end])
				if cc <= Ne
					# Excitatory neuron fired, modify outputs to 2nd i-population
					@simd for dd in nzColsByRowEI[cc]
						weights[cc, dd] = min(weights[cc, dd] + eta_ie * trace_eistdp[dd], jiemax)
					end
				elseif cc > (Ncells-Ni2)
					# 2nd i-population neuron fired, modify inputs from E neurons
					@simd for dd in nzRowsByColEI[cc]
						weights[dd, cc] = clamp(weights[dd, cc] + eta_ie * (trace_eistdp[dd] - Adep_ie), jiemin, jiemax)
					end
				end
			end # End, eiSTDP

			# ___________ vSTDP ___________
			# vSTDP, LTD component
			if spiked[cc] && (t > stdpdelay) && (cc <= Ne)
				@simd for dd in nzColsByRowEE[cc]
					if u_vstdp[dd] > thetaltd
						weights[cc, dd] = max(weights[cc, dd] - altd * (u_vstdp[dd] - thetaltd), jeemin)
					end
				end
			end  # End, LTD

			# vSTDP, LTP component
			if (t > stdpdelay) && (cc <= Ne) && (v[cc] > thetaltp) && (v_vstdp[cc] > thetaltd)
				@simd for dd in nzRowsByColEE[cc]
					weights[dd, cc] = min(weights[dd, cc] + dt * altp * x_vstdp[dd] * (v[cc] - thetaltp) * (v_vstdp[cc] - thetaltd), jeemax)
				end
			end # End LTP
		end  # End, loop over cells
		@. forwardInputsEPrev = forwardInputsE
		@. forwardInputsIPrev = forwardInputsI
		@. forwardInputsI2Prev = forwardInputsI2
		@. spiked = false

		# --- TRACKERS ---
		if mod(tt-1, tracker_dt) == 0
			for ipop = 1:Npop
				pmembers = filter(i->i>0, popmembers[:, ipop])
				for iipop = 1:Npop					
					ppmembers = filter(i->i>0, popmembers[:, iipop])
					weightsEE[iipop, ipop, tracker_indx] = sum(weights[ppmembers, pmembers]) / count(i->i>0, weights[ppmembers, pmembers])
				end
				for cc = 1:Ni
					weightsIE[cc, ipop, tracker_indx] = sum(weights[cc+Ne, pmembers]) / count(i->i>0, weights[cc+Ne, pmembers])
					if cc > (Ni - Ni2)
						weightsEI[(cc-(Ni-Ni2)), ipop, tracker_indx] = sum(weights[pmembers, (cc+Ne)]) / count(i->i>0, weights[pmembers, (cc+Ne)])
					end
				end
			end
			tracker_indx += 1
		end
	end # End, loop over time

	print("\r")
	times = times[:, 1:maximum(ns)]

	return times, weights, popmembers, weightsEE, weightsIE, weightsEI
end
