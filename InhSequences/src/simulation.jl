# this file is part of litwin-kumar_doiron_formation_2014
# Copyright (C) 2014 Ashok Litwin-Kumar
# see README for more information

function simnew(stim::Matrix{Float64}, T::Int64)
	#  Generates new weights and populations with unpotentiated synapses, runs simulation

	@info "Setting up weights"
	@unpack Ne, Ni, Ni2, jee0, jei0, jie, ji2e, jii, jii12, jii2, p, pmembership, Nmaxmembers = InitializationParameters()
	Ncells::Int64 = Ne + Ni

	# --- Set up weights ---
	# Weights are set up so that w[i,j] is weight from presynaptic i to postsynaptic j
	weights::Matrix{Float64} = zeros(Ncells,Ncells)
	weights[1:Ne, 1:Ne] .= jee0
	@. weights[1:Ne, (Ne+1):Ncells] = jie

	@. weights[1:Ne, (Ncells-Ni2+1):Ncells] = ji2e

	@. weights[(Ne+1):Ncells, 1:Ne] = jei0

	@. weights[(Ne+1):Ncells, (Ne+1):Ncells] = jii
	@. weights[(Ne+1):(Ncells-Ni2), (Ncells-Ni2+1):Ncells] = jii12
	@. weights[(Ncells-Ni2+1):Ncells, (Ncells-Ni2+1):Ncells] = jii2
	weights[rand(Ncells, Ncells) .> p] .= 0.

	weights[diagind(weights)] .= 0.

	# --- Populations ---
	Npop::Int64 = maximum(Int, stim[1, :])	# Number of assemblies
	popmembers::Matrix{Int64} = zeros(Int, Nmaxmembers, Npop)	# Contains indexes of neurons for each population
	@simd for pp = 1:Npop
		members::Vector{Int64} = findall(rand(Ne) .< pmembership)
		popmembers[1:length(members), pp] = members
	end
	
	return sim(stim, weights, popmembers, T)
end


function sim(stim::Matrix{Float64}, weights::Matrix{Float64}, popmembers::Matrix{Int64}, T::Int64)
	#  Runs the simulation given weight matrix and populations

	@info "Setting up parameters..."
	@unpack Ne, Ni, Ni2, Nmaxmembers = InitializationParameters()
	@unpack taue, taui, vleake, vleaki, deltathe, C, erev, irev = NeuronalParameters()
	@unpack vth0, ath, tauth, vpeak, vre, taurefrac, aw_adapt, bw_adapt, tauw_adapt = NeuronalParameters()
	@unpack tauerise, tauedecay, tauirise, tauidecay, rex, rix, jex, jix = SynapticParameters()
	@unpack jeemin, jeemax, jeimin, jeimax, jiemin , jiemax, jei2min, jei2max = SynapticParameters()
	@unpack altd, altp, thetaltd, thetaltp, tauu, tauv, taux = PlasticityParameters()
	@unpack tauy, eta, r0, tau_i, mi, alfa, ilamda, tau_ie, eta_ie, Adep_ie = PlasticityParameters()
	@unpack dt, dtnormalize, stdpdelay, Nspikes = SimulationParameters()

	# _________________________--- Simulation ---______________________________
	Ncells::Int64 = Ne + Ni									# Total number of neurons
	times::Matrix{Float64} = zeros(Ncells, Nspikes)			# Times of neuron spikes
	ns::Vector{Int64} = zeros(Int, Ncells)					# Number of spikes per neuron
	forwardInputsE::Vector{Float64} = zeros(Ncells) 		# Summed weight of incoming E spikes
	forwardInputsI::Vector{Float64} = zeros(Ncells)			# Summed weight of incoming I spikes
	forwardInputsEPrev::Vector{Float64} = zeros(Ncells) 	# As above, for previous timestep
	forwardInputsIPrev::Vector{Float64} = zeros(Ncells)		# As above, for I population
	v::Vector{Float64} = zeros(Ncells) 						# Membrane voltage
	nextx::Vector{Float64} = zeros(Ncells) 					# Time of next external excitatory input
	sumwee0::Vector{Float64} = zeros(Ne) 					# Initial summed E weight, for normalization
	Nee::Vector{Int64} = zeros(Int, Ne) 					# Number of E->E inputs, for normalization

	# E➡I₂ synaptic normalization
	sumwie0::Vector{Float64} = zeros(Ni2) 					# Initial summed I2 weight, for normalization
	Nie::Vector{Int64} = zeros(Int, Ni2) 					# Number of E->I₂ inputs, for normalization
	
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
	xerise::Vector{Float64} = zeros(Ncells)		# for E/I currents (difference of exponentials)
	xedecay::Vector{Float64} = zeros(Ncells)
	xirise::Vector{Float64} = zeros(Ncells)
	xidecay::Vector{Float64} = zeros(Ncells)
	t::Float64 = 0.
	tprev::Float64 = 0.
	ipop::Int64 = 0
	sumwee::Float64 = 0.
	sumwie::Float64 = 0.
	ge::Float64 = 0.
	gi::Float64 = 0.
	spiked::BitVector = falses(Ncells)
	indx::Int64 = 0.
	stim_size::Int64 = size(stim)[2]

	expdist = Exponential()

	# Single cell initialization
	for cc = 1:Ncells
		v[cc] = vre + (vth0 - vre) * rand()
		if cc <= Ne
			rx[cc] = rex
			nextx[cc] = rand(expdist) / rx[cc]
			# Calculate normalization parameters
			@simd for dd = 1:Ne
				sumwee0[cc] += weights[dd, cc]
				(weights[dd, cc] > 0.) && (Nee[cc] += 1)
			end
		else
			rx[cc] = rix
			nextx[cc] = rand(expdist) / rx[cc]
			if cc > Ncells-Ni2
				indx = round(Int, cc-Ncells+Ni2)
				@simd for dd = 1:Ne
					sumwie0[indx] += weights[dd, cc]
					(weights[dd, cc] > 0.) && (Nie[indx] += 1)
				end
			end
		end
	end

	# --- TRACKERS ---
	tracker_indx::Int64 = 1
	Npop::Int64 = 20#maximum(Int, stim[1, :])	# Number of assemblies
	tracker = Tracker(T=T, Ni=Ni, Ni2=Ni2, Npop=Npop, tracker_dt=round(Int, T/10_000))
	@unpack tracker_dt , weightsEE, weightsEI, weightsIE = tracker
	tracker_dt /= dt

	@info "Starting simulation"
	# --- Begin main simulation loop ---
	iterSteps = ProgressBar(1:Nsteps)
	# @inbounds @fastmath 
	for tt in iterSteps
		
		t = dt * tt
		tprev = dt * (tt - 1)

		@. forwardInputsE = 0.
		@. forwardInputsI = 0.

		# Check if we have entered or exited a stimulation period
		for ss = 1:stim_size
			if (tprev < stim[2, ss]) && (t >= stim[2, ss])	# Just entered stimulation period
				ipop = round(Int, stim[1, ss])
				@simd for ii = 1:Nmaxmembers
					# (popmembers[ipop, ii] == 0 || popmembers[ipop, ii] == -1) && break
					# rx[popmembers[ipop, ii]] += stim[ss, 4]	# Add external drive
					!(popmembers[ii, ipop] == 0 || popmembers[ii, ipop] == -1) && (rx[popmembers[ii, ipop]] += stim[4, ss])	# Add external drive
				end
			end
			if (tprev < stim[3, ss]) && (t >= stim[3, ss]) 	# Just exited stimulation period
				ipop = round(Int, stim[1, ss])
				@simd for ii = 1:Nmaxmembers
					# (popmembers[ipop, ii] == 0 || popmembers[ipop, ii] == -1) && break
					# rx[popmembers[ipop, ii]] -= stim[ss, 4]	# Subtract external drive
					!(popmembers[ii, ipop] == 0 || popmembers[ii, ipop] == -1) && (rx[popmembers[ii, ipop]] -= stim[ss, 4])	# Subtract external drive
				end
			end
		end  # End loop over stimuli

		# Synaptic normalization
		if mod(tt, inormalize) == 0
			for cc = 1:Ne
				sumwee = 0.
				@simd for dd = 1:Ne
					sumwee += weights[dd, cc]
				end
				@simd for dd = 1:Ne
					if weights[dd, cc] > 0.
						weights[dd, cc] = (weights[dd, cc] / sumwee) * sumwee0[cc]
						if weights[dd, cc] < jeemin
							weights[dd, cc] = jeemin
						elseif weights[dd, cc] > jeemax
							weights[dd, cc] = jeemax
						end
					end
				end
			end
			for cc = (Ncells-Ni2+1):Ncells
				indx = cc - Ncells + Ni2
				sumwie = 0.
				@simd for dd = 1:Ne
					sumwie += weights[dd, cc]
				end
				@simd for dd = 1:Ne
					if weights[dd, cc] > 0.
						weights[dd, cc] = (weights[dd, cc] / sumwie) * sumwie0[indx]
						if weights[dd, cc] < jiemin
							weights[dd, cc] = jiemin
						elseif weights[dd, cc] > jiemax
							weights[dd, cc] = jiemax
						end
					end
				end
			end
		end  # End normalization

		# --- Update single cells ---
		for cc = 1:Ncells
			trace_istdp[cc] -= dt * trace_istdp[cc] / tauy
			trace_eistdp[cc] -= dt * trace_eistdp[cc] / tau_ie
			trace_expDecay[cc] -= dt * trace_expDecay[cc] / tau_i

			while(t > nextx[cc]) 	# External input
				nextx[cc] += rand(expdist) / rx[cc]
				cc <= Ne ? forwardInputsEPrev[cc] += jex : forwardInputsEPrev[cc] += jix
			end

			xerise[cc] += -dt * xerise[cc] / tauerise + forwardInputsEPrev[cc]
			xedecay[cc] += -dt * xedecay[cc] / tauedecay + forwardInputsEPrev[cc]
			xirise[cc] += -dt * xirise[cc] / tauirise + forwardInputsIPrev[cc]
			xidecay[cc] += -dt * xidecay[cc] / tauidecay + forwardInputsIPrev[cc]

			if cc <= Ne
				vth[cc] += dt * (vth0 - vth[cc]) / tauth;
				wadapt[cc] += dt * (aw_adapt * (v[cc] - vleake) - wadapt[cc]) / tauw_adapt;
				u_vstdp[cc] += dt * (v[cc] - u_vstdp[cc]) / tauu;
				v_vstdp[cc] += dt * (v[cc] - v_vstdp[cc]) / tauv;
				x_vstdp[cc] -= dt * x_vstdp[cc] / taux;
			end

			if t > (lastSpike[cc] + taurefrac)  # Not in refractory period

				# Update membrane voltage
				ge = (xedecay[cc] - xerise[cc]) / (tauedecay - tauerise)
				gi = (xidecay[cc] - xirise[cc]) / (tauidecay - tauirise)

				if cc <= Ne  # Excitatory neuron, has adaptation
					v[cc] += dt * ((vleake - v[cc] + deltathe * exp((v[cc] - vth[cc]) / deltathe)) / taue + ge * (erev - v[cc]) / C + gi * (irev - v[cc]) / C - wadapt[cc] / C)
					if v[cc] > vpeak
						spiked[cc] = true
						wadapt[cc] += bw_adapt
					end
				else
					v[cc] += dt * ((vleaki - v[cc]) / taui + ge * (erev - v[cc]) / C + gi * (irev - v[cc]) / C)
					if v[cc] > vth0
						spiked[cc] = true
					end
				end

				if spiked[cc]	# Spike occurred
					v[cc] = vre
					lastSpike[cc] = t
					trace_istdp[cc] += 1.
					trace_eistdp[cc] += 1.
					trace_expDecay[cc] += 1.

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
						else
							forwardInputsI[dd] += weights[cc, dd]
						end
					end
				end  # End if(spiked)
			end  # End, if(not refractory)

			# ________________________ Plasticity Rules ________________________
			# _____ iSTDP₁ _____
			if spiked[cc] && (t > stdpdelay)
				if cc <= Ne
					# Excitatory neuron fired, modify inputs from 1st i-population
					for dd = (Ne+1):(Ncells-Ni2)
						(weights[dd, cc] == 0.) && (continue)
						weights[dd, cc] += eta * trace_istdp[dd]
						(weights[dd, cc] > jeimax) && (weights[dd, cc] = jeimax)
					end
				elseif cc <= (Ncells - Ni2)
					# 1st i-population neuron fired, modify outputs to excitatory neurons
					for dd = 1:Ne
						(weights[cc, dd] == 0.) && (continue)
						weights[cc, dd] += eta * (trace_istdp[dd] - 2. * r0 * tauy)
						if weights[cc, dd] > jeimax
							weights[cc, dd] = jeimax
						elseif weights[cc, dd] < jeimin
							weights[cc, dd] = jeimin
						end
					end
				end
			end	# End, iSTDP₁

			# _____ iSTDP₂ _____
			if spiked[cc] && (t > stdpdelay) && (t < stim[3, end])
				if cc <= Ne
					# Excitatory neuron fired, modify inputs from 2nd i-population
					for dd = (Ncells - Ni2 + 1):Ncells
						(weights[dd, cc] == 0.) && (continue)
						weights[dd, cc] -= ilamda * ((1. - weights[dd, cc])^mi) * trace_expDecay[dd]
						(weights[dd, cc] < jei2min) && (weights[dd, cc] = jei2min)
					end
				elseif cc > (Ncells - Ni2)
					# 2nd i-population neuron fired, modify outputs to excitatory neurons
					for dd = 1:Ne
						(weights[cc, dd] == 0.) && (continue)
						weights[cc, dd] += ilamda * alfa * (weights[cc, dd]^mi) * trace_expDecay[dd]
						(weights[cc, dd] > jei2max) && (weights[cc, dd] = jei2max)
					end
				end
			end # End, iSTDP₂

			#  _____ eiSTDP _____
			if spiked[cc] && (t > stdpdelay) && (t < stim[3, end])
				if cc <= Ne
					# Excitatory neuron fired, modify outputs to 2nd i-population
					for dd = (Ncells-Ni2+1):Ncells
						(weights[cc, dd] == 0.) && (continue)
						weights[cc, dd] += eta_ie * trace_eistdp[dd]
						(weights[cc, dd] > jiemax) && (weights[cc, dd] = jiemax)
					end
				elseif cc > (Ncells-Ni2)
					# 2nd i-population neuron fired, modify inputs from E neurons
					for dd = 1:Ne
						(weights[dd, cc] == 0.) && (continue)
						weights[dd, cc] += eta_ie * (trace_eistdp[dd] - Adep_ie)
						# weights[dd, cc] -= eta_ie * trace_eistdp[dd]
						if weights[dd, cc] > jiemax
							weights[dd, cc] = jiemax
						elseif weights[dd, cc] < jiemin
							weights[dd, cc] = jiemin
						end
					end
				end
			end # End, eiSTDP

			# ___________ vSTDP ___________
			# vSTDP, LTD component
			if spiked[cc] && (t > stdpdelay) && (cc < Ne)
				for dd = 1:Ne 	# Depress weights from cc to dd
					(weights[cc, dd] == 0.) && (continue)
					if u_vstdp[dd] > thetaltd
						weights[cc, dd] -= altd * (u_vstdp[dd] - thetaltd)
						(weights[cc, dd] < jeemin) && (weights[cc, dd] = jeemin)
					end
				end
			end  # End, LTD

			# vSTDP, LTP component
			if (t > stdpdelay) && (cc < Ne) && (v[cc] > thetaltp) && (v_vstdp[cc] > thetaltd)
				for dd = 1:Ne
					(weights[dd, cc] == 0.) && (continue)
					weights[dd, cc] += dt * altp * x_vstdp[dd] * (v[cc] - thetaltp) * (v_vstdp[cc] - thetaltd)
					(weights[dd, cc] > jeemax) && (weights[dd, cc] = jeemax)
				end
			end # End LTP
		end  # End, loop over cells
		@. forwardInputsEPrev = forwardInputsE
		@. forwardInputsIPrev = forwardInputsI
		@. spiked = false

		# --- TRACKERS ---
		if mod(tt-1, tracker_dt) == 0
			for ipop = 1:Npop
				pmembers = filter(i->i>0, popmembers[:, ipop])
				for iipop = 1:Npop					
					ppmembers = filter(i->i>0, popmembers[:, ipop])
					weightsEE[iipop, ipop, tracker_indx] = sum(weights[ppmembers, pmembers])
				end
				for cc = 1:Ni
					weightsIE[cc, ipop, tracker_indx] = sum(weights[cc+Ne, pmembers])
					if cc > (Ni - Ni2)
						weightsEI[(cc-(Ni-Ni2)), ipop, tracker_indx] = sum(weights[pmembers, (cc+Ne)])
					end
				end
			end
			tracker_indx += 1
		end
	end # End, loop over time

	print("\r")
	times = times[:, 1:maximum(ns)]

	return times, ns, Ne, Ncells, T, weights, weightsEE, weightsIE, weightsEI, popmembers
end
