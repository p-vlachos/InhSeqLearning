# this file is part of litwin-kumar_doiron_formation_2014
# Copyright (C) 2014 Ashok Litwin-Kumar
# see README for more information

function simnew(stim::Matrix{Float64})
	#  Generates new weights and populations with unpotentiated synapses, runs simulation

	@info "Setting up weights"
	@unpack Ne, Ni, Ni2, jee0, jei0, jie, jii, jii12, jii2, p, pmembership, Nmaxmembers = InitParams()
	Ncells::Int64 = Ne + Ni

	# --- Set up weights ---
	# Weights are set up so that w[i,j] is weight from presynaptic i to postsynaptic j
	weights::Matrix{Float64} = zeros(Ncells,Ncells)
	weights[1:Ne, 1:Ne] .= jee0
	weights[1:Ne, (Ne+1):Ncells] .= jie

	weights[1:Ne, (Ncells-Ni2+1):Ncells] .= 1.52 #.*= 1.2

	weights[(Ne+1):Ncells, 1:Ne] .= jei0

	weights[(Ne+1):Ncells, (Ne+1):Ncells] .= jii
	weights[(Ne+1):(Ncells-Ni2), (Ncells-Ni2+1):Ncells] .= jii12
	weights[(Ncells-Ni2+1):Ncells, (Ncells-Ni2+1):Ncells] .= jii2
	weights[rand(Ncells, Ncells) .> p] .= 0.

	weights[diagind(weights)] .= 0.

	# --- Populations ---
	Npop::Int64 = maximum(Int, stim[:, 1])	# Number of assemblies
	popmembers::Matrix{Int64} = zeros(Int, Npop, Nmaxmembers)	# Contains indexes of neurons for each population
	for pp = 1:Npop
		members::Vector{Int64} = findall(rand(Ne) .< pmembership)
		popmembers[pp, 1:length(members)] .= members
	end
	
	return sim(stim, weights, popmembers), popmembers
end


function sim(stim::Matrix{Float64}, weights::Matrix{Float64}, popmembers::Matrix{Int64})
	#  Runs the simulation given weight matrix and populations

	@info "Setting up parameters..."
	@unpack Ne, Ni, Ni2, Nmaxmembers = InitParams()
	@unpack taue, taui, vleake, vleaki, deltathe, C, erev, irev = NeuronalParams()
	@unpack vth0, ath, tauth, vpeak, vre, taurefrac, aw_adapt, bw_adapt, tauw_adapt = NeuronalParams()
	@unpack tauerise, tauedecay, tauirise, tauidecay, rex, rix, jex, jix = SynapticParams()
	@unpack jeemin, jeemax, jeimin, jeimax, jiemin , jiemax, jei2min, jei2max = SynapticParams()
	@unpack altd, altp, thetaltd, thetaltp, tauu, tauv, taux = PlasticityParams()
	@unpack tauy, eta, r0, tau_i, mi, alfa, ilamda, tau_ie, eta_ie, Adep_ie = PlasticityParams()
	@unpack dt, dtnormalize, stdpdelay, Nspikes = SimParams()

	# _________________________--- Simulation ---______________________________
	T::Int64 = 6000#1500000#6000#1500000#707000#800000								# Simulation time (ms)
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

	# # NOTE:eiNorm
	# sumwei20::Vector{Float64} = zeros(Ne) 					# Initial summed E weight, for normalization
	# Nei2::Vector{Int64} = zeros(Int, Ne) 					# Number of I₂->E inputs, for normalization

	# # NOTE:eiNorm - all
	# sumwei0::Vector{Float64} = zeros(Ne) 					# Initial summed E weight, for normalization
	# Nei::Vector{Int64} = zeros(Int, Ne) 					# Number of I->E inputs, for normalization

	# NOTE:i2eNorm
	sumwie0::Vector{Float64} = zeros(Ni2) 					# Initial summed I2 weight, for normalization
	Nie::Vector{Int64} = zeros(Int, Ni2) 					# Number of E->I₂ inputs, for normalization
	
	rx::Vector{Float64} = zeros(Ncells) 							# Rate of external input
	vth::Vector{Float64} = vth0 * ones(Ncells) 						# Adaptive threshold
	wadapt::Vector{Float64} = aw_adapt * (vre - vleake) * ones(Ne) 	# Adaptation current
	lastSpike::Vector{Float64} = -100 * ones(Ncells)				# Last time the neuron spiked
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
	t::Float64 = 0
	tprev::Float64 = 0
	ipop::Int64 = 0
	sumwee::Float64 = 0.
	sumwei::Float64 = 0.	# NOTE:eiNorm
	sumwie::Float64 = 0.	# NOTE:ieNorm
	ge::Float64 = 0
	gi::Float64 = 0
	spiked::BitVector = falses(Ncells)

	expdist = Exponential()

	# # Single cell initialization
	# for cc = 1:Ncells
	# 	v[cc] = vre + (vth0 - vre) * rand()
	# 	if cc <= Ne
	# 		rx[cc] = rex
	# 		nextx[cc] = rand(expdist) / rx[cc]
	# 		# Calculate normalization parameters
	# 		for dd = 1:Ne
	# 			sumwee0[cc] += weights[dd, cc]
	# 			if weights[dd, cc] > 0
	# 				Nee[cc] += 1
	# 			end
	# 		end
	# 	else
	# 		rx[cc] = rix
	# 		nextx[cc] = rand(expdist) / rx[cc]
	# 	end
	# end

	# # Single cell initialization (with ei)
	# for cc = 1:Ncells
	# 	v[cc] = vre + (vth0 - vre) * rand()
	# 	if cc <= Ne
	# 		rx[cc] = rex
	# 		nextx[cc] = rand(expdist) / rx[cc]
	# 		# Calculate normalization parameters
	# 		for dd = 1:Ne
	# 			sumwee0[cc] += weights[dd, cc]
	# 			if weights[dd, cc] > 0
	# 				Nee[cc] += 1
	# 			end
	# 		end
	# 		for dd = (Ncells-Ni2+1):Ncells
	# 			sumwei20[cc] += weights[dd, cc]
	# 			if weights[dd, cc] > 0
	# 				Nei2[cc] += 1
	# 			end
	# 		end
	# 	else
	# 		rx[cc] = rix
	# 		nextx[cc] = rand(expdist) / rx[cc]
	# 	end
	# end

	# # Single cell initialization (with ei All)
	# for cc = 1:Ncells
	# 	v[cc] = vre + (vth0 - vre) * rand()
	# 	if cc <= Ne
	# 		rx[cc] = rex
	# 		nextx[cc] = rand(expdist) / rx[cc]
	# 		# Calculate normalization parameters
	# 		for dd = 1:Ne
	# 			sumwee0[cc] += weights[dd, cc]
	# 			if weights[dd, cc] > 0
	# 				Nee[cc] += 1
	# 			end
	# 		end
	# 		for dd = (Ne+1):Ncells
	# 			sumwei0[cc] += weights[dd, cc]
	# 			if weights[dd, cc] > 0
	# 				Nei[cc] += 1
	# 			end
	# 		end
	# 	else
	# 		rx[cc] = rix
	# 		nextx[cc] = rand(expdist) / rx[cc]
	# 	end
	# end

	# Single cell initialization ieNorm	# NOTE::: THIS WAS THE CORRECT
	for cc = 1:Ncells
		v[cc] = vre + (vth0 - vre) * rand()
		if cc <= Ne
			rx[cc] = rex
			nextx[cc] = rand(expdist) / rx[cc]
			# Calculate normalization parameters
			for dd = 1:Ne
				sumwee0[cc] += weights[dd, cc]
				if weights[dd, cc] > 0
					Nee[cc] += 1
				end
			end
		else
			rx[cc] = rix
			nextx[cc] = rand(expdist) / rx[cc]
			if cc > Ncells-Ni2
				for dd = 1:Ne
					indx = round(Int, cc-Ncells+Ni2)
					sumwie0[indx] += weights[dd, cc]
					(weights[dd, cc] > 0) && (Nie[indx] += 1)
				end
			end
		end
	end
	# THIS IS THE NEWLY ADDED THING
	# sumwie0 .= mean(sumwie0)

	# # Single cell initialization ieNorm&eiNorm
	# for cc = 1:Ncells
	# 	v[cc] = vre + (vth0 - vre) * rand()
	# 	if cc <= Ne
	# 		rx[cc] = rex
	# 		nextx[cc] = rand(expdist) / rx[cc]
	# 		# Calculate normalization parameters
	# 		for dd = 1:Ne
	# 			sumwee0[cc] += weights[dd, cc]
	# 			if weights[dd, cc] > 0
	# 				Nee[cc] += 1
	# 			end
	# 		end
	# 		for dd = (Ncells-Ni2+1):Ncells
	# 			sumwei20[cc] += weights[dd, cc]
	# 			if weights[dd, cc] > 0
	# 				Nei2[cc] += 1
	# 			end
	# 		end
	# 	else
	# 		rx[cc] = rix
	# 		nextx[cc] = rand(expdist) / rx[cc]
	# 		if cc > Ncells-Ni2
	# 			for dd = 1:Ne
	# 				indx = round(Int, cc-Ncells+Ni2)
	# 				sumwie0[indx] += weights[dd, cc]
	# 				(weights[dd, cc] > 0) && (Nie[indx] += 1)
	# 			end
	# 		end
	# 	end
	# end

	@info "Starting simulation"
	# --- Begin main simulation loop ---
	iterSteps = ProgressBar(1:Nsteps)
	@inbounds @fastmath for tt in iterSteps
		
		t = dt * tt
		tprev = dt * (tt - 1)

		if mod(tt, Nsteps/100) == 1
			# Kill switch for epileptic spiking activity
			avgHzE = mean(1000 * ns[1:Ne] / (tt*dt))
			avgHzI = mean(1000 * ns[(Ne + 1):Ncells] / (tt*dt))
			avgHzI1 = mean(1000 * ns[(Ne + 1):(Ncells - Ni2)] / (tt*dt))
			avgHzI2 = mean(1000 * ns[(Ncells - Ni2 + 1):Ncells] / (tt*dt))
			if maximum([avgHzE, avgHzI, avgHzI1, avgHzI2]) > 20.
				@info "\r", "Epileptic spiking activity detected. (Firing rate: ", maxfr, ")"
				exit()
			end
		end

		@. forwardInputsE = 0.
		@. forwardInputsI = 0.

		# Check if we have entered or exited a stimulation period
		for ss = 1:size(stim)[1]
			if (tprev < stim[ss, 2]) && (t >= stim[ss, 2])	# Just entered stimulation period
				ipop = round(Int, stim[ss, 1])
				for ii = 1:Nmaxmembers
					(popmembers[ipop, ii] == 0 || popmembers[ipop, ii] == -1) && break
					rx[popmembers[ipop, ii]] += stim[ss, 4]	# Add external drive
				end
			end
			if (tprev < stim[ss, 3]) && (t >= stim[ss, 3]) 	# Just exited stimulation period
				ipop = round(Int, stim[ss, 1])
				for ii = 1:Nmaxmembers
					(popmembers[ipop, ii] == 0 || popmembers[ipop, ii] == -1) && break
					rx[popmembers[ipop, ii]] -= stim[ss, 4]	# Subtract external drive
				end
			end
		end  # End loop over stimuli

		# if mod(tt, inormalize) == 0 # Excitatory synaptic normalization, additive
		# 	for cc = 1:Ne
		# 		sumwee = 0.
		# 		for dd = 1:Ne
		# 			sumwee += weights[dd,cc]
		# 		end
		# 		for dd = 1:Ne
		# 			if weights[dd,cc] > 0.
		# 				weights[dd,cc] -= (sumwee-sumwee0[cc])/Nee[cc]
		# 				if weights[dd,cc] < jeemin
		# 					weights[dd,cc] = jeemin
		# 				elseif weights[dd,cc] > jeemax
		# 					weights[dd,cc] = jeemax
		# 				end
		# 			end
		# 		end
		# 	end
		# end # End normalization

		# # Excitatory synaptic scaling, multiplicative
		# if mod(tt, inormalize) == 0
		# 	for cc = 1:Ne
		# 		sumwee = 0.
		# 		for dd = 1:Ne
		# 			sumwee += weights[dd, cc]
		# 		end
		# 		for dd = 1:Ne
		# 			if weights[dd, cc] > 0.
		# 				weights[dd, cc] = (weights[dd, cc] / sumwee) * sumwee0[cc]
		# 				if weights[dd, cc] < jeemin
		# 					weights[dd, cc] = jeemin
		# 				elseif weights[dd, cc] > jeemax
		# 					weights[dd, cc] = jeemax
		# 				end
		# 			end
		# 		end
		# 	end
		# end  # End normalization

		# # I₂➡E synaptic scaling, multiplicative
		# if mod(tt, inormalize) == 0
		# 	for cc = 1:Ne
		# 		sumwee = 0.
		# 		sumwei = 0.
		# 		for dd = 1:Ne
		# 			sumwee += weights[dd, cc]
		# 		end
		# 		for dd = (Ncells-Ni2+1):Ncells
		# 			sumwei += weights[dd, cc]
		# 		end
		# 		for dd = 1:Ne
		# 			if weights[dd, cc] > 0.
		# 				weights[dd, cc] = (weights[dd, cc] / sumwee) * sumwee0[cc]
		# 				if weights[dd, cc] < jeemin
		# 					weights[dd, cc] = jeemin
		# 				elseif weights[dd, cc] > jeemax
		# 					weights[dd, cc] = jeemax
		# 				end
		# 			end
		# 		end
		# 		for dd = (Ncells-Ni2+1):Ncells
		# 			if weights[dd, cc] > 0.
		# 				weights[dd, cc] = (weights[dd, cc] /sumwei) * sumwei20[cc]
		# 				if weights[dd, cc] < jei2min
		# 					weights[dd, cc] = jei2min
		# 				elseif weights[dd, cc] > jei2max
		# 					weights[dd, cc] = jei2max
		# 				end
		# 			end
		# 		end
		# 	end
		# end  # End normalization

		# # I₂➡E synaptic normalization, additive
		# if mod(tt, inormalize) == 0
		# 	for cc = 1:Ne
		# 		sumwee = 0.
		# 		sumwei = 0.
		# 		for dd = 1:Ne
		# 			sumwee += weights[dd,cc]
		# 		end
		# 		for dd = (Ncells-Ni2+1):Ncells
		# 			sumwei += weights[dd, cc]
		# 		end
		# 		for dd = 1:Ne
		# 			if weights[dd,cc] > 0.
		# 				weights[dd,cc] -= (sumwee-sumwee0[cc])/Nee[cc]
		# 				if weights[dd,cc] < jeemin
		# 					weights[dd,cc] = jeemin
		# 				elseif weights[dd,cc] > jeemax
		# 					weights[dd,cc] = jeemax
		# 				end
		# 			end
		# 		end
		# 		for dd = (Ncells-Ni2+1):Ncells
		# 			if weights[dd, cc] > 0.
		# 				weights[dd, cc] -= (sumwei-sumwei20[cc])/Nei2[cc]
		# 				if weights[dd, cc] < jei2min
		# 					weights[dd, cc] = jei2min
		# 				elseif weights[dd, cc] > jei2max
		# 					weights[dd, cc] = jei2max
		# 				end
		# 			end
		# 		end
		# 	end
		# end # End normalization

		# # E➡E & I➡E synaptic scaling, multiplicative
		# if mod(tt, inormalize) == 0
		# 	for cc = 1:Ne
		# 		sumwee = 0.
		# 		sumwei = 0.
		# 		for dd = 1:Ne
		# 			sumwee += weights[dd, cc]
		# 		end
		# 		for dd = (Ncells-Ni2+1):Ncells
		# 			sumwei += weights[dd, cc]
		# 		end
		# 		for dd = 1:Ne
		# 			if weights[dd, cc] > 0.
		# 				weights[dd, cc] = (weights[dd, cc] / sumwee) * sumwee0[cc]
		# 				if weights[dd, cc] < jeemin
		# 					weights[dd, cc] = jeemin
		# 				elseif weights[dd, cc] > jeemax
		# 					weights[dd, cc] = jeemax
		# 				end
		# 			end
		# 		end
		# 		for dd = (Ne+1):Ncells
		# 			if weights[dd, cc] > 0.
		# 				weights[dd, cc] = (weights[dd, cc] / sumwei) * sumwei0[cc]
		# 				if dd > (Ncells - Ni2 +1)
		# 					if weights[dd, cc] < jei2min
		# 						weights[dd, cc] = jei2min
		# 					elseif weights[dd, cc] > jei2max
		# 						weights[dd, cc] = jei2max
		# 					end
		# 				else
		# 					if weights[dd, cc] < jeimin
		# 						weights[dd, cc] = jeimin
		# 					elseif weights[dd, cc] > jeimax
		# 						weights[dd, cc] = jeimax
		# 					end
		# 				end
		# 			end
		# 		end
		# 	end
		# end  # End normalization


		# Excitatory synaptic scaling (& etoiNorm), multiplicative	NOTE::: THIS WAS THE CORRECT ONE
		if mod(tt, inormalize) == 0
			for cc = 1:Ne
				sumwee = 0.
				for dd = 1:Ne
					sumwee += weights[dd, cc]
				end
				for dd = 1:Ne
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
				for dd = 1:Ne
					sumwie += weights[dd, cc]
				end
				for dd = 1:Ne
					if weights[dd, cc] > 0
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


		# # Excitatory synaptic scaling(&ie/eiNorm), multiplicative
		# if mod(tt, inormalize) == 0
		# 	for cc = 1:Ne
		# 		sumwee = 0.
		# 		sumwei = 0.
		# 		for dd = 1:Ne
		# 			sumwee += weights[dd, cc]
		# 		end
		# 		for dd = (Ncells-Ni2+1):Ncells
		# 			sumwei += weights[dd, cc]
		# 		end
		# 		for dd = 1:Ne
		# 			if weights[dd, cc] > 0.
		# 				weights[dd, cc] = (weights[dd, cc] / sumwee) * sumwee0[cc]
		# 				if weights[dd, cc] < jeemin
		# 					weights[dd, cc] = jeemin
		# 				elseif weights[dd, cc] > jeemax
		# 					weights[dd, cc] = jeemax
		# 				end
		# 			end
		# 		end
		# 		for dd = (Ncells-Ni2+1):Ncells
		# 			if weights[dd, cc] > 0.
		# 				weights[dd, cc] = (weights[dd, cc] /sumwei) * sumwei20[cc]
		# 				if weights[dd, cc] < jei2min
		# 					weights[dd, cc] = jei2min
		# 				elseif weights[dd, cc] > jei2max
		# 					weights[dd, cc] = jei2max
		# 				end
		# 			end
		# 		end
		# 	end
		# 	for cc = (Ncells-Ni2+1):Ncells
		# 		indx = cc - Ncells + Ni2
		# 		sumwie = 0.
		# 		for dd = 1:Ne
		# 			sumwie += weights[dd, cc]
		# 		end
		# 		for dd = 1:Ne
		# 			if weights[dd, cc] > 0
		# 				weights[dd, cc] = (weights[dd, cc] / sumwie) * sumwie0[indx]
		# 				if weights[dd, cc] < jiemin
		# 					weights[dd, cc] = jiemin
		# 				elseif weights[dd, cc] > jiemax
		# 					weights[dd, cc] = jiemax
		# 				end
		# 			end
		# 		end
		# 	end
		# end  # End normalization
	

		# --- Update single cells ---
		@. spiked = 0.

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

				if cc <= Ne  # Excitatory neuron (EIF), has adaptation
					v[cc] += dt * ((vleake - v[cc] + deltathe * exp((v[cc] - vth[cc]) / deltathe)) / taue + ge * (erev - v[cc]) / C + gi * (irev - v[cc]) / C - wadapt[cc] / C);
					if v[cc] > vpeak
						spiked[cc] = true
						wadapt[cc] += bw_adapt
					end
				else
					v[cc] += dt * ((vleaki - v[cc]) / taui + ge * (erev - v[cc]) / C + gi * (irev - v[cc]) / C);
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
						if weights[dd, cc] == 0.
							continue
						end
						weights[dd, cc] += eta * trace_istdp[dd]
						if weights[dd, cc] > jeimax
							weights[dd, cc] = jeimax
						end
					end
				elseif cc <= (Ncells - Ni2)
					# 1st i-population neuron fired, modify outputs to excitatory neurons
					for dd = 1:Ne
						if weights[cc, dd] == 0.
							continue
						end
						weights[cc, dd] += eta * (trace_istdp[dd] - 2 * r0 * tauy)
						if weights[cc, dd] > jeimax
							weights[cc, dd] = jeimax
						elseif weights[cc, dd] < jeimin
							weights[cc, dd] = jeimin
						end
					end
				end
			end	# End, iSTDP₁

			# _____ iSTDP₂ _____
			if spiked[cc] && (t > stdpdelay) && (t < stim[end,3])
				if cc <= Ne
					# Excitatory neuron fired, modify inputs from 2nd i-population
					for dd = (Ncells - Ni2 + 1):Ncells
						if weights[dd, cc] == 0.
							continue
						end
						weights[dd, cc] -= ilamda * ((1 - weights[dd, cc])^mi) * trace_expDecay[dd]
						if weights[dd, cc] < jei2min
							weights[dd, cc] = jei2min
						end
					end
				elseif cc > (Ncells - Ni2)
					# 2nd i-population neuron fired, modify outputs to excitatory neurons
					for dd = 1:Ne
						if weights[cc, dd] == 0.
							continue
						end
						weights[cc, dd] += ilamda * alfa * (weights[cc, dd]^mi) * trace_expDecay[dd]
						if weights[cc, dd] > jei2max
							weights[cc, dd] = jei2max
						end
					end
				end
			end # End, iSTDP₂

			#  _____ eiSTDP _____
			if spiked[cc] && (t > stdpdelay)# && (t < stim[end,3])
				if cc <= Ne
					# Excitatory neuron fired, modify outputs to 2nd i-population
					for dd = (Ncells-Ni2+1):Ncells
						if weights[cc, dd] == 0.
							continue
						end
						weights[cc, dd] += eta_ie * trace_eistdp[dd]
						if weights[cc, dd] > jiemax
							weights[cc, dd] = jiemax
						end
					end
				elseif cc > (Ncells-Ni2)
					# 2nd i-population neuron fired, modify inputs from E neurons
					for dd = 1:Ne
						if weights[dd, cc] == 0.
							continue
						end
						weights[dd, cc] += eta_ie * (trace_eistdp[dd] - Adep_ie)
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
					if weights[cc, dd] == 0.
						continue
					end
					if u_vstdp[dd] > thetaltd
						weights[cc, dd] -= altd * (u_vstdp[dd] - thetaltd)
						if weights[cc, dd] < jeemin
							weights[cc, dd] = jeemin
						end
					end
				end
			end  # End, LTD

			# vSTDP, LTP component
			if (t > stdpdelay) && (cc < Ne) && (v[cc] > thetaltp) && (v_vstdp[cc] > thetaltd)
				for dd = 1:Ne
					if weights[dd, cc] == 0.
						continue
					end
					weights[dd, cc] += dt * altp * x_vstdp[dd] * (v[cc] - thetaltp) * (v_vstdp[cc] - thetaltd)
					if weights[dd, cc] > jeemax
						weights[dd, cc] = jeemax
					end
				end
			end # End LTP
		end  # End, loop over cells
		@. forwardInputsEPrev = copy(forwardInputsE)
		@. forwardInputsIPrev = copy(forwardInputsI)
	end # End, loop over time

	print("\r")
	times = times[:, 1:maximum(ns)]

	return times, ns, Ne, Ncells, T, weights
end