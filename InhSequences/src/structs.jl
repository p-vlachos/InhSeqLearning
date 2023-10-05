@with_kw struct InitializationParameters
	#  Parameters needed to generate weight matrix
	Ne::Int64 = 4000		    # Excitatory no. neurons
	Ni::Int64 = 1000		    # Total Inhibitory no. neurons
	Ni2::Int64 = 500	        # Inhibitory no. neurons (I₂)
	jee0::Float64 = 2.86 	    # Initial E➡E strength (pF)
	jei0::Float64 = 48.7    	# Initial I➡E strength (pF)
	jie::Float64 = 1.27 	    # Initial E➡I1 strength (pF)
	ji2e::Float64 = 1.52 	    # Initial E➡I2 strength (pF)
	jii::Float64 = 16.2 	    # I1➡I1 & I2➡I1 strength (not plastic; pF)
	jii12::Float64 = 24.3	    # I1➡I2 strength (not plastic; pF)
	jii2::Float64 = 32.4	    # I2➡I2 strength (not plastic; pF)
	p::Float64 = 0.2		    # Connection probability
	pmembership::Float64 = .05  # Probability of a neuron to belong to any assembly
	Nmaxmembers::Int64 = 300  	# Maximum number of neurons in a population (to set size of matrix)
end

@with_kw struct NeuronalParameters
    # --- Membrane dynamics parameters ---
    taue::Float64 = 20. 	    # E resting membrane time constant (ms)
	taui::Float64 = 20.  		# I resting membrane time constant (ms)
	vleake::Float64 = -70.  	# E resting potential (mV)
	vleaki::Float64 = -62.  	# I resting potential (mV)
	deltathe::Float64 = 2.  	# EIF slope parameter (mV)
	C::Float64 = 300.  			# Capacitance (pF)
	erev::Float64 = 0.  		# E synapse reversal potential (mV)
	irev::Float64 = -75.  		# I synapse reversal potntial (mV)
	vth0::Float64 = -52.  		# Initial spike voltage threshold (mV)
	ath::Float64 = 10.  		# Increase in threshold post spike (mV)
	tauth::Float64 = 30.  		# Threshold decay timescale (ms)
	vpeak::Float64 = 20.        # Cutoff for voltage.  When crossed, record a spike and reset
	vre::Float64 = -60.  		# Reset potential (mV)
	taurefrac::Float64 = 1.  	# Absolute refractory period (ms)
	aw_adapt::Float64 = 4.  	# Adaptation parameter a (subthreshold adaptation; nS)
	bw_adapt::Float64 = 0.  	# Adaptation parameter b (spike-triggered adaptation; pA) original: .805
	tauw_adapt::Float64 = 150.	# Adaptation timescale (ms)
end

@with_kw struct SynapticParameters
    # --- Synaptic dynamics parameters ---
	tauerise::Float64 = 1. 	# E synapse rise time (ms)
	tauedecay::Float64 = 6.	# E synapse decay time (ms)
	tauirise::Float64 = .5	# I synapse rise time (ms)
	tauidecay::Float64 = 2. # I synapse decay time (ms)
	rex::Float64 = 4.5 		# External input rate to E (khz)
	rix::Float64 = 2.25  	# External input rate to I (khz)
	jex::Float64 = 1.78 	# External to E strength (pF)
	jix::Float64 = 1.27 	# External to I strength (pF)
    # --- Plastic synapses, strength limits ---
	jeemin::Float64 = 1.78  # Minimum EE strength (pF)
	jeemax::Float64 = 21.4  # Maximum EE strength (pF)
	jeimin::Float64 = 48.7 	# Minimum I1➡E strength (pF)
	jeimax::Float64 = 243. 	# Maximum I1➡E strength (pF)
	jiemin::Float64 = .1 	# Minimum E➡I2 strength (pF)
	jiemax::Float64 = 4. 	# Maximum E➡I2 strength (pF)
	jei2min::Float64 = .001	# Minimum I2➡E strength (pF)
	jei2max::Float64 = 243. # Maximum I2➡E strength (pF)
end

@with_kw struct PlasticityParameters
    # --- Voltage based STDP ---
	altd::Float64 = .0008 	    # LTD strength (pA / mV)
	altp::Float64 = .0014 	    # LTP strength (pA / mV^2)
	thetaltd::Float64 = -70.    # LTD voltage threshold (mV)
	thetaltp::Float64 = -49. 	# LTP voltage threshold (mV)
	tauu::Float64 = 10. 		# Timescale for u variable (ms)
	tauv::Float64 = 7.  		# Timescale for v variable (ms)
	taux::Float64 = 15. 		# Timescale for x variable (ms)
	# --- iSTDP₁ ---
	tauy::Float64 = 20. 		# Width of iSTDP₁ curve (ms)
	eta::Float64 = 1. 	    	# iSTDP₁ learning rate (pA)
	r0::Float64 = .003 		    # Target rate (kHz)
	# --- iSTDP₂ ---
	tau_i::Float64 = 100.    	# iSTDP₂ time constant (ms)
	mi::Float64 = 0.			# Weight dependence parameter (zero: no dependence)
	alfa::Float64 = 1.		    # Asymmetry between potentiation and depression (one: symmetric)
	ilamda::Float64 = 1.		# iSTDP₂ learning rate
	# --- eiSTDP ---
	tau_ie::Float64 = 20.		# eiSTDP time constant (ms)
	eta_ie::Float64 = .0015  	# eiSTDP learning rate (pA)
	Adep_ie::Float64 = .12	    # Amplitude of depression (kHz*ms; NOTE: then it has no unit (?))
end

@with_kw struct SimulationParameters
    dt::Float64 = .1            # Integration timestep (ms)
    dtnormalize::Int64 = 20 	# How often to normalize rows of EE weights (ms)
	stdpdelay::Int64 = 10000 	# Time before stdp is activated, to allow transients to die out (ms)
	Nspikes::Int64 = 10000	 	# Maximum number of spikes to record per neuron
end