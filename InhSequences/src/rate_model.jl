function rate_simnew(T::Int64; random_seed::Int64=2817)
	#  Generates new weights and populations for the rate model, runs simulation   
    @info "Setting up weights"
    Random.seed!(random_seed)
	Np::Int64 = 4
    Ne::Int64 = Np + 1
    Npop::Int64 = 2 * Ne
    aux::UnitRange{Int64} = (Npop-Np):(Npop-1)    # Auxiliary variable with indexes of I-assemblies

    # --- Set up weights ---
	weights::Matrix{Float64} = zeros(Npop, Npop)
    
    @. weights[1:Np, 1:Np] = .01                # Weights between E-assemblies
    @. weights[1:Np, Np+1] = .01                # From E-assemblies to E non-members
    @. weights[Np+1, 1:Np] = .01                # From E non-members to E-assemblies
	
    @. weights[1:Np, Npop] = .01                # From E-assemblies to I non-members
    @. weights[Npop, 1:Np] = .01                # From I non-members to E-assemblies

    @. weights[aux, aux] = .01                  # From I-assemblies to I-assemblies
    @. weights[aux, Np+1] = .01                 # From I-assemblies to E non-members
    @. weights[Np+1, aux] = .01                 # From E non-members to I-assemblies

    @. weights[aux, Npop] = .01                 # From I-assemblies to I non-members
    @. weights[Npop, aux] = .01                 # From I non-members to I-assemblies

    weights[Np+1, Npop] = .01                   # From E non-members to I non-members
    weights[Np+1, Np+1] = .01                   # From E non-members to E non-members
    weights[Npop, Np+1] = .01                   # From I non-members to E non-members
    weights[Npop, Npop] = .01                   # From I non-members to I non-members
	

    @. weights[1:Np, aux] = .01                 # From E-assemblies to I-assemblies (default)
    @. weights[aux, 1:Np] = .01                 # From I-assemblies to E-assemblies (default)
	
    # Feedforward structure
    for pop = 1:Np
        weights[pop, (Np+1+pop)] = .04                      # From E-assemblies to I-assemblies (E/I)
        weights[(Np+1+pop), pop] = .01                      # From I-assemblies to E-assemblies (E/I)
        (pop > 1) && (weights[(Ne+pop), pop-1] = .01)       # From I-assemblies to E-assemblies (previous)
        (pop < Np) && (weights[(Ne+pop), pop+1] = .01)      # From I-assemblies to E-assemblies (next)
    end

	return rate_sim(weights, T)
end

function rate_sim(weights::Matrix{Float64}, T::Int64; random_seed::Int64=2817)
    # Runs the simulation given weight matrix and populations
    # NOTE: The structure of vectors is as follows: E-assemblies, E non-members, I-assemblies, I non-members
    # _________________________--- Simulation ---______________________________
    @info "Setting up parameters..."    
    Random.seed!(random_seed)
    dt::Float64 = 1.
    Nsteps::Int64 = round(Int, T/dt)
    Np::Int64 = 4
    Ne::Int64 = Np+1
    Npop::Int64 = 2*Ne

    x::Vector{Float64} = zeros(Npop)        # Currents
    xpre::Vector{Float64} = zeros(Npop)     # Need it to calculate input (for now at least)
    taue::Float64 = 20.                     # Time constants for E
    taui::Float64 = 20.                     # Time constants for I
    
    # rex::Float64 = 0. 					    # Rate of external input to E (kHz)
    # rix::Float64 = 0. 					    # Rate of external input to I (kHz)
    jex::Vector{Float64} = ones(Npops) * .01   # Strength of external input to E
    jix::Vector{Float64} = ones(Npops) * .01   # Strength of external input to I

    # Auxiliary variables
    input::Float64 = 0.
    t::Float64 = 0.
    ipop::Int64 = 0
    stim_strength::Float64 = .01

    gausdist = Normal(0, .00001)

    # --- TRACKERS ---
    tracker_indx::Int64 = 1
    tracker_dt::Float64 = 1. #round(Int, T/1_000)
    Nsteps_tr::Int64 = round(Int, T/tracker_dt)
    tracker_x::Vector{Vector{Float64}} = [zeros(Nsteps_tr) for _ in 1:Npop]
    tracker_dt /= dt
    @info "Starting simulation"
    # --- Begin main simulation loop ---
    iterSteps = ProgressBar(1:Nsteps)
    @inbounds @fastmath for tt in iterSteps

        t = dt * tt

        # Stimulation
        if tt == 200
            jex[1] += stim_strength
        end
        if tt == 300
            jex[1] -= stim_strength
        end

        # --- Update population ---
        for ipop = 1:Npop
            if ipop <= Ne
                # Excitatory population
                input = sum(weights[1:Ne, ipop] .* Φ.(xpre[1:Ne])) - sum(weights[Ne+1:Npop, ipop] .* Φ.(xpre[Ne+1:Npop])) + jex[ipop]
                x[ipop] += dt/taue * (- xpre[ipop] + max(input, 0.)) + rand(gausdist)
            else
                # Inhibitory assembly
                input = sum(weights[1:Ne, ipop] .* Φ.(xpre[1:Ne])) - sum(weights[Ne+1:Npop, ipop] .* Φ.(xpre[Ne+1:Npop])) + jix[ipop]
                x[ipop] += dt/taui * (- xpre[ipop] + max(input, 0.)) + rand(gausdist)
            end
        end
        @. xpre = copy(x)

        # --- TRACKERS ---
        if mod(tt-1, tracker_dt) == 0
            for ipop = 1:Npop
                tracker_x[ipop][tracker_indx] = x[ipop]
            end
            tracker_indx += 1
        end
    end # End, loop over time

    print("\r")

    return weights, tracker_x
end

# Define the sigmoid function Φ
Φ(x::Float64) = 1. / (1. + exp(-x))