module InhSequences
    using Revise
    using UnPack
    using Random
    using Parameters
    using Statistics
    using ProgressBars
    using LinearAlgebra
    using Distributions

    include("structs.jl")
    include("simulation.jl")
    include("utils.jl")

    export simnew, sim
    export make_seq, make_stim, findI2populations, convolveSpikes, getPopulationRates
    export NeuronalParams, SynapticParams, PlasticityParams
end # module InhSequences
