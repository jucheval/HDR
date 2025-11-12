using MeanFieldGraph
using Random
using DrWatson

begin # parameters
    n = 900
    μ = 0.01
    γ = 0.1
    p = 0.8
    model = MarkovChainModel(μ, 1 - γ, p)
    r₊ = 0.5
end
excitatory = MeanFieldGraph.N2excitatory(n, r₊) # encodes the two classes

# simulation
Random.seed!(1)
simulation = rand(model, excitatory, 2000)

safesave("data/MCRE-animation.jld2", Dict("spikes" => simulation.X))