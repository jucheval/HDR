using MeanFieldGraph
using Random
using CairoMakie

function Makie.plot(::Type{MarkovChainModel}, simulation)
    fig = Figure(size=(600, 300))

    ax = Axis(fig[1, 1], xlabel=L"t", ylabel="neuron")
    ax.yticks = [1, findfirst(excitatory .== 0), n]

    scatter!(ax, [p.I for p in findall(transpose(simulation.X))], markersize=5, color=:black)

    fig
end

begin # parameters
    n = 50
    μ = 0.1
    γ = 0.7
    p = 0.8
    model = MarkovChainModel(μ, 1 - γ, p)
    r₊ = 0.7
end
excitatory = MeanFieldGraph.N2excitatory(n, r₊) # encodes the two classes

# simulation
Random.seed!(1)
simulation = rand(model, excitatory, 100)
# plot
fig = Makie.plot(MarkovChainModel, simulation)

# save figure
#save("plots/MCRE-1.png", fig, px_per_unit=2)


begin # parameters
    n = 50
    μ = 0.3
    γ = 0.3
    p = 0.5
    model = MarkovChainModel(μ, 1 - γ, p)
    r₊ = 0.4
end
excitatory = MeanFieldGraph.N2excitatory(n, r₊) # encodes the two classes

# simulation
Random.seed!(1)
simulation = rand(model, excitatory, 100)
# plot
fig = Makie.plot(MarkovChainModel, simulation)

# save figure
#save("plots/MCRE-2.png", fig, px_per_unit=2)