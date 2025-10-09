using MeanFieldGraph

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
data = rand(model, excitatory, 100)

begin # plot
    Random.seed!(1)
    fig = Figure()

    ax = Axis(fig[1, 1], xlabel=L"t", ylabel="neuron")
    ax.yticks = [1, findfirst(excitatory .== 0), n]

    scatter!(ax, [p.I for p in findall(transpose(data.X))])

    fig
end
