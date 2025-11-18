using MeanFieldGraph
using Random
using GLMakie

fig = Figure(size=(400, 500))

begin # parameters
    r₊ = 0.5
    μ = 0.25
    γ = 0.5
    p = 0.5
    model = MarkovChainModel(μ, 1.0 - γ, p)
    r₋ = 1.0 - r₊
    D = γ * p * (r₊ - r₋)
    m = (μ + γ * p * r₋) / (1 - D)
    σ₊, σ₋ = (γ * p * m * (1 - m)) .* ((γ^2 * p^2 * (r₊ - r₋)) / (1 - D^2) .+ (1.0, -1.0))

    sg = SliderGrid(fig[1, 1],
        (label="n", range=50:10:200, startvalue=50),
        (label="t", range=100:100:100000, startvalue=100),
        tellwidth=false)
    n, t = [sl.value for sl in sg.sliders]

    perms = Dict([n => randperm(n) for n in 50:10:200])
end

begin # observables
    excitatory = @lift MeanFieldGraph.N2excitatory($n, r₊)
    data = @lift rand(model, $excitatory, 100000)
end

begin # axis and plot
    σ̂ = @lift MeanFieldGraph.covariance_vector($data[1:$t])
    ax = Axis(fig[2, 1], xlabel=L"\hat{\sigma}^n_t", ylabel="neuron",
        yticks=@lift([1, $n]), xticks=([σ₊, σ₋, 0], [L"\sigma_+", L"\sigma_-", "0"]),
        limits=@lift (0.2 .* (-1.0, 1.0), (1 - 0.05 * $n, 1.05 * $n)))

    scatter!(ax, @lift([($σ̂[k], perms[$n][k]) for k in eachindex($σ̂)]), markersize=5, color=:black)

    vlines!(ax, [σ₊, σ₋], color=[:red, :blue])
end

begin # Adjust layout

end
