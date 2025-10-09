include("../src/coupling_diffusion_SDE.jl")
using CairoMakie

begin # parameters
    sd = 0.3
    offset = 0.2
    diffusion = Diffusion(
        drift=x -> -x,
        standarddeviation=x -> sd
    )
    diffusiontilde = Diffusion(
        drift=x -> -(x - offset),
        standarddeviation=x -> sd
    )
    initial_conditions = [0.0, offset]
end

begin # time grid 
    tmin = 0.
    tmax = 100.
    domain = [tmin, tmax]
    dt = 0.1
end

begin # plot
    Random.seed!(1)
    fig = Figure()
    ax = Axis(fig[1, 1], yticks=(
        [1.5, 2.0, 2.5, 3.5, 4.0, 4.5, 5.5, 6.0, 6.5],
        string.(repeat([-1, 0, 1], 3))))
    for k in 1:3
        ts, Xs, X̃s = coupling(diffusion, diffusiontilde, initial_conditions, domain, dt)
        lines!(ax, ts, Xs .+ 2 * k, alpha=0.8, color=:darkblue)
        lines!(ax, ts, X̃s .+ 2 * k, alpha=0.8, color=:red)
    end
    hlines!(ax, [3, 5], color=:grey)
    fig
end