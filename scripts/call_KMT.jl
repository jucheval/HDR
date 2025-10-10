include("../src/coupling_KMT.jl")
using CairoMakie

begin # time grid
    len_ts = 2^12
    dt = 2^(-8)
    tmax = len_ts * dt
end

begin # simulation and plot
    Random.seed!(2)
    fig = Figure(size=(800, 400))

    ax = Axis(fig[1, 1], xlabel=L"t")
    kwargs = (; alpha=0.7)

    time, X, W = KMT_simulation(len_ts, dt)
    N = counting2point(time, X)
    lines!(ax, time, X .- time, alpha=0.7, color=:darkblue)
    lines!(ax, time, W .- time; color=:orange, kwargs...)

    time, X, W = KMT_simulation(len_ts, dt)
    N = counting2point(time, X)
    lines!(ax, time, X .- time, alpha=0.7, color=:darkblue, linestyle=:dashdot)
    lines!(ax, time, W .- time; color=:orange, linestyle=:dashdot, kwargs...)

    fig
end

# save figure
#save("plots/coupling_KMT.png", fig, px_per_unit=2)