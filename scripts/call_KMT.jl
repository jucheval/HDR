include("../src/coupling_KMT.jl")
using CairoMakie

begin # time grid
    len_ts = 2^12
    dt = 2^(-8)
    tmax = len_ts * dt
end

begin # simulation

end

begin # simulation and plot
    Random.seed!(1)
    fig = Figure()

    ax = Axis(fig[1, 1], xlabel=L"t")
    kwargs = (; alpha=0.7)

    time, X, W = KMT_simulation(len_ts, dt)
    N = counting2point(time, X)
    lines!(ax, time, X .- time, alpha=0.7, color=:red)
    lines!(ax, time, W .- time; color=:darkblue, kwargs...)

    time, X, W = KMT_simulation(len_ts, dt)
    N = counting2point(time, X)
    lines!(ax, time, X .- time, alpha=0.7, color=:red, linestyle=:dashdot)
    lines!(ax, time, W .- time; color=:darkblue, linestyle=:dashdot, kwargs...)

    fig
end