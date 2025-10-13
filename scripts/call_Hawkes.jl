include("../src/simu_Hawkes.jl")
using CairoMakie

begin # parameters 
    μ = [0.8, 1.2]
    α = 4.
    β = [0.5 -0.2; -0.5 0.3]
    pp = MultiHawkesProcess(μ, α, β)
end

begin # time grid
    tmin = 0.
    tmax = 5.
    dt = 0.1
end

begin # simulation
    Random.seed!(1)
    ts, λs, h = rand(pp, tmin, tmax, dt)
    spikes1 = event_times(h)[event_marks(h).==1]
    spikes2 = event_times(h)[event_marks(h).==2]
end

begin # plot
    fig = Figure(size=(800, 400))

    ax = Axis(fig[1, 1], ylabel=L"\lambda", xlabel=L"t", title="Bivariate Hawkes process")

    lines!(ax, ts, λs[1, :], color=:darkblue, label=L"\lambda^1_t")
    lines!(ax, ts, λs[2, :], color=:orange, label=L"\lambda^2_t")
    hlines!(ax, μ[1], color=:darkblue, linestyle=:dot)
    hlines!(ax, μ[2], color=:orange, linestyle=:dot)
    scatter!(ax, spikes1, 0.0 * spikes1, color=:darkblue, markersize=15, label=L"N^1")
    scatter!(ax, spikes2, 0.0 * spikes2, color=:orange, markersize=15, label=L"N^2")

    Legend(fig[1, 2], ax)

    fig
end

# save figure
#save("plots/Hawkes-process.png", fig, px_per_unit=2)