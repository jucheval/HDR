using PointProcesses
using Random
using CairoMakie

begin # parameters
    n = 10
    tmin = 0.0
    tmax = 8π
    ϵ = 1
    Λ(t) = t + ϵ * sin(t)
end

begin # simulation
    Random.seed!(1)
    h = rand(MultivariatePoissonProcess(ones(n)), tmin, tmax)
    h = time_change(h, Λ)
end

begin # plot
    fig = Figure(size=(300, 300))

    ax = Axis(fig[1, 1], ylabel="neuron", xlabel=L"t", title="Raster plot")

    scatter!(ax, event_times(h), event_marks(h), markersize=6, color=:black)

    fig
end

# save figure
#save("plots/raster-plot.png", fig, px_per_unit=2)