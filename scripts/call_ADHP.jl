include("../src/simu_ADHP.jl")
using CairoMakie

# Sepúlveda 2023, Example 2
φ(ξ) = 10.0 * ξ^2 / (ξ^2 + 1.0) + 0.5
φ̅ = 10.5
# Modified
φ(ξ) = 5.0 * ξ / (ξ + 2.0) + 0.5
φ̅ = 5.5

begin # parameters
    adhp = AgeDependentHawkesProcess(
        n=Int(1e4),
        refractoryperiod=1.0,
        decayrate=10.0,
        firingrate=φ,
        firingratebound=φ̅
    )
    initial_condition = rand(Exponential(), adhp.n) .+ 1
end
begin # time domain
    domain = [0., 10.]
    length_ts = 1000
end

# Simulation
h, ξs, activities = simulate(adhp, initial_condition, domain)
begin # plot of the activity
    ids = 1:floor(Int, length(h) / length_ts):length(h)
    ts = event_times(h)[ids]
    ys = activities[ids]
    lines(ts, ys)
end

begin # raster plot
    fig = Figure()

    axscatt = Axis(fig[1, 1], ylabel="neuron",)
    axlines = Axis(fig[1, 1], yticklabelcolor=:red, yaxisposition=:right, ylabel="mean activity", ylabelcolor=:red)
    hidespines!(axlines)
    hidexdecorations!(axlines)
    hideydecorations!(axscatt, label=false)

    ids = findall(event_marks(h) .!= 0)
    scatter!(axscatt, event_times(h)[ids], event_marks(h)[ids], markersize=2, color=:black)

    ids = 1:floor(Int, length(h) / length_ts):length(h)
    ts = event_times(h)[ids]
    ys = activities[ids]
    lines!(axlines, ts, ys, color=:red)

    fig
end

using CairoMakie

fig = Figure()

ax1 = Axis(fig[1, 1], yticklabelcolor=:blue)
ax2 = Axis(fig[1, 1], yticklabelcolor=:red, yaxisposition=:right)
hidespines!(ax2)
hidexdecorations!(ax2)

lines!(ax1, 0 .. 10, sin, color=:blue)
lines!(ax2, 0 .. 10, x -> 100 * cos(x), color=:red)

fig