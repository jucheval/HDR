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

begin # Simulation and plot of the activity
    h, ξs, activities = simulate(adhp, initial_condition, domain)
    ids = 1:floor(Int, length(h) / length_ts):length(h)
    ts = event_times(h)[ids]
    ys = activities[ids]
    lines(ts, ys)
end