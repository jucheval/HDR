include("../src/simu_RDE.jl")
using CairoMakie

# Sepúlveda 2023, Example 2
φ(ξ) = 10.0 * ξ^2 / (ξ^2 + 1.0) + 0.5
# Modified
φ(ξ) = 5.0 * ξ / (ξ + 2.0) + 0.5

begin # Parameters
    srde = StochasticRDE(
        n=Inf,
        refractoryperiod=1.0,
        decayrate=50.0,
        firingrate=φ
    )
end
begin # time grid
    tmin = amin = 0.
    tmax = 10.0
    amax = 3.0
    domains = [
        [tmin, tmax],
        [amin, amax]
    ]
    dt = 1e-3
    length_ts = 1000
end
begin # initial condition
    u₀(a) = pdf(Exponential(), a - 1.)
end

begin # Simulation and plot of the activity
    ts, sol = simulate(srde, u₀, domains, dt; saveat=Int((tmax - tmin) / dt / length_ts))
    activity = sol[:, 1]
    lines(ts, activity)
end