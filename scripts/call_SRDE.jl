include("../src/simu_RDE.jl")

# Periodic solution taken from Sepúlveda 2023, Example 2
begin # Parameters
    φ(ξ) = 10.0 * ξ^2 / (ξ^2 + 1.0) + 0.5
    srde = StochasticRDE(
        n=100.,
        refractoryperiod=1.0,
        decayrate=10.0,
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
    length_ts = 2000
end
begin # initial condition
    u₀(a) = pdf(Exponential(), a - 1.)
end

# simulation
Random.seed!(1)
simulation = simulate(srde, u₀, domains, dt; saveat=Int((tmax - tmin) / dt / length_ts));
# plot
fig = plot(StochasticRDE, simulation);

# save figure
#save("plots/SRDE_1.png", fig, px_per_unit=2)


# Aperiodic solution
begin # Parameters
    φ(ξ) = 5.0 * ξ / (ξ + 2.0) + 0.5
    srde = StochasticRDE(
        n=100.,
        refractoryperiod=1.0,
        decayrate=10.0,
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
    length_ts = 2000
end
begin # initial condition
    u₀(a) = pdf(Exponential(), a - 1.)
end

# simulation
Random.seed!(1)
simulation = simulate(srde, u₀, domains, dt; saveat=Int((tmax - tmin) / dt / length_ts));
fig = plot(StochasticRDE, simulation);

# save figure
#save("plots/SRDE_2.png", fig, px_per_unit=2)
