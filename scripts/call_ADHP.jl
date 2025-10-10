include("../src/simu_ADHP.jl")

# Periodic solution taken from Sepúlveda 2023, Example 2
begin # parameters
    φ(ξ) = 10.0 * ξ^2 / (ξ^2 + 1.0) + 0.5
    φ̅ = 10.5
    adhp = AgeDependentHawkesProcess(
        n=100,
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

# simulation
simulation = simulate(adhp, initial_condition, domain)
# plot
fig = plot(AgeDependentHawkesProcess, simulation);

# save figure
#save("plots/ADHP_1.png", fig, px_per_unit=2)


# Aperiodic solution
begin # parameters
    φ(ξ) = 5.0 * ξ / (ξ + 2.0) + 0.5
    φ̅ = 5.5
    adhp = AgeDependentHawkesProcess(
        n=100,
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

# simulation
simulation = simulate(adhp, initial_condition, domain)
# plot
fig = plot(AgeDependentHawkesProcess, simulation);

# save figure
#save("plots/ADHP_2.png", fig, px_per_unit=2)