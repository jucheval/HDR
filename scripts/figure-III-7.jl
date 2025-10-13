include("../src/simu_SHP.jl")

begin # Coombes and Laing 2011, Figure 1
    begin # parameters
        Nneur = 300
        J(x) = 1 + 0.4 * sin(x)
        w(y, x) = J(y) * exp(-abs(x - y)) / 2
        κ = 1 / 20
        ρ = 0.3
        f(u) = (1 + exp(-(u - ρ) / κ))^(-1)
        α = 1.0
        u₀(x) = 1.0 * (x < -20)
    end
    begin # time grid
        tmin = 0.0
        tmax = 100.0
        domain = [tmin, tmax]
        xmin = -50
        xmax = 50
        shp = SpatialHawkesProcess(
            n=Nneur,
            decayrate=α,
            synapticweight=(y, x) -> (xmax - xmin) * w(y, x),
            positions=range(xmin, xmax, Nneur),
            firingrate=f,
            firingratebound=1.0
        )
    end
end;

# simulation
Random.seed!(1)
simulation = simulate(shp, u₀, domain, length_ts=1000)
# plot
fig = plot(SpatialHawkesProcess, simulation)

# save figure
#save("plots/SHP_1.png", fig, px_per_unit=2)

begin # Agathe Nerine 2025, Figure 2
    begin # parameters
        Nneur = 300
        κ = 1 / 20
        ρ = 1 / 2
        f(u) = (1 + exp(-(u - ρ) / κ))^(-1)
        w(y, x) = 2 * pi * cos(y - x)
        α = 1.0
        u₀(x) = 1.8 * cos(x)
    end
    begin # time grid
        tmin = 0.0
        tmax = 100.0
        xmin = -pi
        xmax = pi
        domains = [tmin, tmax]
        shp = SpatialHawkesProcess(
            n=Nneur,
            decayrate=α,
            synapticweight=(y, x) -> (xmax - xmin) * w(y, x),
            positions=range(xmin, xmax, Nneur),
            firingrate=f,
            firingratebound=1.0
        )
    end
end;

# simulation
Random.seed!(1)
simulation = simulate(shp, u₀, domain, length_ts=1000)
# plot
fig = plot(SpatialHawkesProcess, simulation)

# save figure
#save("plots/SHP_2.png", fig, px_per_unit=2)