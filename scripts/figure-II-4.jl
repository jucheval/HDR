include("../src/simu_NFE.jl")

begin # Coombes and Laing 2011, Figure 1
    begin # parameters
        J(x) = 1 + 0.4 * sin(x)
        w(y, x) = J(y) * exp(-abs(x - y)) / 2
        κ = 1 / 20
        ρ = 0.3
        f(u) = (1 + exp(-(u - ρ) / κ))^(-1)
        α = 1.0
        snfe = StochasticNFE(
            n=Inf,
            decayrate=α,
            synapticweight=w,
            firingrate=f
        )
        u₀(x) = 1.0 * (x < -20)
    end
    begin # time grid
        tmin = 0.0
        tmax = 100.0
        xmin = -50
        xmax = 50
        domains = [
            [tmin, tmax],
            [xmin, xmax]
        ]
        dt = 1e-1
        dx = 1e-1
        length_ts = 1000
    end
end

# simulation
simulation = simulate(snfe, u₀, domains, dt, dx; saveat=convert(Integer, fld(tmax - tmin, length_ts * dt)))
# plot
fig = plot(StochasticNFE, simulation)

# save figure
#save("plots/NFE_1.png", fig, px_per_unit=2)


begin # Agathe Nerine 2025, Figure 2
    begin # parameters
        κ = 1 / 20
        ρ = 1 / 2
        f(u) = (1 + exp(-(u - ρ) / κ))^(-1)
        w(y, x) = 2 * pi * cos(y - x)
        α = 1.0
        snfe = StochasticNFE(
            n=Inf,
            decayrate=α,
            synapticweight=w,
            firingrate=f
        )
        u₀(x) = 1.8 * cos(x)
    end
    begin # time grid
        tmin = 0.0
        tmax = 100.0
        xmin = -pi
        xmax = pi
        domains = [
            [tmin, tmax],
            [xmin, xmax]
        ]
        dt = 1e-1
        dx = 3e-2
        length_ts = 1000
    end
end

# simulation
simulation = simulate(snfe, u₀, domains, dt, dx; saveat=convert(Integer, fld(tmax - tmin, length_ts * dt)))
# plot
fig = plot(StochasticNFE, simulation)

# save figure
#save("plots/NFE_2.png", fig, px_per_unit=2)