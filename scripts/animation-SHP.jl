include("../src/simu_SHP.jl")
using Random
using DrWatson

function square_torusdistance(a::Point2d, b::Point2d)
    xdistance = rem(a[1] - b[1], 1, RoundNearest)
    ydistance = rem(a[2] - b[2], 1, RoundNearest)
    return xdistance^2 + ydistance^2
end

begin # 2D Neural field with 900 neurons which exhibit a wandering bump
    begin # parameters
        sqN = 30
        Nneur = sqN^2
        κ = 0.2
        ϵ = 0.03
        A₁ = 30.0
        A₂ = 5.0
        σ₁ = 0.1
        σ₂ = 0.2
        f(u) = (1 + exp(-(u - κ) / ϵ))^(-1)
        function w(y, x)
            r = square_torusdistance(x, y)

            return A₁ * exp(-r / σ₁^2) - A₂ * exp(-r / σ₂^2)
        end
        α = 1.0
        u₀(x) = 2 * w(x, Point2d(0.5, 0.5))
    end
    begin # time grid
        tmin = 0.0
        tmax = 150.0
        xmin = ymin = 0.0
        xmax = ymax = 1.0
        sdpositions = 1 / sqN / 2
        xs = repeat(range(xmin, xmax, sqN + 1)[begin:end-1], inner=sqN)
        xs += sdpositions * rand(Nneur)
        xs = rem.(xs, 1, RoundDown)
        ys = repeat(range(ymin, ymax, sqN + 1)[begin:end-1], outer=sqN)
        ys += sdpositions * rand(Nneur)
        ys = rem.(ys, 1, RoundDown)
        domain = [tmin, tmax]
        shp = SpatialHawkesProcess(
            n=Nneur,
            decayrate=α,
            synapticweight=(y, x) -> (xmax - xmin) * (ymax - ymin) * w(y, x),
            positions=[Point2d(xs[i], ys[i]) for i in 1:Nneur],
            firingrate=f,
            firingratebound=1.0
        )
    end
end;

# simulation
Random.seed!(1)
simulation = simulate(shp, u₀, domain, length_ts=1000)

safesave("data/SHP-animation.jld2", Dict("history" => simulation[1], "positions" => simulation[3]))