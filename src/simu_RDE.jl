using Random
using Distributions
using Plots

struct StochasticRDE
    n::Float64
    refractoryperiod::Float64
    decayrate::Float64
    firingrate::Function
end

StochasticRDE(; n, refractoryperiod, decayrate, firingrate) = StochasticRDE(n, refractoryperiod, decayrate, firingrate)

function simulate(srde::StochasticRDE, initial_condition::Vector{Float64}, domains::Vector{Vector{Float64}}, dt::Float64; saveat::Int=1)::Tuple{Vector{Float64},Matrix{Float64}}
    da = dt
    t = domains[1][1]
    tmax = domains[1][2]
    amin = domains[2][1]
    amax = domains[2][2]
    n = srde.n
    a₀ = srde.refractoryperiod
    α = srde.decayrate
    φ = srde.firingrate

    arange = amin:da:amax
    i₀ = findfirst(arange .> a₀)
    u = copy(initial_condition)
    ξ = 0

    len_ts = convert(Integer, fld(tmax - tmin, saveat * dt)) + 2
    ts = zeros(Float64, len_ts)
    sol = zeros(Float64, (len_ts, length(u)))
    j = 1
    ts[j] = t
    sol[j, :] = u
    elapsedtime_save = 0
    while (t < tmax)
        t += dt
        elapsedtime_save += 1
        ξ *= exp(-α * dt)
        ξ += α * u[1] * dt
        newactivity = u[end]
        for i in length(u):-1:2
            uφdt = u[i-1] * φ(ξ) * (i - 1 >= i₀) * dt
            increment = max(0,
                uφdt + rand(Normal()) * sqrt((max(uφdt, 0.0)) / n)
            )
            u[i] = u[i-1] - increment
            newactivity += increment
        end
        u[1] = newactivity

        if elapsedtime_save == saveat
            j += 1
            ts[j] = t
            sol[j, :] = u
            elapsedtime_save = 0
        end
    end
    return ts, sol
end