using Random
using Distributions
using CairoMakie

struct StochasticRDE
    n::Float64
    refractoryperiod::Float64
    decayrate::Float64
    firingrate::Function
end

StochasticRDE(; n, refractoryperiod, decayrate, firingrate) = StochasticRDE(n, refractoryperiod, decayrate, firingrate)
getfields(srde::StochasticRDE) = srde.n, srde.refractoryperiod, srde.decayrate, srde.firingrate

function simulate(srde::StochasticRDE, initial_condition::Function, domains::Vector{Vector{Float64}}, dt::Float64; saveat::Int=1)
    # notations
    n_srde, a₀, α, φ = getfields(srde)
    da = dt
    t = tmin = domains[1][1]
    tmax = domains[1][2]
    amin = domains[2][1]
    amax = domains[2][2]

    # age grid and initial conditions
    arange = amin:da:amax
    i₀ = findfirst(arange .> a₀)
    u = initial_condition.(arange)
    u = u / sum(dt * u)
    ξ = 0

    # pre-allocation
    len_ts = convert(Integer, fld(tmax - tmin, saveat * dt)) + 2 + (saveat == 1)
    ts = zeros(Float64, len_ts)
    sol = zeros(Float64, (len_ts, length(u)))

    # initialization
    idt = 1
    ts[idt] = t
    sol[idt, :] = u
    elapsedtime_save = 0
    while (t < tmax)
        t += dt
        elapsedtime_save += 1
        ξ *= exp(-α * dt)
        ξ += α * u[1] * dt
        newactivity = u[end]
        for i in length(u):-1:2
            uφdt = u[i-1] * φ(ξ) * (i - 1 >= i₀) * dt
            increment = max(0, uφdt + rand(Normal()) * sqrt((max(uφdt, 0.0)) / n_srde))
            u[i] = u[i-1] - increment
            newactivity += increment
        end
        u[1] = newactivity

        if elapsedtime_save == saveat
            idt += 1
            ts[idt] = t
            sol[idt, :] = u
            elapsedtime_save = 0
        end
    end
    return ts, arange, sol
end

function Makie.plot(::Type{StochasticRDE}, simulation)
    ts, as, sol = simulation
    fig = Figure(size=(600, 300))

    axleft = Axis(fig[1, 1], xlabel=L"t", ylabel=L"a", title=L"heatmap of $u(t,a)$")
    axright = Axis(fig[1, 3], xlabel=L"t", ylabel=L"u(t,0)")

    hm = heatmap!(axleft, ts, as, sol, colorscale=z -> log(z + 1 / 2))
    Colorbar(fig[1, 2], hm)

    lines!(axright, ts, sol[:, 1], color=:darkblue)

    fig
end