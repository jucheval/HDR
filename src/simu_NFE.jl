using Random
using Distributions

struct StochasticNFE
    n::Float64
    decayrate::Float64
    synapticweight::Function
    firingrate::Function
end

StochasticNFE(; n, decayrate, synapticweight, firingrate) = StochasticNFE(n, decayrate, synapticweight, firingrate)
getfields(snfe::StochasticNFE) = snfe.n, snfe.decayrate, snfe.synapticweight, snfe.firingrate

function simulate(snfe::StochasticNFE, initial_condition::Function, domains::Vector{Vector{Float64}}, dt::Float64, dx::Float64; saveat::Int=1)
    # notations
    n_snfe, α, w, f = getfields(snfe)
    t = tmin = domains[1][1]
    tmax = domains[1][2]
    xmin = domains[2][1]
    xmax = domains[2][2]

    # spatial grid and initial conditions
    xgrid = xmin:dx:xmax
    nx = length(xgrid)
    u = initial_condition.(xgrid)
    u_new = copy(u)

    # pre-allocation
    len_ts = convert(Integer, fld(tmax - tmin, saveat * dt)) + 2 + (saveat == 1)
    ts = zeros(Float64, len_ts)
    sol = zeros(Float64, (len_ts, nx))
    noise = zeros(Float64, nx)

    # initialization
    idt = 1
    ts[idt] = t
    sol[idt, :] = u
    elapsedtime_save = 0
    while (t < tmax)
        t += dt
        elapsedtime_save += 1
        if (n_snfe < Inf)
            noise .= rand(Normal(), nx)
        end
        for (i, x) in enumerate(xgrid)
            increment = 0
            for (j, y) in enumerate(xgrid)
                fdtdx = f(u[j]) * dt * dx
                increment += w(y, x) * (fdtdx + sqrt(fdtdx) * noise[j] / sqrt(n_snfe))
            end
            u_new[i] *= (1 - α * dt)
            u_new[i] += increment
        end
        u .= u_new

        if elapsedtime_save == saveat
            idt += 1
            ts[idt] = t
            sol[idt, :] = u
            elapsedtime_save = 0
        end
    end
    return ts, xgrid, sol
end
