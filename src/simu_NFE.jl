using Random
using Distributions
using CairoMakie

# Space is 1 dimensional
# The spatial distribution ρ is the Lebesgue measure on the spatial domain interval

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
    scalingfactor = xmax - xmin
    eαdt = exp(-α * dt)

    # spatial grid and initial conditions
    xgrid = xmin:dx:xmax
    wtranspose = scalingfactor * transpose(synapticweightmatrix(xgrid, w))
    nx = length(xgrid)
    u = initial_condition.(xgrid)

    # pre-allocation
    len_ts = convert(Integer, fld(tmax - tmin, saveat * dt)) + 2 + (saveat == 1)
    ts = zeros(Float64, len_ts)
    sol = zeros(Float64, (nx, len_ts))
    noise = zeros(Float64, nx)

    # initialization
    idt = 1
    ts[idt] = t
    sol[:, idt] = u
    elapsedtime_save = 0
    while (t < tmax)
        t += dt
        elapsedtime_save += 1
        if (n_snfe < Inf)
            noise .= rand(Normal(), nx)
        end
        fdtdx = dt * f.(u) * dx / scalingfactor
        u = eαdt * u + wtranspose * (fdtdx + sqrt.(fdtdx / n_snfe) .* noise)

        if elapsedtime_save == saveat
            idt += 1
            ts[idt] = t
            sol[:, idt] = u
            elapsedtime_save = 0
        end
    end
    return ts, xgrid, transpose(sol)
end

function synapticweightmatrix(xs, w)
    [w(xs[i], xs[j]) for i in eachindex(xs), j in eachindex(xs)]
end

function Makie.plot(::Type{StochasticNFE}, simulation)
    ts, xs, sol = simulation
    fig = Figure(size=(350, 300))

    ax = Axis(fig[1, 1], xlabel=L"t", ylabel=L"x", title=L"heatmap of $u(t,x)$")

    hm = heatmap!(ax, ts, xs, sol)
    Colorbar(fig[1, 2], hm)

    fig
end
