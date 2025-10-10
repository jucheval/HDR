using PointProcesses
using Distributions
using CairoMakie

struct AgeDependentHawkesProcess
    n::Int
    refractoryperiod::Float64
    decayrate::Float64
    firingrate::Function
    firingratebound::Float64
end

AgeDependentHawkesProcess(; n, refractoryperiod, decayrate, firingrate, firingratebound) = AgeDependentHawkesProcess(n, refractoryperiod, decayrate, firingrate, firingratebound)
getfields(adhp::AgeDependentHawkesProcess) = adhp.n, adhp.refractoryperiod, adhp.decayrate, adhp.firingrate, adhp.firingratebound

include("simu_Poisson.jl") # function to simulate an homogeneous Poisson process

function simulate(adhp::AgeDependentHawkesProcess, initial_condition::Vector{Float64}, domain::Vector{Float64})
    # notations
    Nneur, a₀, α, φ, φ̅ = getfields(adhp)
    tmin = domain[1]
    tmax = domain[2]

    # initial conditions
    ages = copy(initial_condition)
    ξ = 0

    # simulation of a dominating Poisson Process with notations
    dominatingPoisson = rand(MultivariatePoissonProcess(fill(φ̅, Nneur)), tmin, tmax)
    npts = length(dominatingPoisson)
    times = event_times(dominatingPoisson)
    marks = event_marks(dominatingPoisson)

    # pre-allocation
    ξs = zeros(npts)
    activities = zeros(npts)
    for k in 1:npts
        time_step = k == 1 ? times[1] : times[k] - times[k-1]
        ξ *= exp(-α * time_step)
        ages .+= time_step
        pt_mark = marks[k]
        intensity = φ(ξ)
        activities[k] = intensity * sum(ages .> a₀) / Nneur
        if ((ages[pt_mark] <= a₀) || (rand() > intensity / φ̅))
            marks[k] = 0
        else
            ages[pt_mark] = 0
            ξ = ξ + α / Nneur
        end
        ξs[k] = ξ
    end

    h = History(times, marks, tmin, tmax)
    return h, ξs, activities
end

function Makie.plot(::Type{AgeDependentHawkesProcess}, simulation)
    h, _, activities = simulation
    fig = Figure(size=(600, 300))

    axscatt = Axis(fig[1, 1], ylabel="neuron",)
    axlines = Axis(fig[1, 1], yticklabelcolor=:darkblue, yaxisposition=:right, ylabel="mean firing rate", ylabelcolor=:darkblue)
    hidespines!(axlines)
    hidexdecorations!(axlines)
    hideydecorations!(axscatt, label=false, ticklabels=false)

    ids = findall(event_marks(h) .!= 0)
    scatter!(axscatt, event_times(h)[ids], event_marks(h)[ids], markersize=4, color=:black)

    ids = 1:floor(Int, length(h) / length_ts):length(h)
    ts = event_times(h)[ids]
    ys = activities[ids]
    lines!(axlines, ts, ys, color=:darkblue)

    fig
end