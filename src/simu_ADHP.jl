using PointProcesses
using Distributions

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