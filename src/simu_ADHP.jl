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

function simulate(adhp::AgeDependentHawkesProcess, initial_condition::Vector{Float64}, domain::Vector{Float64})
    # notations
    Nneur, a₀, α, φ, φ̅ = getfields(adhp)
    tmin = domain[1]
    tmax = domain[2]

    # initial conditions
    ages = initial_condition
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

# Needed here because the current version of PointProcesses.jl on the general repository outputs an unsorted event_times vector.
function Base.rand(pp::MultivariatePoissonProcess, tmin::Float64, tmax::Float64
)
    mark_dist = mark_distribution(pp)
    N = rand(Poisson(float(ground_intensity(pp) * (tmax - tmin))))
    times = sort!(rand(Uniform(tmin, tmax), N))
    marks = [rand(mark_dist) for n in 1:N]
    return History(; times=times, marks=marks, tmin=tmin, tmax=tmax)
end