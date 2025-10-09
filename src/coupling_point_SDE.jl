using Random
using Distributions
using PointProcesses

struct InhomogeneousPoissonProcess
    intensity::Function
    intensitybound::Float64
end

InhomogeneousPoissonProcess(; intensity, intensitybound) = InhomogeneousPoissonProcess(intensity, intensitybound)
getfields(ipp::InhomogeneousPoissonProcess) = ipp.intensity, ipp.intensitybound

include("simu_Poisson.jl") # function to simulate an homogeneous Poisson process

function coupling(ipp::InhomogeneousPoissonProcess, ipp_tilde::InhomogeneousPoissonProcess, domain::Vector{Float64})
    # notations
    λ, λbound = getfields(ipp)
    λ̃, λ̃bound = getfields(ipp_tilde)
    bound = max(λbound, λ̃bound)
    tmin, tmax = domain

    # simulation of a dominating Poisson Process with notations
    dominatingPoisson = rand(MultivariatePoissonProcess([bound]), tmin, tmax)
    times = event_times(dominatingPoisson)

    N = Float64[]
    Ñ = Float64[]
    for t in times
        z = rand() * bound
        z < λ(t) && push!(N, t)
        z < λ̃(t) && push!(Ñ, t)
    end

    return N, Ñ
end