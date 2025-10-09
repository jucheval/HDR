using Random
using Distributions
using PointProcesses

struct InhomogeneousTimeChange
    inversecumulativeintensity::Function
end

InhomogeneousTimeChange(; inversecumulativeintensity) = InhomogeneousTimeChange(inversecumulativeintensity)
getfields(itc::InhomogeneousTimeChange) = itc.inversecumulativeintensity

include("simu_Poisson.jl") # function to simulate an homogeneous Poisson process

function coupling(itc::InhomogeneousTimeChange, itc_tilde::InhomogeneousTimeChange, domain::Vector{Float64})
    # notations
    Λ⁻ = getfields(itc)
    Λ⁻̃ = getfields(itc_tilde)
    tmin, tmax = domain

    # simulation of a common Poisson Process with notations
    commonPoisson = rand(MultivariatePoissonProcess([1.0]), tmin, tmax)

    N = time_change(commonPoisson, Λ⁻)
    Ñ = time_change(commonPoisson, Λ⁻̃)

    return event_times(N), event_times(Ñ)
end