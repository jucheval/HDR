using PointProcesses
using Distributions
using Random

struct MultiHawkesProcess <: AbstractPointProcess{Int}
    μ::Vector{Float64}
    α::Float64
    β::Matrix{Float64}
end

getfields(pp::MultiHawkesProcess) = pp.μ, pp.α, pp.β

function Base.rand(pp::MultiHawkesProcess, tmin, tmax, dt)
    μ, α, β = getfields(pp)
    h = History(; times=Float64[], marks=Int[], tmin=tmin, tmax=tmax)
    λ = μ
    t = tmin
    Ξ = fill(0, length(μ))
    λs = λ
    ts = [t]
    while t < tmax
        bound = sum(max.(μ, λ))
        τ = bound > 0 ? rand(Exponential(inv(bound))) : typemax(inv(bound))
        if τ > dt
            Ξ *= exp(-α * dt)
            λ = max.(0, μ .+ Ξ)
            λs = hcat(λs, λ)
            t += dt
            push!(ts, t)
        else
            Ξ *= exp(-α * τ)
            λ = max.(0, μ .+ Ξ)
            λs = hcat(λs, λ)
            t += τ
            push!(ts, t)
            U_max = sum(max.(0, λ)) / bound
            U = rand(typeof(U_max))
            if (U < U_max) && (t < tmax)
                m = rand(Distributions.Categorical(λ / sum(λ)))
                push!(h, t, m)
                Ξ += β[:, m]
                λ = max.(0, μ .+ Ξ)
                λs = hcat(λs, λ)
                push!(ts, t)
            end
        end
    end
    return ts, λs, h
end