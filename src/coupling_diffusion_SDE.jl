using Random
using Distributions

struct Diffusion
    drift::Function
    standarddeviation::Function
end

Diffusion(; drift, standarddeviation) = Diffusion(drift, standarddeviation)
getfields(diffusion::Diffusion) = diffusion.drift, diffusion.standarddeviation

function coupling(diffusion::Diffusion, diffusion_tilde::Diffusion, initial_conditions::Vector{Float64}, domain::Vector{Float64}, dt::Float64)
    # notations
    v, σ = getfields(diffusion)
    ṽ, σ̃ = getfields(diffusion_tilde)

    # time grid
    tmin, tmax = domain
    ts = tmin:dt:tmax
    len_ts = length(ts)

    # pre-allocation
    Xs = zeros(Float64, len_ts)
    X̃s = zeros(Float64, len_ts)
    dBs = sqrt(dt) * rand(Normal(), len_ts - 1) # increments of the coupling Brownian motion

    X = Xs[1] = initial_conditions[1]
    X̃ = X̃s[1] = initial_conditions[2]
    for idt in eachindex(dBs)
        X += v(X) * dt + σ(X) * dBs[idt]
        Xs[idt+1] = X
        X̃ += ṽ(X̃) * dt + σ̃(X̃) * dBs[idt]
        X̃s[idt+1] = X̃
    end

    return ts, Xs, X̃s
end