using Random
using Distributions

struct TimeChangeBrownianMotion
    timescaling::Function
end

TimeChangeBrownianMotion(; timescaling) = TimeChangeBrownianMotion(timescaling)
getfields(tcbm::TimeChangeBrownianMotion) = tcbm.timescaling

function coupling(tcbm::TimeChangeBrownianMotion, tcbm_tilde::TimeChangeBrownianMotion, initial_condition::Float64, domain::Vector{Float64}, dt::Float64)
    # notations
    τ = getfields(tcbm)
    τ̃ = getfields(tcbm_tilde)
    sdt = sqrt(dt)

    # time initialization
    tmin, tmax = domain
    t = t̃ = tmin
    B = initial_condition

    ts = [t]
    t̃s = [t̃]
    Bs = [B]
    while min(t, t̃) < tmax
        t += dt / τ(B)
        t̃ += dt / τ̃(B)
        push!(ts, t)
        push!(t̃s, t̃)
        B += sdt * rand(Normal())
        push!(Bs, B)
    end

    return ts, t̃s, Bs
end