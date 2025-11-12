include("../src/simu_ADHP.jl")
using Random
using DrWatson

begin # 900 neurons which exhibits synchronization
    begin # parameters
        φ(ξ) = 10.0 * ξ^2 / (ξ^2 + 1.0) + 0.5
        φ̅ = 10.5
        adhp = AgeDependentHawkesProcess(
            n=900,
            refractoryperiod=1.0,
            decayrate=10.0,
            firingrate=φ,
            firingratebound=φ̅
        )
        initial_condition = rand(Exponential(), adhp.n) .+ 1
    end
    begin # time domain
        domain = [0., 10.]
        length_ts = 1000
    end
end

# simulation
Random.seed!(1)
simulation = simulate(adhp, initial_condition, domain)
h = simulation[1]
ids = findall(event_marks(h) .!= 0)
true_h = History(event_times(h)[ids], event_marks(h)[ids], min_time(h), max_time(h))

safesave("data/ADHP-animation.jld2", Dict("history" => true_h))