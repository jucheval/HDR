include("../src/simu_2CHP.jl")
using Random
using DrWatson

begin # parameters
    nb_pop = 2
    nb_neuron = fill(450, nb_pop)

    f1(x) = x < log(20) ? 10 * exp(x) : 400 / (1 + 400 * exp(-2 * x))
    f2(x) = x < log(20) ? exp(x) : 40 / (1 + 400 * exp(-2 * x))
    intensity_function = [f1, f2]

    c_vec = [-1, 1]
    nu_vec = fill(1, nb_pop)
    eta_vec = [3, 2] .+ 1  # eta+1 as in the paper
    X_init = vcat(fill(-2, eta_vec[1]), fill(2, eta_vec[2]))
end

begin # time grid
    len_ts = 2^17
    dt = 0.5
end

begin # simulation
    Random.seed!(1)
    # pre-allocation
    Pi_list = Vector{Vector{Float64}}(undef, nb_pop)
    len_Pi = zeros(Int, nb_pop)

    multiplyfactor = [1, 1]

    for k in 1:nb_pop
        coupling_KMT = KMT_simulation(len_ts * multiplyfactor[k], dt)
        Z = counting2point(coupling_KMT.time, coupling_KMT.X)
        len_Pi[k] = length(Z)
        Pi_list[k] = Z
    end

    Pi = zeros(maximum(len_Pi), nb_pop)

    for k in 1:nb_pop
        Pi[1:len_Pi[k], k] = Pi_list[k]
    end

    cc = Pi_2_Hawkes_Erlang(nb_pop, nb_neuron, intensity_function, c_vec, nu_vec, eta_vec, X_init, Pi)
    h = History(cc.spike_train, cc.type, 0.0, maximum(cc.spike_train))
end

safesave("data/2CHP-animation.jld2", Dict("history" => h))