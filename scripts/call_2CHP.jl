include("../src/simu_2CHP.jl")
using CairoMakie

begin # parameters
    nb_pop = 2
    nb_neuron = fill(10, nb_pop)

    f1(x) = x < log(20) ? 10 * exp(x) : 400 / (1 + 400 * exp(-2 * x))
    f2(x) = x < log(20) ? exp(x) : 40 / (1 + 400 * exp(-2 * x))
    intensity_function = [f1, f2]

    c_vec = [-1, 1]
    nu_vec = fill(1, nb_pop)
    eta_vec = [3, 2] .+ 1  # eta+1 as in the paper
    X_init = vcat(fill(-2, eta_vec[1]), fill(2, eta_vec[2]))
end

begin # time grid
    len_ts = 2^13
    dt = 0.5
end

begin # simulation
    Random.seed!(1)
    # pre-allocation
    Pi_list = Vector{Vector{Float64}}(undef, nb_pop)
    B_list = Vector{Vector{Float64}}(undef, nb_pop)
    len_Pi = zeros(Int, nb_pop)
    len_B = zeros(Int, nb_pop)

    multiplyfactor = [1, 1]

    for k in 1:nb_pop
        coupling_KMT = KMT_simulation(len_ts * multiplyfactor[k], dt)
        Z = counting2point(coupling_KMT.time, coupling_KMT.X)
        len_Pi[k] = length(Z)
        Pi_list[k] = Z
        B = coupling_KMT.W
        len_B[k] = length(B)
        B_list[k] = B
    end

    Pi = zeros(maximum(len_Pi), nb_pop)
    B = zeros(maximum(len_B), nb_pop + 1)
    time_vec = collect(dt:dt:dt*maximum(len_B))
    B[:, nb_pop+1] = time_vec

    for k in 1:nb_pop
        Pi[1:len_Pi[k], k] = Pi_list[k]
        B[1:len_B[k], k] = B_list[k]
    end

    B_test = zeros(maximum(len_B), nb_pop + 1)
    for k in 1:nb_pop
        B_test[1:len_B[k], k] = cumsum(randn(len_B[k]) * sqrt(dt)) .+ time_vec[1:len_B[k]]
    end
    B_test[:, nb_pop+1] = time_vec

    coupled_counting = Pi_2_Hawkes_Erlang(nb_pop, nb_neuron, intensity_function, c_vec, nu_vec, eta_vec, X_init, Pi)
    coupled_diffusion = B_2_Hawkes_Erlang(nb_pop, nb_neuron, intensity_function, c_vec, nu_vec, eta_vec, X_init, B)
    test_diffusion = B_2_Hawkes_Erlang(nb_pop, nb_neuron, intensity_function, c_vec, nu_vec, eta_vec, X_init, B_test)

    Tfinal = 150
    dt = 0.1
    ODE_traj = MF_limit_ODE(nb_pop, intensity_function, c_vec, nu_vec, eta_vec, X_init, Tfinal, dt)
end

begin # plot
    fig = Figure(size=(800, 480))

    axtop = Axis(fig[1, 1], ylabel=L"U^1")
    axbot = Axis(fig[2, 1], ylabel=L"U^2", xlabel=L"t")
    linkxaxes!(axtop, axbot)
    hidexdecorations!(axtop, grid=false)
    xlims!(axbot, 0, 140)

    kwargs = (; alpha=0.7)
    lines!(axtop, coupled_diffusion.time, coupled_diffusion.X[1, :]; color=:orange, kwargs...)
    lines!(axtop, coupled_counting.spike_train[begin:20:end], coupled_counting.X[1, begin:20:end]; color=:darkblue, kwargs...)
    lines!(axtop, test_diffusion.time, test_diffusion.X[1, :]; color=:orange, linestyle=:dot, kwargs...)
    lines!(axtop, ODE_traj.time, ODE_traj.X[1, :]; color=:black, kwargs...)

    lines!(axbot, coupled_diffusion.time, coupled_diffusion.X[eta_vec[1]+1, :]; color=:orange, kwargs...)
    lines!(axbot, coupled_counting.spike_train[begin:20:end], coupled_counting.X[eta_vec[1]+1, begin:20:end]; color=:darkblue, kwargs...)
    lines!(axbot, test_diffusion.time, test_diffusion.X[eta_vec[1]+1, :]; color=:orange, linestyle=:dot, kwargs...)
    lines!(axbot, ODE_traj.time, ODE_traj.X[eta_vec[1]+1, :]; color=:black, kwargs...)

    fig
end

# save figure
save("plots/coupling-2CHP.png", fig, px_per_unit=2)