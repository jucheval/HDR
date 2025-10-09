include("coupling_KMT.jl")
using ExponentialUtilities

function linear_ODE_matrix(nu_vec, eta_vec)
    dimension = sum(eta_vec)
    A = zeros(dimension, dimension)
    nb_pop = length(nu_vec)
    D = Float64[]
    for k in 1:nb_pop
        append!(D, fill(-nu_vec[k], eta_vec[k]))
    end
    for i in 1:dimension
        A[i, i] = D[i]
    end
    upper = Int[]
    for k in 1:(nb_pop-1)
        append!(upper, fill(1, eta_vec[k] - 1))
        push!(upper, 0)
    end
    append!(upper, fill(1, eta_vec[nb_pop] - 1))
    indx = 1:(dimension-1)
    for (i, val) in enumerate(upper)
        A[indx[i], indx[i]+1] = val
    end
    return A
end

function MF_limit_ODE(nb_pop, intensity_function, c_vec, nu_vec, eta_vec, X_init, Tfinal, time_step)
    t = 0.0
    X = copy(X_init)
    A = linear_ODE_matrix(nu_vec, eta_vec)
    h = time_step
    times = collect(0:h:Tfinal)
    nb_points = length(times)
    intensity_in_time = zeros(nb_pop, nb_points)
    X_in_time = zeros(length(X_init), nb_points)
    lambda = intensity(X, nb_pop, eta_vec, intensity_function)
    for i in 1:nb_points
        t += h
        X = expv(h, A, X)
        X[cumsum(eta_vec)] .+= h .* c_vec[1:nb_pop] .* reverse(lambda[1:nb_pop])
        lambda = intensity(X, nb_pop, eta_vec, intensity_function)
        intensity_in_time[:, i] .= lambda
        X_in_time[:, i] .= X
    end
    return (time=times, intensity=intensity_in_time, X=X_in_time)
end

function intensity(X, nb_pop, eta_vec, intensity_function)
    value = zeros(nb_pop)
    for k in 1:nb_pop
        potential = X[sum(eta_vec[1:(k-1)])+1]
        value[k] = intensity_function[k](potential)
    end
    return value
end

function B_2_Hawkes_Erlang(nb_pop, nb_neuron, intensity_function, c_vec, nu_vec, eta_vec, X_init, B)
    nb_points_B = size(B, 1)
    nb_points = nb_pop * nb_points_B
    t = 0.0
    B_time = B[:, nb_pop+1]
    increments_B = diff(vcat(zeros(1, nb_pop), B[:, 1:nb_pop]), dims=1)
    X = copy(X_init)
    A = linear_ODE_matrix(nu_vec, eta_vec)
    diffusion_in_time = zeros(nb_pop, nb_points)
    intensity_in_time = zeros(nb_pop, nb_points)
    X_in_time = zeros(length(X_init), nb_points)
    time_vec = zeros(nb_points)
    lambda = intensity(X, nb_pop, eta_vec, intensity_function)
    Lambda = zeros(nb_pop)
    current_B = ones(Int, nb_pop)
    i = 1
    while all(current_B .<= nb_points_B) && i <= nb_points
        h_vec = (B_time[current_B] .- Lambda) ./ (nb_neuron .* lambda)
        spiking_pop = argmin(h_vec)
        h = h_vec[spiking_pop]
        t += h
        Lambda .+= h .* nb_neuron .* lambda
        update_B = increments_B[current_B[spiking_pop], spiking_pop]
        influenced_pop = spiking_pop == 1 ? nb_pop : spiking_pop - 1
        X = expv(h, A, X)
        X[sum(eta_vec[1:influenced_pop])] += update_B * c_vec[influenced_pop] / nb_neuron[spiking_pop]
        lambda = intensity(X, nb_pop, eta_vec, intensity_function)
        current_B[spiking_pop] += 1
        index_vecB = current_B .+ nb_points_B .* (0:(nb_pop-1))
        diffusion_in_time[:, i] .= B[index_vecB]
        time_vec[i] = t
        intensity_in_time[:, i] .= lambda
        X_in_time[:, i] .= X
        i += 1
    end
    good_index = 1:(i-1)
    return (
        time=time_vec[good_index],
        diffusion=diffusion_in_time[:, good_index],
        intensity=intensity_in_time[:, good_index],
        X=X_in_time[:, good_index]
    )
end

function Pi_2_Hawkes_Erlang(nb_pop, nb_neuron, intensity_function, c_vec, nu_vec, eta_vec, X_init, Pi, tol=0.01)
    nb_points_vec = [sum(Pi[:, k] .> 0) for k in 1:nb_pop]
    nb_points = sum(nb_points_vec)
    t = 0.0
    X = copy(X_init)
    A = linear_ODE_matrix(nu_vec, eta_vec)
    points = zeros(nb_points)
    type = zeros(Int, nb_points)
    intensity_in_time = zeros(nb_pop, nb_points)
    X_in_time = zeros(length(X_init), nb_points)
    lambda = intensity(X, nb_pop, eta_vec, intensity_function)
    h = minimum(tol ./ (nb_neuron .* lambda))
    Lambda = zeros(nb_pop)
    current_Pi = ones(Int, nb_pop)
    i = 1
    while i <= nb_points
        if any(current_Pi .> nb_points_vec)
            break
        end
        index_vecPi = current_Pi .+ size(Pi, 1) .* (0:(nb_pop-1))
        while all(Lambda .< Pi[index_vecPi])
            if any(Pi[index_vecPi] .- Lambda .< 2 * tol)
                h_vec = (Pi[index_vecPi] .- Lambda) ./ (nb_neuron .* lambda)
                spiking_pop = argmin(h_vec)
                h = h_vec[spiking_pop]
            end
            t += h
            Lambda .+= h .* nb_neuron .* lambda
            X = expv(h, A, X)
            lambda = intensity(X, nb_pop, eta_vec, intensity_function)
        end
        spiking_pop = findfirst(x -> x, Lambda .>= Pi[index_vecPi])
        influenced_pop = spiking_pop == 1 ? nb_pop : spiking_pop - 1
        X[sum(eta_vec[1:influenced_pop])] += c_vec[influenced_pop] / nb_neuron[spiking_pop]
        neuron_number = rand(1:nb_neuron[spiking_pop])
        current_Pi[spiking_pop] += 1
        lambda = intensity(X, nb_pop, eta_vec, intensity_function)
        points[i] = t
        type[i] = sum(nb_neuron[1:(spiking_pop-1)]) + neuron_number
        intensity_in_time[:, i] .= lambda
        X_in_time[:, i] .= X
        i += 1
    end
    good_index = findall(x -> x > 0, type)
    return (
        spike_train=points[good_index],
        type=type[good_index],
        intensity=intensity_in_time[:, good_index],
        X=X_in_time[:, good_index]
    )
end