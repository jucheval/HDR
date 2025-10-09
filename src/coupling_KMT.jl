include("auxiliary_quantiles.jl")

function KMT_coupling_W2X(time, W)
    # Augment time horizon to next power of 2
    N = length(time)
    N = 2^ceil(Int, log2(N))
    Delta_my = time[1]
    time = collect(Delta_my:Delta_my:Delta_my*N)

    normal = (diff([0; W]) .- Delta_my) ./ sqrt(Delta_my)
    normal = vcat(normal, randn(N - length(W)))

    N1 = N - 1
    K = Int(log2(N))
    K1 = K - 1

    V_mat = zeros(K, N1)
    V_tilde_mat = zeros(K, div(N - 2, 2))
    U_mat = zeros(K, N1)
    U_tilde_mat = zeros(K, div(N - 2, 2))

    T_mat = cumsum(normal)

    for j in 0:K1
        term = div(N, 2^j) - 1
        for k in 1:term
            V_mat[j+1, k] = T_mat[(k+1)*2^j] - T_mat[k*2^j]
        end
    end

    for j in 2:K
        term = div(N, 2^(j - 1)) - 1
        for k in 1:term
            V_tilde_mat[j, k] = V_mat[j-1, 2*k] - V_mat[j-1, 2*k+1]
        end
    end

    U_mat[:, 1] .= quant.(V_mat[:, 1], collect(0:K1), Delta=Delta_my)
    X0 = quant(normal[1], 0, Delta=Delta_my)

    for i in 1:K
        U_tilde_mat[i, 1] = condquant(V_tilde_mat[i, 1], U_mat[i, 1], i - 1, Delta=Delta_my)
    end

    M = K - 2
    counter1 = 1
    counter2 = 1
    for m in 0:M
        K_new = K - m
        for p in 1:counter2
            for a in 2:K_new
                U_mat[a-1, 2*counter1] = 0.5 * (U_mat[a, counter1] + U_tilde_mat[a, counter1])
                U_mat[a-1, 2*counter1+1] = 0.5 * (U_mat[a, counter1] - U_tilde_mat[a, counter1])
                if a - 1 != 1
                    U_tilde_mat[a-1, 2*counter1] = condquant(V_tilde_mat[a-1, 2*counter1], U_mat[a-1, 2*counter1], a - 2, Delta=Delta_my)
                    U_tilde_mat[a-1, 2*counter1+1] = condquant(V_tilde_mat[a-1, 2*counter1+1], U_mat[a-1, 2*counter1+1], a - 2, Delta=Delta_my)
                end
            end
            counter1 += 1
        end
        counter2 *= 2
    end

    Xsummands = vcat(X0, U_mat[1, :]) * sqrt(Delta_my) .+ Delta_my
    Wsummands = normal * sqrt(Delta_my) .+ Delta_my
    W_new = cumsum(Wsummands)
    X = cumsum(Xsummands)

    return (time=time, W=W_new, X=X)
end

function KMT_coupling_X2W(time, X; construction=false)
    # Augment time horizon to next power of 2
    N = length(time)
    N = 2^ceil(Int, log2(N))
    Delta_my = time[1]
    time = collect(Delta_my:Delta_my:Delta_my*N)

    poisson = round.(diff([0; X]))
    poisson = vcat(poisson, rand(Poisson(Delta_my), N - length(X)))
    standard_poisson = (poisson .- Delta_my) ./ sqrt(Delta_my)

    if N <= 2
        Wsummands = quant_norm.(standard_poisson, 0, Delta=Delta_my) * sqrt(Delta_my) .+ Delta_my
        W = cumsum(Wsummands)
    else
        N1 = N - 1
        K = Int(log2(N))
        K1 = K - 1
        V_mat = zeros(K, N1)
        V_tilde_mat = zeros(K, div(N - 2, 2))
        U_mat = zeros(K, N1)
        U_tilde_mat = zeros(K, div(N - 2, 2))

        S_mat = cumsum(standard_poisson)

        for j in 0:K1
            term = div(N, 2^j) - 1
            for k in 1:term
                U_mat[j+1, k] = S_mat[(k+1)*2^j] - S_mat[k*2^j]
            end
        end

        for j in 2:K
            term = div(N, 2^(j - 1)) - 1
            for k in 1:term
                U_tilde_mat[j, k] = U_mat[j-1, 2*k] - U_mat[j-1, 2*k+1]
            end
        end

        V_mat[:, 1] .= quant_norm.(U_mat[:, 1], collect(0:K1), Delta=Delta_my)
        Y0 = quant_norm(standard_poisson[1], 0, Delta=Delta_my)

        for i in 1:K
            V_tilde_mat[i, 1] = condquant_norm(U_tilde_mat[i, 1], U_mat[i, 1], i - 1, Delta=Delta_my)
        end

        M = K - 2
        counter1 = 1
        counter2 = 1
        for m in 0:M
            K_new = K - m
            for p in 1:counter2
                for a in 2:K_new
                    V_mat[a-1, 2*counter1] = 0.5 * (V_mat[a, counter1] + V_tilde_mat[a, counter1])
                    V_mat[a-1, 2*counter1+1] = 0.5 * (V_mat[a, counter1] - V_tilde_mat[a, counter1])
                    if a - 1 != 1
                        V_tilde_mat[a-1, 2*counter1] = condquant_norm(U_tilde_mat[a-1, 2*counter1], U_mat[a-1, 2*counter1], a - 2, Delta=Delta_my)
                        V_tilde_mat[a-1, 2*counter1+1] = condquant_norm(U_tilde_mat[a-1, 2*counter1+1], U_mat[a-1, 2*counter1+1], a - 2, Delta=Delta_my)
                    end
                end
                counter1 += 1
            end
            counter2 *= 2
        end

        Wsummands = vcat(Y0, V_mat[1, :]) * sqrt(Delta_my) .+ Delta_my
        W = cumsum(Wsummands)
    end

    if construction
        return (time=time, X=X, W=W, Wsummands=Wsummands)
    end
    return (time=time, X=X, W=W)
end

function counting2point(time, X)
    # Returns the Poisson point process (vector Z) from the discretized Poisson counting process X
    X = round.(Int, X)
    max_X = X[end]
    Z = zeros(max_X)
    jumps = diff(X)
    pos_jumps_index = findall(x -> x > 0.1, jumps)
    i = 1
    for index in pos_jumps_index
        nb_jumps = round(Int, jumps[index])
        j = i + nb_jumps - 1
        Z[i:j] .= sort(rand(Uniform(time[index], time[index+1]), nb_jumps))
        i = j + 1
    end
    return Z
end

function point2counting(time, Z)
    # Returns the counting process X from the Poisson point process Z
    X = 0
    max_X = length(Z)
    nbtime = length(time)
    Xvec = zeros(Int, nbtime)
    for t in 1:nbtime
        while X < max_X && time[t] > Z[X+1]
            X += 1
        end
        Xvec[t] = X
    end
    return Xvec
end

function KMT_simulation(N, Delta_my)
    N = 2^ceil(Int, log2(N))
    normal = randn(N)
    Wsummands = normal * sqrt(Delta_my) .+ Delta_my
    time = collect(Delta_my:Delta_my:Delta_my*N)
    W = cumsum(Wsummands)
    X = KMT_coupling_W2X(time, W).X
    return (time=time, X=X, W=W)
end