using Distributions
using Random

# This code is inspired from the following article :
# A review of the deterministic and diffusion approximations for stochastic chemical reaction networks
# P. Mozgunov, M. Beccuti, A. Horvath, T. Jaki, R. Sirovich, and E. Bibbona
# Reaction Kinetics, Mechanisms and Catalysis, 123(2):289–312, 2018

# The functions defined below correspond to the functions defined on top of page 300

# Gj function
function Gj(t, jei, Delta, lambda)
    mu = 2.0^jei * lambda * Delta
    k = quantile(Poisson(mu), t)
    y = (k - lambda * 2.0^jei * Delta) / sqrt(lambda * Delta)
    return y
end

# Fj function
function Fj(x, jei, Delta, lambda)
    mu = 2.0^jei * lambda * Delta
    y = round(Int, sqrt(lambda * Delta) * x + lambda * 2.0^jei * Delta)
    borne_t = [cdf(Poisson(mu), y - 1), cdf(Poisson(mu), y)]
    t = rand(Uniform(borne_t[1], borne_t[2]))
    return t
end

# qcUtildedatoU_R function
function qcUtildedatoU_R(p, j, par)
    r = -j
    s = 2.0^(-j)
    t = -j * log(2.0)
    n = j
    d = 1
    while p > s
        t += log(n) - log(d)
        n -= 1
        d += 1
        s += exp(t)
        r += 2
    end
    return r
end

# pcUtildedatoU_R function
function pcUtildedatoU_R(t, j, par)
    if t < -j
        @warn "conditional probability less than 0"
    end
    if t > j
        @warn "conditional probability larger than 1"
    end
    succ_i = collect(-j:2:t)
    p_up = min(1, sum(pdf(Poisson(par), (succ_i .+ j) ./ 2.0) .* pdf(Poisson(par), (j .- succ_i) ./ 2.0)) / pdf(Poisson(2.0 * par), j))
    i_max = succ_i[end]
    p_down = max(0, p_up - pdf(Poisson(par), (i_max + j) / 2.0) * pdf(Poisson(par), (j - i_max) / 2.0) / pdf(Poisson(2.0 * par), j))
    return (p_down, p_up)
end

# Gndatoy function
function Gndatoy(t, y; Delta=1.0, lambda=1.0, prm)
    p = t
    if p <= 0
        return 0
    elseif p > 1
        return Inf
    end
    j = round(Int, sqrt(lambda * Delta) * y + lambda * 2.0^prm * Delta)
    par = lambda * 2.0^(prm - 1) * Delta
    z = qcUtildedatoU_R(p, j, par)
    return z / sqrt(lambda * Delta)
end

# Fndatoy function
function Fndatoy(x, y; Delta=1.0, lambda=1.0, prm=0)
    t = round(Int, sqrt(lambda * Delta) * x)
    j = round(Int, sqrt(lambda * Delta) * y + lambda * 2.0^prm * Delta)
    par = lambda * 2.0^(prm - 1) * Delta
    z = pcUtildedatoU_R(t, j, par)
    z_rand = rand(Uniform(z[1], z[2]))
    return z_rand
end

# quant function : G_j \circ \Phi
function quant(x, prm; Delta=1.0, lambda=1.0)
    y = cdf(Normal(), 2.0^(-prm / 2.0) * x)
    y = Gj(y, prm, Delta, lambda)
    return y
end

# quant.norm function
function quant_norm(x, prm; Delta=1.0, lambda=1.0)
    y = Fj(x, prm, Delta, lambda)
    y = 2.0^(prm / 2.0) * quantile(Normal(), y)
    return y
end

# condquant function : condition G_j \circ \Phi
function condquant(x, y, prm; Delta=1.0, lambda=1.0)
    z = cdf(Normal(), 2.0^(-prm / 2.0) * x)
    z = Gndatoy(z, y, Delta=Delta, lambda=lambda, prm=prm)
    return z
end

# condquant.norm function
function condquant_norm(x, y, prm; Delta=1.0, lambda=1.0)
    z = Fndatoy(x, y, Delta=Delta, lambda=lambda, prm=prm)
    z = 2.0^(prm / 2.0) * quantile(Normal(), z)
    return z
end