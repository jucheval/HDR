using GLMakie
using PointProcesses
using Distributions
import Random: seed!

# Compute a 2d Hawkes process and its intensities via thinning
function thinning(h_poisson, μ, α, β, dt)
    tmin = min_time(h_poisson)
    tmax = max_time(h_poisson)
    times_poisson = event_times(h_poisson)
    marks_poisson = event_marks(h_poisson)

    h = History(; times=Float64[], marks=Int[], tmin=tmin, tmax=tmax)

    λ = μ
    t = tmin
    Ξ = fill(0, length(μ))
    λs = λ
    ts = [t]
    i = 1
    while t < tmax
        τ = (i < length(times)) ? times_poisson[i] - t : Inf
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

            j, u = marks[i]
            m = Int(j)
            if u < λ[m]
                push!(h, t, m)
                Ξ += β[:, m]
                λ = max.(0, μ .+ Ξ)

                λs = hcat(λs, λ)
                push!(ts, t)
            end

            i += 1
        end
    end
    return ts, λs, h
end

intensity_bound = 5.0
begin # time grid
    tmin = 0.0
    tmax = 5.0
    dt = 0.1
end

begin # simulation of two 2d Poisson for thinning
    seed!(1)
    N = rand(Poisson(float(2 * intensity_bound * (tmax - tmin))))
    times = sort!(rand(Uniform(tmin, tmax), N))
    marks = [rand(product_distribution(DiscreteUniform(1, 2), Uniform(0.0, intensity_bound))) for n in 1:N]
    h_poisson = History(; times=times, marks=marks, tmin=tmin, tmax=tmax)
end

fig = Figure(size=(1400, 500))

begin # parameters 
    sg = SliderGrid(fig[1, 1],
        (label="μ₁", range=0.5:0.01:1.5, startvalue=0.8),
        (label="μ₂", range=0.5:0.01:1.5, startvalue=1.2),
        (label="α", range=1.0:0.01:5.0, startvalue=4.0),
        (label="β₁₁", range=-1.0:0.01:1.0, startvalue=0.5),
        (label="β₁₂", range=-1.0:0.01:1.0, startvalue=-0.2),
        (label="β₂₁", range=-1.0:0.01:1.0, startvalue=-0.5),
        (label="β₂₂", range=-1.0:0.01:1.0, startvalue=0.3),
        tellheight=false)
    obs_parameters = sg.sliders
    μ = @lift [$(obs_parameters[1].value), $(obs_parameters[2].value)]
    α = obs_parameters[3].value
    β = @lift [$(obs_parameters[4].value) $(obs_parameters[6].value);
        $(obs_parameters[5].value) $(obs_parameters[7].value)]
end

begin # Observables
    obs_thinned = @lift thinning(h_poisson, $μ, $α, $β, dt)
    ts = @lift $obs_thinned[1]
    λs = @lift $obs_thinned[2]
    h = @lift $obs_thinned[3]

    spikes1 = @lift event_times($h)[event_marks($h).==1]
    spikes2 = @lift event_times($h)[event_marks($h).==2]
end

begin # Plot
    ax = Axis(fig[1, 2], ylabel=L"\lambda", xlabel=L"t", title="Bivariate Hawkes process")

    lines!(ax, ts, (@lift $λs[1, :]), color=:darkblue, label=L"\lambda^1_t")
    lines!(ax, ts, (@lift $λs[2, :]), color=:orange, label=L"\lambda^2_t")
    hlines!(ax, (@lift $μ[1]), color=:darkblue, linestyle=:dot)
    hlines!(ax, (@lift $μ[2]), color=:orange, linestyle=:dot)
    scatter!(ax, spikes1, (@lift 0.0 * $spikes1), color=:darkblue, markersize=15, label=L"N^1")
    scatter!(ax, spikes2, (@lift 0.0 * $spikes2), color=:orange, markersize=15, label=L"N^2")

    axislegend(ax)
end

begin # Adjust layout
    colsize!(fig.layout, 1, Fixed(300))
end

fig
