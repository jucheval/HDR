using GLMakie
using PointProcesses
using Distributions
import Random: seed!

begin # time grid
    const tmin = 0.0
    const tmax = 5.0
    ts = tmin:0.01:tmax
end

begin # user-defined Poisson for thinning
    h_poisson = History(;
        times=[0.75, 2.01, 2.87, 3.78, 4.94],
        marks=fill(nothing, 5),
        tmin=tmin, tmax=tmax
    )
end

fig = Figure(size=(700, 400))

begin # parameters 
    sg = SliderGrid(fig[2, 1],
        (label="t", range=tmin:0.01:tmax, startvalue=tmin), tellwidth=false)
    t_animation = sg.sliders[1].value

    Λ1(t) = (t - tmin) * (t <= tmin + (tmax - tmin) / 2) + ((tmax - tmin) / 2 + 0.8 * (t - (tmax - tmin) / 2)) * (tmin + (tmax - tmin) / 2 < t)
    Λ1_inv(τ) = τ * (τ <= (tmax - tmin) / 2) + ((tmax - tmin) / 2 + (τ - (tmax - tmin) / 2) / 0.8) * ((tmax - tmin) / 2 < τ)
    Λ2(t) = 0.8 * (t - tmin) * (0 < t <= tmin + (tmax - tmin) / 2) + (0.8 * (tmax - tmin) / 2 + (t - (tmax - tmin) / 2)) * (tmin + (tmax - tmin) / 2 < t <= tmax)
    Λ2_inv(τ) = τ / 0.8 * (τ <= 0.8 * (tmax - tmin) / 2) + ((tmax - tmin) / 2 + (τ - 0.8 * (tmax - tmin) / 2)) * (0.8 * (tmax - tmin) / 2 < τ)
end

begin # observables
    t_range = @lift tmin:0.01:$t_animation

    visibles1 = @lift h_poisson.times .< Λ1($t_animation)
    visibles2 = @lift h_poisson.times .< Λ2($t_animation)

    spikes1_animated = @lift Λ1_inv.(h_poisson.times[$visibles1])
    spikes2_animated = @lift Λ2_inv.(h_poisson.times[$visibles2])
end

begin # Plot layout parameters
    alpha_intensity = 0.3
    marksize = 13
    arrowspad = 0.2
end

begin # Plot
    ax = Axis(fig[1, 1], ylabel="τ", xlabel="t", title="Time change coupling", xticks=[0, 2.5, 5])

    # vertical line to represent t_animation
    vlines!(ax, t_animation, color=:red, alpha=0.5)

    # transparent lines for future cumulative intensity
    lines!(ax, ts, Λ1.(ts), color=:darkblue, alpha=alpha_intensity)
    lines!(ax, ts, Λ2.(ts), color=:orange, alpha=alpha_intensity)

    # dots for the 1d Poisson on y axis
    scatter!(ax, fill(0.0, length(h_poisson)), h_poisson.times, color=:black, markersize=marksize, label=L"\Pi")

    # solid lines for past cumulative intensity
    lines!(ax, t_range, (@lift Λ1.($t_range)), color=:darkblue, label=L"\lambda^1(t)")
    lines!(ax, t_range, (@lift Λ2.($t_range)), color=:orange, label=L"\lambda^2(t)")

    # dots for the time changed points
    scatter!(ax, spikes1_animated, (@lift 0.0 * $spikes1_animated), color=:darkblue, markersize=marksize, label=L"N^1")
    scatter!(ax, spikes2_animated, (@lift 0.0 * $spikes2_animated), color=:orange, markersize=marksize, label=L"N^2")

    # projection arrows
    for i in 1:4
        lines!(ax, [0.0, Λ1_inv(h_poisson.times[i])], fill(h_poisson.times[i], 2), color=:darkblue, alpha=0.5, visible=@lift $visibles1[i])
        arrows2d!(ax, Point2f(Λ1_inv(h_poisson.times[i]), h_poisson.times[i]), Vec2f(0, -h_poisson.times[i] + arrowspad), minshaftlength=0, color=:darkblue, alpha=0.5, visible=@lift $visibles1[i])
        lines!(ax, [0.0, Λ2_inv(h_poisson.times[i])], fill(h_poisson.times[i], 2), color=:orange, alpha=0.5, visible=@lift $visibles2[i])
        arrows2d!(ax, Point2f(Λ2_inv(h_poisson.times[i]), h_poisson.times[i]), Vec2f(0, -h_poisson.times[i] + arrowspad), minshaftlength=0, color=:orange, alpha=0.5, visible=@lift $visibles2[i])
    end


    axislegend(ax, position=:ct, orientation=:horizontal)
end

begin # Adjust layout

end

fig