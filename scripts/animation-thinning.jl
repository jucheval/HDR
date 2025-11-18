using GLMakie
using PointProcesses
using Distributions
import Random: seed!

intensity_bound = 1.5
begin # time grid
    tmin = 0.0
    tmax = 5.0
end

begin # user-defined Poisson for thinning
    h_poisson = History(;
        times=[0.75, 1.32, 2.01, 2.87, 3.78, 4.54],
        marks=[1.23, 0.43, 0.95, 0.67, 1.18, 0.87],
        tmin=tmin, tmax=tmax
    )
    spikes12 = [1.32, 2.87]
    marks12 = [0.43, 0.67]
    spikes1 = [2.01]
    marks1 = [0.95]
    spikes2 = [4.54]
    marks2 = [0.87]

    # Defines points and vectors to draw arrows
    crosspad = 0.05
    crosses12 = Point2f.(spikes12, marks12 .- crosspad)
    vec12 = Vec2f.(0, -(marks12 .- 2 * crosspad))
    crosses1 = Point2f.(spikes1, marks1 .- crosspad)
    vec1 = Vec2f.(0, -(marks1 .- 2 * crosspad))
    crosses2 = Point2f.(spikes2, marks2 .- crosspad)
    vec2 = Vec2f.(0, -(marks2 .- 2 * crosspad))
end

fig = Figure(size=(700, 400))

begin # parameters 
    sg = SliderGrid(fig[2, 1],
        (label="t", range=tmin:0.01:tmax, startvalue=tmin), tellwidth=false)
    t_animation = sg.sliders[1].value
end

begin # observables
    x_range_1 = @lift [tmin, min($t_animation, tmin + (tmax - tmin) / 2)]
    x_range_2 = @lift [min($t_animation, tmin + (tmax - tmin) / 2), min($t_animation, tmax)]

    spikes12_animated = @lift spikes12[spikes12.<$t_animation]
    spikes1_animated = @lift spikes1[spikes1.<$t_animation]
    spikes2_animated = @lift spikes2[spikes2.<$t_animation]
end

begin # Plot layout parameters
    alpha_intensity = 0.05
    marksize = 13
    ypad = 0.1
end

begin # Plot
    ax = Axis(fig[1, 1], ylabel=L"\lambda", xlabel="t", title="SDE coupling", xticks=[0, 2.5, 5])
    limits!(ax, (nothing, nothing), (0.0 - ypad, intensity_bound + ypad))

    # vertical line to represent t_animation
    vlines!(ax, t_animation, color=:red, alpha=0.5)

    # transparent rectangles for area under the intensity
    poly!(ax, Rect(tmin, 0.0, tmin + (tmax - tmin) / 2, 1.0), color=:darkblue, alpha=alpha_intensity)
    poly!(ax, Rect(tmin + (tmax - tmin) / 2, 0.0, tmin + (tmax - tmin) / 2, 0.8), color=:darkblue, alpha=alpha_intensity)
    poly!(ax, Rect(tmin, 0.0, tmin + (tmax - tmin) / 2, 0.8), color=:orange, alpha=alpha_intensity)
    poly!(ax, Rect(tmin + (tmax - tmin) / 2, 0.0, tmin + (tmax - tmin) / 2, 1.0), color=:orange, alpha=alpha_intensity)

    # transparent lines for future intensity
    lines!(ax, [tmin, tmin + (tmax - tmin) / 2], [1.0, 1.0], color=:darkblue, alpha=7 * alpha_intensity)
    lines!(ax, [tmin + (tmax - tmin) / 2, tmax], [0.8, 0.8], color=:darkblue, alpha=7 * alpha_intensity)
    lines!(ax, [tmin, tmin + (tmax - tmin) / 2], [0.8, 0.8], color=:orange, alpha=7 * alpha_intensity)
    lines!(ax, [tmin + (tmax - tmin) / 2, tmax], [1.0, 1.0], color=:orange, alpha=7 * alpha_intensity)

    # crosses for the 2d Poisson
    scatter!(ax, h_poisson.times, h_poisson.marks, color=:black, marker=:xcross, markersize=marksize, label=L"\Pi")

    # solid lines for past intensity
    lines!(ax, x_range_1, [1.0, 1.0], color=:darkblue, label=L"\lambda^1(t)")
    lines!(ax, x_range_2, [0.8, 0.8], color=:darkblue)
    lines!(ax, x_range_1, [0.8, 0.8], color=:orange, label=L"\lambda^2(t)")
    lines!(ax, x_range_2, [1.0, 1.0], color=:orange)

    # dots for the projection of thinned points
    scatter!(ax, spikes12_animated, (@lift 0.0 * $spikes12_animated), color=:brown, markersize=marksize)
    scatter!(ax, spikes1_animated, (@lift 0.0 * $spikes1_animated), color=:darkblue, markersize=marksize, label=L"N^1")
    scatter!(ax, spikes2_animated, (@lift 0.0 * $spikes2_animated), color=:orange, markersize=marksize, label=L"N^2")

    # projection arrows
    arrows2d!(ax, crosses12[1], vec12[1], minshaftlength=0, color=:brown, visible=@lift spikes12[1] < $t_animation)
    arrows2d!(ax, crosses12[2], vec12[2], minshaftlength=0, color=:brown, visible=@lift spikes12[2] < $t_animation)
    arrows2d!(ax, crosses1, vec1, minshaftlength=0, color=:darkblue, visible=@lift spikes1[1] < $t_animation)
    arrows2d!(ax, crosses2, vec2, minshaftlength=0, color=:orange, visible=@lift spikes2[1] < $t_animation)

    axislegend(ax, position=:ct, orientation=:horizontal)
end

begin # Adjust layout

end

fig
