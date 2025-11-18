using GLMakie
using Distributions

fig = Figure(size=(700, 400))

begin # buttons
    buttongrid = fig[1, 1] = GridLayout(tellwidth=false)
    buttonlabels = ["Independent", "Translation"]
    buttons = buttongrid[1, 1:2] = [Makie.Button(fig, label=buttonlabels[i]) for i in 1:2]
end

begin # observables
    X1 = Observable(rand(Normal(), 4) .+ 0.1)
    X2 = Observable(rand(Normal(), 4) .- 0.1)
end

begin # Plot layout parameters
    marksize = 13
    ypad = 0.1
end

begin # Plot
    ax = Axis(fig[2, 1], ylabel="X", xlabel="realization")
    limits!(ax, (nothing, nothing), (-3.0 - ypad, 3.0 + ypad))

    scatter!(ax, 1:4, X1, label=L"X_1 \sim \mathcal{N}(-0.1,1)", markersize=marksize, color=:darkblue)
    scatter!(ax, 1:4, X2, label=L"X_2 \sim \mathcal{N}(+0.1,1)", markersize=marksize, color=:orange)

    axislegend(ax, position=:ct)
end

# Connect buttons to actions
on(buttons[1].clicks) do _
    X1[] = rand(Normal(), 4) .+ 0.1
    X2[] = rand(Normal(), 4) .- 0.1
end
on(buttons[2].clicks) do _
    X1[] = rand(Normal(), 4) .+ 0.1
    X2[] = X1[] .- 0.2
end

begin # Adjust layout

end

fig
