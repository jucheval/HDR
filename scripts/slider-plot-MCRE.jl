using MeanFieldGraph
using GLMakie

fig = Figure(size=(1300, 500))

begin # parameters
    sg = SliderGrid(fig[1, 1],
        (label="n", range=50:10:100, startvalue=100),
        (label="μ", range=0.0:0.05:1.0, startvalue=0.25),
        (label="γ", range=0.0:0.05:1.0, startvalue=0.5),
        (label="p", range=0.0:0.05:1.0, startvalue=0.5),
        (label="r₊", range=0.0:0.05:1.0, startvalue=0.5),
        tellheight=false)
    n, μ, γ, p, r₊ = [sl.value for sl in sg.sliders]

    warning_text = @lift $μ + $γ > 1 ? "Warning: μ + γ must be less than 1!" : " "
    Label(fig[2, 1], warning_text, tellwidth=false, tellheight=false, color=:red, fontsize=20)
end

begin # observables
    model = @lift MarkovChainModel($μ, 1.0 - $γ, $p)
    excitatory = @lift MeanFieldGraph.N2excitatory($n, $r₊)
    simulation = @lift rand($model, $excitatory, 100)

    yticks = @lift [1, floor(Int, $n * $r₊), $n]
end

begin # axis and plot
    ax = Axis(fig[:, 2], xlabel=L"t", ylabel="neuron", yticks=yticks, limits=@lift ((nothing, nothing), (1 - 0.05 * $n, 1.05 * $n)))

    poly!((@lift Rect(0.0, 1.0, 100, floor(Int, $n * $r₊) - 0.5)), color=:red, alpha=0.5)
    poly!((@lift Rect(0.0, floor(Int, $n * $r₊) + 0.5, 100, $n - floor(Int, $n * $r₊) - 0.5)), color=:blue, alpha=0.5)
    scatter!(ax, (@lift [p.I for p in findall(transpose($simulation.X))]), markersize=5, color=:black)
end

begin # Adjust layout
    colsize!(fig.layout, 1, Fixed(300))
end