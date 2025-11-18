using DrWatson
using DataFrames
using GLMakie

begin # Load data
    data = load("notebooks/assets/quantiles_abserror.jld2")

    labels = ["q25",
        "q50",
        "q75",
        "q25 large time",
        "q50 large time",
        "q75 large time"
    ]

    q25, q50, q75, q25inf, q50inf, q75inf = [data[label] for label in labels]

    estimatornames = [L"\hat{m}", L"\hat{v}", L"\hat{w}", L"\hat{\mu}", L"\hat{\lambda}", L"\hat{p}"]
end


fig = Figure(size=(1200, 600))

# tog = Toggle(fig[2, 2], active=true)
# visibility = tog.active

ax = Axis(fig[:, 1], yscale=log10, ylabel="absolute error", xlabel=L"t")
ylims!(ax, 4.0e-5, 2.7)

for id in 1:6
    lines!(ax, q50[:, 1], q50[:, id+1]; color=id, colormap=:tab10, colorrange=(1, 10), label=estimatornames[id], visible=id == 1)
    band!(ax, q50[:, 1], q25[:, id+1], q75[:, id+1]; color=id, colormap=:tab10, colorrange=(1, 10), label=estimatornames[id], alpha=0.1, visible=id == 1)
    scatter!(ax, q50[end, 1], q50inf[1, id]; color=id, colormap=:tab10, colorrange=(1, 10), label=estimatornames[id], visible=id == 1)
end

Legend(fig[1, 2], ax, merge=true)
