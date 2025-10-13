using CairoMakie

u_in(x) = x
u₀(t) = t
v = 1
u(t, x) = x - v * t > 0 ? u_in(x - v * t) : u₀(t - x / v)

xs = 0:0.1:10
ts = 0:0.1:10

begin # plot
    fig = Figure(size=(300, 300))

    ax = Axis(fig[1, 1], aspect=DataAspect(), xlabel=L"t", ylabel=L"x", title=L"heatmap of $u(t,x)$")

    hm = heatmap!(ax, ts, xs, u, colormap=Reverse(:viridis))
    for k in -10:2:10
        lines!(ax, ts, ts .+ k, color=:grey)
    end
    limits!(ax, 0, 10, 0, 10)

    Colorbar(fig[:, end+1], hm)

    fig
end

# save figure
#save("plots/EDP-characteristics.png", fig, px_per_unit=4)
