using DrWatson
using PointProcesses
using CairoMakie

begin # Load data
    data = load("notebooks/assets/histories_NFE.jld2")

    labels = ["neurons positions and spatial grid for potential SHP",
        "history SHP",
        "history NFE",
        "history SNFE",
        "time grid",
        "spatial grid for potential NFE and SNFE",
        "potential SHP",
        "potential NFE",
        "potential SNFE"
    ]

    xs_shp, h_shp, h_nfe, h_snfe, ts, xs_nfe, U_shp, U_nfe, U_snfe = [data[label] for label in labels]
end

begin # Heatmap layout parameters
    resolution_hm = (400, 380)
    xlabel_hm = L"t"
    ylabel_hm = L"x"
end

# NFE heatmap
fig = Figure(size=resolution_hm);
ax = Axis(fig[1, 1], xlabel=xlabel_hm, ylabel=ylabel_hm)
hm = heatmap!(ax, ts, xs_nfe, U_nfe)
Colorbar(fig[1, 2], hm)
safesave("notebooks/assets/raw/nfe-hm.png", fig)

# SNFE heatmap
fig = Figure(size=resolution_hm);
ax = Axis(fig[1, 1], xlabel=xlabel_hm, ylabel=ylabel_hm)
hm = heatmap!(ax, ts, xs_nfe, U_snfe)
Colorbar(fig[1, 2], hm)
safesave("notebooks/assets/raw/snfe-hm.png", fig)

# SHP heatmap
fig = Figure(size=resolution_hm);
ax = Axis(fig[1, 1], xlabel=xlabel_hm, ylabel=ylabel_hm)
hm = heatmap!(ax, ts, xs_shp, U_shp)
Colorbar(fig[1, 2], hm)
safesave("notebooks/assets/raw/shp-hm.png", fig)

begin # Raster layout parameters
    resolution_raster = (400, 380)
    xlabel_raster = L"t"
    ylabel_raster = L"x"
    marksize = 2
end

# NFE raster
fig = Figure(size=(500, 380));
ax = Axis(fig[1, 1], xlabel=xlabel_raster, ylabel=ylabel_raster, limits=((nothing, nothing), (-pi, pi)))
scatter!(ax, event_times(h_nfe), xs_shp[event_marks(h_nfe)], markersize=marksize, color=:black)
safesave("notebooks/assets/raw/nfe-raster.png", fig)

# SNFE raster
fig = Figure(size=(500, 380));
ax = Axis(fig[1, 1], xlabel=xlabel_raster, ylabel=ylabel_raster, limits=((nothing, nothing), (-pi, pi)))
scatter!(ax, event_times(h_snfe), xs_shp[event_marks(h_snfe)], markersize=marksize, color=:black)
safesave("notebooks/assets/raw/snfe-raster.png", fig)

# SHP raster
fig = Figure(size=(500, 380));
ax = Axis(fig[1, 1], xlabel=xlabel_raster, ylabel=ylabel_raster, limits=((nothing, nothing), (-pi, pi)))
scatter!(ax, event_times(h_shp), xs_shp[event_marks(h_shp)], markersize=marksize, color=:black)
safesave("notebooks/assets/raw/shp-raster.png", fig)