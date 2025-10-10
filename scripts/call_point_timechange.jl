include("../src/coupling_point_timechange.jl")
using CairoMakie

begin # parameters
    high = 1.0
    low = 0.9
    width = 50.
    itc = InhomogeneousTimeChange(
        inversecumulativeintensity=τ ->
            τ / high * (τ <= high * width) +
            (width + (τ - high * width) / low) * (τ > high * width)
    )
    itctilde = InhomogeneousTimeChange(
        inversecumulativeintensity=τ ->
            τ / low * (τ <= low * width) +
            (width + (τ - low * width) / high) * (τ > low * width)
    )
end
tmin = 0.0
τmin = 0.0
tmax = 100.
τmax = high * width + low * (tmax - width)
domain = [τmin, τmax]

begin # plot
    Random.seed!(1)
    fig = Figure(size=(800, 300))

    axleft = Axis(fig[1, 1], xlabel=L"t", ylabel="count")
    axright = Axis(fig[1, 2], xlabel=L"t")
    linkyaxes!(axleft, axright)
    hideydecorations!(axright, grid=false)
    axleft.yticks = 0:25:100
    axright.yticks = 0:25:100

    N, Ñ = coupling(itc, itctilde, domain)
    stairs!(axleft, [0.0; N; 100.0], [0:length(N); length(N)], alpha=0.7, color=:darkblue)
    stairs!(axleft, [0.0; Ñ; 100.0], [0:length(Ñ); length(Ñ)], alpha=0.7, color=:orange)

    N, Ñ = coupling(itc, itctilde, domain)
    stairs!(axright, [0.0; N; 100.0], [0:length(N); length(N)], alpha=0.7, color=:darkblue)
    stairs!(axright, [0.0; Ñ; 100.0], [0:length(Ñ); length(Ñ)], alpha=0.7, color=:orange)

    fig
end

# save figure
#save("plots/coupling-point-timechange.png", fig, px_per_unit=2)