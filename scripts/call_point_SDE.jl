include("../src/coupling_point_SDE.jl")
using CairoMakie

begin # parameters
    high = 1.0
    low = 0.95
    width = 50.
    ipp = InhomogeneousPoissonProcess(
        intensity=t ->
            high * (rem(t, 2 * width, RoundNearest) >= 0) +
            low * (rem(t, 2 * width, RoundNearest) < 0),
        intensitybound=high
    )
    ipptilde = InhomogeneousPoissonProcess(
        intensity=t ->
            low * (rem(t, 2 * width, RoundNearest) >= 0) +
            high * (rem(t, 2 * width, RoundNearest) < 0),
        intensitybound=high
    )
end
tmin = 0.0
tmax = 100.0
domain = [tmin, tmax]

begin # plot
    Random.seed!(1)
    fig = Figure()

    axleft = Axis(fig[1, 1], xlabel=L"t", ylabel="count")
    axright = Axis(fig[1, 2], xlabel=L"t")
    linkyaxes!(axleft, axright)
    hideydecorations!(axright, grid=false)
    axleft.yticks = 0:25:tmax
    axright.yticks = 0:25:tmax

    N, Ñ = coupling(ipp, ipptilde, domain)
    stairs!(axleft, [0.0; N; tmax], [0:length(N); length(N)], alpha=0.7, color=:darkblue)
    stairs!(axleft, [0.0; Ñ; tmax], [0:length(Ñ); length(Ñ)], alpha=0.7, color=:red)

    N, Ñ = coupling(ipp, ipptilde, domain)
    stairs!(axright, [0.0; N; tmax], [0:length(N); length(N)], alpha=0.7, color=:darkblue)
    stairs!(axright, [0.0; Ñ; tmax], [0:length(Ñ); length(Ñ)], alpha=0.7, color=:red)

    fig
end