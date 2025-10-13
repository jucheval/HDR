include("../src/coupling_diffusion_timechange.jl")
using CairoMakie

begin # parameters
    width = 5
    fast = 1.0
    slow = 0.8
    tcbm = TimeChangeBrownianMotion(
        timescaling=x ->
            fast * (rem(x, 2 * width, RoundNearest) >= 0) +
            slow * (rem(x, 2 * width, RoundNearest) < 0)
    )
    tcbmtilde = TimeChangeBrownianMotion(
        timescaling=x ->
            slow * (rem(x, 2 * width, RoundNearest) >= 0) +
            fast * (rem(x, 2 * width, RoundNearest) < 0)
    )
end
tmin = 0.
tmax = 100.
domain = [tmin, tmax]
dt = 0.1

begin # plot
    Random.seed!(1)
    fig = Figure(size=(800, 480))

    axtop = Axis(fig[1, 1])
    axbot = Axis(fig[2, 1], xlabel=L"t")
    linkxaxes!(axtop, axbot)
    hidexdecorations!(axtop, grid=false)
    axtop.yticks = (-10*width):width:(10*width)
    axbot.yticks = (-10*width):width:(10*width)

    ts, t̃s, Bs = coupling(tcbm, tcbmtilde, 0.0, domain, dt)
    lines!(axtop, ts[ts.<tmax], Bs[ts.<tmax], alpha=0.8, color=:darkblue)
    lines!(axtop, t̃s[t̃s.<tmax], Bs[t̃s.<tmax], alpha=0.8, color=:orange)

    ts, t̃s, Bs = coupling(tcbm, tcbmtilde, 0.0, domain, dt)
    lines!(axbot, ts[ts.<tmax], Bs[ts.<tmax], alpha=0.8, color=:darkblue)
    lines!(axbot, t̃s[t̃s.<tmax], Bs[t̃s.<tmax], alpha=0.8, color=:orange)

    fig
end

# save figure
#save("plots/coupling_diffusion_timechange.png", fig, px_per_unit=2)