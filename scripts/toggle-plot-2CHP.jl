using GLMakie
using DrWatson


begin # Load data
    data = load("notebooks/assets/2CHP-coupling.jld2")

    labels = ["2CHP",
        "diffusion coupled",
        "diffusion independent",
        "MF limit"
    ]

    counting, diffusion_coupled, diffusion_indep, ode = [data[label] for label in labels]
end;

fig = Figure(size=(800, 480))

begin # observables
    gl = GridLayout(fig[3, 1], tellwidth=false)
    Label(gl[1, 1], "MC-HP")
    Label(gl[1, 3], "Diffusion coupled")
    Label(gl[1, 5], "Diffusion independent")
    Label(gl[1, 7], "MF limit")
    toggles = [Toggle(gl[1, 2*i], active=(i == 1)) for i in 1:4]

    gl_alpha = GridLayout(fig[3, 2], tellwidth=false)
    Label(gl_alpha[1, 1], "transparency")
    sl_alpha = Slider(gl_alpha[2, 1], range=0:0.01:1, startvalue=0.25)
    alpha_factor = sl_alpha.value
end

begin # axis layout 
    axtop = Axis(fig[1, 1], ylabel=L"U^{1,m'}")
    axbot = Axis(fig[2, 1], ylabel=L"U^{2,m'}", xlabel=L"t")
    linkxaxes!(axtop, axbot)
    hidexdecorations!(axtop, grid=false)
    xlims!(axbot, 0, 140)
end

begin # plots
    # 2CHP lines
    for m in 1:4
        lines!(axtop, counting.spike_train[begin:20:end], counting.X[m, begin:20:end];
            color=:darkblue, label="MC-HP", visible=toggles[1].active, alpha=@lift 1 - $alpha_factor * (m - 1))
    end
    for m in 1:3
        lines!(axbot, counting.spike_train[begin:20:end], counting.X[4+m, begin:20:end];
            color=:darkblue, visible=toggles[1].active, alpha=@lift 1 - $alpha_factor * (m - 1))
    end

    # Coupled diffusion lines
    for m in 1:4
        lines!(axtop, diffusion_coupled.time, diffusion_coupled.X[m, :],
            color=:orange, label="diff. coupled", visible=toggles[2].active, alpha=@lift 1 - $alpha_factor * (m - 1))
    end
    for m in 1:3
        lines!(axbot, diffusion_coupled.time, diffusion_coupled.X[4+m, :],
            color=:orange, visible=toggles[2].active, alpha=@lift 1 - $alpha_factor * (m - 1))
    end

    # Independent diffusion lines
    for m in 1:4
        lines!(axtop, diffusion_indep.time, diffusion_indep.X[m, :],
            color=:red, label="diff. indep.", visible=toggles[3].active, alpha=@lift 1 - $alpha_factor * (m - 1))
    end
    for m in 1:3
        lines!(axbot, diffusion_indep.time, diffusion_indep.X[4+m, :],
            color=:red, visible=toggles[3].active, alpha=@lift 1 - $alpha_factor * (m - 1))
    end

    # MF limit lines
    for m in 1:4
        lines!(axtop, ode.time, ode.X[m, :];
            color=:black, label="ODE", visible=toggles[4].active, alpha=@lift 1 - $alpha_factor * (m - 1))
    end
    for m in 1:3
        lines!(axbot, ode.time, ode.X[4+m, :];
            color=:black, visible=toggles[4].active, alpha=@lift 1 - $alpha_factor * (m - 1))
    end
end

# legend
Legend(fig[1:2, 2], axtop, merge=true)
