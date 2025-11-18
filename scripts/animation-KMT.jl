using GLMakie
using DrWatson

fig = Figure(size=(1000, 700))

begin # time grid
    kmax = 6
    tmax = 2^3
    ts = tmax * (0:2^kmax) / 2^kmax
end

begin # load data
    data = load("notebooks/assets/coupling_KMT.jld2")
    W = data["Brownian"]
    N = data["Poisson"]
    counting_values = [0.0]
    for i in eachindex(N)
        append!(counting_values, -N[i] .+ [i - 1, i])
    end
    push!(counting_values, -tmax + length(N))
end

begin # figure layout parameters
    pad = 0.2
end

begin # observables
    btn = Makie.Button(fig[2, 1], label="Start animation", tellwidth=false)

    ts_black = Observable(Float64[])
    ts_gray = Observable(Float64[])
    Ws = Observable(Float64[])
end

begin # plots
    ax = Axis(fig[1, 1], limits=((0, tmax) .+ pad .* (-1, 1), (minimum(W), maximum(W)) .+ pad .* (-1, 1)),
        title="KMT coupling", xlabel=L"t", ylabel="",
        xgridvisible=false)
    vlines!(ax, ts_black, color=:black)
    vlines!(ax, ts_gray, color=:gray90)

    lines!(ax, [0; repeat(N, inner=2); tmax], counting_values, alpha=0.7, color=:darkblue)
    initial_W_plot = lines!(ax, ts, W; color=:orange)
    lines!(ax, ts_gray, Ws; color=:orange)
end

# Update visualization on each button click
on(btn.clicks) do clicks

    if clicks == 1
        initial_W_plot.visible = false
        global ks = collect(0:(kmax-1))
        global ids_black = Int[]
        global ids_gray = Int[]
    end

    clicks > 2 * kmax && return
    step = div(clicks, 2)

    if isodd(clicks)  # Odd click: new lines in black
        step == 0 && push!(ids_black, 0, 1)
        for k in ks
            append!(ids_black, 2^k .+ 2^(k - step) * (1:2:2^step))
        end
        popfirst!(ks)

        ts_black[] = ts[ids_black.+1]

        # Update button label
        btn.label[] = step == 0 ? "Compute quantile coupling" : "Compute conditional quantile coupling"
    else  # Even click: turn black lines into gray
        append!(ids_gray, ids_black)
        sort!(ids_gray)
        empty!(ids_black)
        empty!(ts_black[])
        ts_gray[] = ts[ids_gray.+1]
        Ws[] = W[ids_gray.+1]

        # Update button label
        btn.label[] = "Add next ticks"
    end

    clicks == 2 * kmax && (btn.label[] = "Coupling complete!")
end
