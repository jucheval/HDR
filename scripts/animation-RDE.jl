using GLMakie
using DrWatson
using PointProcesses

img_RDE = load("notebooks/assets/RDE.png")
img_SRDE = load("notebooks/assets/SRDE.png") # n=100
img_ADHP = load("notebooks/assets/ADHP.png")

h_adhp, h_rde, h_srde = load("notebooks/assets/histories_RDE.jld2", "histories")

fig = Figure(size=(650, 800))

begin # buttons
    buttongrid = fig[1, 1:2] = GridLayout(tellwidth=false)
    buttonlabels = ["Heatmap", "Raster plot"]
    buttons = buttongrid[1, 1:length(buttonlabels)] = [Makie.Button(fig, label=buttonlabels[i]) for i in eachindex(buttonlabels)]
    buttons[1].buttoncolor = Makie.COLOR_ACCENT_DIMMED[]
end;

begin # Grid Layout
    g_hm = fig[2, 1] = GridLayout()
    g_raster = fig[2, 2] = GridLayout()
end

begin # observables

end

begin # Plot layout parameters
    marksize = 5
end

begin # Plot
    axs_hm = [Axis(g_hm[row, 1], aspect=DataAspect()) for row in 1:3]
    hidedecorations!.(axs_hm)
    hidespines!.(axs_hm)
    tightlimits!.(axs_hm)
    image!(axs_hm[1], rotr90(img_ADHP))
    image!(axs_hm[2], rotr90(img_SRDE))
    image!(axs_hm[3], rotr90(img_RDE))

    axs_raster = [Axis(g_raster[row, 1], ylabel="neuron", xlabel=L"t") for row in 1:3]
    tightlimits!.(axs_raster)
    scatter!(axs_raster[1], event_times(h_adhp), event_marks(h_adhp), markersize=marksize, color=:black)
    scatter!(axs_raster[2], event_times(h_srde), event_marks(h_srde), markersize=marksize, color=:black)
    scatter!(axs_raster[3], event_times(h_rde), event_marks(h_rde), markersize=marksize, color=:black)
end

begin # Adjust layout
    colsize!(fig.layout, 2, Fixed(0))

    ax_adhp_raster.blockscene.visible[] = false
    ax_adhp_raster.scene.visible[] = false
    ax_rde_raster.blockscene.visible[] = false
    ax_rde_raster.scene.visible[] = false
    ax_srde_raster.blockscene.visible[] = false
    ax_srde_raster.scene.visible[] = false
end

# Connect buttons to actions
on(buttons[1].clicks) do _

    colsize!(fig.layout, 1, Auto())
    ax_adhp_hm.blockscene.visible[] = true
    ax_adhp_hm.scene.visible[] = true
    ax_rde_hm.blockscene.visible[] = true
    ax_rde_hm.scene.visible[] = true
    ax_srde_hm.blockscene.visible[] = true
    ax_srde_hm.scene.visible[] = true

    colsize!(fig.layout, 2, Fixed(0))
    ax_adhp_raster.blockscene.visible[] = false
    ax_adhp_raster.scene.visible[] = false
    ax_rde_raster.blockscene.visible[] = false
    ax_rde_raster.scene.visible[] = false
    ax_srde_raster.blockscene.visible[] = false
    ax_srde_raster.scene.visible[] = false

    buttons[1].buttoncolor = Makie.COLOR_ACCENT_DIMMED[]
    buttons[2].buttoncolor = RGBf(0.94, 0.94, 0.94)
end
on(buttons[2].clicks) do _

    colsize!(fig.layout, 2, Auto())
    ax_adhp_raster.blockscene.visible[] = true
    ax_adhp_raster.scene.visible[] = true
    ax_rde_raster.blockscene.visible[] = true
    ax_rde_raster.scene.visible[] = true
    ax_srde_raster.blockscene.visible[] = true
    ax_srde_raster.scene.visible[] = true

    colsize!(fig.layout, 1, Fixed(0))
    ax_adhp_hm.blockscene.visible[] = false
    ax_adhp_hm.scene.visible[] = false
    ax_rde_hm.blockscene.visible[] = false
    ax_rde_hm.scene.visible[] = false
    ax_srde_hm.blockscene.visible[] = false
    ax_srde_hm.scene.visible[] = false

    buttons[2].buttoncolor = Makie.COLOR_ACCENT_DIMMED[]
    buttons[1].buttoncolor = RGBf(0.94, 0.94, 0.94)
end


fig



### Transforms history of ADHP into an age distribution to plot a heatmap
function history2sol(h, ts, as)
    times = event_times(h)
    marks = event_marks(h)
    sol = zeros(Float64, (length(ts), length(as)))
    dt = ts[2]

    ages = rand(Exponential(), 100) .+ 1
    tdata = Float64[]
    adata = Float64[]
    for idt in eachindex(ts)
        # append!(tdata, fill(ts[idt], 100))
        # append!(adata, ages)
        # kernel = kde(ages, boundary=(0.0,4.0))
        for ida in eachindex(as)
            previous_a = ida == 1 ? 0.0 : as[ida-1]
            sol[idt, ida] = sum(previous_a .< ages .<= as[ida]) / 100 / as[2]
        end
        # sol[idt,:] .= pdf(kernel, as)

        current_time = ts[idt]
        previous_time = idt == 1 ? 0.0 : ts[idt-1]

        ages .+= dt
        for id_spike in findall(previous_time .< times .<= current_time)
            ages[marks[id_spike]] = current_time - times[id_spike]
        end
    end

    # kernel = kde((tdata, adata), boundary=((0.0,10.0),(0.0,4.0)))
    # sol = pdf(kernel, ts, as)

    return sol
end