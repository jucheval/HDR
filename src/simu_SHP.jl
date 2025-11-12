using PointProcesses
using Distributions
using CairoMakie

struct SpatialHawkesProcess
    n::Int
    decayrate::Float64
    synapticweight::Function
    positions::Vector
    firingrate::Function
    firingratebound::Float64

    function SpatialHawkesProcess(n, decayrate, synapticweight, positions, firingrate, firingratebound)
        length(positions) == n || error("The length of positions must be equal to n")
        new(n, decayrate, synapticweight, positions, firingrate, firingratebound)
    end
end

SpatialHawkesProcess(; n, decayrate, synapticweight, positions, firingrate, firingratebound) = SpatialHawkesProcess(n, decayrate, synapticweight, positions, firingrate, firingratebound)
getfields(shp::SpatialHawkesProcess) = shp.n, shp.decayrate, shp.synapticweight, shp.positions, shp.firingrate, shp.firingratebound

include("simu_Poisson.jl") # function to simulate an homogeneous Poisson process

function simulate(shp::SpatialHawkesProcess, initial_condition::Function, domain::Vector{Float64}; length_ts::Int=1000)
    # notations
    Nneur, α, w, xs, f, f̅ = getfields(shp)
    tmin, tmax = domain
    ts = range(tmin, tmax, length_ts)

    # spatial grid and initial conditions
    wmat = synapticweightmatrix(xs, w)
    Us = initial_condition.(xs)

    # simulation of a dominating Poisson Process with notations
    dominatingPoisson = rand(MultivariatePoissonProcess(fill(f̅, Nneur)), tmin, tmax)
    npts = length(dominatingPoisson)
    times = event_times(dominatingPoisson)
    marks = event_marks(dominatingPoisson)

    # pre-allocation
    Umat = zeros(Float64, (Nneur, length_ts))

    # initialization
    idt = 1
    Umat[:, idt] = Us
    for k in 1:npts
        time_step = k == 1 ? times[1] : times[k] - times[k-1]
        Us .*= exp(-α * time_step)
        pt_mark = marks[k]
        intensity = f(Us[pt_mark])
        if (rand() > intensity / f̅)
            marks[k] = 0
        else
            synaptic_vec = @view wmat[pt_mark, :]
            Us += synaptic_vec * (1 / Nneur)
        end

        if times[k] > ts[idt+1]
            idt += 1
            Umat[:, idt] = Us
        end
    end
    Umat[:, end] = Us

    actual_spikes = marks .!= 0
    h = History(times[actual_spikes], marks[actual_spikes], tmin, tmax)
    return h, ts, xs, transpose(Umat)
end


function synapticweightmatrix(xs, w)
    [w(xs[i], xs[j]) for i in eachindex(xs), j in eachindex(xs)]
end

function Makie.plot(::Type{SpatialHawkesProcess}, simulation)
    h, ts, xs, Umat = simulation
    fig = Figure(size=(600, 300))

    axleft = Axis(fig[1, 1], xlabel=L"t", ylabel=L"x", title=L"heatmap of $U^n(t,x)$")
    axright = Axis(fig[1, 3], xlabel=L"t", ylabel=L"x", title="raster plot")

    hm = heatmap!(axleft, ts, xs, Umat)
    Colorbar(fig[1, 2], hm)

    ids = findall(event_marks(h) .!= 0)
    scatter!(axright, event_times(h)[ids], xs[event_marks(h)[ids]], markersize=1, color=:black)

    fig
end