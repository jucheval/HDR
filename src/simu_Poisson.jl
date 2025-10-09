# Needed here because the current version of PointProcesses.jl on the general repository outputs an unsorted event_times vector.
function Base.rand(pp::MultivariatePoissonProcess, tmin::Float64, tmax::Float64
)
    mark_dist = mark_distribution(pp)
    N = rand(Poisson(float(ground_intensity(pp) * (tmax - tmin))))
    times = sort!(rand(Uniform(tmin, tmax), N))
    marks = [rand(mark_dist) for n in 1:N]
    return History(; times=times, marks=marks, tmin=tmin, tmax=tmax)
end