# HDR

This code base contains the material used to plot most of the figures in my HDR manuscript. Here is the exhaustive list of figures:

- Figure I.1: [julia script file](scripts/figure-I-1.jl)
- Figure II.1: [julia script file](scripts/figure-II-1.jl)
- Figure II.2: [julia script file](scripts/figure-II-2.jl)
- Figure II.3: [R script file](scripts/figure-II-3.R)
- Figure II.4: [julia script file](scripts/figure-II-4.jl)
- Figure II.5: [julia script file](scripts/figure-II-5.jl)
- Figure III.1: [julia script file](scripts/figure-III-1.jl)
- Figure III.2: [julia script file](scripts/figure-III-2.jl)
- Figure III.3: [julia script file](scripts/figure-III-3.jl)
- Figure III.4: [julia script file](scripts/figure-III-4.jl)
- Figure III.5: [julia script file](scripts/figure-III-5.jl)
- Figure III.6: [julia script file](scripts/figure-III-6.jl)
- Figure III.7: [julia script file](scripts/figure-III-7.jl)
- Figure III.8: [julia script file](scripts/figure-III-8.jl)
- Figure III.9: [julia script file](scripts/figure-III-9.jl)
- Figure III.10: [julia script file](scripts/figure-III-10.jl)
- Figure III.11: ask [Anna Melnykova](https://www.amelnykova.com)
- Figure III.12: it is a preliminary figure and should be implemented in a julia script file soon.
- Figure IV.1: see the [MeanFieldGraph.jl github repo](https://github.com/jucheval/MeanFieldGraph.jl)
- Figure IV.2: see the [MeanFieldGraph.jl github repo](https://github.com/jucheval/MeanFieldGraph.jl)
- Figure IV.3: [R script file](scripts/figure-IV-3.R)
- Figure IV.4: it is a preliminary figure and should be implemented in a julia script file soon.
- Figure A.1: [julia script file](scripts/figure-A-1.jl)

## How to reproduce ?

To (locally) reproduce this project, do the following (as advised by [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)):

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:

   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.
