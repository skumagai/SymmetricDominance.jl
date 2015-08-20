using SymmetricDominance
using JLD

sttime = Dates.format(now(), "yyyy-m-d-H-M")
hostname = chomp(readall(`hostname`))
filename = "result.$hostname.$sttime.jld"

params = ModelParameters(
    1000,    # population size
    3,     # Number of loci (1 selected + n neutral loci)
    0.8,   # Fitness of heterozygote (absolute)
    1.0,   # Fitness of homozygote (absolute)
    [0.0, 0.1],   # Inter-loci recommbination rate
    [0.000001, 0.000001, 0.000001] # Per-locus mutation rate per-generation
)

results  = simulate(params, 5000, 30000, 3, 100)

@save filename results
