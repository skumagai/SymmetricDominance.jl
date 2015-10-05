using SymmetricDominance
if length(ARGS) == 1
    startsimulation(ARGS[1])
else
    error("usage: julia symmetric-dominance.jl CONFIG_FILE")
end

