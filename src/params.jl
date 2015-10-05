import YAML

immutable Parameters
    burnin::Int
    tmax::Int
    coalterm::Int
    tclean::Int
    repeat::Int
    nloci::Int
    popsize::Int
    fithet::Float64
    fithom::Float64
    mutrate::Array{Float64, 1}
    recrate::Array{Float64, 1}
end

function loadparams(path)
    conf = YAML.load_file(path)
    Parameters(
        burnin(conf),
        tmax(conf),
        coalterm(conf),
        tclean(conf),
        repeat(conf),
        nloci(conf),
        popsize(conf),
        fithet(conf),
        fithom(conf),
        mutrate(conf),
        recrate(conf)
    )
end

burnin(conf) = Int(conf["burn in"])
tmax(conf) = Int(conf["simulated generations"])
coalterm(conf) = Int(conf["terminate"])
nloci(conf) = Int(conf["number of loci"])
popsize(conf) = Int(conf["population size"])
fithet(conf) = Float64(conf["fitness of heterozygotes"])
fithom(conf) = Float64(conf["fitness of homozygotes"])
mutrate(conf) = convert(Array{Float64, 1}, conf["mutation rates"])
recrate(conf) = convert(Array{Float64, 1}, conf["recombination rates"])

# parameters with default values
function repeat(conf)
    if haskey(conf, "repeat")
        Int(conf["repeat"])
    else
        1
    end
end

function tclean(conf)
    if haskey(conf, "clean up frequency")
        Int(conf["clean up frequency"])
    else
        100
    end
end
