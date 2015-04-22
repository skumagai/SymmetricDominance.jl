module SymmetricDominance

export ModelParameters, simulate

immutable ModelParameters
    popsize::Int
    numberofloci::Int
    heterozygousfitness::Float64
    homozygousfitness::Float64
    recombinationrate::Float64
    mutationrates::Vector{Float64}
end

immutable LineageRecord
    parent::Int
    epoch::Int
end

immutable Gene
    state::Int
    lineage::Int

    Gene(s::Int) = (gene = new(); gene.state = s)
    Gene(s::Int, l::Int) = new(s, l)
end

immutable Organism
    genes::Array{Gene, 2}

    function Organism(nloci::Int, idx0::Int)
        genes = Array(Gene, 2, nloci)
        for locus = 1:nloci, chr = 1:2
            genes[locus, chr] = Gene(idx0, idx0)
            idx0 += 1
        end
        new(genes)
    end
    Organism(l::Array{Gene, 2}) = new(l)
end

Base.length(o::Organism) = size(o.genes, 1)

isidbystate(g1::Gene, g2::Gene) = g1.state == g2.state
# isidbydescent(g1::Gene, g2::Gene) = g1.lineage == g2.lineage

function mutate(g::Gene, midx::Int, ldb::Vector{LineageRecord}, epoch::Int)
    midx += 1
    push!(ldb, LineageRecord(g.lineage, epoch))
    midx, Gene(midx, length(ldb))
end

function getparentids!(ps, n)
    rand!(ps, 1:n)
    while ps[1] == ps[2]
        ps[2] = rand(1:n)
    end
    nothing
end

function hascoalesced(pops::Array{Organism, 2}, ldb::Vector{LineageRecord}, cidx::Int)
    coalesced = true
    n = size(pops, 1)
    nloci = length(pops[1,1])
    for locus = 1:nloci
        ancestor = getancestor(pops[1, cidx].genes[locus, 1], ldb)
        for org = 1:n, chr = 1:2
            if ancestor != getancestor(pops[org, cidx].genes[locus, chr], ldb)
                coalesced = false
                break
            end
        end
        coalesced || break
    end
    coalesced
end

function getancestor(gene::Gene, ldb::Vector{LineageRecord})
    lineage = gene.lineage
    while lineage > 0
        anc = ldb[lineage].parent
        anc == 0 && break
        lineage = anc
    end
    lineage
end

function evolve!(
    pops::Array{Organism, 2},
    ldb::Vector{LineageRecord},
    params::ModelParameters,
    t::Int,
    pidx::Int,
    cidx::Int,
    midx::Int,
    termoncoal::Bool)
    # unpacking parameters
    n = params.popsize
    heterofit = params.heterozygousfitness
    homofit = params.homozygousfitness
    recomb = params.recombinationrate
    muts = params.mutationrates
    nloci = params.numberofloci

    # normalize mutation rates
    maxfit = max(heterofit, homofit)
    heterofit /= maxfit
    homofit /= maxfit

    mutarray = Array(Bool, nloci, 2) # boolean value for each gene if it will be mutated.
    ps = Array(Int, 2) # indices of parents of an offspring.
    parchrs = Array(Int, 2) # a gene from which chromosome is passed on to offspring.

    gen = 1 # current generation
    for gen = 1:t
        for i = 1:n # iterate over offspring
            while true
                getparentids!(ps, n)
                # determine if mutations occur.
                for chr = 1:2, locus = 1:nloci
                    mutarray[locus, chr] = rand() < muts[locus] ? true : false
                end
                # process the first locus, which is under selection. A offspring is homozygous only when it
                # inherits identical-by-state genes from both parents without mutation. Otherwise, the offspring
                # is heterozygous.
                rand!(parchrs, 1:2)
                if isidbystate(pops[ps[1], pidx].genes[1, parchrs[1]], pops[ps[2], pidx].genes[1, parchrs[2]]) &&
                    mutarray[1,1] == mutarray[1,2] == false

                    if rand() > homofit
                        continue
                    end

                else
                    if rand() > heterofit
                        continue
                    end
                end

                for par = 1:2,  locus = 1:nloci
                    if mutarray[locus, par]
                        midx, g = mutate(pops[ps[par], pidx].genes[locus, parchrs[par]], midx, ldb, gen)
                        pops[i, cidx].genes[locus, par] = g
                    else
                        pops[i, cidx].genes[locus, par] = pops[ps[par], pidx].genes[locus, parchrs[par]]
                    end
                    parchrs[par] = rand() < recomb ? 3 - parchrs[par] : parchrs[par]
                end
                break
            end
        end
        pidx, cidx = cidx, pidx
        if termoncoal && hascoalesced(pops, ldb, cidx)
            println("Info: All lineages share a common ancestor at generation ", gen)
            break
        end
    end
    pidx, cidx, midx
end

function initialize!(pops::Array{Organism, 2}, ldb::Vector{LineageRecord}, params::ModelParameters)
    n = params.popsize
    nloci = params.numberofloci
    # Initialize a parental population. Genes are distinct.
    for i = 1:n
        pops[i, 1] = Organism(nloci, 2 * nloci * (i - 1) + 1)
    end
    # Initialize an offspring population. All organisms are just placeholders, as such values don't matter.
    for i = 1:n
        pops[i, 2] = Organism(nloci, 0)
    end
    append!(ldb, [LineageRecord(0, 0) for _ in 1:(2 * nloci * n)])
    nothing
end

function recalibratelineages!(pops::Array{Organism, 2}, ldb::Vector{LineageRecord}, params::ModelParameters, pidx::Int)
    empty!(ldb)
    lidx = 1
    for i = 1:params.popsize # iterate over a parental population
        for locus = 1:params.numberofloci, chr = 1:2
            g = pops[i, pidx].genes[locus, chr]
            pops[i, pidx].genes[locus, chr] = Gene(g.state, lidx)
            push!(ldb, LineageRecord(0, 0))
            lidx += 1
        end
    end
    nothing
end

function simulate(params::ModelParameters, burnin::Int, t::Int)
    ldb = Array(LineageRecord, 0) # a simple database of lineage, its parent, and epoch of the lineage.
    pops = Array(Organism, params.popsize, 2) # two populations, parental and offspring populations, stored as a 2-d array.
    midx = 2 * params.popsize * params.numberofloci # state of the gene, which have the last mutation.

    # Initialization
    # All genes are distinct.
    initialize!(pops, ldb, params)

    # Burnin
    # Execute the exact-same sequence as main-loop of evolution and throws out lineage information afterwords.
    # This loop runs exacctly "burnin" generations regardless of the presence of coalescence.
    pidx, cidx, midx = evolve!(pops, ldb, params, burnin, 1, 2, midx, false)
    # pidx, cidx: indices of parental and offspring population in "pops".

    # # Reset lineage information
    recalibratelineages!(pops, ldb, params, pidx)

    # Main loop of evolution
    # This loop terminates upon the first coalescence or after "t" generations.
    _, cidx, _ = evolve!(pops, ldb, params, t, pidx, cidx, midx, true)
    pops[:, cidx], ldb
end

end
