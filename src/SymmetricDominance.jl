module SymmetricDominance

export Population,
       ModelParameters,

       simulate,
       tmrca,
       nsegsites,
       distances,
       spectra

using ForwardPopGenSimulations

immutable ModelParameters
    popsize::Int
    numberofloci::Int
    heterozygousfitness::Float64
    homozygousfitness::Float64
    recombinationrates::Vector{Float64}
    mutationrates::Vector{Float64}
end

function evolve!(
    core::BasicData,
    parpop::Population,
    params::ModelParameters,
    termon::Int,
    tclean::Int)

    info("evolution started on ", now(), ".")

    # unpacking parameters
    n = params.popsize
    heterofit = params.heterozygousfitness
    homofit = params.homozygousfitness
    # The last element is just there as a placeholder, and its value does not affect runs.
    recombs = [params.recombinationrates; 0.0]
    muts = [params.mutationrates params.mutationrates]
    nloci = params.numberofloci

    # normalize mutation rates
    maxfit = max(heterofit, homofit)
    heterofit /= maxfit
    homofit /= maxfit

    mutarray = Array{Bool}(nloci, 2) # boolean value for each gene if it will be mutated.
    ps = Array{Int}(2) # indices of parents of an offspring.
    parchrs = Array{Int}(2) # a gene from which chromosome is passed on to offspring.

    ncoals = [0 for _ = 1:nloci]
    lastcoals = [0 for _ = 1:nloci]

    chpop = deepcopy(parpop)

    gids = Array{Int}(2 * n)
    nnewids = 2 * n * nloci

    gen = 0 # current generation
    for gen = core
        for i = 1:n # iterate over offspring
            while true
                selectparents!(ps, n)
                # determine if mutations occur.
                selectmutatedsites!(mutarray, muts)
                # process the first locus, which is under selection. A offspring is homozygous only when it
                # inherits identical-by-state genes from both parents without mutation. Otherwise, the offspring
                # is heterozygous.
                rand!(parchrs, 1:2)
                if isidbystate(db(core), parpop[ps[1], 1, parchrs[1]], parpop[ps[2], 1, parchrs[2]]) &&
                    mutarray[1,1] == mutarray[1,2] == false

                    rand() > homofit && continue
                else
                    rand() > heterofit && continue
                end

                for par = 1:2,  locus = 1:nloci
                    pid = parpop[ps[par], locus, parchrs[par]]
                    if mutarray[locus, par]
                        chpop[i, locus, par] = transmit!(db(core), gen, pid, state=nextstate!(core))
                    else
                        chpop[i, locus, par] = transmit!(db(core), gen, pid)
                    end
                    parchrs[par] = rand() < recombs[locus] ? 3 - parchrs[par] : parchrs[par]
                end
                break
            end
        end
        parpop, chpop = chpop, parpop
        for locus in 1:nloci
            anc = mrca(db(core), locus, parpop)
            if anc.epoch > lastcoals[locus]
                lastcoals[locus] = anc.epoch
                ncoals[locus] += 1
                info("turn over ", ncoals[locus], " of locus ", locus, " on generation ", gen)
            end
        end
        if termon == minimum(ncoals)
            break
        end
        gen % tclean == 0 && clean!(db(core), nnewids * gen + 1, nnewids * (gen + 1))
    end
    clean!(db(core), nnewids * gen + 1, nnewids * (gen + 1))

    if termon == minimum(ncoals)
        info("evolution terminated by reaching the min # of turn-overs in gen ", gen, " on ", now(), ".")
    else
        info("evolution terminated by reaching max gen ", gen, " on ", now(), ".")
    end

    parpop, core
end

function simulate(params::ModelParameters, burnin::Int, t::Int, termon::Int, tclean::Int, rep=1)
    info("process started on ", now(), ".")
    core = BasicData()

    # This is a parental population, a population of offspring is created within evolve! function.
    pop = Population(params.popsize, params.numberofloci, 2)

    # Burnin
    # Execute the exact-same sequence as main-loop of evolution and throws out lineage information afterwords.
    # This loop runs exacctly "burnin" generations regardless of the presence of coalescence.
    initialize!(core, pop) # all genes are distinct
    settmax!(core, burnin)
    pop, core = evolve!(core, pop, params, -1, tclean)

    datastore = Array{Tuple{typeof(pop), typeof(db(core)), typeof(time(core))}}(0)

    # Main loop of evolution
    # This loop terminates upon the first coalescence or after "t" generations.
    for _ = 1:rep
        pop = deepcopy(pop)
        core = reinitialize!(core, pop)
        settmax!(core, t)
        pop, core = evolve!(core, pop, params, termon, tclean)
        push!(datastore, (pop, db(core), time(core)))
    end
    info("process terminated on ", now(), ".")
    datastore
end

end
