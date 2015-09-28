import SymmetricDominance
using Base.Test

const fpgs = ForwardPopGenSimulations
const sd = SymmetricDominance

gdb = fpgs.GeneDB()

type StateCounter
    state::Int
end

nextstate!(s::StateCounter) = (s.state += 1; s.state)

sc = StateCounter(0)

for _ = 1:2
    sd.insert!(gdb, sd.GeneRecord(1, nextstate!(sc)))
end
for i = 1:2, _ = 1:2
    sd.insert!(gdb, sd.GeneRecord(2, gdb[i], state=nextstate!(sc)))
end
for i = 3:6, _ = 1:2
    sd.insert!(gdb, sd.GeneRecord(3, gdb[i], state=nextstate!(sc)))
end
sd.insert!(gdb, sd.GeneRecord(4, gdb[1]))

for i = 1:15
    gdb[i].state += 10
end

pop = sd.Population(2, 1, 2)
pop[1, 1, 1] = 7
pop[2, 1, 1] = 8
pop[1, 1, 2] = 9
pop[2, 1, 2] = 10

sdata3 = Array(Int, 2, 1, 2)
sdata3[1, 1, 1] = 7
sdata3[2, 1, 1] = 8
sdata3[1, 1, 2] = 9
sdata3[2, 1, 2] = 10
@test sd.toarray(gdb, pop, :id) == sdata3
sdata3[1, 1, 1] = 17
sdata3[2, 1, 1] = 18
sdata3[1, 1, 2] = 19
sdata3[2, 1, 2] = 20
@test sd.toarray(gdb, pop, :state) == sdata3
sdata3[1, 1, 1] = 3
sdata3[2, 1, 1] = 3
sdata3[1, 1, 2] = 4
sdata3[2, 1, 2] = 4
@test sd.toarray(gdb, pop, :parent) == map(x -> gdb[x], sdata3)
sdata3[1, 1, 1] = 3
sdata3[2, 1, 1] = 3
sdata3[1, 1, 2] = 3
sdata3[2, 1, 2] = 3
@test sd.toarray(gdb, pop, :epoch) == sdata3

const acs, gcs, hcs = sd.counts(gdb, pop)
@test acs == [Dict{Int,Int}(17=>1, 18=>1, 19=>1, 20=>1)]
@test gcs == [Dict{Tuple{Int,Int}, Int}((17, 19)=>1, (18, 20)=>1)]
@test hcs == Dict{Tuple{Int}, Int}((17,)=>1, (18,)=>1, (19,)=>1, (20,)=>1)

const afs, gfs, hfs = sd.spectra(gdb, pop)
@test afs == [Dict{Int,Int}(1=>4)]
@test gfs == [Dict{Int,Int}(1=>2)]
@test hfs == Dict{Int,Int}(1=>4)
