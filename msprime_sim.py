import math
import msprime
import argparse

parser = argparse.ArgumentParser(description = 'Simulating pseudo sequences by chromosomes')
parser.add_argument('-c', '--CHR', type = int, default = 1, help = 'Chromosomes number')
parser.add_argument('-m', '--migration_percent', type = int, default = 5, help = 'Gene flow percent')
parser.add_argument('-i', '--iteration', type = int, default = 1, help ='Iteration and seed')
args = parser.parse_args()

#eff. pop. sizes
O = 1000
U = 1000
V1 = 1000
V2 = 1000
X2 = 1000
X1 = 1000
ABC1 = 1000
ABC2 = 1000
BC1 = 1000
BC2 = 1000
AB1 = 1000
AB2 = 1000
C1 = 1000
C2 = 1000
W1 = 1000
W2 = 1000
B1 = 1000
B2 = 1000
A1 = 1000
A2 = 1000

gen = 29
pops = [
	msprime.PopulationConfiguration(initial_size = O), #0
	msprime.PopulationConfiguration(initial_size = U), #1
	msprime.PopulationConfiguration(initial_size = V1), #2
	msprime.PopulationConfiguration(initial_size = V2), #3
	msprime.PopulationConfiguration(initial_size = X2), #4
	msprime.PopulationConfiguration(initial_size = X1), #5
	msprime.PopulationConfiguration(initial_size = ABC1), #6
	msprime.PopulationConfiguration(initial_size = ABC2), #7
	msprime.PopulationConfiguration(initial_size = BC1), #8
	msprime.PopulationConfiguration(initial_size = BC2), #9
	msprime.PopulationConfiguration(initial_size = AB1), #10
	msprime.PopulationConfiguration(initial_size = AB2), #11
	msprime.PopulationConfiguration(initial_size = C1), #12
	msprime.PopulationConfiguration(initial_size = C2), #13
	msprime.PopulationConfiguration(initial_size = W1), #14
	msprime.PopulationConfiguration(initial_size = W2), #15
	msprime.PopulationConfiguration(initial_size = B1), #16
	msprime.PopulationConfiguration(initial_size = B2), #17
	msprime.PopulationConfiguration(initial_size = A1), #18
	msprime.PopulationConfiguration(initial_size = A2) #19
]

#ind. dates
samples = [
	msprime.Sample(0, 0/gen), #O
	msprime.Sample(1, 0/gen), #U
	msprime.Sample(2, 0/gen), #V1
	msprime.Sample(3, 0/gen), #V2
	msprime.Sample(4, 0/gen), #X2
	msprime.Sample(5, 0/gen), #X1
	msprime.Sample(6, 0/gen), #ABC1
	msprime.Sample(7, 0/gen), #ABC2
	msprime.Sample(8, 0/gen), #BC1
	msprime.Sample(9, 0/gen), #BC2
	msprime.Sample(10, 0/gen), #AB1
	msprime.Sample(11, 0/gen), #AB2
	msprime.Sample(12, 0/gen), #C1
	msprime.Sample(13, 0/gen), #C2
	msprime.Sample(14, 0/gen), #W1
	msprime.Sample(15, 0/gen), #W2
	msprime.Sample(16, 0/gen), #B1
	msprime.Sample(17, 0/gen), #B2
	msprime.Sample(18, 0/gen), #A1
	msprime.Sample(19, 0/gen)  #A2
]

#pop. split dates
T_R = 100
T_U = 200
T_V2 = 300
T_V1 = 400
T_V2 = 500
T_ABC1 = 600
T_ABC2 = 700
T_BC1 = 800
T_BC2 = 900
T_AB2 = 1000
T_AB1 = 1100
T_AB2 = 1200

events = [
	msprime.MassMigration(time = T_R/gen, source = 1, destination = 0, proportion = 1),
	msprime.MassMigration(time = T_U/gen, source = 5, destination = 1, proportion = 1),
	msprime.MassMigration(time = T_V2/gen, source = 4, destination = 5, proportion = 1),
	msprime.MassMigration(time = T_V1/gen, source = 2, destination = 1, proportion = 1),
	msprime.MassMigration(time = T_V2/gen, source = 6, destination = 5, proportion = 1),
	msprime.MassMigration(time = T_ABC1/gen, source = 4, destination = 2, proportion = 1),
	msprime.MassMigration(time = T_ABC2/gen, source = 8, destination = 6, proportion = 1),
	msprime.MassMigration(time = T_BC1/gen, source = 3, destination = 2, proportion = 1),
	msprime.MassMigration(time = T_BC2/gen, source = 7, destination = 6, proportion = 1),
	msprime.MassMigration(time = T_AB2/gen, source = 1, destination = 8, proportion = 1),
	msprime.MassMigration(time = T_AB1/gen, source = 3, destination = 4, proportion = 1),
	msprime.MassMigration(time = T_AB2/gen, source = 7, destination = 8, proportion = 1)
]

tree_sequence = msprime.simulate(population_configurations = pops,
                                 length = 1e+09,
                                 samples = samples,
                                 demographic_events = events,
                                 recombination_rate = 1e-09,
                                 mutation_rate = 1e-09,
                                 random_seed = args.iteration)

with open('msprimesim.vcf', 'w') as f:
  tree_sequence.write_vcf(f, ploidy = 2)
