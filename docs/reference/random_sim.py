import numpy
import math
import msprime
import multiprocess

indnam = ["D_1", "C_1", "B_1", "E_1", "A_1"]

samples = [msprime.SampleSet(num_samples=i, population=j, time=k) for i, j, k in zip( (1, 1, 1, 1, 1), ("D", "C", "B", "E", "A"), (4000, 1000, 1000, 1000, 0) )]
nhap = int(len(indnam) * 2)

demo = msprime.Demography()
demo.add_population(initial_size = 1000, name = "R")
demo.add_population(initial_size = 1000, name = "D")
demo.add_population(initial_size = 1000, name = "Rr")
demo.add_population(initial_size = 1000, name = "Rrl")
demo.add_population(initial_size = 1000, name = "C")
demo.add_population(initial_size = 1000, name = "Rrr")
demo.add_population(initial_size = 1000, name = "B")
demo.add_population(initial_size = 1000, name = "E")
demo.add_population(initial_size = 1000, name = "Rl")
demo.add_population(initial_size = 1000, name = "admix")
demo.add_population(initial_size = 1000, name = "A")
demo.add_population_split(time=1000, ancestral="admix", derived=["A"])
demo.add_population_split(time=5000, ancestral="R", derived=["D","Rl"])
demo.add_population_split(time=4000, ancestral="Rl", derived=["Rr"])
demo.add_population_split(time=3000, ancestral="Rr", derived=["Rrl","Rrr"])
demo.add_population_split(time=2000, ancestral="Rrl", derived=["C"])
demo.add_population_split(time=2000, ancestral="Rrr", derived=["B","E"])
demo.add_admixture(time=1000, derived="admix", ancestral=["Rrl","Rl"], proportions=[0.5,0.5])
demo.sort_events()

r_chrom = 2e-08
r_break = math.log(2)
chrom_positions = [0, 1000]
map_positions = [0, 1000]
rates = [r_chrom]
rate_map = msprime.RateMap(position=map_positions, rate=rates)

ts = msprime.sim_ancestry(samples=samples, model=[msprime.DiscreteTimeWrightFisher(duration=25), msprime.StandardCoalescent()], demography=demo, recombination_rate=rate_map)
tree_sequence = msprime.sim_mutations(ts, rate=1.25e-08, model='binary')
hap_gt = tree_sequence.genotype_matrix()
gt = hap_gt[:, range(0,nhap,2)] + hap_gt[:, range(1,nhap,2)]
nsnps = gt.shape[0]
ts_chroms = numpy.searchsorted(numpy.array(chrom_positions[1:]), tree_sequence.tables.sites.position, side='right') + 1

numpy.savetxt('random_sim.geno', gt, '%d', '')

with open('random_sim.snp', 'w') as f:
  for i, j in enumerate(tree_sequence.tables.sites.position):
    chrm = int(ts_chroms[i])
    if chrm == 1:
      pos = int(j) 
    else:
      pos = int(j - chrom_positions[chrm-1])
    nsnps_chr = numpy.sum(ts_chroms == chrm)
    pos_chr = numpy.sum(ts_chroms[:i] == chrm)
    bytes = f.write('rs'+str(i+1)+'\t' + str(chrm) + '\t' + str(pos_chr*50/float(nsnps_chr-1)) + '\t' + str(pos) + '\tA\tC\n')

with open('random_sim.ind', 'w') as f:
  for i in range(len(indnam)):
    bytes = f.write(indnam[i] + '\tU\t' + indnam[i].rstrip("1234567890").rstrip("_") + '\n')

