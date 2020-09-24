import math
import numpy
import msprime
gen = 29
pops = [
	msprime.PopulationConfiguration(initial_size = 100), #0 R
	msprime.PopulationConfiguration(initial_size = 100), #1 Rl
	msprime.PopulationConfiguration(initial_size = 100), #2 Denisova.DG
	msprime.PopulationConfiguration(initial_size = 100), #3 Rlx
	msprime.PopulationConfiguration(initial_size = 100), #4 admix1za
	msprime.PopulationConfiguration(initial_size = 100), #5 Rll
	msprime.PopulationConfiguration(initial_size = 100), #6 Chimp.REF
	msprime.PopulationConfiguration(initial_size = 100), #7 Rllr
	msprime.PopulationConfiguration(initial_size = 100), #8 Rllrl
	msprime.PopulationConfiguration(initial_size = 100), #9 admix1a
	msprime.PopulationConfiguration(initial_size = 100), #10 Vindija.DG
	msprime.PopulationConfiguration(initial_size = 100), #11 Altai_Neanderthal.DG
	msprime.PopulationConfiguration(initial_size = 100), #12 Rlr
	msprime.PopulationConfiguration(initial_size = 100), #13 Mbuti.DG
	msprime.PopulationConfiguration(initial_size = 100), #14 Rlrx
	msprime.PopulationConfiguration(initial_size = 100), #15 admix1b
	msprime.PopulationConfiguration(initial_size = 100), #16 admix1
	msprime.PopulationConfiguration(initial_size = 100), #17 Switzerland_Bichon.SG
	msprime.PopulationConfiguration(initial_size = 100), #18 Russia_Ust_Ishim.DG
	msprime.PopulationConfiguration(initial_size = 100), #19 admix1zb
	msprime.PopulationConfiguration(initial_size = 100) #20 admix1z
]

indnam = ["Denisova.DG_1", "Denisova.DG_2", "Denisova.DG_3", "Denisova.DG_4", "Denisova.DG_5", "Denisova.DG_6", "Denisova.DG_7", "Denisova.DG_8", "Denisova.DG_9", "Denisova.DG_10", "Chimp.REF_1", "Chimp.REF_2", "Chimp.REF_3", "Chimp.REF_4", "Chimp.REF_5", "Chimp.REF_6", "Chimp.REF_7", "Chimp.REF_8", "Chimp.REF_9", "Chimp.REF_10", "Vindija.DG_1", "Vindija.DG_2", "Vindija.DG_3", "Vindija.DG_4", "Vindija.DG_5", "Vindija.DG_6", "Vindija.DG_7", "Vindija.DG_8", "Vindija.DG_9", "Vindija.DG_10", "Altai_Neanderthal.DG_1", "Altai_Neanderthal.DG_2", "Altai_Neanderthal.DG_3", "Altai_Neanderthal.DG_4", "Altai_Neanderthal.DG_5", "Altai_Neanderthal.DG_6", "Altai_Neanderthal.DG_7", "Altai_Neanderthal.DG_8", "Altai_Neanderthal.DG_9", "Altai_Neanderthal.DG_10", "Mbuti.DG_1", "Mbuti.DG_2", "Mbuti.DG_3", "Mbuti.DG_4", "Mbuti.DG_5", "Mbuti.DG_6", "Mbuti.DG_7", "Mbuti.DG_8", "Mbuti.DG_9", "Mbuti.DG_10", "Switzerland_Bichon.SG_1", "Switzerland_Bichon.SG_2", "Switzerland_Bichon.SG_3", "Switzerland_Bichon.SG_4", "Switzerland_Bichon.SG_5", "Switzerland_Bichon.SG_6", "Switzerland_Bichon.SG_7", "Switzerland_Bichon.SG_8", "Switzerland_Bichon.SG_9", "Switzerland_Bichon.SG_10", "Russia_Ust_Ishim.DG_1", "Russia_Ust_Ishim.DG_2", "Russia_Ust_Ishim.DG_3", "Russia_Ust_Ishim.DG_4", "Russia_Ust_Ishim.DG_5", "Russia_Ust_Ishim.DG_6", "Russia_Ust_Ishim.DG_7", "Russia_Ust_Ishim.DG_8", "Russia_Ust_Ishim.DG_9", "Russia_Ust_Ishim.DG_10"]
samples = [
	msprime.Sample(2, 0/gen), # Denisova.DG_1
	msprime.Sample(2, 0/gen), # Denisova.DG_1
	msprime.Sample(2, 0/gen), # Denisova.DG_2
	msprime.Sample(2, 0/gen), # Denisova.DG_2
	msprime.Sample(2, 0/gen), # Denisova.DG_3
	msprime.Sample(2, 0/gen), # Denisova.DG_3
	msprime.Sample(2, 0/gen), # Denisova.DG_4
	msprime.Sample(2, 0/gen), # Denisova.DG_4
	msprime.Sample(2, 0/gen), # Denisova.DG_5
	msprime.Sample(2, 0/gen), # Denisova.DG_5
	msprime.Sample(2, 0/gen), # Denisova.DG_6
	msprime.Sample(2, 0/gen), # Denisova.DG_6
	msprime.Sample(2, 0/gen), # Denisova.DG_7
	msprime.Sample(2, 0/gen), # Denisova.DG_7
	msprime.Sample(2, 0/gen), # Denisova.DG_8
	msprime.Sample(2, 0/gen), # Denisova.DG_8
	msprime.Sample(2, 0/gen), # Denisova.DG_9
	msprime.Sample(2, 0/gen), # Denisova.DG_9
	msprime.Sample(2, 0/gen), # Denisova.DG_10
	msprime.Sample(2, 0/gen), # Denisova.DG_10
	msprime.Sample(6, 0/gen), # Chimp.REF_1
	msprime.Sample(6, 0/gen), # Chimp.REF_1
	msprime.Sample(6, 0/gen), # Chimp.REF_2
	msprime.Sample(6, 0/gen), # Chimp.REF_2
	msprime.Sample(6, 0/gen), # Chimp.REF_3
	msprime.Sample(6, 0/gen), # Chimp.REF_3
	msprime.Sample(6, 0/gen), # Chimp.REF_4
	msprime.Sample(6, 0/gen), # Chimp.REF_4
	msprime.Sample(6, 0/gen), # Chimp.REF_5
	msprime.Sample(6, 0/gen), # Chimp.REF_5
	msprime.Sample(6, 0/gen), # Chimp.REF_6
	msprime.Sample(6, 0/gen), # Chimp.REF_6
	msprime.Sample(6, 0/gen), # Chimp.REF_7
	msprime.Sample(6, 0/gen), # Chimp.REF_7
	msprime.Sample(6, 0/gen), # Chimp.REF_8
	msprime.Sample(6, 0/gen), # Chimp.REF_8
	msprime.Sample(6, 0/gen), # Chimp.REF_9
	msprime.Sample(6, 0/gen), # Chimp.REF_9
	msprime.Sample(6, 0/gen), # Chimp.REF_10
	msprime.Sample(6, 0/gen), # Chimp.REF_10
	msprime.Sample(10, 0/gen), # Vindija.DG_1
	msprime.Sample(10, 0/gen), # Vindija.DG_1
	msprime.Sample(10, 0/gen), # Vindija.DG_2
	msprime.Sample(10, 0/gen), # Vindija.DG_2
	msprime.Sample(10, 0/gen), # Vindija.DG_3
	msprime.Sample(10, 0/gen), # Vindija.DG_3
	msprime.Sample(10, 0/gen), # Vindija.DG_4
	msprime.Sample(10, 0/gen), # Vindija.DG_4
	msprime.Sample(10, 0/gen), # Vindija.DG_5
	msprime.Sample(10, 0/gen), # Vindija.DG_5
	msprime.Sample(10, 0/gen), # Vindija.DG_6
	msprime.Sample(10, 0/gen), # Vindija.DG_6
	msprime.Sample(10, 0/gen), # Vindija.DG_7
	msprime.Sample(10, 0/gen), # Vindija.DG_7
	msprime.Sample(10, 0/gen), # Vindija.DG_8
	msprime.Sample(10, 0/gen), # Vindija.DG_8
	msprime.Sample(10, 0/gen), # Vindija.DG_9
	msprime.Sample(10, 0/gen), # Vindija.DG_9
	msprime.Sample(10, 0/gen), # Vindija.DG_10
	msprime.Sample(10, 0/gen), # Vindija.DG_10
	msprime.Sample(11, 0/gen), # Altai_Neanderthal.DG_1
	msprime.Sample(11, 0/gen), # Altai_Neanderthal.DG_1
	msprime.Sample(11, 0/gen), # Altai_Neanderthal.DG_2
	msprime.Sample(11, 0/gen), # Altai_Neanderthal.DG_2
	msprime.Sample(11, 0/gen), # Altai_Neanderthal.DG_3
	msprime.Sample(11, 0/gen), # Altai_Neanderthal.DG_3
	msprime.Sample(11, 0/gen), # Altai_Neanderthal.DG_4
	msprime.Sample(11, 0/gen), # Altai_Neanderthal.DG_4
	msprime.Sample(11, 0/gen), # Altai_Neanderthal.DG_5
	msprime.Sample(11, 0/gen), # Altai_Neanderthal.DG_5
	msprime.Sample(11, 0/gen), # Altai_Neanderthal.DG_6
	msprime.Sample(11, 0/gen), # Altai_Neanderthal.DG_6
	msprime.Sample(11, 0/gen), # Altai_Neanderthal.DG_7
	msprime.Sample(11, 0/gen), # Altai_Neanderthal.DG_7
	msprime.Sample(11, 0/gen), # Altai_Neanderthal.DG_8
	msprime.Sample(11, 0/gen), # Altai_Neanderthal.DG_8
	msprime.Sample(11, 0/gen), # Altai_Neanderthal.DG_9
	msprime.Sample(11, 0/gen), # Altai_Neanderthal.DG_9
	msprime.Sample(11, 0/gen), # Altai_Neanderthal.DG_10
	msprime.Sample(11, 0/gen), # Altai_Neanderthal.DG_10
	msprime.Sample(13, 0/gen), # Mbuti.DG_1
	msprime.Sample(13, 0/gen), # Mbuti.DG_1
	msprime.Sample(13, 0/gen), # Mbuti.DG_2
	msprime.Sample(13, 0/gen), # Mbuti.DG_2
	msprime.Sample(13, 0/gen), # Mbuti.DG_3
	msprime.Sample(13, 0/gen), # Mbuti.DG_3
	msprime.Sample(13, 0/gen), # Mbuti.DG_4
	msprime.Sample(13, 0/gen), # Mbuti.DG_4
	msprime.Sample(13, 0/gen), # Mbuti.DG_5
	msprime.Sample(13, 0/gen), # Mbuti.DG_5
	msprime.Sample(13, 0/gen), # Mbuti.DG_6
	msprime.Sample(13, 0/gen), # Mbuti.DG_6
	msprime.Sample(13, 0/gen), # Mbuti.DG_7
	msprime.Sample(13, 0/gen), # Mbuti.DG_7
	msprime.Sample(13, 0/gen), # Mbuti.DG_8
	msprime.Sample(13, 0/gen), # Mbuti.DG_8
	msprime.Sample(13, 0/gen), # Mbuti.DG_9
	msprime.Sample(13, 0/gen), # Mbuti.DG_9
	msprime.Sample(13, 0/gen), # Mbuti.DG_10
	msprime.Sample(13, 0/gen), # Mbuti.DG_10
	msprime.Sample(17, 0/gen), # Switzerland_Bichon.SG_1
	msprime.Sample(17, 0/gen), # Switzerland_Bichon.SG_1
	msprime.Sample(17, 0/gen), # Switzerland_Bichon.SG_2
	msprime.Sample(17, 0/gen), # Switzerland_Bichon.SG_2
	msprime.Sample(17, 0/gen), # Switzerland_Bichon.SG_3
	msprime.Sample(17, 0/gen), # Switzerland_Bichon.SG_3
	msprime.Sample(17, 0/gen), # Switzerland_Bichon.SG_4
	msprime.Sample(17, 0/gen), # Switzerland_Bichon.SG_4
	msprime.Sample(17, 0/gen), # Switzerland_Bichon.SG_5
	msprime.Sample(17, 0/gen), # Switzerland_Bichon.SG_5
	msprime.Sample(17, 0/gen), # Switzerland_Bichon.SG_6
	msprime.Sample(17, 0/gen), # Switzerland_Bichon.SG_6
	msprime.Sample(17, 0/gen), # Switzerland_Bichon.SG_7
	msprime.Sample(17, 0/gen), # Switzerland_Bichon.SG_7
	msprime.Sample(17, 0/gen), # Switzerland_Bichon.SG_8
	msprime.Sample(17, 0/gen), # Switzerland_Bichon.SG_8
	msprime.Sample(17, 0/gen), # Switzerland_Bichon.SG_9
	msprime.Sample(17, 0/gen), # Switzerland_Bichon.SG_9
	msprime.Sample(17, 0/gen), # Switzerland_Bichon.SG_10
	msprime.Sample(17, 0/gen), # Switzerland_Bichon.SG_10
	msprime.Sample(18, 0/gen), # Russia_Ust_Ishim.DG_1
	msprime.Sample(18, 0/gen), # Russia_Ust_Ishim.DG_1
	msprime.Sample(18, 0/gen), # Russia_Ust_Ishim.DG_2
	msprime.Sample(18, 0/gen), # Russia_Ust_Ishim.DG_2
	msprime.Sample(18, 0/gen), # Russia_Ust_Ishim.DG_3
	msprime.Sample(18, 0/gen), # Russia_Ust_Ishim.DG_3
	msprime.Sample(18, 0/gen), # Russia_Ust_Ishim.DG_4
	msprime.Sample(18, 0/gen), # Russia_Ust_Ishim.DG_4
	msprime.Sample(18, 0/gen), # Russia_Ust_Ishim.DG_5
	msprime.Sample(18, 0/gen), # Russia_Ust_Ishim.DG_5
	msprime.Sample(18, 0/gen), # Russia_Ust_Ishim.DG_6
	msprime.Sample(18, 0/gen), # Russia_Ust_Ishim.DG_6
	msprime.Sample(18, 0/gen), # Russia_Ust_Ishim.DG_7
	msprime.Sample(18, 0/gen), # Russia_Ust_Ishim.DG_7
	msprime.Sample(18, 0/gen), # Russia_Ust_Ishim.DG_8
	msprime.Sample(18, 0/gen), # Russia_Ust_Ishim.DG_8
	msprime.Sample(18, 0/gen), # Russia_Ust_Ishim.DG_9
	msprime.Sample(18, 0/gen), # Russia_Ust_Ishim.DG_9
	msprime.Sample(18, 0/gen), # Russia_Ust_Ishim.DG_10
	msprime.Sample(18, 0/gen)  # Russia_Ust_Ishim.DG_10
]

events = [
	msprime.MassMigration(time = 0/gen, source = 17, destination = 16, proportion = 1),
	msprime.MassMigration(time = 10000/gen, source = 10, destination = 8, proportion = 1),
	msprime.MassMigration(time = 10000/gen, source = 11, destination = 8, proportion = 1),
	msprime.MassMigration(time = 10000/gen, source = 16, destination = 9, proportion = 0.5),
	msprime.MassMigration(time = 10000/gen, source = 16, destination = 15, proportion = 1),
	msprime.MassMigration(time = 20000/gen, source = 8, destination = 7, proportion = 1),
	msprime.MassMigration(time = 20000/gen, source = 9, destination = 7, proportion = 1),
	msprime.MassMigration(time = 30000/gen, source = 6, destination = 5, proportion = 1),
	msprime.MassMigration(time = 30000/gen, source = 7, destination = 5, proportion = 1),
	msprime.MassMigration(time = 40000/gen, source = 5, destination = 20, proportion = 1),
	msprime.MassMigration(time = 50000/gen, source = 20, destination = 4, proportion = 0.5),
	msprime.MassMigration(time = 50000/gen, source = 20, destination = 19, proportion = 1),
	msprime.MassMigration(time = 60000/gen, source = 18, destination = 14, proportion = 1),
	msprime.MassMigration(time = 60000/gen, source = 19, destination = 14, proportion = 1),
	msprime.MassMigration(time = 70000/gen, source = 13, destination = 12, proportion = 1),
	msprime.MassMigration(time = 70000/gen, source = 14, destination = 12, proportion = 1),
	msprime.MassMigration(time = 80000/gen, source = 12, destination = 3, proportion = 1),
	msprime.MassMigration(time = 80000/gen, source = 15, destination = 3, proportion = 1),
	msprime.MassMigration(time = 90000/gen, source = 3, destination = 1, proportion = 1),
	msprime.MassMigration(time = 90000/gen, source = 4, destination = 1, proportion = 1),
	msprime.MassMigration(time = 1e+05/gen, source = 1, destination = 0, proportion = 1),
	msprime.MassMigration(time = 1e+05/gen, source = 2, destination = 0, proportion = 1)
]

gt = numpy.zeros((int(1e+06),len(samples)/2), 'int')
for i in range(int(1e+06)):
	tree_sequence = msprime.simulate(population_configurations = pops, samples = samples, demographic_events = events, mutation_rate = 0.001)
	if tree_sequence.genotype_matrix().shape[0] > 0:
		gt[i,:] = (tree_sequence.genotype_matrix()[0,range(0,len(samples),2)] + tree_sequence.genotype_matrix()[0,range(1,len(samples),2)])
	else:
		gt[i,:] = numpy.zeros((1, len(samples)/2))

numpy.savetxt('msprime_sim.geno', gt, '%d', '')

with open('msprime_sim.snp', 'w') as f:
	for i in range(int(1e+06)):
		bytes = f.write('. 1 ' + str(i*50/float(999999)) + ' ' + str(i*100) + ' A C\n')

with open('msprime_sim.ind', 'w') as f:
	for i in range(len(indnam)):
		bytes = f.write(indnam[i] + ' U ' + indnam[i].rstrip("1234567890").rstrip("_") + '\n')
