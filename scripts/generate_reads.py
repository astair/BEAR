#!/usr/bin/env python


from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys, csv, io, random, decimal, argparse


parser = argparse.ArgumentParser(description='Generate uniform-length single or paired-end metagenomic reads.')
parser.add_argument('-r', metavar='<reference_fasta>', dest="ref", help="Multi-FASTA file containing genomic sequences from which reads will be sampled.")
parser.add_argument('-a', metavar='<abundance_file>', dest="abund", help="Tab-delimited abundance file with an abundance value for each corre- sponding  sequence in <reference fasta>")
parser.add_argument('-o', metavar='<output_file>', dest="output", help="Name for output file containing simulated uniform-length reads")
parser.add_argument('-t', metavar='<maximum_reads>', type=int, dest="total", default='100', help="The maximum number of reads to sample from all s")
parser.add_argument('-l', metavar='<longest_read>', type=int, dest="length", default='100', help="The length, in bp, of the longest possible read to simulate")
parser.add_argument('-i', metavar='<insert_mean_length>', type=int, dest="insert", default="0", help="Average length of insert for paired-end reads.")
parser.add_argument('-s', metavar='<insert_stddev>', type=int, dest="stddev", default="0", help="Standard deviation of insert length for paired-end reads" )
args = parser.parse_args()


#Reference meta database file (FASTA)
f1 = open(args.ref);

#abundance file (tab-delimited .txt)
f2 = open(args.abund);

total_reads = args.total

max_read_length = args.length

insert_avg = args.insert
insert_stddev = args.stddev

if(insert_avg):
	f4 = open(args.output + '.1.fasta', 'w')
	f5 = open(args.output + '.2.fasta', 'w')
else:
	f4 = open(args.output + '.fasta', 'w')

frags = []

div_file = csv.reader(f2, delimiter='\t')
species = []
abundance = []

freqs = []
comp = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'N', 'R':'Y', 'Y':'R'}
for row in div_file:
	species.append(row[0])
	abundance.append(decimal.Decimal(row[1]))

# Normalize abundances by sequence length (all scaffolds)
# --> abjusted abundance A_adj = (A * len_genome) / sum(A_i * gen_genome_i)
all_seqs = {}
lengths = [0] * len(species)
for seq in SeqIO.parse(f1, 'fasta'):
	all_seqs[seq.description] = seq.seq
	for sp_num in range(len(species)):
		if (species[sp_num] in seq.description):
			lengths[sp_num] += decimal.Decimal(len(seq.seq))
adj_lengths = [abundance[i] * lengths[i] for i in range(len(species))]
adj_abundance = [adj_lengths[i] / sum(adj_lengths) for i in range(len(species))]

# Here's to double checking
# print('\t'.join([str(float(i)) for i in lengths]))
# print('\t'.join([str(float(i)) for i in abundance]))
# print('\t'.join([str(float(i)) for i in adj_lengths]))
# print('\t'.join([str(float(i)) for i in adj_abundance]))

# Generate the reads
for (ID, seq) in all_seqs.items():
	sp_num = -1
	for spec in species:
		if spec in ID:
			sp_num = species.index(spec)
	if sp_num == -1:
		continue

	# Sequence abundance for each scaffold
	seq_div = decimal.Decimal(adj_abundance[sp_num] * len(seq) / lengths[sp_num])
	coverage = max(1, int(seq_div * total_reads))

	# Genreate reads w/ appropriate coverage
	limit = len(seq)
	for j in range(0, coverage):
		rand = random.random()
		rand_length = 0

		if ((insert_avg != 0) & (insert_stddev != 0)):
			cur_insert = int(random.gauss(insert_avg, insert_stddev))
			if (limit > (max_read_length * 2 + cur_insert)):
				start1 = random.randint(0, limit-(2*max_read_length + cur_insert))
				end1 = start1 + max_read_length
				start2 = end1 + cur_insert
				end2 = start2 + max_read_length
			else:
				start1 = 0
				end1 = limit
				start2 = 0
				end2 = limit
				
			read1 = seq[start1:end1]
			read2 = ''.join([comp[b] for b in seq[end2:start2:-1]])
			f4.write(">%s|#%d/1\n" % (ID, j+1))
			f4.write("%s\n" % read1)
			f5.write(">%s|#%d/2\n" % (ID, j+1))
			f5.write("%s\n" % read2)

		else:
			if (limit > max_read_length) :
				start = random.randint(0, limit-max_read_length)
				end = start+max_read_length
			else:
				start = 0
				end = limit
			read = seq[start:end]
			if random.random() < 0.5:
				#reverse orientation
				read = ''.join([comp[b] for b in seq[end:start:-1]])
				f4.write(">%s|#%d\n" % (ID, j+1))
				f4.write("%s\n" % read)
				
			else:
				#forward orientation
				f4.write(">%s|#%d\n" % (ID, j+1))
				f4.write("%s\n" % read)


f1.close()
f2.close()
f4.close()
