### Necessary Libraries

import sys,os
from operator import itemgetter
from Bio import SeqIO

### Functions

def BLAST(seq_file,blast_out): # blast reads of user-provided sequence file against pangenome
	os.system('time blastn -query '+seq_file+' -db pan_genome_reference -perc_identity 80 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qcov" -out '+blast_out+' -num_threads 8 -qcov_hsp_perc 100')

def process_blast_results(blast_out):
	blast_genes = {}
		
	for i in open(blast_out): # get pangenome genes with blast hits from blast file. 
		mapped_read,pan_gene = i.split('\t')[0],i.split('\t')[1]
		if pan_gene in blast_genes:
			blast_genes[pan_gene] = blast_genes[pan_gene] | set([mapped_read])
		else:
			blast_genes[pan_gene] = set([mapped_read])
		
	for k,v in blast_genes.items(): # Get number of unique reads that map to pangenome genes
		blast_genes[k] = len(v)
	
	return blast_genes
	
def Isolate_find(blast_genes,group2gene):

	isolate_vote = {}
	
	for h,i in enumerate(open('gene_presence_absence.Rtab')): # get number of reads that map to SRA isolates in pangenome
		if h == 0:
			isolates = [z for z in i.strip().split('\t')[1:]]
		else:
			i = i.strip().split('\t')
			if group2gene[i[0]] in blast_genes:
				for k,j in enumerate(i[1:]):
					if j == '1':
						if isolates[k] in isolate_vote:
							isolate_vote[isolates[k]] += 1
						else:
							isolate_vote[isolates[k]] = 1
	
	return isolate_vote
	
def writeout(isolate_results_file,isolate_vote):
	o = open(isolate_results_file,'w')

	for k,v in sorted(isolate_vote.items(),reverse=True,key=lambda x:x[1]):
		o.write(str(k)+'\t'+str(v)+'\n')

	o.close()


### Calls and Processes
	
seq_file = sys.argv[1]
blast_out = seq_file[:seq_file.find('.')]+'_blast.txt'
isolate_results_file = seq_file[:seq_file.find('.')]+'_isolate_counts.txt'

print "Processing pangenome file to get group and gene ids..."
group2gene = {i.strip().split(' ')[1]:i.strip().split(' ')[0].replace('>','') for i in open('pan_genome_reference.fa') if i.startswith('>')}
print "Group and gene ids obtained!!!"

print "Running blast. Blasting reads of sequence file against pangenome..."
BLAST(seq_file,blast_out)
print "Blast results obtained!!!"

print "Processing blast results..."
blast_genes = process_blast_results(blast_out)
print "Blast results processed!!!"

print "Running isolate finder program..."
isolate_vote = Isolate_find(blast_genes,group2gene)
print "Isolate finder program finished!!!"

print "Writing results to isolate_results.txt ..."
writeout(isolate_results_file,isolate_vote)
print "Results written to file!!!  Analysis finished!!!  Thank you!!!"

	
	