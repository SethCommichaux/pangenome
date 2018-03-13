import sys

try:
	sys.argv[1]
except IndexError:
	print
	print "Please add arguments: python get_divergence.py pangenome_file all_proteins_clustered_file"
	print
	sys.exit(1)

pangenome = {i.split("(")[0].replace('>',''):'' for i in open(sys.argv[1]) if i.startswith('>')}


furthest = {}
tmp = ''
divergence_max = 100.00

for i in open(sys.argv[2]):
	if i.startswith('>'):
		furthest[tmp] = divergence_max
		tmp = i.strip()
		divergence_max = 100.00
	else:
		i = i.strip().split('\t')[1]
		sra_gene = i.split('>')[1].split('...')[0]
		if sra_gene in pangenome:
			tmp = sra_gene
			if i.endswith('*'):
				continue
			furthest[sra_gene] = min(divergence_max,float(i.split('at ')[1].split('%')[0]))
		else:
			if i.endswith('*'):
				continue
			divergence_max = min(divergence_max,float(i.split('at ')[1].split('%')[0]))

o = open("max_divergence_by_gene_of_pangenome.txt",'w')

for k,v in sorted(furthest.items(),key=lambda x:x[1]):
	o.write(k+'\t'+str(v)+'\n')

o.close()

### was problem child: ERR251847_08357