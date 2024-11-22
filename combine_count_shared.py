# In this script, generate the counts.tsv file for mds. 
import os

FILES_DIR = 'results/default_reference/' #here, update path 

all_counts = []

filenames = ["4376", "4377", "4378", "4379", "4380", "4381", "4382", "4383", "4384", "4385", "4386", "4387", "4388", "4389", "4390", "4391", "4392", "4393", "4394", "4395", "4396", "4397", "4398", "4399", "4400", "4401", "4402", "4403", "4404", "4405", "4406", "4407", "4408", "4409", "4410", "4411"]

for file in filenames:
    with open(os.path.join(FILES_DIR, file + '_default_featurecounts.tsv'), 'r') as f:
        count = []
        for line in f:
            if not line.startswith('#'):
                if not line.startswith('Geneid'):
                    count.append(line.rstrip('\n').split('\t')[6])
        all_counts.append(count)

n_genes = len(all_counts[0])
print(n_genes)

outfile = open('results/default_reference/counts_all_samples.tsv', 'w') #here, update path

for i in range(n_genes):
    out = []
    for count in all_counts:
        out.append(count[i])
    print('\t'.join(out), file=outfile)