# this script computes gene length

import argparse

parser = argparse.ArgumentParser(description="Compute gene length for TPM calculation.")
parser.add_argument("--input",required=True,help="Input. This is the output of the file convert_gtf_to_bed.py /gpfs/home6/tphung/gencode_v46/chr1_info.txt")
parser.add_argument("--unique_genes",required=True,help="Input. This is the unique genes /gpfs/home6/tphung/gencode_v46/chr1_info_uniqgenes.txt")
parser.add_argument("--output_bed",required=True,help="Input the path to the output bed file")

args = parser.parse_args()

chr_dict = {}
with open(args.unique_genes, "r") as f:
    for line in f: 
        if not line.startswith("gene_id"):
            chr_dict[line.rstrip("\n").split()[0]] = []

with open(args.input, "r") as f:
    for line in f:
        items = line.rstrip("\n").split("\t")
        if items[0] != "chr":
            gene = items[4]
            if gene in chr_dict:
                chr_dict[gene].append(int(items[2]))
                chr_dict[gene].append(int(items[3]))
            else: 
                print(gene, " is not in the unique list. Check!")

outfile = open(args.output_bed, "w")
header = ["gene", "start", "end", "length"]
print("\t".join(header), file=outfile)
for k, v in chr_dict.items():
    out = [k, str(min(v)), str(max(v)), str(max(v)-min(v)+1)]
    print("\t".join(out), file=outfile)

