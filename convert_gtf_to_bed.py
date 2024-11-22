import argparse

parser = argparse.ArgumentParser(description="Convert GTF file to bed file format.")
parser.add_argument("--input_gtf",required=True,help="Input the path to the GTF file")
parser.add_argument("--output_bed",required=True,help="Input the path to the output bed file")
parser.add_argument("--chromosome",required=True,help="Input the chromosome name. For example: chrX")

args = parser.parse_args()


outfile = open(args.output_bed, "w")

header = ["chr", "type", "start", "end", "gene_id", "gene_type", "gene_name"]
print ('\t'.join(header), file=outfile)

with open(args.input_gtf, "r") as f:
    for line in f:
        # if line.startswith(args.chromosome): #bug here because it will match chr1 and chr11
        items = line.rstrip('\n').split('\t')
        if items[0] == args.chromosome:
            # gene_name_list = items[8].split('; ')
            # for i in gene_name_list:
            #     if i.startswith('gene_name'):
            #         gene_name = i.split(" ")[1]
            # out = [items[0], str(int(items[3])-1), items[4], gene_name.strip("\"")]
            if items[2] == "exon":
                out = [items[0], items[2], items[3], items[4], "NA", "NA", "NA"]
                info_list = items[8].split("; ")
                for i in info_list:
                    if i.startswith('gene_id'):
                        out[4] = i.split(" ")[1].strip("\"")
                    if i.startswith('gene_type'):
                        out[5] = i.split(" ")[1].strip("\"")
                    if i.startswith('gene_name'):
                        out[6] = i.split(" ")[1].strip("\"")
                print ('\t'.join(out), file=outfile)
