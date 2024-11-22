# In this script, we are updating the config json file (ipsc_config.jsob) in order to add the read group

import json
from collections import defaultdict
import os
import gzip

config_json = "ipsc_config_original.json" #here, need to update this to match with the name of the json config file
rg_config_json = "ipsc_config_rg.json" #here, need to update this to match with the name of the json config file


with open(config_json) as json_data:
    data = json.load(json_data)

    for sample in data["rna_samples"]:
        fq_path = os.path.join(data["out_directory"], "trimmed_fastqs_rna")
        read_group_info = {}
        fq_1 = os.path.join(data["out_directory"], "trimmed_fastqs_rna", sample + "_trimmed_read1.fastq.gz")
        fq_2 = os.path.join(data["out_directory"], "trimmed_fastqs_rna", sample + "_trimmed_read2.fastq.gz")

        # find pu
        with gzip.open(fq_1, 'rt', encoding="utf8", errors='ignore') as f:
            first_line = f.readline() #@LH00371:58:22KCTMLT3:6:1101:22867:3323 1:N:0:CACCTGTT
            items = first_line.split(':')
            pu = items[2] + '.' + items[3]

        read_group_info[sample] = {
            'fq_path': fq_path,
            'fq_1': fq_1,
            'fq_2': fq_2,
            'ID': sample,
            'SM': sample,
            'LB': sample,
            'PU': pu,
            'PL': 'Illumina'
        }
        data.update(read_group_info)

with open(rg_config_json, 'w') as outfile:
    json.dump(data, outfile, indent=4)
